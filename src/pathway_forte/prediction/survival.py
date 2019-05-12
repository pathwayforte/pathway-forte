# -*- coding: utf-8 -*-

"""Survival analysis module.

Uses `scikit-survival <https://github.com/sebp/scikit-survival>`_, which can
be installed with ``pip install scikit-survival`` but is imported with
``import sksurv``.
"""

import logging
import os

import numpy as np
import pandas as pd
from sklearn.model_selection import GridSearchCV, KFold
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.metrics import concordance_index_censored
from tqdm import tqdm

from pathway_forte.constants import (
    CANCER_DATA_SETS, CLINICAL_DATA, NORMAL_EXPRESSION_SAMPLES, PATHWAY_RESOURCES, RESULTS,
)
from pathway_forte.utils import get_num_samples

logger = logging.getLogger(__name__)


def prepare_ssgsea_data_for_survival_analysis(
        enrichment_score_path: str,
        clinical_data_path: str,
        normal_sample_size: int,
):
    """Prepare data for input into survival analysis.

    :param enrichment_score_path: ssgsea normalized enrichment scores file
    :param clinical_data_path: dataFrame of survival status and time to death if death occurred
    :param normal_sample_size: sample size
    :return: dataFrame of pathway scores for each sample, array of survival status and time to death info
    """
    # Read csv files
    clinical_data_df = pd.read_csv(clinical_data_path, sep='\t')
    enrichment_score = pd.read_csv(enrichment_score_path, sep='\t', header=0)

    # Reset index to sample IDs
    clinical_data_df.set_index('Sample ID', inplace=True)

    # Transpose dataFrame to arrange columns as pathways and rows as genes
    enrichment_score_df = enrichment_score.transpose()

    # Set column index to the first row in the dataframe
    enrichment_score_df.columns = enrichment_score_df.iloc[0]

    # Remove the first row because it is already set as column index
    enrichment_score_df = enrichment_score_df.drop("Term|NES")

    # Reset index values
    filtered_ssgsea_df = enrichment_score_df.reset_index()

    # Filter out controls from dataFrame to include cases only
    filtered_ssgsea_df = filtered_ssgsea_df.truncate(before=normal_sample_size)

    # Reset index values
    filtered_ssgsea_df.reset_index(drop=True, inplace=True)

    survival_info = []

    clinical_data_patient_ids = clinical_data_df.index

    for index, ssgsea_row in filtered_ssgsea_df.iterrows():
        # Patient id is only the first 15 characters
        patient_id = ssgsea_row['index'][0:15]

        # Remove patient data if either patient enrichment scores or their clinical data is missing
        if patient_id not in clinical_data_patient_ids:
            filtered_ssgsea_df.drop(index, inplace=True)
            continue

        row = clinical_data_df.loc[patient_id]

        # Remove patient data if their survival information is missing
        if pd.isna(row["Overall Survival (Months)"]):
            filtered_ssgsea_df.drop(index, inplace=True)
            continue

        # Get survival status (i.e. living or deceased) and survival time
        survival_info.append(
            (row["Overall Survival Status"] == "LIVING",
             row["Overall Survival (Months)"])
        )

    dt = np.dtype([('status', 'bool'), ('days_to_death', 'float')])
    event_time_struc_array = np.array(survival_info, dtype=dt)

    # Remove extra column
    del filtered_ssgsea_df['index']

    return filtered_ssgsea_df, event_time_struc_array


def survival_data_to_csv(clinical_data: str, dataset):
    """Get survival data from clinical data file."""

    # Read clinical meta data file
    clinical_data_df = pd.read_csv(clinical_data, sep='\t')

    # Get relevant columns for survival information
    _columns = ['Days to Last Followup', 'Overall Survival (Months)', 'Overall Survival Status']
    clinical_data_df = clinical_data_df[_columns]

    # Convert survival months to days
    clinical_data_df['Survival (Days)'] = round(clinical_data_df['Overall Survival (Months)'] * 30.42, 2)

    # If patient is living, replace their overall survival time in months to NaN
    clinical_data_df.loc[clinical_data_df['Overall Survival Status'] == 'LIVING', 'Overall Survival (Months)'] = np.NaN

    # Rearrange columns
    cols = clinical_data_df.columns.tolist()
    cols = ['Days to Last Followup', 'Overall Survival (Months)', 'Survival (Days)', 'Overall Survival Status']

    clinical_data_df = clinical_data_df[cols]

    clinical_data_df.to_csv('{}_survival_data.tsv'.format(dataset), sep='\t')


def train_survival_model(
        x,
        y,
        *,
        outer_cv_splits,
        inner_cv_splits,
        param_grid,
):
    """Train survival model with ssGSEA normalized enrichment scores (NES) from TCGA expression data and cBioPortal
    survival data on patient survival status and survival times.

    :param pandas.core.frame.DataFrame x: dataFrame of ssGSEA NES where controls are filtered out, as are
     patients with missing enrichment scores or survival data
    :param numpy.ndarray y: Structured array A where binary survival status is first field and survival time is
     second field.
    :param int outer_cv_splits: number of folds to split data in train/test sets in outer cross validation loop
    :param int inner_cv_splits: number of folds to split data in train/test sets in inner cross validation loop
    :param dict param_grid: parameter types and values to try in grid search
    :return: concordance scores
    """
    concordance_scores = []

    kf = KFold(n_splits=outer_cv_splits, shuffle=True)
    inner_cv = KFold(n_splits=inner_cv_splits)

    iterator = tqdm(kf.split(x, y))

    # Iterator for each CV step in the outer loop
    for i, (train_indexes, test_indexes) in enumerate(iterator):
        # Slice main data frame to get the training and test data for this CV step
        x_train = x.iloc[train_indexes]
        x_test = x.iloc[test_indexes]
        y_train = np.asarray([y[train_index] for train_index in train_indexes])
        y_test = np.asarray([y[test_index] for test_index in test_indexes])

        # Instantiate Cox’s proportional hazard’s regression model with elastic net penalty
        coxnet = CoxnetSurvivalAnalysis()

        # Tune hyper-parameters (e.g., L1-ratio) of the estimator using grid search (Inner loop in the nested-CV)
        gcv = GridSearchCV(estimator=coxnet, param_grid=param_grid, cv=inner_cv, return_train_score=True)

        # Run grid search on training data
        gcv.fit(x_train, y_train)

        # Extract best model from the grid
        coxnet = gcv.best_estimator_

        # predict y using the best model from the grid
        prediction = coxnet.predict(x_test)

        # Evaluate the performance of the model during grid search using Harrell's concordance index
        # Note that the main data frame is sliced to use only the test data for this CV step
        cindex, concordant, discordant, tied_risk, tied_time = concordance_index_censored(
            [y[test_index]['status'] for test_index in test_indexes],  # The status array for test set
            [y[test_index]['days_to_death'] for test_index in test_indexes],  # The days to death for test set
            prediction,  # Prediction scores
        )

        # print C-Index and best parameter found in the grid search
        print('best c-index: {}'.format(cindex))
        print('best parameter: {}'.format(gcv.best_params_))

        concordance_scores.append({
            "c-index": cindex,
            "number of concordant pairs": concordant,
            "number of discordant pairs": discordant,
            "tied_risk": tied_risk,
            "tied_time": tied_time,
            "l1-ratio": gcv.best_estimator_.l1_ratio,
            "split": i,
        })

    # avg_c_index = np.average([
    #     iter_result["c-index"]
    #     for iter_result in concordance_scores
    # ])

    # print('Avg C-Index {}'.format(avg_c_index))
    print(concordance_scores)
    # return avg_c_index, concordance_scores
    return concordance_scores


def run_survival_all_datasets(
        path: str,
        outer_cv_splits: int,
        inner_cv_splits: int,
        param_grid,
):
    """Run survival analysis using one ssGSEA file (pathway database) in all cancer datasets."""
    it = _help_run_survival_all_datasets(
        directory=path,
        outer_cv_splits=outer_cv_splits,
        inner_cv_splits=inner_cv_splits,
        param_grid=param_grid,
    )
    return {
        cancer_data_set: result
        for cancer_data_set, result in it
    }


def _help_run_survival_all_datasets(
        directory: str,
        outer_cv_splits: int,
        inner_cv_splits: int,
        param_grid,
):
    # For each cancer dataset analyzed in the dataset
    for file in os.listdir(directory):
        # Skip all files that do not start with tsv
        if not file.endswith('tsv'):
            continue

        # If it is one of the msig file, a different file name coding applies for the ssGSEA files
        # -> (<database>_<msig>_<cancer_dataset>.tsv)
        if 'msig' in file:
            file_name = file[:-4].split('_')  # Remove '.tsv' and split by underscore
            pathway_resource = '{}_{}'.format(file_name[0], file_name[1])
            cancer_data_set = file_name[2]

        # Normal convention for ssGSEA files:
        # -> <database>__<cancer_dataset>.tsv
        else:
            pathway_resource, cancer_data_set = file.strip('.tsv').split('_')

            # Check that is a valid file
            if pathway_resource not in PATHWAY_RESOURCES or cancer_data_set not in CANCER_DATA_SETS:
                logger.warning(f'Skipping file {file}')
                continue

        logger.info('{}-{}'.format(cancer_data_set, pathway_resource))

        # Prepare data
        ssgsea_df, event_time_array = prepare_ssgsea_data_for_survival_analysis(
            os.path.join(directory, file),
            CLINICAL_DATA.format(cancer_data_set, cancer_data_set),
            get_num_samples(NORMAL_EXPRESSION_SAMPLES.format(cancer_data_set))
        )

        # Train model and store results to export it as a pickle
        result = train_survival_model(
            ssgsea_df,
            event_time_array,
            outer_cv_splits=outer_cv_splits,
            inner_cv_splits=inner_cv_splits,
            param_grid=param_grid,
        )

        yield cancer_data_set, result


def run_survival_on_dataset(
        ssgsea_file,
        outer_cv_splits: int,
        inner_cv_splits: int,
        param_grid,
):
    if 'msig' in ssgsea_file:
        file_name = ssgsea_file.strip('.tsv').split('_')
        pathway_resource = '{}_{}'.format(file_name[0], file_name[1])
        cancer_data_set = file_name[2]

    # Normal convention for ssGSEA files:
    # -> <database>__<cancer_dataset>.tsv
    else:
        pathway_resource, cancer_data_set = ssgsea_file.split('_')

        pathway_resource = pathway_resource.split('/')[-1]

        cancer_data_set = cancer_data_set.replace('.tsv', '')

        # Check that is a valid file
        if pathway_resource not in PATHWAY_RESOURCES or cancer_data_set not in CANCER_DATA_SETS:
            logger.warning('Skipping file {}'.format(ssgsea_file))

    logger.info('{}-{}'.format(cancer_data_set, pathway_resource))

    # Prepare data
    ssgsea_df, event_time_array = prepare_ssgsea_data_for_survival_analysis(
        ssgsea_file,
        CLINICAL_DATA.format(cancer_data_set, cancer_data_set),
        get_num_samples(NORMAL_EXPRESSION_SAMPLES.format(cancer_data_set))
    )

    # Train model and store results to export it as a pickle
    return train_survival_model(
        ssgsea_df,
        event_time_array,
        outer_cv_splits=outer_cv_splits,
        inner_cv_splits=inner_cv_splits,
        param_grid=param_grid,
    )


def iterator_input_folders(
        main_directory,
        outer_cv_splits,
        inner_cv_splits,
        param_grid,
):
    """Iterate through all database folder in main folder."""
    return {
        pathway_resource: run_survival_all_datasets(
            os.path.join(RESULTS, main_directory, pathway_resource),
            outer_cv_splits,
            inner_cv_splits,
            param_grid,
        )
        for pathway_resource in PATHWAY_RESOURCES
    }
