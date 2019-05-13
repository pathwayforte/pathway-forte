# -*- coding: utf-8 -*-

""""""

import logging
import os
import json

from pathway_forte.constants import CLASSIFIER_RESULTS, SSGSEA
from pathway_forte.prediction.binary import get_l1_ratios
from pathway_forte.prediction.survival import run_survival_all_datasets

__all__ = [
    'do_survival_prediction',
]

logger = logging.getLogger(__name__)


def do_survival_prediction(
        data,
        *,
        outer_cv_splits,
        inner_cv_splits,
):
    # List of parameters to use in grid search
    param_grid = {
        'l1_ratio': get_l1_ratios(),
    }
    logger.info(f'Parameter grid: {param_grid}')

    # Run survival analysis in all datasets available
    results = run_survival_all_datasets(
        os.path.join(SSGSEA, data),
        outer_cv_splits=outer_cv_splits,
        inner_cv_splits=inner_cv_splits,
        param_grid=param_grid,
    )
    logger.info(results)

    results_path = os.path.join(CLASSIFIER_RESULTS, f'survival_results_{data}.json')

    # Export prediction metrics as a JSON
    with open(results_path, 'w') as file:
        json.dump(results, file, sort_keys=True, indent=2)

    logger.info(f'Done with survival analysis. Results are exported in {results_path}')
