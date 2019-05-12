PathwayForte |build| |docs| |coverage| |zenodo|
===============================================
A Python package for benchmarking pathway databases with functional enrichment and prediction methods
tasks.

Installation |pypi_version| |python_versions| |pypi_license|
------------------------------------------------------------
``pathway_forte`` can be installed from `PyPI <https://pypi.org/project/pathway-forte>`_
with the following command in your terminal:

.. code-block:: sh

    $ python3 -m pip install pathway_forte

The latest code can be installed from `GitHub <https://github.com/pathwayforte/pathway-forte>`_
with:

.. code-block:: sh

    $ python3 -m pip install git+https://github.com/pathwayforte/pathway-forte.git

For developers, the code can be installed with:

.. code-block:: sh

    $ git clone https://github.com/pathwayforte/pathway-forte.git
    $ cd pathway-forte
    $ python3 -m pip install -e .

Main Commands
-------------

The table below lists the main commands of PathwayForte.

+------------+--------------------------------+
| Command    | Action                         |
+============+================================+
| datasets   | Lists of Cancer Datasets       |
+------------+--------------------------------+
| export     | Export Gene Sets using ComPath |
+------------+--------------------------------+
| ora        | List of ORA Analyses           |
+------------+--------------------------------+
| fcs        | List of FCS Analyses           |
+------------+--------------------------------+
| prediction | List of Prediction Methods     |
+------------+--------------------------------+

More details of some commands:

1. **ora**. Over-Representation Analysis (e.g., one-tailed hyper-geometric test).

2. **fcs**. Functional Class Score Analaysis such as GSEA and ssGSEA using
   `GSEAPy <https://github.com/ostrokach/gseapy>`_.

3. **prediction**. Prediction methods include training elastic nets for binary classification, training SVMs for
   multi-classification tasks, or survival analysis.
   
Prediction Methods
------------------
``pathway_forte`` enables three classification methods (i.e., binary classification, training SVMs for multi-classification tasks, or survival analysis) using individualized pathway activity scores. The scores can be calculated from any pathway with a variety of tools (see [1]_) using any pathway database that enables to export its gene sets.

1. **binary**. Trains an elastic net model for a binary classification task (e.g., tumor vs. normal patients). The training is conducted using a nested cross validation approach (the number of cross validation in both loops can be selected). The model used can be easily changed since most of the models in `scikit-learn <https://scikit-learn.org/>`_ (the machine learning library used by this package) required the same input.

2. **subtype**. Trains a SVM model for a multi-class classification task (e.g., predict tumor subtypes). The training is conducted using a nested cross validation approach (the number of cross validation in both loops can be selected). Similarly as the previous classification task, other models can quickly be implemented.

3. **survival**. Trains a Cox's proportional hazard's model with elastic net penalty. The training is conducted using a nested cross validation approach with a grid search in the inner loop. This analysis requires pathway activity scores, patient classes and lifetime patient information.
   

References
----------

.. [1] Lim, S., *et al.* (2018). `Comprehensive and critical evaluation of individualized pathway activity measurement tools on pan-cancer data <https://doi.org/10.1093/bib/bby097>`_. * Briefings in bioinformatics*, bby125.
    

.. |build| image:: https://travis-ci.com/pathwayforte/pathway-forte.svg?branch=master
    :target: https://travis-ci.com/pathwayforte/pathway-forte
    :alt: Build Status

.. |docs| image:: http://readthedocs.org/projects/pathwayforte/badge/?version=latest
    :target: https://pathwayforte.readthedocs.io/en/latest/
    :alt: Documentation Status

.. |coverage| image:: https://codecov.io/gh/pathwayforte/pathway-forte/coverage.svg?branch=master
    :target: https://codecov.io/gh/pathwayforte/pathway-forte?branch=master
    :alt: Coverage Status

.. |python_versions| image:: https://img.shields.io/pypi/pyversions/pathway_forte.svg
    :target: https://pypi.org/project/pathway-forte
    :alt: Stable Supported Python Versions

.. |pypi_version| image:: https://img.shields.io/pypi/v/pathway_forte.svg
    :target: https://pypi.org/project/pathway-forte
    :alt: Current version on PyPI

.. |pypi_license| image:: https://img.shields.io/pypi/l/pathway_forte.svg
    :target: https://github.com/pathwayforte/pathway-forte/blob/master/LICENSE
    :alt: Apache-2.0

.. |zenodo| image:: https://zenodo.org/badge/178654585.svg
    :target: https://zenodo.org/badge/latestdoi/178654585
