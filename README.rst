PathwayForte |build| |docs| |zenodo|
====================================
A Python package for benchmarking pathway databases in functional enrichment and prediction methods
tasks.


Installation |pypi_version| |python_versions| |pypi_license|
------------------------------------------------------------
1. ``pathway_forte`` can be installed with the following commands:

.. code-block:: sh

    $ python3 -m pip install pathway_forte

2. or in editable mode with:

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

1. **ora**. Over-Representation Analysis (e.g., one-tail hyper-geometric test).

2. **.fcs**.. Functional Class Score Analaysis such as GSEA and ssGSEA using `GSEAPy <https://github.com/ostrokach/gseapy>`_.

3. **.prediction**.. Prediction methods include training elastic nets for binary classification, training SVMs for multi-classification tasks, or survival analysis.


.. |build| image:: https://travis-ci.com/pathwayforte/pathway-forte.svg?branch=master
    :target: https://travis-ci.com/pathwayforte/pathway-forte
    :alt: Build Status

.. |docs| image:: http://readthedocs.org/projects/pathwayforte/badge/?version=latest
    :target: https://pathwayforte.readthedocs.io/en/latest/
    :alt: Documentation Status

.. |coverage| image:: https://codecov.io/gh/pathwayforte/pathway-forte/coverage.svg?branch=master
    :target: https://codecov.io/gh/pathwayforte/pathway-forte?branch=master
    :alt: Coverage Status

.. |climate| image:: https://codeclimate.com/github/pathwayforte/pathway-forte/badges/gpa.svg
    :target: https://codeclimate.com/github/pathwayforte/pathway-forte
    :alt: Code Climate

.. |python_versions| image:: https://img.shields.io/pypi/pyversions/pathway_forte.svg
    :alt: Stable Supported Python Versions

.. |pypi_version| image:: https://img.shields.io/pypi/v/pathway_forte.svg
    :alt: Current version on PyPI

.. |pypi_license| image:: https://img.shields.io/pypi/l/pathway_forte.svg
    :alt: Apache-2.0

.. |zenodo| image:: https://zenodo.org/badge/178654585.svg
    :target: https://zenodo.org/badge/latestdoi/178654585



