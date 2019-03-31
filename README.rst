PathwayForte |build| |docs|
===========================
A Python package for comparing the effect of pathway database choice in functional enrichment and classification
methods.

Structure
---------
- src: Source code of the Python package
- docs: Python package Documentation
- data: Directories containing gmt files, TCGA gene expression data, enrichment scores and test data used in the paper
  (in case proprocessing has been already conducted)
- R: Scripts to download and handle TCGA gene expression data

Installation
------------
1. ``pathway_forte`` can be installed with the following commands:

.. code-block:: sh

    $ python3 -m pip install git+https://github.com/pathwayforte/pathway-forte.git

2. or in editable mode with:

.. code-block:: sh

    $ git clone https://github.com/pathwayforte/pathway-forte.git
    $ cd pathway-forte
    $ python3 -m pip install -e .


.. |build| image:: https://travis-ci.com/pathwayforte/pathway-forte.svg?branch=master
    :target: https://travis-ci.com/pathwayforte/pathway-forte
    :alt: Build Status

.. |docs| image:: http://readthedocs.org/projects/pathwayforte/badge/?version=latest
    :target: https://pathwayforte.readthedocs.io/en/latest/
    :alt: Documentation Status
