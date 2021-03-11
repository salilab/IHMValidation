# Validation pipeline for integrative and hybrid models

## Project objective 
- `1.` Develop a validation pipeline (including data validation, model validation, fit of input to model, fit of data not used for modeling and uncertainty of the model) for assessing IHM structures deposited to [PDB-Dev](https://pdb-dev.wwpdb.org/index.html)

## List of files and directories:
- `docs` documentation for all classes and functions (sphinx)

- `src`  relevant source code for execution

- `static`  relevant css, js files for static and dynamic HTML reports along with output files from validation [htmls,images,pdfs,json,results,supplementary]

- `templates`  all HTML templates for static HTML reports, PDF files and supplementary table

- `templates_for_flask`  all HTML templates for dynamic HTML reports

- `tests`  tests for classes

- `example_sas`  example script to create validation reports for SAS data

- `example_imp_models`  example script to create validation reports for IMP models

- `example_summary_table`  example script to create validation reports for summary table 


## List of classes:

- `WriteReport`  class to write dictionary for jinja2, HTML and PDF outputs

- `plots`  quality at glance plots from all different analysis

- `get_excluded_volume` class to calculate excluded volume and other relevant statistics   

- `get_molprobity_information` class to get data from molprobity, used only for atomistic models

- `sas_validation` class to perform validation of models built using SAS datasets

- `sas_validation_plots` plots based on SAS validation analysis

- `cx_validation`  class to perform validation of models built using CX-MS datasets  

- `cx_validation_plots`  plots based on CX-MS validation analysis  

- `em_validation`  class to perform validation of models built using EM datasets

## Test site
- [`test web page`](https://modbase.compbio.ucsf.edu/pdbdev-test/home/) 

## Initial setup 

This initial setup is performed once.

Create and activate a Python3.6 virtual environment.

.. code-block:: RST

    $ python3 -m venv .venv
    $ source .venv/bin/activate

.. code-block:: RST

    $ pip3 install -r dependencies.txt

Create a local environment file and add the relevant variables.

.. code-block:: RST

    $ touch .env
    $ nano .env

The variables to add to the .env file can be seen below (fill in the quotations with paths to the relevant values).

.. code-block:: RST

    ATSAS=""
    MOLPROBITY=""
    PDFKIT=""

## Information

_Author(s)_: Sai J. Ganesan


