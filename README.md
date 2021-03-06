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

Create and activate a Python3.8 virtual environment.


    $ python3 -m venv .venv
    $ source .venv/bin/activate


    $ pip3 install -r dependencies.txt

Install the following packages based on your OS.

- [`ATSAS`](https://www.embl-hamburg.de/biosaxs/download.html) 

- [`Molprobity`](https://github.com/rlabduke/MolProbity) 

- [`wkhtmltopdf`](https://wkhtmltopdf.org/) 

Create a local environment file and add the relevant variables.

    $ touch .env
    $ nano/vi .env

The variables to add to the .env file can be seen below (fill in the quotations with paths to the relevant values).

    ATSAS=""
    Molprobity_ramalyze=""
    Molprobity_molprobity=""
    Molprobity_clashscore=""
    Molprobity_rotalyze=""
    wkhtmltopdf=""

## Common errors in the installation process

1. One common error, depending on your OS and webdriver is from bokeh/selenium. This error is usually displayed as:

- `RuntimeError: Neither firefox and geckodriver nor a variant of chromium browser and chromedriver are available on system PATH. You can install the former with 'conda install -c conda-forge firefox geckodriver'.` 

This error originates from converting htmls to svgs. Please install/update your webdriver. You can do this by adding pre-installed binaries to path variable or install packages using the suggested conda command. 

2. Another potential error could arise from having another env variable file or not having the file in the same directory. This error is displayed as:

- `ValueError: UndefinedValueError('{} not found. Declare it as envvar or define a default value.'.format(option)).` 

A solution to this error is using Autoconfig of decouple library to add the path to your .env file. See this [`stackoverflow post`](https://stackoverflow.com/questions/43570838/how-do-you-use-python-decouple-to-load-a-env-file-outside-the-expected-paths) for specific details.


## Information

_Author(s)_: Sai J. Ganesan


