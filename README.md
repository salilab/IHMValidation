# Validation pipeline for integrative and hybrid models

## Project objective 
- `1.` Develop a validation pipeline (including data validation, model validation, fit of input to model, fit of data not used for modeling and uncertainty of the model) for assessing IHM structures deposited to [PDB-IHM](https://pdb-ihm.org/)

## List of files and directories:
- `docs` documentation for all classes and functions (sphinx)

- `src`  relevant source code for execution

- `static`  relevant css, js files for static and dynamic HTML reports along with output files from validation [htmls,images,pdfs,json,results,supplementary]

- `templates`  all HTML templates for static HTML reports, PDF files and supplementary table

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

## Docker image

Rather than installing and configuring all dependencies (below), you can
build a Docker (or podman) container by

1. Downloading the [`ATSAS`](https://www.embl-hamburg.de/biosaxs/download.html)    CentOS 7 RPM and placing it in the `docker` subdirectory (this cannot be
   redistributed by us; you must sign up at their web site for an academic
   license).
2. Building the image with `docker build -t ihm-validation docker`

The resulting image has all dependencies in the default PATH and this
repository available in the `/IHMValidation` directory; no further configuration
should be necessary.

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

Few pointers:
1. ATSAS variable should contain the path to datcmp functionality, example : `ATSAS-3.0.3-1/bin/datcmp`

2. Molprobity variables should point to respective functionalities, example : `build/bin/molprobity.ramalyze`

3. wkhtmltopdf variable should point to the binary

## Common errors in the installation process

1. One common error, depending on your OS and webdriver is from bokeh/selenium. This error is usually displayed as:

- `RuntimeError: Neither firefox and geckodriver nor a variant of chromium browser and chromedriver are available on system PATH. You can install the former with 'conda install -c conda-forge firefox geckodriver'.` 

This error originates from converting htmls to svgs. Please install/update your webdriver. You can do this by adding pre-installed binaries to path variable or install packages using the suggested conda command. 

To add pre-installed binaries (firefox and geckodriver), find the path of the binaries. Please try using conda first, add paths to binaries only if you are unable to use conda. You can do that in two steps:

- You can use the command `which firefox` or `which geckodriver` to get path to respective binaries. If you don't have the pre-installed binaries, install fireforx/geckodriver using brew and then locate binaries. You can also download geckodriver from the [`mozilla github page`](https://github.com/mozilla/geckodriver/releases). 

- You should then open the webdriver.py file in bokeh and add the path to appropriate functions. Open the webdriver.py file using the following path `.venv/lib/python3.8/site-packages/bokeh/io/webdriver.py`. Edit `create_firefox_webdriver` function by changing the variables for `firefox` to the binary path (delete which firefox) and change `geckodriver` to it's pre-installed binary path (delete which geckodriver). 

2. Another potential error could arise from having another env variable file or not having the file in the same directory. This error is displayed as:

- `ValueError: UndefinedValueError('{} not found. Declare it as envvar or define a default value.'.format(option)).` 

A solution to this error is using Autoconfig of decouple library to add the path to your .env file. See this [`stackoverflow post`](https://stackoverflow.com/questions/43570838/how-do-you-use-python-decouple-to-load-a-env-file-outside-the-expected-paths) for specific details.

3. An error could arise from not being able to access the executable, even though the path is found. This error can occur with the ATSAS package and is displayed as: 

- `PermissionError: [Errno 13] Permission denied: '/ATSAS-3.0.3-1'` 

A solution to this is to open your .venv/bin/activate file and add in the above six variables at the top as 6 lines of code using the format `export KEY=VALUE`. See this [`stackoverflow post`](https://stackoverflow.com/questions/9554087/setting-an-environment-variable-in-virtualenv) for specific details.

## Running an example

After the initial setup, you can start executing the scripts to generate validation reports. Here are the steps: 

- Go to the `example` directory.
- Command to execute: `python Execute.py -f PDBDEV_00000009.cif`
- If new software was used to build the structure, update reference.csv, located in the templates directory with appropriate software name, PubMed link, and citation.

The input to the `Execute.py` script is a PDBDEV file in cif format. The output includes directories and files that are listed below:

- `Directory  ../Validation/PDBDEV_00000009`
- `Directory  ../Validation/PDBDEV_00000009/images` 
- `Directory  ../Validation/PDBDEV_00000009/pdf  `
- `Directory  ../Validation/PDBDEV_00000009/json `
- `Directory  ../Validation/PDBDEV_00000009/supplementary  `
- `Directory  ../Validation/PDBDEV_00000009/htmls  `
- `Directory  ../Validation/PDBDEV_00000009/csv `

Here's the description of all the directories:

- `Validation` is the main head directory and is located one step above the example directory [you can see that in the way this repo is structured]
- `PDBDEV_00000009` is the entry directory with relevant files 
- `PDBDEV_00000009/images` contains all the images generated for this entry
- `PDBDEV_00000009/pdf` contains all the pdf files generated for this entry, including pdf version of the validation report 
- `PDBDEV_00000009/json` contains the validation report in a json format as key-value pairs 
- `PDBDEV_00000009/supplementary` contains the summary table in a PDF format
- `PDBDEV_00000009/htmls` contains corresponding html pages 
- `PDBDEV_00000009/csv` contains detailed molprobity tables for download 

## Transferring files to the server

The `Validation` folder that is generated needs to be transferred to the server.
If individual entries are being evaluated, move/copy the entry directory from the local validation directory into the server's validation directory.

## Information

_Author(s)_: Sai J. Ganesan


