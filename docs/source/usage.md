# Usage

## PDB-IHM Validation System

The easiest way to run PBD-IHM Validation pipeline it to use standalone validation server [https://validate.pdb-ihm.org](https://validate.pdb-ihm.org). Register an account using [Globus](https://app.globus.org/groups/99da042e-64a6-11ea-ad5f-0ef992ed7ca1/about), upload your [IHMCIF](https://pubmed.ncbi.nlm.nih.gov/38508301/) file, and in 5-115 minutes you will get both full and summary validation reports in the PDF format.

## Local deployment

```{include} ../../singularity/README.md
```

### Running the pipeline

A typical IHMValidation run looks like this:

```
singularity run --pid ihmv.sif --cache-root cache --output-root . -f input/8ZZU.cif
```
where
* `--pid` ensures that in case of error or termination all subprocesses will be killed;
* `ihmv.sif` name of the apptainer image with the IHMValidation pipeline;
* `--cache-root cache` path to a cache directory. Some validation steps, like MolProbity assesment
or assessment of model fit to 3DEM density can be resource and time consuming, thus we recommend
preserving these intermidiate results in a dedicated directory. Otherwise you can point to `/tmp` or current directory `.`;
* `--output-root .` path to the directory where output will be created. `.` means the current directory. The actual output directory will have the same stem as the input file; 
* `-f 8ZZU.cif` path to the input file in the [IHMCIF](https://pubmed.ncbi.nlm.nih.gov/38508301/) format. Entry [8ZZU](https://files.rcsb.org/view/8ZZU.cif) is provided only as an example.

The output will have a following structure:
* `8ZZU` - directory with the results;
* `8ZZU/8ZZU_full_validation.pdf` - full validation report in the PDF format;  
* `8ZZU/8ZZU_summary_validation.pdf` - summary table in the PDF format;
* `8ZZU/8ZZU_html.tar.gz` - compressed archive with the full validation report in HTML format. **NB:** If you're planning to review HTML version of the report, you should add `--html-mode local` flag to the pipeline call, this will replace links to [PDB-IHM](https://pdb-ihm.org) resources with local copies.

To get a full list of available options you can call the pipeline with a `-h` flag:
```
singularity run --pid ihmv.sif -h
```

## Interpreting the reports

If you're not an expert in structral biology or integrative modeling, we recommend checking a detailed [user guide](https://pdb-ihm.org/validation_help.html) with descriptions of all metrics and plots presented in the report. 
