# IHMValidation

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) ![version](https://img.shields.io/github/v/release/salilab/IHMValidation) [![docs](https://app.readthedocs.org/projects/ihmvalidation/badge/?version=latest&style=flat-default)](https://ihmvalidation.readthedocs.io/en/latest/)

IHMValidation is a software pipeline that follows guidelines and recommendations from IHM TaskForce [(Berman et al. 2019)](https://pubmed.ncbi.nlm.nih.gov/31780431/) for assesment of integrative biomolecular structures. The current version of the PDB-IHM validation report consists of six sections: (i) overview; (ii) model details; (iii) data quality assessments; (iv) local geometry assessments (i.e., model quality); (v) fit of the model to the data used to generate it; and (vi) fit of the model to the data used for validation. The sixth category, fit to the data used to validate the model is under development. Data quality assessments and fit to the data used for modeling sections are dependent on the different types of experimental data used in integrative modeling. We are breaking down this section based on experimental data type and addressing each method separately. The current version of the validation report is focused on validating models built using Small Angle Scattering (SAS), Chemical Crosslinking Mass Spectrometry (crosslinking-MS), and 3D Electron Microscopy (3DEM) data and is based on the model and data validation guidelines published by the [wwPDB SAS validation task force](https://www.wwpdb.org/task/sas) [(Trewhella et al., 2017)](https://pubmed.ncbi.nlm.nih.gov/28876235/), crosslinking-MS [(Leitner et al., 2020)](https://pubmed.ncbi.nlm.nih.gov/33065067/), and 3DEM [(Kleywegt et al., 2024)](https://pubmed.ncbi.nlm.nih.gov/38358351/) communities.

Validation of models built using Förster Resonance Energy Transfer (FRET) is under development and will be included in subsequent versions of the validation report.

Please refer to the [documentation](https://ihmvalidation.readthedocs.io/en/latest/) for details on running the pipeline and interpreting the validation reports. 

**References**

1.  Berman, Helen M., Paul D. Adams, Alexandre Bonvin, Stephen K. Burley, Bridget Carragher, Wah Chiu, Frank DiMaio, et al. 2019. "Federating Structural Models and Data: Outcomes from A Workshop on Archiving Integrative Structures." Structure 27 (12): 1745–59.
    
2.  Trewhella, Jill, Anthony P. Duff, Dominique Durand, Frank Gabel, J. Mitchell Guss, Wayne A. Hendrickson, Greg L. Hura, et al. 2017. "2017 Publication Guidelines for Structural Modelling of Small-Angle Scattering Data from Biomolecules in Solution: An Update." Acta Crystallographica. Section D, Structural Biology 73 (Pt 9): 710–28.
    
3.  Leitner, Alexander, Alexandre Bonvin, et al. 2020. "Toward Increased Reliability, Transparency, and Accessibility in Cross-linking Mass Spectrometry." Structure 28 (11): 1259-1268.
    
4.  Kleywegt, Gerard J., Paul D. Adams, Sarah J. Butcher, Catherine L. Lawson, et al. 2024. "Community recommendations on cryoEM data archiving and validation." IUCrJ, 11, 140-151.
