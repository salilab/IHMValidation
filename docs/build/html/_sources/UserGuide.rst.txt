User Guide
**********
User guide to the wwPDB integrative structure supplementary table

Table of contents
=================
        1. `Model composition`_
                a. `Entry composition`_
                b. `Datasets used for modeling`_
        2. `Representation`_
                a. `Atomic structural coverage`_
                b. `Rigid bodies`_
                c. `Flexible units`_
                d. `Interface units`_
                e. Resolution_
        3. Restraints_
                a. `Physical restraints`_
                b. `Experimental information`_
        4. Validation_
                a. `Sampling validation`_
                b. `Clustering algorithm`_
                c. `Clustering feature`_
                d. `Number of ensembles`_
                e. `Number of models in ensembles`_
                f. `Model precision`_
                g. `Quality of data`_
                h. `Model Quality`_
                i. `Assessment of atomic segments`_
                j. `Excluded volume satisfaction`_
                k. `Fit of the model to information used to compute it`_
                l. `Fit of the model to information not used to compute it`_
        5. `Methodology and software`_
                a. `Method name`_
                b. `Method details`_
                c. Software_

Introduction
=============
*The wwPDB integrative structure (IM) validation reports are prepared according to the recommendations of the wwPDB Integrative/Hybrid Methods (IHM) Task Force (cite). The report summarises the quality of the structure and includes measures of the quality of the data on which the structures were based, the standard criteria for assessing atomic models, the fit of a model to information used to compute it, the fit of a model to information not used to compute it, and uncertainty in the model. The summary table is an overview of the validation report.*

Definitions
===========
.. _`Model composition`:

1. *Model composition*: The total number of unique molecules that are present in the entry, and information used for modeling.

.. _`Entry composition`:

        a. *Entry composition*: List of unique molecules that are present in the entry.

.. _`Datasets used for modeling`:

        b. *Datasets used for modeling*: List of input datasets used for modeling.

.. _`Representation`:

2. *Representation*: Representation of modeled structure.

.. _`Atomic structural coverage`:

        a.*Atomic structural coverage*: Percentage of modeled structure or residues for which atomic structures are available. These structures can include X-ray, NMR, EM, and other comparative models.

.. _`Rigid bodies`:

        b. *Rigid bodies*: A rigid body consists of multiple coarse-grained (CG) beads or atomic residues. In a rigid body, the beads (or residues) have their relative distances constrained during conformational sampling.

.. _`Flexible units`:

        c. *Flexible units*: Flexible units consist of strings of beads that are restrained by the sequence connectivity.

.. _`Interface units`:

        d. *Interface units*: An automatic definition based on identified interface for each model. Applicable to models built with HADDOCK.

.. _Resolution:

        e. *Resolution*: Resolution of segments of modeled structure.

.. _Restraints:

3. *Restraints*: A set of restraints used to compute modeled structure.

.. _`Physical restraints`:

        a. *Physical restraints*: A list of restraints derived from physical principles to compute modeled structure.

.. _`Experimental information`:
        
        b. *Experimental information*: A list of restraints derived from experimental datasets to compute modeled structure.

.. _Validation:

4. *Validation*: Assessment of models based on validation criteria set by IHM task force (cite).

.. _`Sampling validation`:

        a. *Sampling validation*: Validation metrics used to assess sampling convergence for stochastic sampling. Sampling precision is defined as the largest allowed Root-mean-square deviation (RMSD) between the cluster centroid and a model within any cluster in the finest clustering for which each sample contributes structures proportionally to its size (considering both the significance and magnitude of the difference) and for which a sufficient proportion of all structures occur in sufficiently large clusters.See (cite) for more details.

.. _`Clustering algorithm`:

        b. *Clustering algorithm*: Clustering algorithm used to analyze resulting solution. See (cite) for more details.
 
.. _`Clustering feature`:

        c. *Clustering feature*: Feature or reaction co-ordinate used to cluster solution.

.. _`Number of ensembles`:
        
        d. *Number of ensembles*: Number of solutions or ensembles of modeled structure.

.. _`Number of models in ensembles`:

        e. *Number of models in ensemble(s)*: Number of structures in the solution ensemble(s).        

.. _`Model precision`:

        f. *Model precision*: Measurement of variation among the models in the ensemble upon a global least-squares superposition. See (cite) for more details.

.. _`Quality of data`:

        g. *Quality of data* : Assessment of data on which modeled structures are based. See user guide on IM validation report.

.. _`Model Quality`:

        h. *Model quality* : Assessment of modeled structures based on physical principles

.. _`Assessment of atomic segments`:

        i. *Assessment of atomic segments* : Assessment of atomic segments in the integrative structure. Quality statistics in this section are calculated using standard compilations of covalent geometry parameters (Engh & Huber, 2001; Parkinson et al., 1996), tools in MolProbity (Chen et al., 2010), Validation-pack (Feng et al.) and the wwPDB chemical component dictionary (CCD). See user guide on PDB validation report.

.. _`Excluded volume satisfaction`:

        j. *Excluded volume satisfaction* : Assessment of excluded volume satisfaction of CG beads in the modeled structure. Excluded volume between two beads not connected in sequence are satisfied if the distance between them is greater than that of the sum of their radii.

.. _`Fit of the model to information used to compute it`:

        k. *Fit of the model to information used to compute it* : Assessment of modeled structure based on data used for modeling. See user guide on IM validation report.

.. _`Fit of the model to information not used to compute it`:

        l. *Fit of the model to information not used to compute it* : Assessment of modeled structure based on data not used for modeling. See user guide on IM validation report.


.. _`Methodology and software`:

5. *Methodology and software*: List of methods on which modeled structures are based and software used to obtain structures.

.. _`Method name`:

        a. *Method name* : Name(s) of method(s) used to generate modeled structures.

.. _`Method details`:
 
        b. *Method details* : Details of method(s) used to generate modeled structures.

.. _Software:
        
        c. *Software* : Software used to compute modeled structure, also includes scripts used to generate and analyze models.



*Updated on October 25th, 2019*


