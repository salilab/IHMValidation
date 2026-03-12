# -*- coding: utf-8 -*-
#
# report.py - Generation PDF and HTML reports
#
# Copyright (C) 2019-2025 Arthur Zalevsky, Sai Ganesan, Benjamin M. Webb, Brinda Vallat
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
Generation of PDF and HTML reports
"""

import os
from pathlib import Path
import logging
from mmcif_io import GetInputInformation
import excludedvolume
import molprobity
import get_plots, sas, sas_plots
import utility
import pickle
import json
from multiprocessing import Manager
from collections import Counter
import numpy as np
import cx
import precision
import em

REPORT_VERSION = '3.2'

class WriteReport(object):
    report_version = REPORT_VERSION
    def __init__(self,
                 mmcif_file: str,
                 db: str='.',
                 cache: str='.',
                 nocache: bool=False,
                 enable_sas: bool=False,
                 enable_cx: bool=False,
                 enable_em: bool=False,
                 enable_prism: bool=False,
                 ):

        self.mmcif_file = mmcif_file
        self.cache = cache
        self.nocache = nocache
        self.enable_sas = enable_sas
        self.enable_cx = enable_cx
        self.enable_em = enable_em
        self.enable_prism = enable_prism

        # Read input
        self.input = GetInputInformation(mmcif_file=self.mmcif_file,
                                         cache=self.cache,
                                         nocache=self.nocache)
        self.system = self.input.system
        self.encoding = self.input.encoding

    def clean(self) -> None:
        '''cleanup'''
        pass

    def run_entry_composition(self, Template_Dict: dict) -> dict:
        '''
        get entry composition, relies on IHM library
        '''
        # check for ensembles
        if self.input.get_ensembles():
            ensemble_info = utility.dict_to_JSlist(self.input.get_ensembles())
        else:
            ensemble_info = None
        # from here on, we just fill the template dict with terms
        # we draw these terms from mmcif using python-ihm library
        # the terms are typically drawn out as dictionaries and then
        # converted to list of lists. These list of lists get fed into javascript
        # to print out tables. Why Sai?--because JS is annoying and it is just easier
        # to construct tables with lists than any other data struc
        # utility module has functions to check outputs from python-ihm library and convert to JS friendly format
        Template_Dict['report_version'] = self.report_version
        Template_Dict['python_ihm_version'] = utility.get_python_ihm_version()
        Template_Dict['ensemble_info'] = ensemble_info
        Template_Dict['sphere'] = self.input.check_sphere()
        Template_Dict['num_ensembles'] = self.input.check_ensembles()
        RB, flex, RB_nos, all_nos = self.input.get_RB_flex_dict()
        Template_Dict['Rigid_Body'] = RB_nos
        Template_Dict['Flexible_Unit'] = all_nos-RB_nos
        Template_Dict['RB_list'] = utility.dict_to_JSlist_rows(RB, flex)
        Template_Dict['RB'] = utility.get_RB(
            utility.dict_to_JSlist_rows(RB, flex))
        Template_Dict['flex'] = utility.get_flex(
            utility.dict_to_JSlist_rows(RB, flex))
        Template_Dict['ID'] = self.input.get_id()
        Template_Dict['ID_f'] = self.input.get_file_id()
        Template_Dict['PDB_ID'] = self.input.get_pdb_id()
        Template_Dict['PDBx_ID'] = self.input.get_pdbx_id()
        Template_Dict['PDBDEV_ID'] = self.input.get_pdb_dev_id()
        Template_Dict['ranked_id_list'] = self.input.get_ranked_id_list()
        Template_Dict['Molecule'] = self.input.get_struc_title()
        Template_Dict['Authors'] = self.input.get_authors()
        title, authors = self.input.get_primary_citation_info()
        Template_Dict['deposition_date'] = self.input.deposition_date
        # Template_Dict['Citation_Title'] = title
        # Template_Dict['Citation_Authors'] = authors
        Template_Dict['Entry_list'] = utility.dict_to_JSlist(
            self.input.get_composition())
        # Template_Dict['number_of_molecules'] = self.input.get_number_of_models()
        Template_Dict['number_of_models'] = self.input.get_number_of_models()
        Template_Dict['model_names'] = self.input.get_model_names()
        Template_Dict['number_of_software'] = self.input.get_software_length()
        Template_Dict['soft_list'] = utility.dict_to_JSlist(
            self.input.get_software_comp())
        Template_Dict['references'] = list(self.input.ref_cit.values())
        Template_Dict['references'].sort()
        Template_Dict['number_of_datasets'] = self.input.get_dataset_length()
        Template_Dict['cx_present'] = self.input.has_crosslinking_ms_dataset
        Template_Dict['sas_present'] = self.input.has_sas_dataset
        Template_Dict['em_present'] = self.input.has_em_dataset
        Template_Dict['Data'] = [i.upper() for i in list(set(self.input.get_dataset_comp(
        )['Dataset type']).difference({'Experimental model', 'Comparative model'}))]
        Template_Dict['Datasets_list'] = utility.dict_to_JSlist(
            self.input.get_dataset_comp())
        Template_Dict['Datasets_summary'] = utility.get_datasets_summary(self.system)
        Template_Dict['Unique_dataset'] = utility.get_unique_datasets(
            self.input.get_dataset_comp())
        Template_Dict['Protocols_number'] = self.input.get_protocol_number()
        Template_Dict['Sampling_list'] = utility.dict_to_JSlist(
            self.input.get_sampling())
        Template_Dict['num_chains'] = int(len(self.input.get_composition(
        )['Chain ID']))/int(len(list(Counter(self.input.get_composition()['Model ID']).keys())))
        Template_Dict['ChainL'] = self.input.get_composition()['Chain ID [auth]']
        Template_Dict['number_of_fits'] = 0
        Template_Dict['MAXPLOTS'] = get_plots.MAXPLOTS
        Template_Dict['rep_info'] = self.input.get_representation_info()
        return Template_Dict

    def run_model_quality(self, Template_Dict: dict, csvDirName: str, htmlDirName: str) -> (dict, dict, dict, dict, dict):
        '''
        get excluded volume for multiscale models
        get MolProbity info for atomic models
        exception: models with DNA--we need a way to assess models with DNA
        '''

        Template_Dict['disclaimer'] = 0
        Template_Dict['NumModels'] = self.input.num_models
        Template_Dict['atomic'] = False
        Template_Dict['assess_atomic_segments'] = None
        molprobity_dict = None
        Template_Dict['cg'] = False
        Template_Dict['assess_excluded_volume'] = None
        exv_data = None

        if self.input.atomic:
            Template_Dict['atomic'] = True
            # if there are no spheres, wed have atoms, so go ahead and set exv to 0/none
            # global clashscore; global rama; global sidechain;
            I_mp = molprobity.GetMolprobityInformation(mmcif_file=self.mmcif_file,
                                                       system=self.system,
                                                       encoding=self.encoding,
                                                       cache=self.cache,
                                                       nocache=self.nocache)
            molprobity_raw_data = I_mp.get_mp_data()
            Template_Dict['molprobity_models'] = I_mp.models
            Template_Dict['molprobity_data'] = I_mp.summarize_mp_data()
            Template_Dict['assess_atomic_segments'] = I_mp.get_summary_table_stats()
            molprobity_dict = I_mp.get_mq_plot_data()
            Template_Dict['molprobity_version'] = I_mp.molprobity_version

        # Run excluded volume for CG models or as a fall-back
        if self.input.cg or (Template_Dict['atomic'] and Template_Dict['molprobity_data'] is None):
            Template_Dict['cg'] = True
            # if there are no spheres, wed have atoms, so go ahead and set exv to 0/none
            logging.info("Getting excluded volume satisfaction")
            I_ev = excludedvolume.GetExcludedVolume(mmcif_file=self.mmcif_file,
                                                    system=self.system,
                                                    encoding=self.encoding,
                                                    cache=self.cache,
                                                    nocache=self.nocache)
            exv_data = I_ev.get_excluded_volume()
            viol_percent = np.asarray(exv_data['Excluded Volume Satisfaction (%)'], dtype=float)

            # let's update template dict with appropriate terms
            Template_Dict['excluded_volume_models'] = exv_data['Model ID']
            Template_Dict['excluded_volume'] = utility.dict_to_JSlist(exv_data)
            r_ = utility.format_range(viol_percent)
            Template_Dict['assess_excluded_volume'] = f'Satisfaction: {r_}%'

        # # we now set the disclaimer tag to see if there are issues while calculating exc vol
        # Template_Dict['disclaimer'] = 0
        # if exv_data:
        #     satisfaction = set(exv_data['Excluded Volume Satisfaction (%)'])
        #     violation = set(exv_data['Number of violations'])
        #     if len(satisfaction) == 1 and len(violation) == 1 and satisfaction == {'0.0'} and violation == {'0.0'}:
        #         Template_Dict['disclaimer'] = 1

        # Model has to be either one or both
        assert Template_Dict['atomic'] or Template_Dict['cg']

        return Template_Dict, molprobity_dict, exv_data

    def run_sas_validation(self, Template_Dict: dict) -> (dict, dict, dict):
        '''
        get sas validation information from SASCIF or JSON files
        '''
        # we start by checking if sas dataset was used to build model
        Template_Dict['sas'] = False

        if self.input.has_sas_dataset:
            Template_Dict['sas'] = True
            I_sas = sas.SasValidation(mmcif_file=self.mmcif_file,
                                      system=self.system,
                                      encoding=self.encoding,
                                      cache=self.cache,
                                      nocache=self.nocache)
            Template_Dict['freesas_version'] = I_sas.get_freesas_version()
            Template_Dict['p_val'] = utility.dict_to_JSlist(I_sas.get_pvals())
            Template_Dict['sasdb_code'] = I_sas.get_sas_ids()
            Template_Dict['sasdb_code_html'] = I_sas.get_sasbdb_ids()
            Template_Dict['sasdb_sascif'] = I_sas.check_sascif_dicts()
            Template_Dict['parameters_volume'] = utility.dict_to_JSlist(
                I_sas.get_parameters_vol_many())
            Template_Dict['parameters_mw'] = utility.dict_to_JSlist(
                I_sas.get_parameters_mw_many())
            Template_Dict['pddf_info'] = utility.dict_to_JSlist(
                I_sas.get_pddf_info())
            Template_Dict['number_of_fits'] = I_sas.get_total_number_of_fits()
            Template_Dict['rg_table'] = utility.dict_to_JSlist(
                I_sas.get_rg_table_many())
            Template_Dict['sasdb_code_fits'] = I_sas.get_sasdb_code_fits()
            Template_Dict['sas_data_quality'] = utility.get_rg_data(
                I_sas.get_rg_for_plot())
            Template_Dict['sas_fits_stats'] = utility.get_rg_data_fits(
                I_sas.get_fits_for_plot())

            sas_data = I_sas.get_rg_for_plot()
            sas_fit = I_sas.get_fits_for_plot()
        # if there are no sas datasets used to build the model, we set appropriate keys
        else:

            sas_data = {}
            sas_fit = {}
            Template_Dict['sasdb_sascif'] = []
        return Template_Dict, sas_data, sas_fit

    def run_sas_validation_plots(self, Template_Dict: dict, imageDirName: str):
        '''
        get sas validation information from SASCIF or JSON files
        '''
        # again, we start by checking for sas datasets
        if self.input.has_sas_dataset:
            Template_Dict['sas'] = ["True"]
            # I_sas = sas.SasValidation(self.mmcif_file)
            # create all relevant plots
            # try:
            I_sas_plt = sas_plots.SasValidationPlots(mmcif_file=self.mmcif_file,
                                                     system=self.system,
                                                     encoding=self.encoding,
                                                     cache=self.cache,
                                                     nocache=self.nocache,
                                                     imageDirName=imageDirName,
                                                     )
            I_sas_plt.plot_multiple()
            # I_sas.get_pofr_errors()
            I_sas_plt.plot_pf()
            I_sas_plt.plot_Guinier()
            if Template_Dict['number_of_fits'] > 0:
                I_sas_plt.plot_fits()
            # exception occurs if sascif not present
            # except (TypeError, KeyError, ValueError):
            #     pass

    def run_cx_validation(self, Template_Dict: dict) -> (dict):
        '''
        get cx validation information from mmcif files
        NOTE: this function is incomplete
        it currently evaluates satisfaction from mmcif files
        and not the enetire ensemble
        '''
        # if crosslinking-MS dataset was used to build the model, then evaluate satisfaction
        Template_Dict['cx'] = False
        Template_Dict['cx_stats'] = None
        Template_Dict['cx_ertypes'] = None
        Template_Dict['cx_num_of_restraints'] = None
        Template_Dict['cx_num_of_restraint_groups'] = None
        Template_Dict['cx_stats_per_model'] = None
        Template_Dict['cx_data_quality'] = None
        # Template_Dict['pyhmmer_version'] = None
        output = (Template_Dict, None, None)


        if self.input.has_crosslinking_ms_dataset:
            Template_Dict['cx'] = True
            I_cx = cx.CxValidation(mmcif_file=self.mmcif_file,
                                   system=self.system,
                                   encoding=self.encoding,
                                   cache=self.cache,
                                   nocache=self.nocache)
            self.I_cx = I_cx

            raw_data = None
            raw_ertypes = None
            ertypes, nr, nrg = None, None, None
            stats = None

            ertypes, data = I_cx.get_cx_data()

            if ertypes is not None:
                ertypes = I_cx.get_ertypes_df_html()
                nr = I_cx.get_number_of_restraints()
                nrg = I_cx.get_number_of_restraint_groups()
            if data is not None:
                stats = I_cx.get_stats_per_model_group()

            Template_Dict['cx_ertypes'] = ertypes
            Template_Dict['cx_num_of_restraints'] = nr
            Template_Dict['cx_num_of_restraint_groups'] = nrg
            Template_Dict['cx_stats'] = stats
            Template_Dict['cx_stats_per_model'] = I_cx.get_per_model_satifaction_rates()

            cx_data_quality = I_cx.validate_all_pride_data()
            Template_Dict['cx_data_quality'] = cx_data_quality
            if len(cx_data_quality) > 0:
                Template_Dict['pyhmmer_version'] = I_cx.get_pyhmmer_version()

            output = (Template_Dict, raw_data, raw_ertypes)

        return output


    def run_cx_validation_plots(self, Template_Dict: dict, imageDirName: str) -> None:
        '''
        create validation plots for cx datasets
        NOTE: this function is incomplete, the plots are also ugly
        and need to be refined
        '''
        if bool(Template_Dict['cx']):
            if Template_Dict['cx_stats'] is not None:

                html_fn, json_fn, svgs_fn = self.I_cx.plot_distograms_per_model_group(imageDirName)
                with open(json_fn, 'r') as f:
                    plot = json.dumps(json.load(f))
                Template_Dict['cx_distograms_plot_json'] = plot
                Template_Dict['cx_distograms_plots_svg'] = svgs_fn

                html_fn, json_fn, svgs_fn = self.I_cx.plot_satisfaction_per_ensemble(imageDirName)
                with open(json_fn, 'r') as f:
                    plot = json.dumps(json.load(f))
                Template_Dict['cx_satisfaction_plot_json'] = plot
                Template_Dict['cx_satisfaction_plots_svg'] = svgs_fn

    def run_quality_glance(self,
                           molprobity_dict: dict, exv_data: dict,
                           sas_data_quality: dict, sas_fit: dict,
                           cx_data_quality: dict, cx_fit: dict,
                           em_data_quality: dict, em_fit: dict,
                           imageDirName: str) -> dict:
        '''
        get quality at glance image; will be updated as validation report is updated
        '''
        I_plt = get_plots.Plots(mmcif_file=self.mmcif_file,
                                system=self.system,
                                encoding=self.encoding,
                                cache=self.cache,
                                nocache=self.nocache,
                                imageDirName=imageDirName,
                                )
        glance_plots = I_plt.plot_quality_at_glance(
            molprobity_dict, exv_data,
            sas_data_quality, sas_fit,
            cx_data_quality, cx_fit,
            em_data_quality, em_fit,
        )
        return glance_plots

    def run_supplementary_table(self,
                                Template_Dict,
                                location='N/A',
                                physics=['Information about physical principles was not provided'],
                                method_details='N/A',
                                sampling_validation=None,
                                validation_input=['-'],
                                cross_validation='N/A',
                                Data_quality=['-'],
                                clustering=None,
                                resolution='N/A'):
        '''
        get supplementary table, will be updated as validation report is updated
        '''
        # this again uses python-ihm to fill in template dict
        # the output from ihm is modified to fit appropriate formats
        if (self.input.get_ensembles() is not None) and (utility.all_same(self.input.get_ensembles()['Clustering method'])):
            Template_Dict['clustering'] = self.input.get_ensembles()[
                'Clustering method'][0]
        elif self.input.get_ensembles() is not None:
            Template_Dict['clustering'] = ', '.join(
                self.input.get_ensembles()['Clustering method'])
        else:
            Template_Dict['clustering'] = None
        Template_Dict['location'] = location
        Template_Dict['complex_name'] = self.input.get_struc_title()
        Template_Dict['Subunits'] = utility.get_subunits(
            self.input.get_composition())
        Template_Dict['datasets'] = utility.get_datasets(self.input.get_dataset_details(
        )) if self.input.get_dataset_details() is not None else 'Not provided or used'
        Template_Dict['physics'] = physics
        Template_Dict['software'] = utility.get_software(
            self.input.get_software_comp())
        Template_Dict['struc'] = self.input.get_atomic_coverage()
        Template_Dict['method_info'] = self.input.get_sampling()
        # TODO: Remove because deprecated
        # Template_Dict['method'] = utility.get_method_name(
        #     self.input.get_sampling())
        # Template_Dict['protocol_name'] = self.input.get_protocol_name()
        # Template_Dict['method_type'] = utility.get_method_type(
        #     self.input.get_sampling())
        # Template_Dict['method_details'] = utility.get_method_details(self.input.get_sampling())
        Template_Dict['models'] = ', '.join(self.input.get_ensembles(
        )['Number of models']) if self.input.get_ensembles() is not None else 'Not applicable'
        Template_Dict['sampling_validation'] = sampling_validation
        Template_Dict['feature'] = self.input.get_ensembles(
        )['Clustering feature'][0] if self.input.get_ensembles() is not None else 'Not applicable'
        Template_Dict['cross_validation'] = cross_validation
        model_precision = []

        edata = self.input.get_ensembles()
        if edata is not None:
            for p in edata['Cluster precision']:
                if p is None:
                    l = utility.NA
                else:
                    l = f'{p:.2f}, Å'
                model_precision.append(l)
        else:
            model_precision.append(utility.NA)

        Template_Dict['model_precision'] = model_precision

        Template_Dict['restraint_info'] = utility.get_restraints_info(self.input.get_restraints(
        )) if self.input.get_restraints() is not None else 'Not provided or used'

        dq = []

        if 'sas_data_quality' in Template_Dict and Template_Dict['sas_data_quality'] is not None:
            dq.extend(Template_Dict['sas_data_quality'])

        if 'cx_data_quality' in Template_Dict and Template_Dict['cx_data_quality'] is not None:
            for data_ in Template_Dict['cx_data_quality']:
                try:
                    r_ = float(data_['stats']['entry']['matched_pct'])
                    dq.append(f'{data_["pride_id"]}: {r_:.2f}% of crosslinks found in the data.')
                except (ValueError, TypeError, KeyError) as e:
                    pass

            for data_ in Template_Dict['cx_data_quality']:
                try:
                    r_ = float(data_['stats']['ms']['mapped_entities_pct'])
                    dq.append(f'{data_["pride_id"]}: {r_:.2f}% of crosslinks from the data were used for modeling.')
                except (ValueError, TypeError, KeyError) as e:
                    pass

        if 'em_data_quality' in Template_Dict and Template_Dict['em_data_quality'] is not None:
            for data_ in Template_Dict['em_data_quality']:
                try:
                    r_ = float(data_["data_stats"]["resolution"])
                    dq.append(f'{data_["emdbid"]}: resolution is {r_:.2f} Å')
                except (ValueError, TypeError, KeyError) as e:
                    dq.append(f'{data_["emdbid"]}: resolution is {utility.NA}')

        if len(dq) == 0:
            dq.append('Data quality has not been assessed')

        Template_Dict['Data_quality'] = dq

        fq = []
        if 'sas_fits_stats' in Template_Dict and Template_Dict['sas_fits_stats'] is not None:
            fq.extend(Template_Dict['sas_fits_stats'])

        if 'cx_stats_per_model' in Template_Dict and Template_Dict['cx_stats_per_model'] is not None:
            if Template_Dict['cx_stats_per_model']:
                r_ = utility.format_range(Template_Dict['cx_stats_per_model'])
                fq.append(f'Satisfaction of crosslinks: {r_}%')

        if 'em_data_quality' in Template_Dict and Template_Dict['em_data_quality'] is not None:
            qscores = []
            for data_ in Template_Dict['em_data_quality']:
                for mid, fit_stats in data_['fit_stats'].items():
                    try:
                        q = float(fit_stats['q_score']['average'])
                        qscores.append(q)
                    except (ValueError, TypeError, KeyError):
                        pass
            if len(qscores) > 0:
                r_ = utility.format_range(qscores)
                fq.append(f'3DEM q-score(s): {r_}')

            aiscores = []
            for data_ in Template_Dict['em_data_quality']:
                for mid, fit_stats in data_['fit_stats'].items():
                    try:
                        q = float(fit_stats['ai_score']['average'])
                        aiscores.append(q)
                    except (ValueError, TypeError, KeyError):
                        pass
            if len(aiscores) > 0:
                r_ = utility.format_range(aiscores)
                fq.append(f'3DEM atom inclusion score(s): {r_}')

        if len(fq) == 0:
            fq.append('Fit of model to information used to compute it has not been determined')
        Template_Dict['validation_input'] = fq

        scale = utility.pretty_print_representations(self.input.get_representation_details())
#        Template_Dict['clustering'] = clustering
        Template_Dict['summary_scale'] = scale
        Template_Dict['summary_entities'] = utility.summarize_entities(Template_Dict['rep_info'])
        Template_Dict['summary_segments'] = utility.summarize_segments(Template_Dict['rep_info'])

        return Template_Dict


    def run_prism(self, Template_Dict: dict, imageDirName: str) -> dict:
        I_p = precision.PRISM(mmcif_file=self.mmcif_file,
                              system=self.system,
                              encoding=self.encoding,
                              cache=self.cache,
                              nocache=self.nocache)
        Template_Dict['prism_data'] = I_p.get_data()
        Template_Dict['prism_plots'] = I_p.get_plots(imageDirName)

        if len(Template_Dict['prism_plots']) > 0:
            Template_Dict['prism_version'] = I_p.prism_version
            Template_Dict['pymol_version'] = I_p.pymol_version

    def run_em_validation(self, Template_Dict: dict, imageDirName: str) -> (dict):
        '''
        3DEM validation
        '''
        # if 3DEM dataset was used to build the model
        Template_Dict['em'] = False
        Template_Dict['em_stats'] = None
        Template_Dict['em_data_quality'] = None
        output = (Template_Dict, None, None)

        if self.enable_em and self.input.has_em_dataset:
            Template_Dict['em'] = True
            I_em = em.EMValidation(mmcif_file=self.mmcif_file,
                                   system=self.system,
                                   encoding=self.encoding,
                                   cache=self.cache,
                                   nocache=self.nocache)
            self.I_em = I_em

            em_data_quality = I_em.validate_all_emdb_data(imageDirName)
            Template_Dict['em_data_quality'] = em_data_quality
            if len(em_data_quality) > 0:
                Template_Dict['va_version'] = I_em.get_va_version()
                # Check if we have fit to data
                for dataset in em_data_quality:
                    if dataset['fit_stats'] is not None:
                        Template_Dict['chimera_version'] = I_em.get_chimera_version()
                        Template_Dict['chimerax_version'] = I_em.get_chimerax_version()
                        Template_Dict['mapq_version'] = I_em.get_mapq_version()
                        break

            output = (Template_Dict, None, None)
