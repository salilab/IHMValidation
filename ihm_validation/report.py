###################################
# Script :
# 1) Contains class to write dictionary
# for jinja2/write HTML files/PDF files
# and supplementary table
#
# ganesans - Salilab - UCSF
# ganesans@salilab.org
###################################

import os
from pathlib import Path
import logging
from mmcif_io import GetInputInformation
import excludedvolume
import molprobity
import get_plots, sas, sas_plots
import cx
import utility
import pickle
import json
from multiprocessing import Manager
from collections import Counter
import numpy as np

class WriteReport(object):
    def __init__(self, mmcif_file, db, driver, cache, nocache=False):
        self.mmcif_file = mmcif_file
        self.db = db
        self.input = GetInputInformation(self.mmcif_file)
        # Webdriver for figures
        self.driver = driver
        self.cache = cache
        self.nocache = nocache


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
        Template_Dict['ID_w'] = self.input.get_id().split()
        Template_Dict['ID_T'] = self.input.get_id()[0:6]+'_' + \
            self.input.get_id()[6:]
        Template_Dict['ID_MP'] = [str(self.input.get_id()[0:6]+'_' +
                                      self.input.get_id()[6:])]
        Template_Dict['ID_R'] = (
            self.input.get_id()[0:6]+'_'+self.input.get_id()[6:]).split()
        Template_Dict['Molecule'] = self.input.get_struc_title()
        Template_Dict['Authors'] = self.input.get_authors()
        title, authors = self.input.get_primary_citation_info()
        Template_Dict['Citation_Title'] = title
        Template_Dict['Citation_Authors'] = authors
        Template_Dict['Entry_list'] = utility.dict_to_JSlist(
            self.input.get_composition())
        Template_Dict['number_of_molecules'] = self.input.get_number_of_models()
        Template_Dict['number_of_models'] = self.input.get_number_of_models()
        Template_Dict['model_names'] = self.input.get_model_names()
        Template_Dict['number_of_software'] = self.input.get_software_length()
        Template_Dict['soft_list'] = utility.dict_to_JSlist(
            self.input.get_software_comp())
        Template_Dict['references'] = list(self.input.ref_cit.values())
        Template_Dict['references'].sort()
        Template_Dict['number_of_datasets'] = self.input.get_dataset_length()
        Template_Dict['Data'] = [i.upper() for i in list(set(self.input.get_dataset_comp(
        )['Dataset type']).difference({'Experimental model', 'Comparative model'}))]
        Template_Dict['Datasets_list'] = utility.dict_to_JSlist(
            self.input.get_dataset_comp())
        Template_Dict['Unique_dataset'] = utility.get_unique_datasets(
            self.input.get_dataset_comp())
        Template_Dict['Protocols_number'] = self.input.get_protocol_number()
        Template_Dict['Sampling_list'] = utility.dict_to_JSlist(
            self.input.get_sampling())
        Template_Dict['num_chains'] = int(len(self.input.get_composition(
        )['Chain ID']))/int(len(list(Counter(self.input.get_composition()['Model ID']).keys())))
        Template_Dict['ChainL'] = self.input.get_composition()['Chain ID [auth]']
        Template_Dict['ChainLMap'] = self.input.get_auth_label_map()
        Template_Dict['number_of_fits'] = 0
        return Template_Dict

    def check_mmcif(self, Template_Dict: dict) -> None:
        '''
        test function to check rewrite_mmcif function without any hassle.
        not useful for processing, useful only for testing.
        '''
        if self.input.check_sphere() < 1:
            self.input.rewrite_mmcif()
            # I_mp = molprobity.GetMolprobityInformation('test.cif')
            logging.info("File rewritten...")
            logging.info("Molprobity analysis is being calculated...")

    def check_disclaimer_warning(self, exv_data: dict, Template_Dict: dict) -> dict:
        '''
        set a disclaimer if model quality can not be determined
        this func is not used anymore, needs to be deleted
        '''
        # we set the disclaimer as false
        Template_Dict['disclaimer'] = False
        if exv_data:
            # check for exc vol exceptions
            satisfaction = set(exv_data['Excluded Volume Satisfaction'])
            violation = set(exv_data['Number of violations'])
            if len(satisfaction) == 1 and len(violation) == 1 and satisfaction == {'0.0'} and violation == {'0.0'}:
                Template_Dict['disclaimer'] = True
        return Template_Dict

    def run_model_quality(self, Template_Dict: dict, csvDirName: str, htmlDirName: str) -> (dict, dict, dict, dict, dict):
        '''
        get excluded volume for multiscale models
        get molprobity info for atomic models
        exception: models with DNA--we need a way to assess models with DNA
        '''
        Template_Dict['molprobity_version'] = None

        if self.input.check_sphere() < 1:
            # if there are no spheres, wed have atoms, so go ahead and set exv to 0/none
            # global clashscore; global rama; global sidechain;
            exv_data = None
            I_mp = molprobity.GetMolprobityInformation(self.mmcif_file,
                                                       cache = self.cache)
            Template_Dict['molprobity_version'] = I_mp.get_version()
            key = Path(self.mmcif_file).stem
            filename = str(Path(
                self.cache, key + '_temp_mp.pickle'))
            # check if molprobity for this entry has already been detetmined
            if os.path.exists(filename) and not self.nocache:
                d_mp = {}
                logging.info("Molprobity analysis file already exists...\n...assuming clashscores, \
                        Ramachandran and rotamer outliers have already been calculated")
                with open(filename, 'rb') as fp:
                    d_mp['molprobity'] = pickle.load(fp)

                f_rota = str(Path(
                    self.cache, key + '_temp_rota.pickle'))
                with open(f_rota, 'rb') as fp:
                    d_mp['rota'] = pickle.load(fp)

                f_rama = str(Path(
                    self.cache, key + '_temp_rama.pickle'))
                with open(f_rama, 'rb') as fp:
                    d_mp['rama'] = pickle.load(fp)

                f_clash = str(Path(
                    self.cache, key + '_temp_clash.pickle'))
                with open(f_clash, 'rb') as fp:
                    d_mp['clash'] = pickle.load(fp)

            else:
                # if molprobity for these entries have not yet been determined, go ahead and set them up to run
                # we rewrite all files into a format that is suitable for molprobity
                I_mp = molprobity.GetMolprobityInformation(self.mmcif_file, cache=self.cache)
                logging.info("Molprobity analysis is being calculated...")
                try:
                    manager = Manager()
                    d_mp = manager.dict()
                    # utility.runInParallel(
                    I_mp.run_clashscore(d_mp)
                    I_mp.run_ramalyze(d_mp)
                    I_mp.run_rotalyze(d_mp)
                    I_mp.run_molprobity(d_mp)
                    # Cleanup
                    I_mp.cleanup()

                # if by any chance the rewrite doens't help, and we are unable to run molprobity,
                # step out and print an error
                except (TypeError, KeyError, ValueError):
                    logging.error("Molprobity cannot be calculated...")

            # at this stage, we should have all our dictionary terms
            if d_mp:
                # get total number of bond and angle outliers (this is an approx number)
                bond_nos, angle_nos = I_mp.process_molprobity(
                    d_mp['molprobity'])

                # if we do have bond outliers, print all the ones, we will write this into a csv file
                if bond_nos:
                    bonddict, total_bonds = I_mp.process_bonds_list(
                        d_mp['molprobity'], Template_Dict['ChainL'], Template_Dict['ChainLMap'])
                    Template_Dict['total_bonds'] = total_bonds
                    Template_Dict['molp_b_csv'] = utility.dict_to_JSlist(
                        bonddict)
                else:
                    Template_Dict['total_bonds'] = 1

                # if we have angle outliers, print all the ones
                # angle outliers are a little tricky as molprobity outputs angle outliers multiple times
                # in different parts of the file
                angledict, total_angles = I_mp.process_angles_list(
                    d_mp['molprobity'], Template_Dict['ChainL'], Template_Dict['ChainLMap'])
                Template_Dict['total_angles'] = total_angles
                if angle_nos:
                    Template_Dict['molp_a_csv'] = utility.dict_to_JSlist(
                        I_mp.add_angles_outliers(angle_nos, angledict, Template_Dict['ChainL'], Template_Dict['ChainLMap']))
                else:
                    Template_Dict['total_angles'] = 1
                    Template_Dict['molp_a_csv'] = utility.dict_to_JSlist(
                        angledict)

                # we compute total number of bond and angle outliers
                try:
                    Template_Dict['angle'] = len(Template_Dict['molp_a_csv'])-1
                except KeyError:
                    Template_Dict['angle'] = 0
                    Template_Dict['molp_a_csv'] = [[]]
                try:
                    Template_Dict['bond'] = len(Template_Dict['molp_b_csv'])-1
                except KeyError:
                    Template_Dict['bond'] = 0
                    Template_Dict['molp_b_csv'] = [[]]

                # we write summary tables for bonds and angles to be printed to HTML and PDF reports
                Template_Dict['molp_b'] = utility.dict_to_JSlist(
                    I_mp.bond_summary_table(Template_Dict['molp_b_csv']))
                Template_Dict['molp_a'] = utility.dict_to_JSlist(
                    I_mp.angle_summary_table(Template_Dict['molp_a_csv']))

                # we write all the tables to csv and html files
                I_mp.write_table_csv(
                    Template_Dict['molp_a_csv'], csvDirName, table_filename='angle_outliers.csv')
                I_mp.write_table_csv(
                    Template_Dict['molp_b_csv'], csvDirName, table_filename='bond_outliers.csv')
                I_mp.write_table_html(
                    Template_Dict['molp_a_csv'], htmlDirName, table_filename='angle_outliers.html')
                I_mp.write_table_html(
                    Template_Dict['molp_b_csv'], htmlDirName, table_filename='bond_outliers.html')
                # print (Template_Dict['bond'],sum(I_mp.bond_summary_table(Template_Dict['molp_b_csv'])['Number of outliers']))
                # print (Template_Dict['angle'],sum(I_mp.angle_summary_table(Template_Dict['molp_a_csv'])['Number of outliers']))

                # after anle and bond outliers, we move onto processing rotamers, ramachandran outliers, and clashlists
                Template_Dict['rotascore'] = utility.dict_to_JSlist(
                    I_mp.rota_summary_table(I_mp.process_rota(d_mp['rota'])))
                Template_Dict['rotalist'] = utility.dict_to_JSlist(
                    I_mp.rota_detailed_table(I_mp.process_rota(d_mp['rota']), Template_Dict['ChainL'], Template_Dict['ChainLMap']))
                Template_Dict['ramascore'] = utility.dict_to_JSlist(
                    I_mp.rama_summary_table(I_mp.process_rama(d_mp['rama'])))
                Template_Dict['ramalist'] = utility.dict_to_JSlist(
                    I_mp.rama_detailed_table(I_mp.process_rama(d_mp['rama']), Template_Dict['ChainL'], Template_Dict['ChainLMap']))
                clashscores, Template_Dict['tot'] = I_mp.clash_summary_table(
                    d_mp['clash'])
                Template_Dict['clashscore_list'] = utility.dict_to_JSlist(
                    clashscores)
                Template_Dict['clashlist'] = utility.dict_to_JSlist(I_mp.clash_detailed_table(
                    d_mp['clash'], Template_Dict['ChainLMap']))
                Template_Dict['assess_excluded_volume'] = 'Not applicable'
                molprobity_dict = I_mp.get_data_for_quality_at_glance(Template_Dict['clashscore_list'],
                                                                      Template_Dict['rotascore'], Template_Dict['ramascore'])
                Template_Dict['assess_atomic_segments'] = utility.mp_readable_format(
                    molprobity_dict)
                Template_Dict['NumModels'] = len(molprobity_dict['Names'])

        else:
            Template_Dict['bond'] = Template_Dict['angle'] = 0
            Template_Dict['total_bonds'] = Template_Dict['total_angles'] = 1
            # if there are spheres, wed have coarse grained beadas, so we go ahead and calculate excluded volume
            # set the appropriate flag for assessing atomic segments
            Template_Dict['assess_atomic_segments'] = None
            # check if exv has already been evaluated
            file = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..',
                                                'Validation', 'results', str(Template_Dict['ID'])+'exv.txt'))
            if os.path.exists(file) and not self.nocache:
                logging.info("Excluded volume file already exists...")
                with open(file, 'r') as inf:
                    line = [ln.strip().replace('[', '').replace(']', '').replace('"', '').
                            replace(' ', '').split(',')[1:] for ln in inf.readlines()]
                exv_data = {
                    'Models': line[0], 'Excluded Volume Satisfaction (%)':
                    line[1], 'Number of violations': line[2]}
            else:
                logging.info("Excluded volume is being calculated...")
                I_ev = excludedvolume.GetExcludedVolume(self.mmcif_file, cache=self.cache)
                model_dict = I_ev.get_all_spheres()
                exv_data = I_ev.run_exc_vol_parallel(model_dict)

            Template_Dict['NumModels'] = len(exv_data['Models'])
            viol_percent = np.asarray(exv_data['Excluded Volume Satisfaction (%)'], dtype=float)
            min_viol_percent = min(viol_percent)
            max_viol_percent = max(viol_percent)

            # let's update template dict with appropriate terms
            Template_Dict['excluded_volume'] = utility.dict_to_JSlist(exv_data)
            # print (Template_Dict['excluded_volume'])
            # Template_Dict['assess_excluded_volume'] = utility.exv_readable_format(
            #    exv_data)
            Template_Dict['assess_excluded_volume'] = f'Satisfaction: {min_viol_percent:.2f}-{max_viol_percent:.2f}%'
            molprobity_dict = None

        # we now set the disclaimer tag to see if there are issues while calculating exc vol
        Template_Dict['disclaimer'] = 0
        if exv_data:
            satisfaction = set(exv_data['Excluded Volume Satisfaction (%)'])
            violation = set(exv_data['Number of violations'])
            if len(satisfaction) == 1 and len(violation) == 1 and satisfaction == {'0.0'} and violation == {'0.0'}:
                Template_Dict['disclaimer'] = 1
        return Template_Dict, molprobity_dict, exv_data

    def run_sas_validation(self, Template_Dict: dict) -> (dict, dict, dict):
        '''
        get sas validation information from SASCIF or JSON files
        '''
        # we start by checking if sas dataset was used to build model
        if self.input.check_for_sas(self.input.get_dataset_comp()):
            Template_Dict['sas'] = ["True"]
            I_sas = sas.SasValidation(self.mmcif_file, self.db)
            Template_Dict['atsas_version'] = I_sas.get_atsas_version()
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
            Template_Dict['Data_quality'] = utility.get_rg_data(
                I_sas.get_rg_for_plot())
            Template_Dict['validation_input'] = utility.get_rg_data_fits(
                I_sas.get_fits_for_plot())
            # we check if model fits have been deposited
            if len(Template_Dict['validation_input']) < 1:
                Template_Dict['validation_input'] = [
                    'Fit of model to data has not been deposited']
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
        if self.input.check_for_sas(self.input.get_dataset_comp()):
            Template_Dict['sas'] = ["True"]
            # I_sas = sas.SasValidation(self.mmcif_file)
            # create all relevant plots
            # try:
            I_sas_plt = sas_plots.SasValidationPlots(
                self.mmcif_file, imageDirName, self.driver)
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
        # if cx-ms dataset was used to build the model, then evaluate satisfaction
        Template_Dict['cx'] = False
        Template_Dict['cx_stats'] = None
        Template_Dict['cx_ertypes'] = None
        Template_Dict['cx_num_of_restraints'] = None
        Template_Dict['cx_num_of_restraint_groups'] = None
        Template_Dict['cx_stats_per_model'] = None
        output = (Template_Dict, None, None)

        if self.input.check_for_cx(self.input.get_dataset_comp()):
            Template_Dict['cx'] = True
            I_cx = cx.CxValidation(self.mmcif_file)
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

                html_fn, json_fn = self.I_cx.plot_distograms_per_model_group(imageDirName)
                with open(json_fn, 'r') as f:
                    plot = json.dumps(json.load(f))
                Template_Dict['cx_distograms_plot_json'] = plot

                html_fn, json_fn = self.I_cx.plot_satisfaction_per_ensemble(imageDirName)
                with open(json_fn, 'r') as f:
                    plot = json.dumps(json.load(f))
                Template_Dict['cx_satisfaction_plot_json'] = plot

    def run_quality_glance(self, molprobity_dict: dict, exv_data: dict,
                           sas_data: dict, sas_fit: dict,
                           cx_fit: dict, imageDirName: str):
        '''
        get quality at glance image; will be updated as validation report is updated
        '''
        I_plt = get_plots.Plots(self.mmcif_file, imageDirName, driver=self.driver)
        I_plt.plot_quality_at_glance(
            molprobity_dict, exv_data, sas_data, sas_fit, cx_fit)

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
        Template_Dict['PDB_ID'] = self.input.get_id()
        Template_Dict['Subunits'] = utility.get_subunits(
            self.input.get_composition())
        Template_Dict['datasets'] = utility.get_datasets(self.input.get_dataset_details(
        )) if self.input.get_dataset_details() is not None else 'Not provided or used'
        Template_Dict['physics'] = physics
        Template_Dict['software'] = utility.get_software(
            self.input.get_software_comp())
        Template_Dict['struc'] = self.input.get_atomic_coverage()
        Template_Dict['method'] = utility.get_method_name(
            self.input.get_sampling())
#        Template_Dict['protocol_name'] = self.input.get_protocol_name()
        Template_Dict['method_info'] = self.input.get_sampling()
        Template_Dict['method_type'] = utility.get_method_type(
            self.input.get_sampling())
        # Template_Dict['method_details'] = utility.get_method_details(self.input.get_sampling())
        Template_Dict['models'] = ', '.join(self.input.get_ensembles(
        )['Number of models']) if self.input.get_ensembles() is not None else 'Not applicable'
        Template_Dict['sampling_validation'] = sampling_validation
        Template_Dict['feature'] = self.input.get_ensembles(
        )['Clustering feature'][0] if self.input.get_ensembles() is not None else 'Not applicable'
        Template_Dict['cross_validation'] = cross_validation
        Template_Dict['model_precision'] = ', '.join([f'{i}, Å' for i in self.input.get_ensembles(
        )['Cluster precision']]) if self.input.get_ensembles() is not None else \
            'Model precision can not be calculated with one structure'
        Template_Dict['restraint_info'] = utility.get_restraints_info(self.input.get_restraints(
        )) if self.input.get_restraints() is not None else 'Not provided or used'
        if 'Data_quality' not in list(Template_Dict.keys()):
            Template_Dict['Data_quality'] = ['Data quality has not been assessed']
        # if 'validation_input' not in list(Template_Dict.keys()):
        #     Template_Dict['validation_input'] = validation_input

        validation_input = []
        if 'cx_stats_per_model' in Template_Dict:
            if Template_Dict['cx_stats_per_model']:
                min_cx = min(Template_Dict['cx_stats_per_model'])
                max_cx = max(Template_Dict['cx_stats_per_model'])
                validation_input.append(f'Satisfaction of crosslinks: {min_cx:.2f}-{max_cx:.2f}%')

        if len(validation_input) == 0:
            validation_input.append('Fit of model to information used to compute it has not been determined')
        Template_Dict['validation_input'] = validation_input

        resolution = pretty_print_representations(self.input.get_representation_details())
#        Template_Dict['clustering'] = clustering
        Template_Dict['resolution'] = resolution
        return Template_Dict


def pretty_print_representations(reprs: dict) -> str:
    out = ''
    if reprs['atomic'] is True:
        out += 'Atomic'

    if reprs['coarse-grained'] is True:
        if out != '':
            out += '; '

        out += 'Coarse-grained: '
        levels = [str(x) for x in reprs['coarse-grain_levels']]
        out += ', '.join(levels)
        out += ' residue(s) per bead'

    return out
