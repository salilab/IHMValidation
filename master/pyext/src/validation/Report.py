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
import validation
from validation import excludedvolume, GetInputInformation
from validation import molprobity
from validation import get_plots, sas
from validation import cx
from validation import utility
import pickle
from multiprocessing import Manager
from collections import Counter


class WriteReport(object):
    def __init__(self, mmcif_file):
        self.mmcif_file = mmcif_file
        self.Input = GetInputInformation(self.mmcif_file)

    def run_entry_composition(self, Template_Dict: dict) -> dict:
        '''
        get entry composition, relies on IHM library
        '''
        if self.Input.get_ensembles():
            ensemble_info = utility.dict_to_JSlist(self.Input.get_ensembles())
        else:
            ensemble_info = None
        Template_Dict['ensemble_info'] = ensemble_info
        Template_Dict['sphere'] = self.Input.check_sphere()
        Template_Dict['num_ensembles'] = self.Input.check_ensembles()
        RB, flex, RB_nos, all_nos = self.Input.get_RB_flex_dict()
        Template_Dict['Rigid_Body'] = RB_nos
        Template_Dict['Flexible_Unit'] = all_nos-RB_nos
        Template_Dict['RB_list'] = utility.dict_to_JSlist_rows(RB, flex)
        Template_Dict['RB'] = utility.get_RB(
            utility.dict_to_JSlist_rows(RB, flex))
        Template_Dict['flex'] = utility.get_flex(
            utility.dict_to_JSlist_rows(RB, flex))
        Template_Dict['ID'] = self.Input.get_id()
        Template_Dict['ID_w'] = self.Input.get_id().split()
        Template_Dict['ID_T'] = self.Input.get_id()[0:6]+'_'+self.Input.get_id()[6:]
        Template_Dict['ID_R'] = (
            self.Input.get_id()[0:6]+'_'+self.Input.get_id()[6:]).split()
        Template_Dict['Molecule'] = self.Input.get_struc_title()
        Template_Dict['Title'] = self.Input.get_title()
        Template_Dict['Authors'] = self.Input.get_authors()
        Template_Dict['Entry_list'] = utility.dict_to_JSlist(
            self.Input.get_composition())
        Template_Dict['number_of_molecules'] = self.Input.get_number_of_models()
        Template_Dict['model_names'] = self.Input.get_model_names()
        Template_Dict['number_of_software'] = self.Input.get_software_length()
        Template_Dict['soft_list'] = utility.dict_to_JSlist(
            self.Input.get_software_comp())
        Template_Dict['number_of_datasets'] = self.Input.get_dataset_length()
        Template_Dict['Data'] = [i.upper() for i in list(set(self.Input.get_dataset_comp(
        )['Dataset type']).difference({'Experimental model', 'Comparative model'}))]
        Template_Dict['Datasets_list'] = utility.dict_to_JSlist(
            self.Input.get_dataset_comp())
        Template_Dict['Protocols_number'] = self.Input.get_protocol_number()
        Template_Dict['Sampling_list'] = utility.dict_to_JSlist(
            self.Input.get_sampling())
        Template_Dict['num_chains'] = int(len(self.Input.get_composition(
        )['Chain ID']))/int(len(list(Counter(self.Input.get_composition()['Model ID']).keys())))
        return Template_Dict

    def run_model_quality(self, Template_Dict: dict) -> (dict, dict, dict, dict, dict):
        '''
        get excluded volume for multiscale models
        get molprobity info for atomic models
        exception: models with DNA--we need a way to assess models with DNA
        '''
        if self.Input.check_sphere() < 1:
            # global clashscore; global rama; global sidechain;
            exv_data = None
            I_mp = molprobity.GetMolprobityInformation(self.mmcif_file)
            if I_mp.check_for_molprobity():
                filename = os.path.abspath(os.path.join(
                    os.getcwd(), 'static/results/', str(Template_Dict['ID'])+'_temp_mp.txt'))
                if os.path.exists(filename):
                    d_mp = {}
                    print("Molprobity analysis file already exists...\n...assuming clashscores, \
                        Ramachandran and rotamer outliers have already been calculated")
                    with open(filename, 'rb') as fp:
                        d_mp['molprobity'] = pickle.load(fp)
                    f_rota = os.path.abspath(os.path.join(
                        os.getcwd(), 'static/results/', str(Template_Dict['ID'])+'_temp_rota.txt'))
                    with open(f_rota, 'rb') as fp:
                        d_mp['rota'] = pickle.load(fp)
                    f_rama = os.path.abspath(os.path.join(
                        os.getcwd(), 'static/results/', str(Template_Dict['ID'])+'_temp_rama.txt'))
                    with open(f_rama, 'rb') as fp:
                        d_mp['rama'] = pickle.load(fp)
                    f_clash = os.path.abspath(os.path.join(
                        os.getcwd(), 'static/results/', str(Template_Dict['ID'])+'_temp_clash.txt'))
                    with open(f_clash, 'rb') as fp:
                        d_mp['clash'] = pickle.load(fp)
                else:
                    print("Molprobity analysis is being calculated...")
                    manager = Manager()
                    d_mp = manager.dict()
                    utility.runInParallel(I_mp.run_clashscore(d_mp), I_mp.run_ramalyze(
                        d_mp), I_mp.run_rotalyze(d_mp), I_mp.run_molprobity(d_mp))
                a, b = I_mp.process_molprobity(d_mp['molprobity'])
                Template_Dict['bond'] = len(a)
                Template_Dict['angle'] = len(b)
                global clashscore
                global rama
                global sidechain
                clashscore, rama, sidechain = I_mp.get_data_for_quality_at_glance(
                    d_mp['molprobity'])
                Template_Dict['molp_b'] = utility.dict_to_JSlist(
                    I_mp.molprobity_detailed_table_bonds(a))
                Template_Dict['molp_a'] = utility.dict_to_JSlist(
                    I_mp.molprobity_detailed_table_angles(b))
                Template_Dict['rotascore'] = utility.dict_to_JSlist(
                    I_mp.rota_summary_table(I_mp.process_rota(d_mp['rota'])))
                Template_Dict['rotalist'] = utility.dict_to_JSlist(
                    I_mp.rota_detailed_table(I_mp.process_rota(d_mp['rota'])))
                Template_Dict['ramascore'] = utility.dict_to_JSlist(
                    I_mp.rama_summary_table(I_mp.process_rama(d_mp['rama'])))
                Template_Dict['ramalist'] = utility.dict_to_JSlist(
                    I_mp.rama_detailed_table(I_mp.process_rama(d_mp['rama'])))
                clashscores, Template_Dict['tot'] = I_mp.clash_summary_table(
                    d_mp['clash'])
                Template_Dict['clashscore_list'] = utility.dict_to_JSlist(
                    clashscores)
                Template_Dict['clashlist'] = I_mp.clash_detailed_table(
                    d_mp['clash'])
                Template_Dict['assess_atomic_segments'] = 'Clashscore: ' + str(
                    clashscore) + ', Ramachandran outliers: ' + str(rama) + '% '+',\
                     Sidechain outliers: '+str(sidechain)+'%'
                Template_Dict['assess_excluded_volume'] = ['Not applicable']
            else:
                if not I_mp.check_for_molprobity():
                    self.Input.rewrite_mmcif()
                    I_mp = molprobity.get_molprobity_information('test.cif')
                    print("file rewritten")
                if I_mp.check_for_molprobity():
                    print("Molprobity analysis is being calculated...")
                    manager = Manager()
                    d_mp = manager.dict()
                    try:
                        utility.runInParallel(I_mp.run_clashscore(d_mp), I_mp.run_ramalyze(
                            d_mp), I_mp.run_rotalyze(d_mp), I_mp.run_molprobity(d_mp))
                        a, b = I_mp.process_molprobity(d_mp['molprobity'])
                        Template_Dict['bond'] = len(a)
                        Template_Dict['angle'] = len(b)
                        clashscore, rama, sidechain = I_mp.get_data_for_quality_at_glance(
                            d_mp['molprobity'])
                        Template_Dict['molp_b'] = utility.dict_to_JSlist(
                            I_mp.molprobity_detailed_table_bonds(a))
                        Template_Dict['molp_a'] = utility.dict_to_JSlist(
                            I_mp.molprobity_detailed_table_angles(b))
                        Template_Dict['rotascore'] = utility.dict_to_JSlist(
                            I_mp.rota_summary_table(I_mp.process_rota(d_mp['rota'])))
                        Template_Dict['rotalist'] = utility.dict_to_JSlist(
                            I_mp.rota_detailed_table(I_mp.process_rota(d_mp['rota'])))
                        Template_Dict['ramascore'] = utility.dict_to_JSlist(
                            I_mp.rama_summary_table(I_mp.process_rama(d_mp['rama'])))
                        Template_Dict['ramalist'] = utility.dict_to_JSlist(
                            I_mp.rama_detailed_table(I_mp.process_rama(d_mp['rama'])))
                        clashscores, Template_Dict['tot'] = I_mp.clash_summary_table(
                            d_mp['clash'])
                        Template_Dict['clashscore_list'] = utility.dict_to_JSlist(
                            clashscores)
                        Template_Dict['clashlist'] = I_mp.clash_detailed_table(
                            d_mp['clash'])
                        Template_Dict['assess_atomic_segments'] = 'Clashscore: ' + str(
                            clashscore) + ', Ramachandran outliers: ' + str(rama) + '% '+', \
                            Sidechain outliers: '+str(sidechain)+'%'
                        Template_Dict['assess_excluded_volume'] = [
                            'Not applicable']
                    except:
                        print("Molprobity cannot be calculated...")
                        clashscore = None
                        rama = None
                        sidechain = None
        else:
            Template_Dict['assess_atomic_segments'] = 'Not applicable'
            file = os.getcwd()+'Output/results/' + \
                str(Template_Dict['ID'])+'exv.txt'
            if os.path.exists(file):
                print("Excluded volume file already exists...")
                with open(file, 'r+') as inf:
                    line = [ln.replace('[', '').replace(']', '').replace(
                        ',', '').split() for ln in inf.readlines()]
                exv_data = {
                    'Models': line[0], 'Excluded Volume Satisfaction (%)':
                    line[1], 'Number of violations': line[2]}
            else:
                print("Excluded volume is being calculated...")
                I_ev = excludedvolume.GetExcludedVolume(self.mmcif_file)
                model_dict = I_ev.get_all_spheres()
                exv_data = I_ev.run_exc_vol_parallel(model_dict)

            Template_Dict['excluded_volume'] = utility.dict_to_JSlist(exv_data)
            Template_Dict['assess_excluded_volume'] = utility.exv_readable_format(
                exv_data)
            clashscore = None
            rama = None
            sidechain = None
        return Template_Dict, clashscore, rama, sidechain, exv_data

    def run_sas_validation(self, Template_Dict: dict) -> (dict, dict, dict):
        '''
        get sas validation information from SASCIF or JSON files
        '''
        if self.Input.check_for_sas(self.Input.get_dataset_comp()):
            Template_Dict['sas'] = ["True"]
            I_sas = sas.SasValidation(self.mmcif_file)
            Template_Dict['p_val'] = utility.dict_to_JSlist(I_sas.get_pvals())
            Template_Dict['sasdb_code'] = I_sas.get_SASBDB_code()
            try:
                Template_Dict['parameters_volume'] = utility.dict_to_JSlist(
                    I_sas.get_parameters_vol_many())
            except:
                Template_Dict['parameters_volume'] = utility.dict_to_JSlist(
                    I_sas.get_parameters_vol_many_dep())
            try:
                Template_Dict['parameters_mw'] = utility.dict_to_JSlist(
                    I_sas.get_parameters_mw_many())
            except:
                Template_Dict['parameters_mw'] = utility.dict_to_JSlist(
                    I_sas.get_parameters_mw_many_dep())
            Template_Dict['pddf_info'] = utility.dict_to_JSlist(
                I_sas.get_pddf_info())
            Template_Dict['number_of_fits'] = I_sas.get_total_fits()
            Template_Dict['chi_table'] = utility.dict_to_JSlist(
                I_sas.get_chi_table())
            Template_Dict['rg_table'] = utility.dict_to_JSlist(
                I_sas.get_rg_table_many())
            Template_Dict['sasdb_code_fits'] = I_sas.get_sasdb_code_fits()
            Template_Dict['Data_quality'] = utility.get_rg_data(
                I_sas.get_rg_for_plot())
            Template_Dict['validation_input'] = utility.get_rg_data_fits(
                I_sas.get_fits_for_plot())
            if len(Template_Dict['validation_input']) < 1:
                Template_Dict['validation_input'] = [
                    'Fit of model to data has not been deposited']
            sas_data = I_sas.get_rg_for_plot()
            sas_fit = I_sas.get_fits_for_plot()
        else:
            sas_data = {}
            sas_fit = {}
        return Template_Dict, sas_data, sas_fit

    def run_sas_validation_plots(self, Template_Dict: dict):
        '''
        get sas validation information from SASCIF or JSON files
        '''
        if self.Input.check_for_sas(self.Input.get_dataset_comp()):
            Template_Dict['sas'] = ["True"]
            I_sas = sas.SasValidation(self.mmcif_file)
            try:
                I_sas_plt = validation.sas_plots.SasValidationPlots(
                    self.mmcif_file)
                I_sas.modify_intensity()
                I_sas.get_pofr_errors()
                I_sas_plt.plot_multiple()
                I_sas_plt.plot_pf()
                I_sas_plt.plot_Guinier()
                if Template_Dict['number_of_fits'] > 0:
                    I_sas_plt.plot_fits()
            except:
                pass

    def run_cx_validation(self, Template_Dict: dict) -> (dict, dict):
        if self.Input.check_for_cx(self.Input.get_dataset_comp()):
            Template_Dict['cx'] = ["True"]
            I_cx = cx.CxValidation(self.mmcif_file)
            # xl_df = I_cx.get_xl_data()
            model_df = I_cx.get_df_for_models()
            cx_fit = I_cx.get_violation_plot(model_df)
            for key, value in cx_fit.items():
                Template_Dict['Cx Fit '+str(key)] = value
            try:
                Template_Dict['validation_input'].extend(
                    utility.get_cx_data_fits(cx_fit))
            except:
                Template_Dict['validation_input'] = utility.get_cx_data_fits(
                    cx_fit)
        else:
            cx_fit = dict()

        return cx_fit, Template_Dict

    def run_cx_validation_plots(self, Template_Dict: dict):
        if self.Input.check_for_cx(self.Input.get_dataset_comp()):
            Template_Dict['cx'] = ["True"]
            cx_plt = validation.cx_plots.CxValidationPlots(self.mmcif_file)
            cx_plt.make_gridplot_intra()
            cx_plt.make_gridplot_struc()
            cx_plt.plot_distributions()

    def run_quality_glance(self, clashscore: dict, rama: dict,
                           sidechain: dict, exv_data: dict,
                           sas_data: dict, sas_fit: dict,
                           cx_fit: dict):
        '''
        get quality at glance image; will be updated as validation report is updated
        '''
        I_plt = get_plots.Plots(self.mmcif_file)
        I_plt.plot_quality_at_glance(
            clashscore, rama, sidechain, exv_data, sas_data, sas_fit, cx_fit)

    def run_supplementary_table(self,
                                Template_Dict,
                                location='N/A',
                                physics='Information about physical principles was not provided',
                                method_details='N/A',
                                sampling_validation='N/A',
                                validation_input=['-'],
                                cross_validation='N/A',
                                Data_quality=['-'],
                                clustering='N/A',
                                resolution='N/A'):
        '''
        get supplementary table, will be updated as validation report is updated
        '''
        if (self.Input.get_ensembles() is not None) and (utility.all_same(self.Input.get_ensembles()['Clustering method'])):
            Template_Dict['clustering'] = self.Input.get_ensembles()[
                'Clustering method'][0]
        elif self.Input.get_ensembles() is not None:
            Template_Dict['clustering'] = ', '.join(
                self.Input.get_ensembles()['Clustering method'])
        else:
            Template_Dict['clustering'] = 'Not applicable'
        Template_Dict['location'] = location
        Template_Dict['complex_name'] = self.Input.get_struc_title().lower()
        Template_Dict['PDB_ID'] = self.Input.get_id()
        Template_Dict['Subunits'] = utility.get_subunits(
            self.Input.get_composition())
        Template_Dict['datasets'] = utility.get_datasets(self.Input.get_dataset_details(
        )) if self.Input.get_dataset_details() is not None else 'Not provided or used'
        Template_Dict['physics'] = physics
        Template_Dict['software'] = utility.get_software(
            self.Input.get_software_comp()) + location
        Template_Dict['struc'] = self.Input.get_atomic_coverage()
        Template_Dict['method'] = utility.get_method_name(
            self.Input.get_sampling())
        Template_Dict['method_type'] = utility.get_method_type(
            self.Input.get_sampling())
        Template_Dict['method_details'] = method_details
        Template_Dict['models'] = ', '.join(self.Input.get_ensembles(
        )['Number of models']) if self.Input.get_ensembles() is not None else 'Not applicable'
        Template_Dict['sampling_validation'] = sampling_validation
        Template_Dict['feature'] = self.Input.get_ensembles(
        )['Clustering feature'][0] if self.Input.get_ensembles() is not None else 'Not applicable'
        Template_Dict['cross_validation'] = cross_validation
        Template_Dict['model_precision'] = ', '.join([i+'&#8491' for i in self.Input.get_ensembles(
        )['Cluster precision']]) if self.Input.get_ensembles() is not None else \
            'Model precision can not be calculated with one structure'
        Template_Dict['restraint_info'] = utility.get_restraints_info(self.Input.get_restraints(
        )) if self.Input.get_restraints() is not None else 'Not provided or used'
        if 'Data_quality' not in list(Template_Dict.keys()):
            Template_Dict['Data_quality'] = Data_quality
        if 'validation_input' not in list(Template_Dict.keys()):
            Template_Dict['validation_input'] = validation_input
        Template_Dict['clustering'] = clustering
        Template_Dict['resolution'] = resolution
        return Template_Dict
