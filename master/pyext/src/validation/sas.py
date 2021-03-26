###################################
# Script :
# 1) Contains class to validate models
# built using SAS datasets
#
# ganesans - Salilab - UCSF
# ganesans@salilab.org
###################################
import os
import numpy as np
import pandas as pd
import requests
import json
from sklearn.linear_model import LinearRegression
from decimal import Decimal
from validation import GetInputInformation
from subprocess import run
from decouple import config
import operator


class SasValidation(GetInputInformation):
    def __init__(self, mmcif_file):
        super().__init__(mmcif_file)
        self.ID = str(GetInputInformation.get_id(self))
        self.nos = GetInputInformation.get_number_of_models(self)
        self.dataset = GetInputInformation.get_dataset_comp(self)
        self.imagepath = '../static/images/'
        self.saslink = 'https://www.sasbdb.org/media/sascif/sascif_files/'
        self.sasentry = 'https://www.sasbdb.org/rest-api/entry/summary/'

    def get_SASBDB_code(self) -> list:
        '''
        function to get all SASBDB codes used in the model,
        returns a list of SASBDB codes
        '''
        SAS_db_codes = []
        for indx, datatype in enumerate(self.dataset['Dataset type']):
            if 'SAS' in str(datatype):
                SAS_db_codes.append(self.dataset['Data access code'][indx])
        return SAS_db_codes

    def clean_SASBDB_code(self) -> list:
        '''
        function to clean SASBDB list of codes
        as some might have 'None' or can be repetitive
        '''
        codes = list(set(self.get_SASBDB_code()))
        cleaned_code = [i for i in codes if i != 'None']
        return cleaned_code

    def get_data_from_SASBDB(self) -> dict:
        '''
        get data from JSON
        '''
        data_dic = {}
        for code in self.get_SASBDB_code():
            if 'None' not in str(code):
                url_f = self.sasentry+code+'.json'
                response = requests.get(url_f, data={'key': 'value'})
                if response.status_code != 200:
                    print(
                        "Error....unable to fetch data from SASBDB, please check the entry ID")
                data_dic[code] = response.json()
                with open(code+'.json', 'w') as f:
                    formatted_data = json.dumps(
                        response.json(), indent=4, sort_keys=True)
                    f.write(formatted_data)
        return data_dic

    def get_sascif_file(self):
        '''
        get data from SASCIF files
        '''
        for code in self.get_SASBDB_code():
            if 'None' not in str(code):
                url_f = self.saslink+code+'.sascif'
                response = requests.get(url_f)
                if response.status_code != 200:
                    print(
                        "Error....unable to fetch data from SASBDB, please check the entry ID")
                with open(code+'.sascif', 'w') as f:
                    f.write(response.text)

    def get_all_sascif(self, sasbdb) -> list:
        '''
        get a list of all lines in a SASCIF file
        '''
        if 'None' not in str(sasbdb):
            file = open(sasbdb+'.sascif', 'r')
            all_lines = [data.strip().split()
                         for indx, data in enumerate(file.readlines())]
        return all_lines

    def get_intensities(self) -> dict:
        '''
        get intensity data from SASCIF file
        if SASCIF file is not present, this information will not be present/used in the report
        JSON file typically has raw data
        '''
        self.get_sascif_file()
        Int_dict = {}
        for code in self.clean_SASBDB_code():
            # Int_dict={}
            all_lines = self.get_all_sascif(code)
            data = {}
            for indx, sascifline in enumerate(all_lines):
                if len(sascifline) < 2 and len(sascifline) > 0 and 'scan_intensity' in sascifline[0]:
                    data[(sascifline[0].split('.')[1])] = []
                if len(sascifline) > 2 and len(all_lines[indx-1]) > 0 and 'scan_intensity' in all_lines[indx-1][0]:
                    for indx_sub, sascifline_sub in enumerate(all_lines[indx:]):
                        if len(sascifline_sub) > 2 and '#' not in sascifline_sub and 'sas' not in sascifline_sub[0]:
                            for num, key in enumerate(list(data.keys())):
                                data[key].append(sascifline_sub[num])
                        else:
                            break
            I_df = pd.DataFrame(list(data.values()), index=list(data.keys())).T
            I_df_re = I_df[['momentum_transfer',
                            'intensity', 'intensity_su_counting']]
            I_df_re.rename(columns={
                           'momentum_transfer': 'Q', 'intensity': 'I', 'intensity_su_counting': 'E'}, inplace=True)
            Int_dict[code] = I_df_re
        return Int_dict

    def modify_intensity(self) -> dict:
        '''
        modify intensity data to calcualte errors and log values
        '''
        Int_dict = self.get_intensities()
        Int_dict_modify = {}
        rg_and_io = self.get_rg_and_io()
        for key, val in Int_dict.items():
            Rg = rg_and_io[key][0]
            IO = rg_and_io[key][1]
            dim_num = Rg*Rg/IO
            I_df = val.astype({'Q': float, 'I': float, 'E': float})
            I_df = I_df[I_df['I']-I_df['E'] > 0]
            I_df['Q'] = I_df['Q']*10
            I_df['err_x'] = I_df.apply(
                lambda row: (row['Q'], row['Q']), axis=1)
            I_df['err_y'] = I_df.apply(lambda row: (
                np.log(row['I']-row['E']), np.log(row['I']+row['E'])), axis=1)
            I_df['logI'] = np.log(I_df['I'])
            I_df['logQ'] = np.log(I_df['Q'])
            I_df['logX'] = I_df.apply(lambda row: (
                row['logQ'], row['logQ']), axis=1)
            I_df['Ky'] = I_df['Q']*I_df['Q']*I_df['I']*dim_num
            I_df['Kx'] = I_df['Q']*Rg
            I_df['Px'] = I_df['Q']**4
            I_df['Px'].round(3)
            I_df['Py'] = I_df['Px']*I_df['I']
            Int_dict_modify[key] = I_df
        return Int_dict_modify

    def modify_intensity_dep(self) -> dict:
        '''
        depreciated function to get intensities from JSON/raw data
        '''
        Int_dict = self.get_intensities()
        Int_dict_modify = {}
        for key, val in Int_dict.items():
            I_df = val.astype({'Q': float, 'I': float, 'E': float})
            I_df.head()
            I_df = I_df[I_df['I']-I_df['E'] > 0]
            I_df['Q'] = I_df['Q']*10
            I_df['err_x'] = I_df.apply(
                lambda row: (row['Q'], row['Q']), axis=1)
            I_df['err_y'] = I_df.apply(lambda row: (
                np.log(row['I']-row['E']), np.log(row['I']+row['E'])), axis=1)
            I_df['logI'] = np.log(I_df['I'])
            I_df['logQ'] = np.log(I_df['Q'])
            I_df['logX'] = I_df.apply(lambda row: (
                row['logQ'], row['logQ']), axis=1)
            I_df['Ky'] = I_df['Q']*I_df['Q']*I_df['I']
            I_df['Px'] = I_df['Q']**4
            I_df['Px'].round(3)
            I_df['Py'] = I_df['Px']*I_df['I']
            Int_dict_modify[key] = I_df
        return Int_dict_modify

    def get_rg_for_plot(self) -> dict:
        '''
        get Rg values from SASCIF file, if unavailabel, get it from JSON
        '''
        self.get_sascif_file()
        Rg_dict = {}
        for code in self.clean_SASBDB_code():
            rg = {}
            all_lines = self.get_all_sascif(code)
            for indx, sascifline in enumerate(all_lines):
                if len(sascifline) < 3 and len(sascifline) > 0 and 'sas_result.Rg_from_PR' in sascifline[0] and 'sas_result.Rg_from_PR_' not in sascifline[0]:
                    rg[sascifline[0].split('.')[1]] = float(sascifline[1])
                if len(sascifline) < 3 and len(sascifline) > 0 and 'sas_result.Rg_from_Guinier' in sascifline[0] and 'sas_result.Rg_from_Guinier_' not in sascifline[0]:
                    rg[sascifline[0].split('.')[1]] = float(sascifline[1])
            Rg_dict[code] = list(rg.values())

        if len(list(rg.values())) < 1:
            data_dic = self.get_data_from_SASBDB()
            for key, val in data_dic.items():
                Rg_dict[key] = []
                Rg_dict[key].append(round(float(val['guinier_rg']), 2))
                Rg_dict[key].append(round(float(val['pddf_rg']), 2))
        return Rg_dict

    def get_rg_and_io(self) -> dict:
        '''
        get rg information from SASCIF file
        '''
        self.get_sascif_file()
        rg_and_io = {}
        for code in self.clean_SASBDB_code():
            all_lines = self.get_all_sascif(code)
            for indx, sascifline in enumerate(all_lines):
                if len(sascifline) < 3 and len(sascifline) > 0 and 'sas_result.Rg_from_PR' in sascifline[0] and 'sas_result.Rg_from_PR_' not in sascifline[0]:
                    rg = float(sascifline[1])
                if len(sascifline) < 3 and len(sascifline) > 0 and '_sas_result.I0_from_PR' in sascifline[0] and '_sas_result.I0_from_PR_' not in sascifline[0]:
                    io = float(sascifline[1])
            rg_and_io[code] = (rg, io)
        return rg_and_io

    def get_rg_table_many(self) -> dict:
        '''
        get rg information from multiple SASCIF files
        '''
        data_dic = self.get_data_from_SASBDB()
        rg_table = {'SASDB ID': [], 'Rg': [],
                    'Rg error': [], 'MW': [], 'MW error': []}
        for key, val in data_dic.items():
            rg_table['Rg'].append(
                str(round(float(val['guinier_rg']), 2)) + ' nm')
            try:
                rg_table['Rg error'].append(str(round(float(val['guinier_rg_error']), 2)) + ' nm')
            except (TypeError, KeyError, ValueError):
                rg_table['Rg error'].append('N/A')
            try:
                rg_table['MW'].append(str(round(float(val['guinier_i0_mw']), 2)) + ' nm')
            except (TypeError, KeyError, ValueError):
                rg_table['MW'].append('N/A')
            try:
                rg_table['MW error'].append(str(round(float(val['guinier_i0_mw_error']), 2))  + ' nm')
            except (TypeError, KeyError, ValueError):
                rg_table['MW error'].append('N/A')
            rg_table['SASDB ID'].append(key)
        return rg_table

    def get_fits_for_plot(self) -> dict:
        '''
        get chi-squared values from SASCIF files
        '''
        self.get_sascif_file()
        fit_dict = {}
        for code in self.clean_SASBDB_code():
            fits = []
            all_lines = self.get_all_sascif(code)
            for indx, sascifline in enumerate(all_lines):
                if (len(sascifline) < 3) and (len(sascifline) > 0) and ('sas_model_fitting_details.chi_square' in sascifline[0]) and (float(sascifline[1]) > 0.00000):
                    fits.append(round(float(sascifline[1]), 2))
            if len(fits) > 0:
                fit_dict[code] = fits
        return fit_dict

    def get_pofr(self) -> dict:
        '''
        get pair-dist distribution from SASCIF files
        '''
        pofr_dict = {}
        for code in self.clean_SASBDB_code():
            all_lines = self.get_all_sascif(code)
            data = {}
            for indx, sascifline in enumerate(all_lines):
                if len(sascifline) < 2 and len(sascifline) > 0 and 'sas_p_of_R.' in sascifline[0]:
                    data[(sascifline[0].split('.')[1])] = []
                if len(sascifline) > 2 and len(all_lines[indx-1]) > 0 and 'sas_p_of_R.' in all_lines[indx-1][0]:
                    for subindx, subval in enumerate(all_lines[indx:]):
                        if len(subval) > 2 and '#' not in subval and 'sas' not in subval[0]:
                            for num, key in enumerate(list(data.keys())):
                                data[key].append(subval[num])
                        else:
                            break
            pdf = pd.DataFrame(list(data.values()), index=list(data.keys())).T
            pdf_re = pdf[['r', 'P', 'P_error']]
            pdf_re.rename(columns={'r': 'R', 'P': 'P',
                                   'P_error': 'E'}, inplace=True)
            pofr_dict[code] = pdf_re
        return pofr_dict

    def get_pvals(self) -> dict:
        '''
        get p-values from ATSAS
        '''
        data_dic = self.get_data_from_SASBDB()
        num_of_fits = self.get_number_of_fits()
        pval_table = {'SASDB ID': [], 'Model': [], 'p-value': []}
        for key, val in data_dic.items():
            num = num_of_fits[key]
            if num > 0:
                for fitnum in range(0, num):
                    pval_table['SASDB ID'].append(key)
                    pval_table['Model'].append(fitnum+1)
                    target_url = val['fits'][fitnum]['fit_data']
                    fit = requests.get(target_url)
                    if fit.status_code != 200:
                        print(
                            "Error....unable to fetch data from SASBDB, please check the entry ID")
                    fname = key+str(fitnum)+'fit.csv'
                    with open(fname, 'w') as f:
                        f.write(fit.text)
                    f_df = pd.read_csv(fname, skiprows=3, delim_whitespace=True, names=[
                                       'Q', 'Ie', 'Ib', 'E'])
                    if abs(f_df.iloc[22, 2]-f_df.iloc[22, 1]) > abs(f_df.iloc[22, 3]-f_df.iloc[22, 1]):
                        f_df.rename(
                            columns={'Q': 'Q', 'Ie': 'Ie', 'Ib': 'E', 'E': 'Ib'}, inplace=True)
                    fit_1 = f_df[['Q', 'Ie']]
                    fit_1.to_csv('fit1.csv', header=False, index=False)
                    fit_2 = f_df[['Q', 'Ib']]
                    fit_2.to_csv('fit2.csv', header=False, index=False)
                    f1 = open('pval.txt', 'w+')
                    with f1 as outfile:
                        run([config('ATSAS'), 'fit1.csv',
                             'fit2.csv'], stdout=outfile)
                    f2 = open('pval.txt', 'r')
                    all_lines = [j.strip().split()
                                 for i, j in enumerate(f2.readlines())]
                    p_val = [all_lines[i+1][4]
                             for i, j in enumerate(all_lines) if 'adj' in j][0]
                    pval_table['p-value'].append('%.2E' % Decimal(p_val))
            else:
                pval_table['SASDB ID'].append(key)
                pval_table['Model'].append('N/A')
                pval_table['p-value'].append('N/A')
        return pval_table

    def get_pofr_ext(self) -> dict:
        '''
        get pair-distance details from SASCIF files
        '''
        pofr_dict = {}
        for code in self.clean_SASBDB_code():
            all_lines = self.get_all_sascif(code)
            data = {}
            for indx, sascifline in enumerate(all_lines):
                if len(sascifline) < 2 and len(sascifline) > 0 and '_sas_p_of_R_extrapolated_intensity.' in sascifline[0]:
                    data[(sascifline[0].split('.')[1])] = []
                if len(sascifline) > 2 and len(all_lines[indx-1]) > 0 and '_sas_p_of_R_extrapolated_intensity.' in all_lines[indx-1][0]:
                    for subindx, subval in enumerate(all_lines[indx:]):
                        if len(subval) > 2 and '#' not in subval and 'sas' not in subval[0]:
                            for num, key in enumerate(list(data.keys())):
                                data[key].append(subval[num])
                        else:
                            break
            pdf = pd.DataFrame(list(data.values()), index=list(data.keys())).T
            pdf_re = pdf[['momentum_transfer', 'intensity_reg']]
            pdf_re.rename(columns={'momentum_transfer': 'Q',
                                   'intensity_reg': 'I'}, inplace=True)
            pdf_re = pdf_re.astype({'Q': float, 'I': float})
            pdf_re['Q'] = pdf_re['Q']*10
            pdf_re['logI'] = np.log(pdf_re['I'])
            pofr_dict[code] = pdf_re
        return pofr_dict

    def get_pofr_errors(self) -> dict:
        '''
        get pair-distance details and errors from JSON files
        '''
        pofr_dict = self.get_pofr_ext()
        Int_dict = self.modify_intensity()
        compiled_dict = {}
        for code in self.clean_SASBDB_code():
            I_df = Int_dict[code]
            I_df_dict = dict(zip(I_df.Q, I_df.I))
            I_df_err_dict = dict(zip(I_df.Q, I_df.E))
            p_df = pofr_dict[code]
            p_df_dict = dict(zip(p_df.Q, p_df.I))
            errors = []
            for Q, I in p_df_dict.items():
                data_Q = self.findMinDiff(list(I_df_dict.keys()), Q)
                if data_Q != 9999:
                    data_I = I_df_dict[data_Q]
                    delta_I = (I-data_I)
                    if I_df_err_dict[data_Q] != 0:
                        wt_delta_I = delta_I/I_df_err_dict[data_Q]
                    else:
                        wt_delta_I = 0
                    errors.append([Q, delta_I, wt_delta_I])
            errors_df = pd.DataFrame(errors, columns=['Q', 'R', 'WR'])
            compiled_dict[code] = errors_df
        return compiled_dict

    def findMinDiff(self, listn: list, num: int) -> int:
        '''
        quick min diff operation for calculating errors
        '''
        list_sub = [(i, abs(j-num)) for i, j in enumerate(listn)]
        list_sort = sorted(list_sub, key=lambda x: operator.itemgetter(1))
        if list_sort[0][1] < 0.00001:
            return listn[list_sort[0][0]]
        else:
            return 9999

    def get_Guinier_data(self) -> (dict, dict):
        '''
        get Guinier plot data from JSON files
        '''
        Int_dict = self.get_intensities()
        data_dic = self.get_data_from_SASBDB()
        Guinier_dict = {}
        Guinier_score = {}
        for key, val in Int_dict.items():
            G_df = val.astype({'Q': float, 'I': float, 'E': float})
            G_df['logI'] = np.log(G_df['I'])
            rg = float(data_dic[key]['pddf_rg'])
            # dmax = float(data_dic[key]['pddf_dmax'])
            # index_low = int(data_dic[key]['guinier_point_first'])
            # index_high = int(data_dic[key]['guinier_point_last'])
            # q_min = math.pi/(dmax*10)
            q_max = 1.3/(rg*10)
            G_df_range = G_df[G_df['Q'] < q_max].copy()
            G_df_range['Q'] = G_df['Q']*10
            G_df_range['Q2'] = G_df_range['Q']**2
            X = G_df_range[['Q2']].values
            y = G_df_range['logI'].values
            regression = LinearRegression(fit_intercept=True)
            regression.fit(X, y)
            G_df_range['y_pred'] = regression.predict(X)
            G_df_range['res'] = y-regression.predict(X)
            G_df_range['Q2A'] = G_df_range['Q2']*100
            score = '%.2f' % regression.score(X, y)
            Guinier_score[key] = score
            Guinier_dict[key] = G_df_range
        return Guinier_score, Guinier_dict

    def get_parameters_vol_many_dep(self) -> dict:
        '''
        get volume parameters from JSON files
        '''
        data_dic = self.get_data_from_SASBDB()
        parameter_table = {'SASDB ID': [], 'Estimated volume': [
        ], 'Estimated volume method': [], 'Porod volume': []}
        for key, val in data_dic.items():
            try:
                parameter_table['Estimated volume'].append(
                        val['estimated_volume'])
                parameter_table['Estimated volume method'].append(
                        val['estimated_volume_method'])
            except (TypeError, KeyError, ValueError):
                parameter_table['Estimated volume'].append('N/A')
                parameter_table['Estimated volume method'].append('N/A')
            try:
                parameter_table['Porod volume'].append(
                    val['porod_volume']+' nm\u00b3')
            except (TypeError, KeyError, ValueError):
                parameter_table['Porod volume'].append('N/A')
            parameter_table['SASDB ID'].append(key)
        return parameter_table

    def get_parameters_vol_many(self) -> dict:
        '''
        get volume details from SASCIF files
        '''
        self.get_sascif_file()
        parameter_table = {'SASDB ID': [], 'Estimated Volume': [], 'Porod Volume': [], 'Specific Volume': [],
                           'Sample Contrast': [], 'Sample Concentration': []}
        for code in self.clean_SASBDB_code():
            parameter_table['SASDB ID'].append(code)
            all_lines = self.get_all_sascif(code)
            for indx, sascifline in enumerate(all_lines):
                if 0 < len(sascifline) < 3 and '_sas_sample.specimen_concentration' in sascifline[0]:
                    if len(sascifline[1]) > 1:
                        parameter_table['Sample Concentration'].append(
                            str(round(float(sascifline[1]), 2))+' mg/ml')
                    else:
                        parameter_table['Sample Concentration'].append('N/A')
                if 0 < len(sascifline) < 3 and '_sas_sample.contrast' in sascifline[0]:
                    if len(sascifline[1]) > 1:
                        parameter_table['Sample Contrast'].append(
                            str(round(float(sascifline[1]), 2)))
                    else:
                        parameter_table['Sample Contrast'].append('N/A')
                if 0 < len(sascifline) < 3 and '_sas_sample.specific_vol' in sascifline[0]:
                    if len(sascifline[1]) > 1:
                        parameter_table['Specific Volume'].append(
                            str(round(float(sascifline[1]), 2))+' nm\u00b3')
                    else:
                        parameter_table['Specific Volume'].append('N/A')
                if 0 < len(sascifline) < 3 and '_sas_result.Porod_volume' in sascifline[0] and '_sas_result.Porod_volume_error' not in sascifline[0]:
                    if len(sascifline[1]) > 1:
                        parameter_table['Porod Volume'].append(
                            str(round(float(sascifline[1]), 2))+' nm\u00b3')
                    else:
                        parameter_table['Porod Volume'].append('N/A')
                if 0 < len(sascifline) < 3 and '_sas_result.estimated_volume' in sascifline[0] and '_sas_result.estimated_volume_error' not in sascifline[0]:
                    if len(sascifline[1]) > 1:
                        parameter_table['Estimated Volume'].append(
                            str(round(float(sascifline[1]), 2))+' nm\u00b3')
                    else:
                        parameter_table['Estimated Volume'].append('N/A')
        return parameter_table

    def get_parameters_mw_many(self) -> dict:
        '''
        get MW details from SASCIF files
        '''
        self.get_sascif_file()
        parameter_table = {'SASDB ID': [], 'Chemical composition MW': [
        ], 'Standard MW': [], 'Porod Volume/MW': []}
        for code in self.clean_SASBDB_code():
            parameter_table['SASDB ID'].append(code)
            all_lines = self.get_all_sascif(code)
            for indx, sascifline in enumerate(all_lines):
                if len(sascifline) < 3 and len(sascifline) > 0 and '_sas_result.experimental_MW' in sascifline[0] and '_sas_result.experimental_MW_error' not in sascifline[0]:
                    if len(sascifline[1]) > 1:
                        parameter_table['Chemical composition MW'].append(
                            str(round(float(sascifline[1]), 2))+' kDa')
                    else:
                        parameter_table['Chemical composition MW'].append(
                            'N/A')
                if len(sascifline) < 3 and len(sascifline) > 0 and '_sas_result.MW_standard' in sascifline[0] and '_sas_result.MW_standard_error' not in sascifline[0]:
                    if len(sascifline[1]) > 1:
                        parameter_table['Standard MW'].append(
                            str(round(float(sascifline[1]), 2))+' kDa')
                    else:
                        parameter_table['Standard MW'].append('N/A')
                if len(sascifline) < 3 and len(sascifline) > 0 and '_sas_result.MW_Porod' in sascifline[0] and '_sas_result.MW_Porod_error' not in sascifline[0]:
                    if len(sascifline[1]) > 1:
                        Porod_MW = round(float(sascifline[1]), 2)
                    else:
                        Porod_MW = 0
                if len(sascifline) < 3 and len(sascifline) > 0 and '_sas_result.Porod_volume' in sascifline[0] and '_sas_result.Porod_volume_error' not in sascifline[0]:
                    if len(sascifline[1]) > 1 and Porod_MW > 0:
                        Porod_V = round(float(sascifline[1]), 2)/Porod_MW
                        parameter_table['Porod Volume/MW'].append(
                            str(round(Porod_V, 2))+' nm \u00b3/kDa')
                    else:
                        parameter_table['Porod Volume/MW'].append('N/A')

        return parameter_table

    def get_parameters_mw_many_dep(self) -> dict:
        '''
        depreciated function on getting MW from JSON
        '''
        data_dic = self.get_data_from_SASBDB()
        # MW based on chemical composition
        parameter_table = {'SASDB ID': [], 'Sequence MW': [],
                           'Experimental MW': [], 'Porod MW': []}
        for key, val in data_dic.items():
            try:
                parameter_table['Experimental MW'].append(
                    list(data_dic.values())[0]['experimental_mw']+' kDa')
            except (TypeError, KeyError, ValueError):
                parameter_table['Experimental MW'].append('N/A')
            try:
                parameter_table['Porod MW'].append(
                    list(data_dic.values())[0]['porod_mw']+' kDa')
            except (TypeError, KeyError, ValueError):
                parameter_table['Porod MW'].append('N/A')
            try:
                parameter_table['Sequence MW'].append(list(data_dic.values())[
                                                      0]['experiment']['sample']['molecule'][0]['total_mw']+' kDa')
            except (TypeError, KeyError, ValueError):
                parameter_table['Sequence MW'].append('N/A')
            parameter_table['SASDB ID'].append(key)
        return parameter_table

    def get_pddf(self) -> dict:
        '''
        get p(r) data from JSON
        '''
        # data_dic = self.get_data_from_SASBDB()
        pofr_dic = self.get_pofr()
        pddf_dic = {}
        for key, val in pofr_dic.items():
            pd_df = val.astype({'P': float, 'R': float, 'E': float})
            pd_df['R'] = pd_df['R']/10
            pd_df['err_x'] = pd_df.apply(
                lambda row: (row['R'], row['R']), axis=1)
            pd_df['err_y'] = pd_df.apply(lambda row: (
                row['P']-row['E'], row['P']+row['E']), axis=1)
            pddf_dic[key] = pd_df
        return pddf_dic

    def get_pddf_info(self) -> dict:
        '''
        get p(r) related info from JSON
        '''
        data_dic = self.get_data_from_SASBDB()
        pddf_info = {'SASDB ID': [], 'Software used': [],
                     'Dmax': [], 'Dmax error': [], 'Rg': [], 'Rg error': []}
        for key, val in data_dic.items():
            pddf_info['Software used'].append(str(val['pddf_software']))
            try:
                pddf_info['Dmax'].append(str(val['pddf_dmax'])+' nm')
            except (TypeError, KeyError, ValueError):
                pddf_info['Dmax'].append('N/A')
            try:
                pddf_info['Rg'].append(str(val['pddf_rg'])+' nm')
            except (TypeError, KeyError, ValueError):
                pddf_info['Rg'].append('N/A')
            try:
                pddf_info['Dmax error'].append(
                        str(val['pddf_dmax_error'])+' nm')
            except (TypeError, KeyError, ValueError):
                pddf_info['Dmax error'].append('N/A')
            try:
                pddf_info['Rg error'].append(str(val['pddf_rg_error'])+' nm')
            except (TypeError, KeyError, ValueError):
                pddf_info['Rg error'].append('N/A')
            pddf_info['SASDB ID'].append(key)
        return pddf_info

    def get_number_of_fits(self) -> dict:
        '''
        get number of fits from JSON, deprecated
        '''
        data_dic = self.get_data_from_SASBDB()
        num_of_fits = {}
        for key, val in data_dic.items():
            num_of_fits[key] = len(val['fits'])
        return num_of_fits

    def get_chi_table(self) -> dict:
        '''
        get chi value from JSON, deprecated
        '''
        data_dic = self.get_data_from_SASBDB()
        chi_table = {'SASDB ID': [], 'Model': [], '\u03C7\u00b2': []}
        for key, val in data_dic.items():
            numoffits = self.get_number_of_fits()[key]
            if numoffits > 0:
                for fitnum in range(0, numoffits):
                    count = fitnum+1
                    chi_table['SASDB ID'].append(key)
                    chi_table['Model'].append(str(count))
                    chi_value = val['fits'][fitnum]['chi_square_value']
                    chi_value_round = round(chi_value, 2)
                    chi_table['\u03C7'+'\u00b2'].append(chi_value_round)
            else:
                chi_table['SASDB ID'].append(key)
                chi_table['Model'].append('N/A')
                chi_table['\u03C7'+'\u00b2'].append('N/A')
        return chi_table

    def get_sasdb_code_fits(self) -> list:
        '''
        get number of fits per SASBDB ID
        '''
        fit_dict = self.get_number_of_fits()
        return list(fit_dict.values())

    def get_fit_data(self) -> dict:
        '''
        get fit information to make plots, from JSON
        '''
        data_dic = self.get_data_from_SASBDB()
        num_of_fits = self.get_number_of_fits()
        data_fit = {}
        for key, val in data_dic.items():
            num = num_of_fits[key]
            fits = {}
            if num > 0:
                for fitnum in range(0, num):
                    target_url = val['fits'][fitnum]['fit_data']
                    fit = requests.get(target_url)
                    if fit.status_code != 200:
                        print(
                            "Error....unable to fetch data from SASBDB, please check the entry ID")

                    fname = key+str(fitnum)+'fit.csv'
                    with open(fname, 'w') as f:
                        f.write(fit.text)

                    f_df = pd.read_csv(fname, skiprows=3, delim_whitespace=True, names=[
                                       'Q', 'Ie', 'Ib', 'E'])
                    if abs(f_df.iloc[22, 2]-f_df.iloc[22, 1]) > abs(f_df.iloc[22, 3]-f_df.iloc[22, 1]):
                        f_df.rename(
                            columns={'Q': 'Q', 'Ie': 'Ie', 'Ib': 'E', 'E': 'Ib'}, inplace=True)
                    f_df['logIe'] = np.log(f_df['Ie'])
                    f_df['logIb'] = np.log(f_df['Ib'])
                    f_df['r'] = f_df['Ie']-f_df['Ib']

                    if f_df['E'].isnull().values.any():
                        f_df['rsigma'] = 0
                    else:
                        f_df['rsigma'] = f_df['r']/f_df['E']
                    # f_df['rsigma']=f_df['r']/f_df['E']
                    f_df['logr'] = f_df['logIe']-f_df['logIb']
                    f_df['r2a'] = (f_df['Ib']-f_df['Ie'].mean())**2
                    f_df['r2b'] = (f_df['Ie']-f_df['Ie'].mean())**2
                    fits[fitnum] = (self.get_fit_r2(f_df), f_df)
            else:
                fits[0] = (0, pd.DataFrame())
                # data_fit=None
            data_fit[key] = fits
        return data_fit

    def get_fit_r2(self, df: pd.DataFrame) -> int:
        rsquared = df['r2a'].sum()/df['r2b'].sum()
        return round(rsquared, 2)

    def get_total_fits(self) -> int:
        '''
        get number of fits
        '''
        data_dic = self.get_data_from_SASBDB()
        num_of_fits = 0
        for key, val in data_dic.items():
            num_of_fits += len(val['fits'])
        return num_of_fits

    def get_fit_image(self):
        '''
        get fit image from fit, deprecated
        '''
        data_dic = self.get_data_from_SASBDB()
        num_of_fits = self.get_number_of_fits()
        # data_fit = {}
        for key, val in data_dic.items():
            num = num_of_fits[key]
            if num > 0:
                for fitnum in range(0, num):
                    target_url = val['fits'][fitnum]['models'][0]['model_plot']
                    fitdata = requests.get(target_url)
                    if fitdata.status_code != 200:
                        print(
                            "Error....unable to fetch data from SASBDB, please check the entry ID")
                    # dirname = os.path.dirname(os.path.abspath(__file__))
                    filename = os.path.abspath(os.path.join(
                        os.getcwd(), self.imagepath, self.ID+key+str(fitnum)+'fit.png'))
                    with open(filename, 'wb') as f:
                        f.write(fitdata.content)
