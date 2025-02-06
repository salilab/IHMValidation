###################################
# Script :
# 1) Contains class to validate models
# built using SAS datasets
#
# ganesans - Salilab - UCSF
# ganesans@salilab.org
###################################
import os
import re
import subprocess
import numpy as np
import pandas as pd
import requests
import json
from sklearn.linear_model import LinearRegression
from decimal import Decimal
from mmcif_io import GetInputInformation
from subprocess import run
import operator
import logging
from io import StringIO

import sys
import pkgutil
from pathlib import Path

saspath = pkgutil.get_loader('sasciftools').path
saspath = str(Path(saspath).parent)
sys.path.insert(0, saspath)
from mmCif import mmcifIO
import utility
import ihm

class SasValidation(GetInputInformation):
    def __init__(self, mmcif_file, db='.'):
        super().__init__(mmcif_file)
        self.version = self.get_atsas_version()
        self.dataset = GetInputInformation.get_dataset_comp(self)
        self.imagepath = '../static/images/'
        self.saslink = 'https://sasbdb.org/media/sascif/sascif_files/'
        # self.sasentry = 'https://sasbdb.org/rest-api/entry/summary/'
        self.db = db
        self.sasbdb_ids = self.get_sasbdb_ids()
        self.sascif_dicts = self.get_sascif_dicts()
        self.intensities = self.get_intensities()
        self.intensities = self.modify_intensity()
        # self.data_dic = self.get_data_from_SASBDB()

    def get_atsas_version(self, tool: str = 'datcmp') -> str:
        """ Get ATSAS version """
        line = subprocess.check_output(
            [tool, '--version'],
            stderr=subprocess.STDOUT, text=True).strip()

        q = re.search(f'^{tool}, ATSAS (?P<version>.*)\nCopyright', line)

        version = q.group("version")

        return version

    def get_sas_ids(self) -> list:
        '''
        function to get all SASBDB codes used in the model,
        returns a list of SASBDB codes
        '''
        sas_ids = []
        for indx, datatype in enumerate(self.dataset['Database name']):
            if str(datatype) == 'SASBDB':
                sas_id = self.dataset['Data access code'][indx]
                sas_ids.append(sas_id)
        return sas_ids


    def get_sasbdb_ids(self) -> list:
        '''
        function to get all SASBDB codes used in the model,
        returns a list of SASBDB codes
        '''
        sasbdb_ids = []
        for dataset in self.system.orphan_datasets:
            if isinstance(dataset.location, ihm.location.SASBDBLocation):
                try:
                    code = dataset.location.access_code
                except AttributeError as e:
                    logging.error('Missing SASBDB accession code')
                    logging.error(e)
                else:
                    sasbdb_ids.append(code)

        return sasbdb_ids

    def get_sascif_dicts(self):
        sascif_dicts = {}
        sasCIFIn = mmcifIO.CifFileReader()

        for code in self.sasbdb_ids:
            sascif_fn = self.get_sascif_file(code, self.db)
            if sascif_fn is not None:
                sascif_dicts[code] = sasCIFIn.read(sascif_fn)

        return sascif_dicts

    def check_sascif_dicts(self):
        return(['True' for x in self.sascif_dicts.keys()])

    def get_sascif_file(self, code, output_dir='.'):
        '''
        get data from SASCIF files
        '''

        url = f'{self.saslink}{code}.sascif'
        fn = Path(output_dir, f'{code}.sascif')

        # Check if we already requested the data
        if os.path.isfile(fn):
            logging.info(f'Found {fn} in cache!')
        elif not os.path.isfile(fn):
            response = requests.get(url)
            response.encoding = 'ascii'
            if response.status_code != 200:
                logging.error(
                    f"Unable to fetch data for {code} from SASBDB, "
                    "please check the entry ID")
                fn = None
            else:
                with open(fn, 'w') as f:
                    f.write(response.text)

        return fn

    def get_intensities(self) -> dict:
        '''
        get intensity data from SASCIF file
        '''
        ints = {}
        for code in self.sascif_dicts.keys():
            sascif = self.sascif_dicts[code]

            main = f'{code}_MAIN'

            if main not in sascif:
                raise(KeyError(f'Missing MAIN dataset for {code}'))

            data = sascif[main]['_sas_scan_intensity']

            try:
                I_df = pd.DataFrame({
                    'id': np.array(data['id'], dtype=int),
                    'Q': np.array(data['momentum_transfer'], dtype=float),
                    'I': np.array(data['intensity'], dtype=float),
                    'E': np.array(data['intensity_su_counting'], dtype=float),
                })
            except:
                raise(ValueError("Can't parse intensities from {code}"))

            ints[code] = I_df
        return ints

    def modify_intensity(self) -> dict:
        '''
        modify intensity data to calcualte errors and log values
        '''
        Int_dict_modify = {}
        rg_and_io = self.get_rg_and_io()

        for code, I_df in self.intensities.items():
            data = self.sascif_dicts[code][f'{code}_MAIN']
            unit = data['_sas_scan']['unit']
            unitm = self.get_scan_unit_mult(unit)

            Rg = rg_and_io[code][0]
            IO = rg_and_io[code][1]
            dim_num = Rg * Rg / IO
            I_df = I_df[I_df['I']-I_df['E'] > 0]
            I_df['Q'] = I_df['Q'] * unitm
            I_df['err_x'] = I_df.apply(
                lambda row: (row['Q'], row['Q']), axis=1)
            I_df['err_y'] = I_df.apply(
                lambda row: (
                    np.log(row['I'] - row['E']),
                    np.log(row['I'] + row['E'])),
                axis=1)
            I_df['logI'] = np.log(I_df['I'])
            I_df['logQ'] = np.log(I_df['Q'])
            I_df['logX'] = I_df.apply(lambda row: (
                row['logQ'], row['logQ']), axis=1)
            I_df['Ky'] = I_df['Q'] * I_df['Q'] * I_df['I'] * dim_num
            I_df['Kx'] = I_df['Q'] * Rg
            I_df['Px'] = I_df['Q']**4
            I_df['Px'].round(3)
            I_df['Py'] = I_df['Px']*I_df['I']
            Int_dict_modify[code] = I_df
        return Int_dict_modify

    def get_rg_for_plot(self) -> dict:
        '''
        get Rg values from SASCIF file, if unavailabel, get it from JSON
        '''
        Rg_dict = {}
        for code in self.sascif_dicts.keys():
            sascif = self.sascif_dicts[code]
            main = f'{code}_MAIN'
            data = sascif[main]['_sas_result']

            rgs = []
            rg = round(float(data['Rg_from_Guinier']), 2)
            rgs.append(rg)
            rg = round(float(data['Rg_from_PR']), 2)
            rgs.append(rg)
            Rg_dict[code] = rgs

        return Rg_dict

    def get_rg_and_io(self) -> dict:
        '''
        get rg information from SASCIF file
        '''
        rg_and_io = {}
        for code in self.sascif_dicts.keys():
            sascif = self.sascif_dicts[code]
            main = f'{code}_MAIN'
            data = sascif[main]['_sas_result']
            rg = float(data['Rg_from_PR'])
            io = float(data['I0_from_PR'])
            rg_and_io[code] = (rg, io)
        return rg_and_io

    def get_rg_table_many(self) -> dict:
        '''
        get rg information from multiple SASCIF files
        '''
        rg_table = {'SASDB ID': [], 'R<sub>g</sub>': [],
                    'R<sub>g</sub> error': [], 'MW': [], 'MW error': []}
        for code in self.sascif_dicts.keys():
            sascif = self.sascif_dicts[code]
            main = f'{code}_MAIN'
            data = sascif[main]['_sas_result']

            rg_table['SASDB ID'].append(code)

            try:
                val = f'{float(data["Rg_from_Guinier"]):.2f} nm'
            except ValueError:
                val = utility.NA
            rg_table['R<sub>g</sub>'].append(val)

            try:
                val = f'{float(data["Rg_from_Guinier_error"]):.2f} nm'
            except ValueError:
                val = utility.NA
            rg_table['R<sub>g</sub> error'].append(val)

            try:
                val = f'{float(data["MW_standard"]):.1f} kDa'
            except ValueError:
                val = utility.NA
            rg_table['MW'].append(val)

            try:
                val = f'{float(data["MW_standard_error"]):.1f} kDa'
            except ValueError:
                val = utility.NA
            rg_table['MW error'].append(val)

        return rg_table

    def get_fits_for_plot(self) -> dict:
        '''
        get chi-squared values from SASCIF files
        '''
        fit_dict = {}
        for code in self.sascif_dicts.keys():
            fits = []
            sascif = self.sascif_dicts[code]

            for k, data in sascif.items():
                if re.search('FIT', k):
                    chisq = round(float(data['_sas_model_fitting_details']['chi_square']), 2)
                    fits.append(chisq)
            if len(fits) > 0:
                fit_dict[code] = fits
        return fit_dict

    def get_pofr(self) -> dict:
        '''
        get pair-dist distribution from SASCIF files
        '''
        pofr_dict = {}
        for code in self.sascif_dicts.keys():
            data = {}
            sascif = self.sascif_dicts[code]
            main = f'{code}_MAIN'
            data = sascif[main]['_sas_p_of_R']

            try:
                I_df = pd.DataFrame({
                    'id': np.array(data['id'], dtype=int),
                    'ordinal': np.array(data['ordinal'], dtype=int),
                    'R': np.array(data['r'], dtype=float),
                    'P': np.array(data['P'], dtype=float),
                    'E': np.array(data['P_error'], dtype=float),
                })
            except:
                raise(ValueError("Can't parse sas_p_of_R from {code}"))

            pofr_dict[code] = I_df
        return pofr_dict

    def get_pvals(self) -> dict:
        '''
        get p-values from ATSAS
        '''
        num_of_fits = self.get_number_of_fits()
        pval_table = {'SASDB ID': [], 'Model': [], 'χ²': [], 'p-value': []}


        for code in self.sascif_dicts.keys():
            f1fn = f'{code}_fit1.dat'
            f2fn = f'{code}_fit2.dat'
            f3fn = f'{code}_pval.dat'
            sascif = self.sascif_dicts[code]
            main = f'{code}_MAIN'
            data = sascif[main]

            refX = np.array(data['_sas_scan_intensity']['momentum_transfer'], dtype=float)
            refY = np.array(data['_sas_scan_intensity']['intensity'], dtype=float)
            refS = np.array(data['_sas_scan_intensity']['intensity_su_counting'], dtype=float)

            fits = []

            c = 0
            for k, v in sascif.items():

                if re.search('FIT', k):
                    c += 1
                    chisq = round(float(v['_sas_model_fitting_details']['chi_square']), 2)

                    pval_table['SASDB ID'].append(code)
                    pval_table['Model'].append(c)

                    fitX = np.array(v['_sas_model_fitting']['momentum_transfer'], dtype=float)
                    fit_refY = np.array(v['_sas_model_fitting']['intensity'], dtype=float)
                    fitY = np.array(v['_sas_model_fitting']['fit'], dtype=float)

                    fit_1 = pd.DataFrame({
                        'Q': fitX,
                        'Ie': fit_refY
                    })

                    fit_2 = pd.DataFrame({
                        'Q': fitX,
                        'Ib': fitY
                    })

                    fit_1.to_csv(f1fn, header=False, index=False)
                    fit_2.to_csv(f2fn, header=False, index=False)
                    with open(f3fn, 'w+') as f:
                       run(['datcmp', f1fn, f2fn], stdout=f)
                    with open(f3fn, 'r') as f:
                      all_lines = [j.strip().split()
                                 for i, j in enumerate(f.readlines())]
                    p_val = [all_lines[i+1][4]
                             for i, j in enumerate(all_lines) if 'adj' in j][0]

                    for fn in [f1fn, f2fn, f3fn]:
                        os.remove(fn)

                    pval_table['p-value'].append('%.2E' % Decimal(p_val))
                    pval_table['χ²'].append('%.2f' % chisq)

            if c == 0:
                pval_table['SASDB ID'].append(code)
                pval_table['Model'].append(utility.NA)
                pval_table['p-value'].append(utility.NA)
        return pval_table

    def get_pofr_ext(self) -> dict:
        '''
        get pair-distance details from SASCIF files
        '''
        pofr_dict = {}
        for code in self.sascif_dicts.keys():
            sascif = self.sascif_dicts[code]
            main = f'{code}_MAIN'
            data = self.sascif_dicts[code][f'{code}_MAIN']
            unit = data['_sas_scan']['unit']
            unitm = self.get_scan_unit_mult(unit)

            data = sascif[main]['_sas_p_of_R_extrapolated_intensity']

            pdf_re = pd.DataFrame({
                'Q': np.array(data['momentum_transfer'], dtype=float),
                'I': np.array(data['intensity_reg'], dtype=float)
            })
            pdf_re['Q'] = pdf_re['Q'] * unitm
            pdf_re['logI'] = np.log(pdf_re['I'])
            pofr_dict[code] = pdf_re
        return pofr_dict

    def get_pofr_errors(self) -> dict:
        '''
        get pair-distance details and errors from JSON files
        '''
        pofr_dict = self.get_pofr_ext()
        Int_dict = self.intensities
        compiled_dict = {}
        for code in self.sascif_dicts.keys():
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
        list_sort = sorted(list_sub, key=operator.itemgetter(1))
        if list_sort[0][1] < 0.00001:
            return listn[list_sort[0][0]]
        else:
            return 9999

    @staticmethod
    def get_scan_unit_mult(unit) -> int:
        return {'1/A': 10.0, '1/nm': 1.0}[unit]

    def get_Guinier_data(self) -> (dict, dict):
        '''
        get Guinier plot data from JSON files
        '''
        Int_dict = self.intensities
        Guinier_dict = {}
        Guinier_score = {}

        for code, val in Int_dict.items():
            data = self.sascif_dicts[code][f'{code}_MAIN']
            unit = data['_sas_scan']['unit']
            unitm = self.get_scan_unit_mult(unit)
            rg = float(data['_sas_result']['Rg_from_PR'])

            G_df = val.astype({'Q': float, 'I': float, 'E': float})
            G_df['logI'] = np.log(G_df['I'])
            # dmax = float(data_dic[code]['pddf_dmax'])
            # index_low = int(data_dic[code]['guinier_point_first'])
            # index_high = int(data_dic[code]['guinier_point_last'])
            # q_min = math.pi/(dmax*10)
            q_max = 1.3 / rg
            G_df_range = G_df[G_df['Q'] < q_max].copy()
            G_df_range['Q'] = G_df['Q']
            G_df_range['Q2'] = G_df_range['Q']**2
            X = G_df_range[['Q2']].values
            y = G_df_range['logI'].values
            regression = LinearRegression(fit_intercept=True)
            regression.fit(X, y)
            G_df_range['y_pred'] = regression.predict(X)
            G_df_range['res'] = y-regression.predict(X)
            # G_df_range['Q2A'] = G_df_range['Q2']
            score = '%.2f' % regression.score(X, y)
            Guinier_score[code] = score
            Guinier_dict[code] = G_df_range
        return Guinier_score, Guinier_dict

    def get_parameters_vol_many(self) -> dict:
        '''
        get volume details from SASCIF files
        '''
        parameter_table = {'SASDB ID': [], 'Estimated Volume': [], 'Porod Volume': [], 'Specific Volume': [],
                           'Sample Contrast': [], 'Sample Concentration': []}

        for code in self.sascif_dicts.keys():
            sascif = self.sascif_dicts[code]
            main = f'{code}_MAIN'

            parameter_table['SASDB ID'].append(code)

            data = sascif[main]['_sas_sample']

            try:
                val = f'{float(data["specimen_concentration"]):.2f}  mg/mL'
            except ValueError:
                val = utility.NA
            parameter_table['Sample Concentration'].append(val)

            try:
                val = f'{float(data["contrast"]):.2f}'
            except ValueError:
                val = utility.NA
            parameter_table['Sample Contrast'].append(val)

            try:
                val = f'{float(data["specific_vol"]):.2f} nm\u00b3'
            except ValueError:
                val = utility.NA
            parameter_table['Specific Volume'].append(val)

            data = sascif[main]['_sas_result']

            try:
                val = f'{float(data["Porod_volume"]):.2f} nm\u00b3'
            except ValueError:
                val = utility.NA
            parameter_table['Porod Volume'].append(val)

            try:
                val = f'{float(data["estimated_volume"]):.2f} nm\u00b3'
            except ValueError:
                val = utility.NA
            parameter_table['Estimated Volume'].append(val)

        return parameter_table

    def get_parameters_mw_many(self) -> dict:
        '''
        get MW details from SASCIF files
        '''
        parameter_table = {'SASDB ID': [], 'Chemical composition MW': [
        ], 'Standard MW': [], 'Porod Volume/MW': []}

        for code in self.sascif_dicts.keys():
            sascif = self.sascif_dicts[code]
            main = f'{code}_MAIN'

            data = sascif[main]['_sas_result']

            parameter_table['SASDB ID'].append(code)

            try:
                val = f'{float(data["experimental_MW"]):.1f} kDa'
            except ValueError:
                val = utility.NA
            parameter_table['Chemical composition MW'].append(val)

            try:
                val = f'{float(data["MW_standard"]):.1f} kDa'
            except ValueError:
                val = utility.NA
            parameter_table['Standard MW'].append(val)

            try:
                Porod_MW = round(float(data['MW_Porod']), 2)
                Porod_V = round(float(data['Porod_volume']), 2)
                val = f'{(Porod_V / Porod_MW):.2f} nm\u00b3/kDa'
            except ValueError:
                val = utility.NA
            parameter_table['Porod Volume/MW'].append(val)

        return parameter_table

    def get_pddf(self) -> dict:
        '''
        get p(r) data from JSON
        '''
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
        pddf_info = {'SASDB ID': [], 'Software used': [],
                     'D<sub>max</sub>': [], 'D<sub>max</sub> error': [], 'R<sub>g</sub>': [], 'R<sub>g</sub> error': []}

        for code in self.sascif_dicts.keys():
            sascif = self.sascif_dicts[code]
            main = f'{code}_MAIN'
            data = sascif[main]['_sas_p_of_R_details']

            pddf_info['Software used'].append(str(data['software_p_of_R']))

            data = sascif[main]['_sas_result']

            try:
                val = f'{float(data["D_max"]):.3f} nm'
            except ValueError:
                val = utility.NA
            pddf_info['D<sub>max</sub>'].append(val)

            try:
                val = f'{float(data["Rg_from_PR"]):.3f} nm'
            except ValueError:
                val = utility.NA
            pddf_info['R<sub>g</sub>'].append(val)

            try:
                val = f'{float(data["Dmax_error"]):.3f} nm'
            except ValueError:
                val = utility.NA
            pddf_info['D<sub>max</sub> error'].append(val)

            try:
                val = f'{float(data["Rg_from_PR_error"]):.3f} nm'
            except ValueError:
                val = utility.NA
            pddf_info['R<sub>g</sub> error'].append(val)

            pddf_info['SASDB ID'].append(code)
        return pddf_info

    def get_number_of_fits(self) -> dict:
        '''
        get number of fits
        '''
        num_of_fits = {}
        for code in self.sascif_dicts.keys():
            c = 0
            sascif = self.sascif_dicts[code]

            for k in sascif.keys():
                if re.search('FIT', k):
                    c += 1

            num_of_fits[code] = c
        return num_of_fits

    def get_total_number_of_fits(self) -> int:
        c = 0
        for k, v in self.get_number_of_fits().items():
            c += v
        return c

    def get_sasdb_code_fits(self) -> list:
        '''
        get asumber of fits per SASBDB ID
        '''
        fit_dict = self.get_number_of_fits()
        return list(fit_dict.values())

    def get_fit_data(self) -> dict:
        '''
        get fit information to make plots
        '''
        num_of_fits = self.get_number_of_fits()
        data_fit = {}

        for code in self.sascif_dicts.keys():
            sascif = self.sascif_dicts[code]
            main = f'{code}_MAIN'
            data = sascif[main]

            refX = np.array(data['_sas_scan_intensity']['momentum_transfer'], dtype=float)
            refY = np.array(data['_sas_scan_intensity']['intensity'], dtype=float)
            refS = np.array(data['_sas_scan_intensity']['intensity_su_counting'], dtype=float)

            num = num_of_fits[code]
            fits = {}

            c = 0
            for k, data in sascif.items():

                if re.search('FIT', k):
                    fitX = np.array(data['_sas_model_fitting']['momentum_transfer'], dtype=float)
                    fit_refY = np.array(data['_sas_model_fitting']['intensity'], dtype=float)
                    fitY = np.array(data['_sas_model_fitting']['fit'], dtype=float)
                    fitS = np.array([refS[np.argmin(np.abs(refX - x))] for x in fitX], dtype=float)

                    f_df = pd.DataFrame({
                        'Q': fitX,
                        'Ie': fit_refY,
                        'Ib': fitY,
                        'E': fitS,
                    })

                    f_df['logIe'] = np.log(f_df['Ie'])
                    f_df['logIb'] = np.log(f_df['Ib'])
                    f_df['r'] = f_df['Ie']-f_df['Ib']

                    if f_df['E'].isnull().values.any():
                        f_df['rsigma'] = 0
                    else:
                        f_df['rsigma'] = f_df['r']/f_df['E']

                    f_df['logr'] = f_df['logIe']-f_df['logIb']
                    f_df['r2a'] = (f_df['Ib']-f_df['Ie'].mean())**2
                    f_df['r2b'] = (f_df['Ie']-f_df['Ie'].mean())**2
                    fits[c] = (self.get_fit_r2(f_df), f_df)
                    c += 1

            if c == 0:
                fits[0] = (0, pd.DataFrame())

            data_fit[code] = fits
        return data_fit

    def get_fit_r2(self, df: pd.DataFrame) -> int:
        rsquared = df['r2a'].sum()/df['r2b'].sum()
        return round(rsquared, 2)
