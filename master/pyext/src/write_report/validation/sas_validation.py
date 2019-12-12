import pandas as pd
import glob
import sys,os,math
import numpy as np
import pandas as pd
import validation
import ihm
import ihm.reader
import re,pickle,requests,json
import multiprocessing as mp

class sas_validation(validation.get_input_information):
    def __init__(self,mmcif_file):
        super().__init__(mmcif_file)
        self.ID=str(validation.get_input_information.get_id(self))
        self.nos=validation.get_input_information.get_number_of_models(self)
        self.dataset=validation.get_input_information.get_dataset_comp(self) 
    
    def get_SASBDB_code(self):
        SAS_db_codes=[]
        for i,j in enumerate(self.dataset['Dataset type']):
            if 'SAS' in str(j):
                SAS_db_codes.append(self.dataset['Data access code'][i])
        return SAS_db_codes
    
    def get_data_from_SASBDB(self):
        url='https://www.sasbdb.org/rest-api/entry/summary/'
        data_dic={}
        for i in self.get_SASBDB_code():
            if 'None' not in str(i):
                url_f=url+i+'.json'
                print ('fetching data from: %s'%(url_f))
                response=requests.get(url_f, data={'key':'value'})
                if response.status_code==200:
                    print ("fetched data from SASBDB")
                else:
                    print ("Error....unable to fetch data from SASBDB, please check the entry ID")
                data_dic[i]=response.json()
                print (data_dic)
                with open (i+'.json', 'w') as f:
                    formatted_data=json.dumps(response.json(), indent = 4, sort_keys=True)
                    f.write(formatted_data)
        return data_dic

    def get_intensities(self):
        data_dic=self.get_data_from_SASBDB()
        target_url=list(data_dic.values())[0]['intensities_data']
        intensities = requests.get(target_url)
        if intensities.status_code==200:
            print ("fetched data from SASBDB")
        else:
            print ("Error....unable to fetch data from SASBDB, please check the entry ID")
        with open ('intensities.csv','w') as f:
            f.write(intensities.text)
        I_df=pd.read_csv('intensities.csv', skiprows=4,delim_whitespace=True, names=['Q','I','E'])
        I_df=I_df.astype({'Q':float,'I':float,'E':float})
        I_df['err_x']=I_df.apply(lambda row: (row['Q'],row['Q']), axis=1)
        I_df['err_y']=I_df.apply(lambda row: (np.log(row['I']-row['E']),np.log(row['I']+row['E'])),axis=1)
        I_df['logI']=np.log(I_df['I'])
        I_df['logQ']=np.log(I_df['Q'])
        I_df['logX']=I_df.apply(lambda row: (row['logQ'],row['logQ']), axis=1)
        I_df['Ky']=I_df['Q']*I_df['Q']*I_df['I']
        I_df['Px']=I_df['Q']**4
        I_df['Py']=I_df['Px']*I_df['I']
        #print (I_df.head())
        return I_df
    
    def get_experiment_description(self):
        data_dic=self.get_data_from_SASBDB()
        description=['Data description: '+ str(list(data_dic.values())[0]['experiment_description'])]
        return description

    def get_parameters_vol(self):
        data_dic=self.get_data_from_SASBDB()
        print (list(data_dic.values())[0].keys())
        parameter_table={'Estimated volume':[],'Estimated volume method':[],'Porod volume':[]}
        parameter_table['Estimated volume'].append(list(data_dic.values())[0]['estimated_volume'])
        parameter_table['Estimated volume method'].append(list(data_dic.values())[0]['estimated_volume_method'])
        print (list(data_dic.values())[0]['porod_volume'])
        if list(data_dic.values())[0]['porod_volume'] is None:
            parameter_table['Porod volume'].append(list(data_dic.values())[0]['porod_volume'])
        else:
            parameter_table['Porod volume'].append(list(data_dic.values())[0]['porod_volume']+' nm\u00b3')
        return parameter_table

    def get_parameters_mw(self):
        data_dic=self.get_data_from_SASBDB()
        print (list(data_dic.values())[0].keys())
        parameter_table={'Molecule MW':[],'Experimental MW':[],'Porod MW':[],'Guinier MW':[]}
        parameter_table['Experimental MW'].append(list(data_dic.values())[0]['experimental_mw'])
        parameter_table['Guinier MW'].append(list(data_dic.values())[0]['guinier_i0_mw'])
        parameter_table['Porod MW'].append(list(data_dic.values())[0]['porod_mw'])
        parameter_table['Molecule MW'].append(list(data_dic.values())[0]['experiment']['sample']['molecule'][0]['total_mw'])
        return parameter_table

    def get_pddf(self):
        data_dic=self.get_data_from_SASBDB()
        target_url=list(data_dic.values())[0]['pddf_data']
        pddf = requests.get(target_url)
        if pddf.status_code==200:
            print ("fetched data from SASBDB")
        else:
            print ("Error....unable to fetch data from SASBDB, please check the entry ID")
        with open ('pddf.csv','w') as f:
            f.write(pddf.text)
        f=open('pddf.csv');lines=f.readlines()
        for i,j in enumerate(lines):
            if 'Real Space Data' in j:
                list_f=lines[i:]
                print (''.join(list_f),file=open('pddf_output.csv','w'))
        pd_df=pd.read_csv('pddf_output.csv', skiprows=7,delim_whitespace=True, names=['R','P','E'])
        pd_df=pd_df.astype({'P':float,'R':float,'E':float})
        pd_df['err_x']=pd_df.apply(lambda row: (row['R'],row['R']), axis=1)
        pd_df['err_y']=pd_df.apply(lambda row: (row['P']-row['E'],row['P']+row['E']),axis=1)
        return pd_df

    def get_pddf_software(self):
        data_dic=self.get_data_from_SASBDB()
        return str(list(data_dic.values())[0]['pddf_software'])

    def get_pddf_dmax(self):
        data_dic=self.get_data_from_SASBDB()
        return str(list(data_dic.values())[0]['pddf_dmax'])

    def get_pddf_rg(self):
        data_dic=self.get_data_from_SASBDB()
        return str(list(data_dic.values())[0]['pddf_rg'])

    def get_number_of_fits(self):
        data_dic=self.get_data_from_SASBDB()
        number= len(list(data_dic.values())[0]['fits'])
        return number

    def get_chi_table(self):
        data_dic=self.get_data_from_SASBDB()
        number= self.get_number_of_fits()
        chi_table={'Model':[],'\u03C7\u00b2':[]}
        if number>0:
            for i in range(0,number):
                count=i+1
                chi_table['Model'].append(str(count))
                chi_table['\u03C7\u00b2'].append(list(data_dic.values())[0]['fits'][i]['chi_square_value'])
        return chi_table

    def get_rg_table(self):
        data_dic=self.get_data_from_SASBDB()
        rg_table={'Rg from Guinier analysis':[],'Rg from P(r) plot':[]}
        rg_table['Rg from Guinier analysis'].append(list(data_dic.values())[0]['guinier_rg'])
        rg_table['Rg from P(r) plot'].append(list(data_dic.values())[0]['pddf_rg']) 
        return rg_table

