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
            url=url+i+'.json'
            print ('fetching data from: %s'%(url))
            response=requests.get(url, data={'key':'value'})
            if response.status_code==200:
                print ("fetched data from SASBDB")
            else:
                print ("Error....unable to fetch data from SASBDB, please check the entry ID")
            data_dic[i]=response.json()
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
        parameter_table['Porod volume'].append(list(data_dic.values())[0]['porod_volume']+' nm\u00b3')
        return parameter_table

