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
            data=response.json()
            data_dic[i]=response.json
        return data_dic
