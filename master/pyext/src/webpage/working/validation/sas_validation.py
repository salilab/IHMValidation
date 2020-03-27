import pandas as pd
import glob
import sys,os,math
import numpy as np
import pandas as pd
import validation
import ihm
import ihm.reader
import math
import re,pickle,requests,json
import multiprocessing as mp
from sklearn.linear_model import LinearRegression
from decimal import Decimal

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
                #print ('fetching data from: %s'%(url_f))
                response=requests.get(url_f, data={'key':'value'})
                if response.status_code!=200:
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
        if intensities.status_code != 200:
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
        '''
        min_l=[I_df['Q'].min(),I_df['I'].min(),I_df['logI'].min(),I_df['logQ'].min()]
        min_abs=[abs(i) for i in min_l]
        print ("min df value",min_abs)
        if min(min_abs)<0.01:
            print ("min value is less than 0.01, rounding intensity data to three decimal places...")
            I_df=I_df.round(decimals=3)
        elif min(min_abs)<0.001:
            print ("min value is less than 0.001, rounding intensity data to four decimal places...")
            I_df=I_df.round(decimals=4)
        else:
            print ("rounding intensity data two decimal places...")
            I_df=I_df.round(decimals=2)
        #print (I_df.head())
        '''
        return I_df

    def get_Guinier_data(self):
        data_dic=self.get_data_from_SASBDB()
        target_url=list(data_dic.values())[0]['intensities_data']
        intensities = requests.get(target_url)
        if intensities.status_code != 200:
            print ("Error....unable to fetch data from SASBDB, please check the entry ID")
        with open ('intensities.csv','w') as f:
            f.write(intensities.text)
        G_df=pd.read_csv('intensities.csv', skiprows=4,delim_whitespace=True, names=['Q','I','E'])
        G_df=G_df.astype({'Q':float,'I':float,'E':float})
        G_df['logI']=np.log(G_df['I'])
        rg=float(list(data_dic.values())[0]['pddf_rg'])
        dmax=float(list(data_dic.values())[0]['pddf_dmax'])
        q_min=math.pi/(dmax*10)
        q_max=1.3/(rg*10)
        G_df_range=G_df[G_df['Q']<q_max].copy()
        G_df_range['Q2']=G_df_range['Q']**2
        X=G_df_range[['Q2']].values
        y=G_df_range['logI'].values
        regression=LinearRegression(fit_intercept=True)
        regression.fit(X,y)
        G_df_range['y_pred']=regression.predict(X)
        G_df_range['res']=y-regression.predict(X)
        G_df_range['Q2A']=G_df_range['Q2']*100
        score='%.2f' %regression.score(X,y)
        '''
        min_list=[G_df_range['Q2'].min(),G_df_range['logI'].min(),G_df_range['y_pred'].min(),G_df_range['res'].min()]
        min_abs=[abs(i) for i in min_list]
        if min(min_abs)<0.01:
            print ("min value is less than 0.01, rounding Guinier data to three decimal places...")
            G_df_range=G_df_range.round(decimals=3)
        elif min(min_abs)<0.001:
            print ("min value is less than 0.001, rounding Guinier data to four decimal places...")
            G_df_range=G_df_range.round(decimals=4)
        else:
            print ("rounding Guinier data two decimal places...")
            G_df_range=G_df_range.round(decimals=2)        
        '''
        return score,G_df_range
    
    def get_experiment_description(self):
        data_dic=self.get_data_from_SASBDB()
        description=['Data description: '+ str(list(data_dic.values())[0]['experiment_description'])]
        return description

    def get_parameters_vol(self):
        data_dic=self.get_data_from_SASBDB()
        parameter_table={'Estimated volume':[],'Estimated volume method':[],'Porod volume':[]}
        parameter_table['Estimated volume'].append(list(data_dic.values())[0]['estimated_volume'])
        parameter_table['Estimated volume method'].append(list(data_dic.values())[0]['estimated_volume_method'])
        if list(data_dic.values())[0]['porod_volume'] is None:
            parameter_table['Porod volume'].append(list(data_dic.values())[0]['porod_volume'])
        else:
            parameter_table['Porod volume'].append(list(data_dic.values())[0]['porod_volume']+' nm\u00b3')
        return parameter_table

    def get_parameters_mw(self):
        data_dic=self.get_data_from_SASBDB()
        parameter_table={'Molecule MW':[],'Experimental MW':[],'Porod MW':[],'Guinier MW':[]}
        parameter_table['Experimental MW'].append(list(data_dic.values())[0]['experimental_mw']+' kDa')
        parameter_table['Guinier MW'].append(str(list(data_dic.values())[0]['guinier_i0_mw'])+' kDa')
        parameter_table['Porod MW'].append(list(data_dic.values())[0]['porod_mw']+' kDa')
        parameter_table['Molecule MW'].append(list(data_dic.values())[0]['experiment']['sample']['molecule'][0]['total_mw']+' kDa')
        return parameter_table

    def get_pddf(self):
        data_dic=self.get_data_from_SASBDB()
        target_url=list(data_dic.values())[0]['pddf_data']
        pddf = requests.get(target_url)
        if pddf.status_code != 200:
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
        '''
        min_list=[pd_df['R'].min(),pd_df['P'].min()]
        min_abs=[abs(i) for i in min_list]
        if min(min_abs)<0.01:
            print ("min value is less than 0.01, rounding pair distance data to three decimal places...")
            pd_df=pd_df.round(decimals=3)
        elif min(min_abs)<0.001:
            print ("min value is less than 0.001, rounding pair distance data to four decimal places...")
            pd_df=pd_df.round(decimals=4)
        else:
            print ("rounding pair distance data two decimal places...")
            pd_df=pd_df.round(decimals=2)
        '''
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
                chi_value=list(data_dic.values())[0]['fits'][i]['chi_square_value']
                chi_value_round=round(chi_value,2)
                print (chi_value,chi_value_round,"chi value after rounding")
                chi_table['\u03C7'+'\u00b2'].append(round(list(data_dic.values())[0]['fits'][i]['chi_square_value'],2))
        return chi_table

    def get_rg_table(self):
        data_dic=self.get_data_from_SASBDB()
        rg_table={'Rg from Guinier analysis':[],'Rg from P(r) plot':[]}
        rg_table['Rg from Guinier analysis'].append(str(round(float(list(data_dic.values())[0]['guinier_rg']),2))+ ' nm')
        rg_table['Rg from P(r) plot'].append(str(round(float(list(data_dic.values())[0]['pddf_rg']),2))+' nm') 
        return rg_table

    def get_fit_data(self):
        data_dic=self.get_data_from_SASBDB()
        number=self.get_number_of_fits()
        if number>0:
            target_url=list(data_dic.values())[0]['fits'][0]['fit_data']
            fit = requests.get(target_url)
            if fit.status_code !=200:
                print ("Error....unable to fetch data from SASBDB, please check the entry ID")
            with open ('fit.csv','w') as f:
                f.write(fit.text)
            f_df=pd.read_csv('fit.csv', skiprows=3,delim_whitespace=True, names=['Q','Ie','Ib','E'])
            f_df['logIe']=np.log(f_df['Ie'])
            f_df['logIb']=np.log(f_df['Ib'])
            f_df['r']=f_df['Ie']-f_df['Ib']
            f_df['logr']=f_df['logIe']-f_df['logIb']
            f_df['r2a']=(f_df['Ib']-f_df['Ie'].mean())**2
            f_df['r2b']=(f_df['Ie']-f_df['Ie'].mean())**2
            return f_df
        else:
            return None

    def get_fit_r2(self):
        number=self.get_number_of_fits()
        if number>0:
            f_df=self.get_fit_data()
            r2=f_df['r2a'].sum()/f_df['r2b'].sum()
            return round(r2,2)
        else:
            return None

    def get_fit_image(self):
        data_dic=self.get_data_from_SASBDB()
        target_url=list(data_dic.values())[0]['fits'][0]['models'][0]['model_plot']
        fit = requests.get(target_url)
        if fit.status_code !=200:
            print ("Error....unable to fetch data from SASBDB, please check the entry ID")
        dirname=os.path.dirname(os.path.abspath(__file__))
        filename = os.path.abspath(os.path.join(dirname, '../static/images/',self.ID+'fit.png'))
        with open (filename,'wb') as f:
            f.write(fit.content)

