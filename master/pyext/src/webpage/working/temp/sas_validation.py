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
        #print (SAS_db_codes)
        return SAS_db_codes
    
    def get_data_from_SASBDB(self):
        url='https://www.sasbdb.org/rest-api/entry/summary/'
        data_dic={}
        for i in self.get_SASBDB_code():
            if 'None' not in str(i):
                url_f=url+i+'.json'
                #print ('fetching data from: %s'%(url_f))
                response=requests.get(url_f, data={'key':'value'});
                if response.status_code!=200:
                    print ("Error....unable to fetch data from SASBDB, please check the entry ID")
                data_dic[i]=response.json();
                with open (i+'.json', 'w') as f:
                    formatted_data=json.dumps(response.json(), indent = 4, sort_keys=True);
                    f.write(formatted_data);
        return data_dic

    def get_intensities(self):
        data_dic=self.get_data_from_SASBDB()
        Int_dict={}
        for key,val in data_dic.items():
            target_url=val['intensities_data']
            intensities = requests.get(target_url);
            if intensities.status_code != 200:
                print ("Error....unable to fetch data from SASBDB, please check the entry ID")
            with open (key+'intensities.csv','w') as f:
                f.write(intensities.text) 
            #f=open(key+'intensities.csv');
            #lines=f.readlines()
            #fo=open(key+'intensities2.csv','w')
            #for i,j in enumerate(lines):
            #    if 'REMARK' not in j:
            #        fo.write(j)
            I_df=pd.read_csv(key+'intensities.csv', skiprows=4,delim_whitespace=True, names=['Q','I','E','NAN1','NAN2'])
            I_df=I_df.astype({'Q':float,'I':float,'E':float})
            I_df['err_x']=I_df.apply(lambda row: (row['Q'],row['Q']), axis=1)
            I_df['err_y']=I_df.apply(lambda row: (np.log(row['I']-row['E']),np.log(row['I']+row['E'])),axis=1)
            I_df['logI']=np.log(I_df['I'])
            I_df['logQ']=np.log(I_df['Q'])
            I_df['logX']=I_df.apply(lambda row: (row['logQ'],row['logQ']), axis=1)
            I_df['Ky']=I_df['Q']*I_df['Q']*I_df['I']
            I_df['Px']=I_df['Q']**4
            I_df['Px'].round(3)
            I_df['Py']=I_df['Px']*I_df['I']
            Int_dict[key]=I_df
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
        return Int_dict

    def get_Guinier_data(self):
        data_dic=self.get_data_from_SASBDB()
        Gui_dict={};Gui_score={}
        for key,val in data_dic.items():
            target_url=val['intensities_data']
            intensities = requests.get(target_url)
            if intensities.status_code != 200:
                print ("Error....unable to fetch data from SASBDB, please check the entry ID")
            with open (key+'intensities.csv','w') as f:
                f.write(intensities.text)
            f=open(key+'intensities.csv');
            lines=f.readlines()
            print (key, "intensity fetched")
            #fo=open(key+'intensities2.csv','w')
            #for i,j in enumerate(lines):
            #   if 'REMARK' not in j:
            #        print (key,j)
            #        fo.write(j)
            G_df=pd.read_csv(key+'intensities.csv', skiprows=4,delim_whitespace=True, names=['Q','I','E','NAN1','NAN2'])
            print (G_df.head())
            G_df=G_df.astype({'Q':float,'I':float,'E':float})
            print (G_df.head())
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
            Gui_score[key]=score
            Gui_dict[key]=G_df_range
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
        return Gui_score,Gui_dict
    
    def get_experiment_description(self):
        data_dic=self.get_data_from_SASBDB()
        description=['Data description: '+ str(list(data_dic.values())[0]['experiment_description'])]
        return description


    def get_parameters_vol_many(self):
        data_dic=self.get_data_from_SASBDB()
        parameter_table={'SASDB ID':[],'Estimated volume':[],'Estimated volume method':[],'Porod volume':[]}
        for key,val in data_dic.items():
            try:
                parameter_table['Estimated volume'].append(val['estimated_volume'])
                parameter_table['Estimated volume method'].append(val['estimated_volume_method'])
            except:
                parameter_table['Estimated volume'].append('None')
                parameter_table['Estimated volume method'].append('None')
            try:
                parameter_table['Porod volume'].append(val['porod_volume']+' nm\u00b3')
            except:
                parameter_table['Porod volume'].append('None')
            parameter_table['SASDB ID'].append(key)
        return parameter_table


    def get_parameters_mw_many(self):
        data_dic=self.get_data_from_SASBDB()
        parameter_table={'SASDB ID':[],'Molecule MW':[],'Experimental MW':[],'Porod MW':[]}        
        for key,val in data_dic.items():
            try:
                parameter_table['Experimental MW'].append(list(data_dic.values())[0]['experimental_mw']+' kDa')
            except:
                parameter_tabel['Experimental MW'].append('None')
            try:
                parameter_table['Porod MW'].append(list(data_dic.values())[0]['porod_mw']+' kDa')
            except:
                parameter_table['Porod MW'].append('None')
            try:
                parameter_table['Molecule MW'].append(list(data_dic.values())[0]['experiment']['sample']['molecule'][0]['total_mw']+' kDa')
            except:
                parameter_table['Molecule MW'].append('None')
            parameter_table['SASDB ID'].append(key)
        return parameter_table

    def get_pddf(self):
        data_dic=self.get_data_from_SASBDB()
        pddf_dic={}
        for key,val in data_dic.items():
            target_url=val['pddf_data']
            pddf = requests.get(target_url);
            if pddf.status_code != 200:
                print ("Error....unable to fetch data from SASBDB, please check the entry ID")
            with open (key+'pddf.csv','w') as f:
                f.write(pddf.text)
            f=open(key+'pddf.csv');
            lines=f.readlines()
            fo=key+'pddf_output.csv'
            for i,j in enumerate(lines):
                if 'Distance distribution  function of particle' in j:
                    list_f=lines[i:-2]
                    print (''.join(list_f),file=open(fo,'w'))
            pd_df=pd.read_csv(fo, skiprows=5,delim_whitespace=True, names=['R','P','E'])
            pd_df=pd_df.astype({'P':float,'R':float,'E':float})
            pd_df['R']=pd_df['R']/10;
            pd_df['err_x']=pd_df.apply(lambda row: (row['R'],row['R']), axis=1)
            pd_df['err_y']=pd_df.apply(lambda row: (row['P']-row['E'],row['P']+row['E']),axis=1)
            pddf_dic[key]=pd_df
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
        return pddf_dic

    def get_pddf_info(self):
        data_dic=self.get_data_from_SASBDB()
        pddf_info={'SASDB ID':[],'Software used':[],'Dmax':[],'Dmax error':[],'Rg':[],'Rg error':[]}
        for key,val in data_dic.items():
            pddf_info['Software used'].append(str(val['pddf_software']))
            pddf_info['Dmax'].append(str(val['pddf_dmax'])+' nm')
            pddf_info['Rg'].append(str(val['pddf_rg'])+' nm')
            try:
                pddf_info['Dmax error'].append(str(val['pddf_dmax_error'])+' nm')
            except:
                pddf_info['Dmax error'].append('None')
            try:
                pddf_info['Rg error'].append(str(val['pddf_rg_error'])+' nm')
            except:
                pddf_info['Rg error'].append('None')

            pddf_info['SASDB ID'].append(key)
        return pddf_info

    def get_pddf_dmax(self):
        data_dic=self.get_data_from_SASBDB()
        return str(list(data_dic.values())[0]['pddf_dmax'])

    def get_pddf_rg(self):
        data_dic=self.get_data_from_SASBDB()
        return str(list(data_dic.values())[0]['pddf_rg'])

    def get_number_of_fits(self):
        data_dic=self.get_data_from_SASBDB()
        num_of_fits={}
        for key,val in data_dic.items():
            num_of_fits[key]=len(val['fits'])
        return num_of_fits

    def get_total_fits(self):
        data_dic=self.get_data_from_SASBDB()
        num=0
        for key,val in data_dic.items():
            num += len(val['fits'])
        return num


    def get_chi_table(self):
        data_dic=self.get_data_from_SASBDB()
        chi_table={'SASDB ID':[],'Model':[],'\u03C7\u00b2':[]}        
        for key,val in data_dic.items():
            print (self.get_number_of_fits())
            num=self.get_number_of_fits()[key]
            if num>0:
                for i in range(0,num):
                    count=i+1
                    chi_table['SASDB ID'].append(key)            
                    chi_table['Model'].append(str(count))
                    chi_value=val['fits'][i]['chi_square_value']
                    chi_value_round=round(chi_value,2)
                    print (chi_value,chi_value_round,"chi value after rounding")
                    chi_table['\u03C7'+'\u00b2'].append(chi_value_round)
            else:
                chi_table['SASDB ID'].append(key)
                chi_table['Model'].append('None')
                chi_table['\u03C7'+'\u00b2'].append('None')
        return chi_table

    def get_sasdb_code_fits(self):
        fit_dict=self.get_number_of_fits()
        print (list(fit_dict.values()))
        return list(fit_dict.values())


    def get_rg_table_many(self):
        data_dic=self.get_data_from_SASBDB()
        rg_table={'SASDB ID':[],'Rg':[],'Rg error':[],'MW':[],'MW error':[]}        
        for key,val in data_dic.items():
            rg_table['Rg'].append(str(round(float(val['guinier_rg']),2))+ ' nm')
            try:
                rg_table['Rg error'].append(val['guinier_rg_error']+ ' nm')
            except:
                rg_table['Rg error'].append('None')
            try:
                rg_table['MW'].append(val['guinier_i0_mw']+ ' nm')
            except:
                rg_table['MW'].append('None')
            try:
                rg_table['MW error'].append(val['guinier_i0_mw_error']+ ' nm')
            except:
                rg_table['MW error'].append('None')
            rg_table['SASDB ID'].append(key)
        return rg_table

    def get_fit_data(self):
        data_dic=self.get_data_from_SASBDB()
        num_of_fits=self.get_number_of_fits()
        data_fit={}
        for key,val in data_dic.items():
            num=num_of_fits[key]
            fits={}
            if num>0:
                for i in range(0,num):
                    target_url=val['fits'][i]['fit_data']      
                    fit = requests.get(target_url);
                    if fit.status_code !=200:
                        print ("Error....unable to fetch data from SASBDB, please check the entry ID")
                    fname=key+str(i)+'fit.csv'
                    print (fname)
                    with open (fname,'w') as f:
                        f.write(fit.text)
                    f_df=pd.read_csv(fname, skiprows=3,delim_whitespace=True, names=['Q','Ie','Ib','E'])
                    if abs(f_df.iloc[22,2]-f_df.iloc[22,1])>abs(f_df.iloc[22,3]-f_df.iloc[22,1]):
                        f_df.rename(columns={'Q':'Q','Ie':'Ie','Ib':'E','E':'Ib'},inplace=True)
                    print (f_df.head())
                    f_df['logIe']=np.log(f_df['Ie'])
                    f_df['logIb']=np.log(f_df['Ib'])
                    f_df['r']=f_df['Ie']-f_df['Ib']
                    f_df['logr']=f_df['logIe']-f_df['logIb']
                    f_df['r2a']=(f_df['Ib']-f_df['Ie'].mean())**2
                    f_df['r2b']=(f_df['Ie']-f_df['Ie'].mean())**2
                    fits[i]=(self.get_fit_r2(f_df),f_df)
            else:
                fits[0]=(0,pd.DataFrame())
                #data_fit=None
            data_fit[key]=fits 
        return data_fit

    def get_fit_r2(self,df):
        r2=df['r2a'].sum()/df['r2b'].sum()
        return round(r2,2)

    def get_fit_image(self):
        data_dic=self.get_data_from_SASBDB()
        num_of_fits=self.get_number_of_fits()
        data_fit={}
        for key,val in data_dic.items():
            num=num_of_fits[key]
            if num>0:
                for i in range(0,num):    
                    target_url=val['fits'][i]['models'][0]['model_plot']
                    fit = requests.get(target_url);
                    if fit.status_code !=200:
                        print ("Error....unable to fetch data from SASBDB, please check the entry ID")
                    dirname=os.path.dirname(os.path.abspath(__file__))
                    filename = os.path.abspath(os.path.join(dirname, '../static/images/',self.ID+key+str(i)+'fit.png'))
                    with open (filename,'wb') as f:
                        f.write(fit.content)

