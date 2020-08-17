###################################
# Script : 
# 1) Contains class to validate models
# built using SAS datasets
#
# ganesans - Salilab - UCSF
# ganesans@salilab.org
###################################
import pandas as pd
import sys,os,math
import numpy as np
import pandas as pd
import re,pickle,requests,json
from sklearn.linear_model import LinearRegression
from decimal import Decimal
from validation import get_input_information 
from subprocess import run, call, PIPE

class sas_validation(get_input_information):
    def __init__(self,mmcif_file):
        super().__init__(mmcif_file)
        self.ID=str(get_input_information.get_id(self))
        self.nos=get_input_information.get_number_of_models(self)
        self.dataset=get_input_information.get_dataset_comp(self) 
    
    def get_SASBDB_code(self):
        '''
        function to get all SASBDB codes used in the model,
        returns a list of SASBDB codes
        '''
        SAS_db_codes=[]
        for i,j in enumerate(self.dataset['Dataset type']):
            if 'SAS' in str(j):
                SAS_db_codes.append(self.dataset['Data access code'][i])
        #print (SAS_db_codes)
        return SAS_db_codes

    def clean_SASBDB_code(self):
        '''
        function to clean SASBDB list of codes
        as some might have 'None' or can be repetitive 
        '''
        codes=list(set(self.get_SASBDB_code()))
        cleaned_code=[i for i in codes if i != 'None']
        return cleaned_code
        
    def get_data_from_SASBDB(self):
        '''
        get data from JSON
        '''
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

    def get_sascif_file(self):
        '''
        get data from SASCIF files
        '''
        url='https://www.sasbdb.org/media/sascif/sascif_files/'
        for i in self.get_SASBDB_code():
            if 'None' not in str(i):
                url_f=url+i+'.sascif'
                response=requests.get(url_f);
                if response.status_code!=200:
                    print ("Error....unable to fetch data from SASBDB, please check the entry ID")
                with open (i+'.sascif', 'w') as f:
                    f.write(response.text);

    def get_all_sascif(self,sasbdb):
        '''
        get a list of all lines in a SASCIF file
        '''
        if 'None' not in str(sasbdb):
            file=open(sasbdb+'.sascif','r')
            all_lines=[]
            for i,j in enumerate(file.readlines()):
                all_lines.append(j.strip().split())
        return all_lines

    def get_intensities(self):
        '''
        get intensity data from SASCIF file
        if SASCIF file is not present, this information will not be present/used in the report
        JSON file typically has raw data 
        '''
        self.get_sascif_file()
        Int_dict={}
        for code in self.clean_SASBDB_code():
            #Int_dict={}
            all_lines=self.get_all_sascif(code)
            data={}
            for m,n in enumerate(all_lines):
                if len(n)<2 and len(n)>0 and 'scan_intensity' in n[0]:
                    data[(n[0].split('.')[1])]=[]
                if len(n)>2 and len(all_lines[m-1])>0 and'scan_intensity' in all_lines[m-1][0]:
                    for i,j in enumerate(all_lines[m:]):
                        if len(j)>2 and '#' not in j and 'sas' not in j[0] :  
                            for num,key in enumerate(list(data.keys())):
                                data[key].append(j[num])
                        else:
                            break        
            I_df=pd.DataFrame(list(data.values()),index=list(data.keys())).T
            I_df_re=I_df[['momentum_transfer','intensity','intensity_su_counting']]
            I_df_re.rename(columns={'momentum_transfer':'Q','intensity':'I','intensity_su_counting':'E'},inplace=True)
            Int_dict[code]=I_df_re
        return Int_dict

    def modify_intensity(self):
        '''
        modify intensity data to calcualte errors and log values
        '''
        Int_dict=self.get_intensities()
        Int_dict_modify={}
        rg_and_io=self.get_rg_and_io()
        for key,val in Int_dict.items():
            Rg=rg_and_io[key][0]
            IO=rg_and_io[key][1]
            dim_num=Rg*Rg/IO
            I_df=val.astype({'Q':float,'I':float,'E':float})
            I_df.head()
            I_df=I_df[I_df['I']-I_df['E']>0]
            I_df['Q']=I_df['Q']*10
            I_df['err_x']=I_df.apply(lambda row: (row['Q'],row['Q']), axis=1) 
            I_df['err_y']=I_df.apply(lambda row: (np.log(row['I']-row['E']),np.log(row['I']+row['E'])),axis=1)
            I_df['logI']=np.log(I_df['I'])
            I_df['logQ']=np.log(I_df['Q'])
            I_df['logX']=I_df.apply(lambda row: (row['logQ'],row['logQ']), axis=1)
            I_df['Ky']=I_df['Q']*I_df['Q']*I_df['I']*dim_num
            I_df['Kx']=I_df['Q']*Rg
            I_df['Px']=I_df['Q']**4
            I_df['Px'].round(3)
            I_df['Py']=I_df['Px']*I_df['I']
            Int_dict_modify[key]=I_df
        return Int_dict_modify

    def modify_intensity_dep(self):
        '''
        depreciated function to get intensities from JSON/raw data
        '''
        Int_dict=self.get_intensities()
        Int_dict_modify={}
        for key,val in Int_dict.items():
            I_df=val.astype({'Q':float,'I':float,'E':float})
            I_df.head()
            I_df=I_df[I_df['I']-I_df['E']>0]
            I_df['Q']=I_df['Q']*10
            I_df['err_x']=I_df.apply(lambda row: (row['Q'],row['Q']), axis=1) 
            I_df['err_y']=I_df.apply(lambda row: (np.log(row['I']-row['E']),np.log(row['I']+row['E'])),axis=1)
            I_df['logI']=np.log(I_df['I'])
            I_df['logQ']=np.log(I_df['Q'])
            I_df['logX']=I_df.apply(lambda row: (row['logQ'],row['logQ']), axis=1)
            I_df['Ky']=I_df['Q']*I_df['Q']*I_df['I']
            I_df['Px']=I_df['Q']**4
            I_df['Px'].round(3)
            I_df['Py']=I_df['Px']*I_df['I']
            Int_dict_modify[key]=I_df
        return Int_dict_modify

    def get_rg_for_plot(self):
        '''
        get Rg values from SASCIF file, if unavailabel, get it from JSON
        '''
        self.get_sascif_file()
        Rg_dict={}
        for code in self.clean_SASBDB_code():
            rg={};
            all_lines=self.get_all_sascif(code)
            for m,n in enumerate(all_lines):
                if len(n)<3 and len(n)>0 and 'sas_result.Rg_from_PR' in n[0] and'sas_result.Rg_from_PR_' not in n[0] :
                    rg[n[0].split('.')[1]]=float(n[1])
                if len(n)<3 and len(n)>0 and 'sas_result.Rg_from_Guinier' in n[0] and'sas_result.Rg_from_Guinier_' not in n[0] :
                    rg[n[0].split('.')[1]]=float(n[1])
            Rg_dict[code]=list(rg.values())
        if len(list(rg.values()))<1:
            data_dic=self.get_data_from_SASBDB()
            for key,val in data_dic.items():
                Rg_dict[key]=[]
                Rg_dict[key].append(round(float(val['guinier_rg']),2))
                Rg_dict[key].append(round(float(val['pddf_rg']),2))
        return Rg_dict

    def get_rg_and_io(self):
        '''
        get rg information from SASCIF file
        '''
        self.get_sascif_file()
        rg_and_io={}
        for code in self.clean_SASBDB_code():
            all_lines=self.get_all_sascif(code)
            for m,n in enumerate(all_lines):
                if len(n)<3 and len(n)>0 and 'sas_result.Rg_from_PR' in n[0] and'sas_result.Rg_from_PR_' not in n[0] :
                    rg=float(n[1])
                if len(n)<3 and len(n)>0 and '_sas_result.I0_from_PR' in n[0] and'_sas_result.I0_from_PR_' not in n[0] :
                    io=float(n[1])
            rg_and_io[code]=(rg,io)
        return rg_and_io

    def get_rg_table_many(self):
        '''
        get rg information from multiple SASCIF files
        '''
        data_dic=self.get_data_from_SASBDB()
        rg_table={'SASDB ID':[],'Rg':[],'Rg error':[],'MW':[],'MW error':[]}        
        for key,val in data_dic.items():
            rg_table['Rg'].append(str(round(float(val['guinier_rg']),2))+ ' nm')
            try:
                rg_table['Rg error'].append(val['guinier_rg_error']+ ' nm')
            except:
                rg_table['Rg error'].append('N/A')
            try:
                rg_table['MW'].append(val['guinier_i0_mw']+ ' nm')
            except:
                rg_table['MW'].append('N/A')
            try:
                rg_table['MW error'].append(val['guinier_i0_mw_error']+ ' nm')
            except:
                rg_table['MW error'].append('N/A')
            rg_table['SASDB ID'].append(key)
        return rg_table


    def get_fits_for_plot(self):
        '''
        get chi-squared values from SASCIF files
        '''
        self.get_sascif_file()
        fit_dict={}
        for code in self.clean_SASBDB_code():
            fits=[];
            all_lines=self.get_all_sascif(code)
            for m,n in enumerate(all_lines):
                #if len(n)<3 and len(n)>0 and 'sas_model_fitting_details.id' in n[0]:
                #    fitid=int(n[1])
                if (len(n)<3) and (len(n)>0) and ('sas_model_fitting_details.chi_square' in n[0]) and (float(n[1])>0.00000):
                    fits.append(round(float(n[1]),2))
            if len(fits)>0:
                fit_dict[code]=fits
        print (fit_dict)
        return fit_dict

    def get_pofr(self):
        '''
        get pair-dist distribution from SASCIF files
        '''
        pofr_dict={}
        for code in self.clean_SASBDB_code():
            #pofr_dict={}
            all_lines=self.get_all_sascif(code)
            data={}
            for m,n in enumerate(all_lines):
                if len(n)<2 and len(n)>0 and 'sas_p_of_R.' in n[0]:
                    data[(n[0].split('.')[1])]=[]
                if len(n)>2 and len(all_lines[m-1])>0 and 'sas_p_of_R.' in all_lines[m-1][0]:
                    for i,j in enumerate(all_lines[m:]):
                        if len(j)>2 and '#' not in j and 'sas' not in j[0] :
                            for num,key in enumerate(list(data.keys())):
                                data[key].append(j[num])
                        else:
                            break
            pdf=pd.DataFrame(list(data.values()),index=list(data.keys())).T
            pdf_re=pdf[['r','P','P_error']]
            pdf_re.rename(columns={'r':'R','P':'P','P_error':'E'},inplace=True)
            pofr_dict[code]=pdf_re
        return pofr_dict

    def get_pvals(self):
        '''
        get p-values from ATSAS 
        '''
        data_dic=self.get_data_from_SASBDB()
        num_of_fits=self.get_number_of_fits()
        pval_table={'SASDB ID':[],'Model':[],'p-value':[]}        
        for key,val in data_dic.items():
            num=num_of_fits[key]
            if num>0:
                for i in range(0,num):
                    pval_table['SASDB ID'].append(key)
                    pval_table['Model'].append(i+1)
                    target_url=val['fits'][i]['fit_data']      
                    fit = requests.get(target_url);
                    if fit.status_code !=200:
                        print ("Error....unable to fetch data from SASBDB, please check the entry ID")
                    fname=key+str(i)+'fit.csv'
                    with open (fname,'w') as f:
                        f.write(fit.text)
                    f_df=pd.read_csv(fname, skiprows=3,delim_whitespace=True, names=['Q','Ie','Ib','E'])
                    if abs(f_df.iloc[22,2]-f_df.iloc[22,1])>abs(f_df.iloc[22,3]-f_df.iloc[22,1]):
                        f_df.rename(columns={'Q':'Q','Ie':'Ie','Ib':'E','E':'Ib'},inplace=True)
                    fit_1=f_df[['Q','Ie']]
                    fit_1.to_csv('fit1.csv',header=False,index=False)
                    fit_2=f_df[['Q','Ib']]
                    fit_2.to_csv('fit2.csv',header=False,index=False)
                    f1=open('pval.txt','w+')
                    with f1 as outfile:
                        run([r"/Users/saijananiganesan/Applications/ATSAS-3.0.1-1/bin/datcmp",'fit1.csv','fit2.csv'],stdout=outfile)
                    f2=open('pval.txt','r')
                    all_lines=[j.strip().split() for i,j in enumerate(f2.readlines())]
                    p_val=[all_lines[i+1][4] for i,j in enumerate(all_lines) if 'adj' in j][0]
                    pval_table['p-value'].append('%.2E' % Decimal(p_val))
            else:
                pval_table['SASDB ID'].append(key)
                pval_table['Model'].append('N/A')
                pval_table['p-value'].append('N/A')
        print (pval_table)
        return pval_table

    def get_pofr_ext(self):
        '''
        get pair-distance details from SASCIF files
        '''
        pofr_dict={}
        for code in self.clean_SASBDB_code():
            #pofr_dict={}
            all_lines=self.get_all_sascif(code)
            data={}
            for m,n in enumerate(all_lines):
                if len(n)<2 and len(n)>0 and '_sas_p_of_R_extrapolated_intensity.' in n[0]:
                    data[(n[0].split('.')[1])]=[]
                if len(n)>2 and len(all_lines[m-1])>0 and '_sas_p_of_R_extrapolated_intensity.' in all_lines[m-1][0]:
                    for i,j in enumerate(all_lines[m:]):
                        if len(j)>2 and '#' not in j and 'sas' not in j[0] :
                            for num,key in enumerate(list(data.keys())):
                                data[key].append(j[num])
                        else:
                            break
            pdf=pd.DataFrame(list(data.values()),index=list(data.keys())).T
            pdf_re=pdf[['momentum_transfer','intensity_reg']]
            pdf_re.rename(columns={'momentum_transfer':'Q','intensity_reg':'I'},inplace=True)
            pdf_re=pdf_re.astype({'Q':float,'I':float})
            pdf_re['Q']=pdf_re['Q']*10
            pdf_re['logI']=np.log(pdf_re['I'])
            pofr_dict[code]=pdf_re
            #print (pdf_re.head())
        return pofr_dict

    def get_pofr_errors(self):
        '''
        get pair-distance details and errors from JSON files
        '''
        pofr_dict=self.get_pofr_ext()
        Int_dict=self.modify_intensity()
        compiled_dict={}
        for code in self.clean_SASBDB_code():
            I_df=Int_dict[code]
            I_df_dict= dict(zip(I_df.Q, I_df.I))
            I_df_err_dict= dict(zip(I_df.Q, I_df.E))
            p_df=pofr_dict[code]
            p_df_dict=dict(zip(p_df.Q, p_df.I))
            errors=[]
            for Q,I in p_df_dict.items():
                data_Q=self.findMinDiff(list(I_df_dict.keys()),Q)
                if data_Q != 9999:
                    data_I=I_df_dict[data_Q]
                    delta_I=(I-data_I)
                    if I_df_err_dict[data_Q] != 0:
                        wt_delta_I=delta_I/I_df_err_dict[data_Q]
                    else:
                        wt_delta_I=0
                    errors.append([Q,delta_I,wt_delta_I])
            errors_df=pd.DataFrame(errors,columns=['Q','R','WR'])
            print (errors_df.head())
            compiled_dict[code]=errors_df
        return compiled_dict

    def findMinDiff(self,list1, num1): 
        '''
        quick min diff operation for calculating errors
        '''
        list_sub=[(i,abs(j-num1)) for i,j in enumerate(list1)]
        list_sort=sorted(list_sub, key=lambda x: x[1])
        if list_sort[0][1]<0.00001:
        #print (list1[list_sort[0][0]])
            return list1[list_sort[0][0]]
        else:
            return 9999  

    def get_Guinier_data(self):
        '''
        get Guinier plot data from JSON files
        '''
        Int_dict=self.get_intensities()
        data_dic=self.get_data_from_SASBDB()
        Guinier_dict={};Guinier_score={}
        for key,val in Int_dict.items():
            G_df=val.astype({'Q':float,'I':float,'E':float})
            #G_df['Q']=G_df['Q']*10
            G_df['logI']=np.log(G_df['I'])
            rg=float(data_dic[key]['pddf_rg'])
            dmax=float(data_dic[key]['pddf_dmax'])
            index_low=int(data_dic[key]['guinier_point_first'])
            index_high=int(data_dic[key]['guinier_point_last'])
            q_min=math.pi/(dmax*10)
            q_max=1.3/(rg*10)
            #G_df_range=G_df.iloc[[index_low,index_high],:] 
            G_df_range=G_df[G_df['Q']<q_max].copy()
            G_df_range['Q']=G_df['Q']*10
            G_df_range['Q2']=G_df_range['Q']**2
            X=G_df_range[['Q2']].values
            y=G_df_range['logI'].values
            regression=LinearRegression(fit_intercept=True)
            regression.fit(X,y)
            G_df_range['y_pred']=regression.predict(X)
            G_df_range['res']=y-regression.predict(X)
            G_df_range['Q2A']=G_df_range['Q2']*100
            score='%.2f' %regression.score(X,y)
            Guinier_score[key]=score
            Guinier_dict[key]=G_df_range
        return Guinier_score,Guinier_dict

    def get_parameters_vol_many_dep(self):
        '''
        get volume parameters from JSON files 
        '''
        data_dic=self.get_data_from_SASBDB()
        parameter_table={'SASDB ID':[],'Estimated volume':[],'Estimated volume method':[],'Porod volume':[]}
        for key,val in data_dic.items():
            try:
                if parameter_table['Estimated volume'] is None:
                    parameter_table['Estimated volume'].append('N/A')
                else:
                    parameter_table['Estimated volume'].append(val['estimated_volume'])
                if parameter_table['Estimated volume method'] is None:
                    parameter_table['Estimated volume method'].append('N/A')
                else:
                    parameter_table['Estimated volume method'].append(val['estimated_volume_method'])
            except:
                parameter_table['Estimated volume'].append('N/A')
                parameter_table['Estimated volume method'].append('N/A')
            try:
                parameter_table['Porod volume'].append(val['porod_volume']+' nm\u00b3')
            except:
                parameter_table['Porod volume'].append('N/A')
            parameter_table['SASDB ID'].append(key)
        return parameter_table

    def get_parameters_vol_many(self):
        '''
        get volume details from SASCIF files
        '''
        self.get_sascif_file()
        parameter_table={'SASDB ID':[],'Estimated Volume':[],'Porod Volume':[],'Specific Volume':[],
        'Sample Contrast':[],'Sample Concentration':[]}   
        for code in self.clean_SASBDB_code():
            parameter_table['SASDB ID'].append(code)
            all_lines=self.get_all_sascif(code)
            for m,n in enumerate(all_lines):
                if len(n)<3 and len(n)>0 and '_sas_sample.specimen_concentration' in n[0]:
                    if len(n[1])>1:
                        parameter_table['Sample Concentration'].append(str(round(float(n[1]),2))+' mg/ml')
                    else:
                        parameter_table['Sample Concentration'].append('N/A')
                if len(n)<3 and len(n)>0 and '_sas_sample.contrast' in n[0] :
                    if len(n[1])>1:
                        parameter_table['Sample Contrast'].append(str(round(float(n[1]),2)))
                    else:
                        parameter_table['Sample Contrast'].append('N/A')
                if len(n)<3 and len(n)>0 and '_sas_sample.specific_vol' in n[0]:
                    if len(n[1])>1:
                        parameter_table['Specific Volume'].append(str(round(float(n[1]),2))+' nm\u00b3')
                    else:
                        parameter_table['Specific Volume'].append('N/A')
                if len(n)<3 and len(n)>0 and '_sas_result.Porod_volume' in n[0] and '_sas_result.Porod_volume_error' not in n[0]:
                    if len(n[1])>1:
                        parameter_table['Porod Volume'].append(str(round(float(n[1]),2))+' nm\u00b3')
                    else:
                        parameter_table['Porod Volume'].append('N/A')
                if len(n)<3 and len(n)>0 and '_sas_result.estimated_volume' in n[0] and '_sas_result.estimated_volume_error' not in n[0]:
                    if len(n[1])>1:
                        parameter_table['Estimated Volume'].append(str(round(float(n[1]),2))+' nm\u00b3')
                    else:
                        parameter_table['Estimated Volume'].append('N/A')

        return parameter_table


    def get_parameters_mw_many(self):
        '''
        get MW details from SASCIF files
        '''
        self.get_sascif_file()
        parameter_table={'SASDB ID':[],'Chemical composition MW':[],'Standard MW':[],'Porod Volume/MW':[]}   
        for code in self.clean_SASBDB_code():
            parameter_table['SASDB ID'].append(code)
            all_lines=self.get_all_sascif(code)
            for m,n in enumerate(all_lines):
                if len(n)<3 and len(n)>0 and '_sas_result.experimental_MW' in n[0] and '_sas_result.experimental_MW_error' not in n[0]:
                    if len(n[1])>1:
                        parameter_table['Chemical composition MW'].append(str(round(float(n[1]),2))+' kDa')
                    else:
                        parameter_table['Chemical composition MW'].append('N/A')
                if len(n)<3 and len(n)>0 and '_sas_result.MW_standard' in n[0] and '_sas_result.MW_standard_error' not in n[0]:
                    if len(n[1])>1:
                        parameter_table['Standard MW'].append(str(round(float(n[1]),2))+' kDa')
                    else:
                        parameter_table['Standard MW'].append('N/A')
                if len(n)<3 and len(n)>0 and '_sas_result.MW_Porod' in n[0] and '_sas_result.MW_Porod_error' not in n[0]:
                    if len(n[1])>1:
                        Porod_MW=round(float(n[1]),2)
                    else:
                        Porod_MW=0
                if len(n)<3 and len(n)>0 and '_sas_result.Porod_volume' in n[0] and '_sas_result.Porod_volume_error' not in n[0]:
                    if len(n[1])>1 and Porod_MW>0:
                        Porod_V=round(float(n[1]),2)/Porod_MW
                        parameter_table['Porod Volume/MW'].append(str(round(Porod_V,2))+' nm \u00b3/kDa')
                    else:
                        parameter_table['Porod Volume/MW'].append('N/A')

        return parameter_table

    def get_parameters_mw_many_dep(self):
        '''
        depreciated function on getting MW from JSON
        '''
        data_dic=self.get_data_from_SASBDB()
        #MW based on chemical composition
        parameter_table={'SASDB ID':[],'Sequence MW':[],'Experimental MW':[],'Porod MW':[]}        
        for key,val in data_dic.items():
            try:
                parameter_table['Experimental MW'].append(list(data_dic.values())[0]['experimental_mw']+' kDa')
            except:
                parameter_tabel['Experimental MW'].append('N/A')
            try:
                parameter_table['Porod MW'].append(list(data_dic.values())[0]['porod_mw']+' kDa')
            except:
                parameter_table['Porod MW'].append('N/A')
            try:
                parameter_table['Sequence MW'].append(list(data_dic.values())[0]['experiment']['sample']['molecule'][0]['total_mw']+' kDa')
            except:
                parameter_table['Sequence MW'].append('N/A')
            parameter_table['SASDB ID'].append(key)
        return parameter_table
    
    def get_pddf(self):
        '''
        get p(r) data from JSON
        '''
        data_dic=self.get_data_from_SASBDB()
        pofr_dic=self.get_pofr()
        pddf_dic={}
        for key,val in pofr_dic.items():
            pd_df=val.astype({'P':float,'R':float,'E':float})
            pd_df['R']=pd_df['R']/10;
            pd_df['err_x']=pd_df.apply(lambda row: (row['R'],row['R']), axis=1)
            pd_df['err_y']=pd_df.apply(lambda row: (row['P']-row['E'],row['P']+row['E']),axis=1)
            pddf_dic[key]=pd_df
            print (pd_df.head())
        return pddf_dic

    def get_pddf_info(self):
        '''
        get p(r) related info from JSON
        '''
        data_dic=self.get_data_from_SASBDB()
        pddf_info={'SASDB ID':[],'Software used':[],'Dmax':[],'Dmax error':[],'Rg':[],'Rg error':[]}
        for key,val in data_dic.items():
            pddf_info['Software used'].append(str(val['pddf_software']))
            try:
                pddf_info['Dmax'].append(str(val['pddf_dmax'])+' nm')
            except: 
                pddf_info['Dmax'].append('N/A')

            try:
                pddf_info['Rg'].append(str(val['pddf_rg'])+' nm')
            except:
                pddf_info['Rg'].append('N/A')
            try:
                if val['pddf_dmax_error'] is None:
                    pddf_info['Dmax error'].append('N/A')
                else:
                    pddf_info['Dmax error'].append(str(val['pddf_dmax_error'])+' nm')

            except:
                pddf_info['Dmax error'].append('N/A')
            try:
                pddf_info['Rg error'].append(str(val['pddf_rg_error'])+' nm')
            except:
                pddf_info['Rg error'].append('N/A')

            pddf_info['SASDB ID'].append(key)
        return pddf_info

    def get_number_of_fits(self):
        '''
        get number of fits from JSON, deprecated
        '''
        data_dic=self.get_data_from_SASBDB()
        num_of_fits={}
        for key,val in data_dic.items():
            num_of_fits[key]=len(val['fits'])
        return num_of_fits

    def get_chi_table(self):
        '''
        get chi value from JSON, deprecated
        '''
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
                chi_table['Model'].append('N/A')
                chi_table['\u03C7'+'\u00b2'].append('N/A')
        return chi_table

    def get_sasdb_code_fits(self):
        '''
        get number of fits per SASBDB ID
        '''
        fit_dict=self.get_number_of_fits()
        print (list(fit_dict.values()))
        return list(fit_dict.values())

    def get_fit_data(self):
        '''
        get fit information to make plots, from JSON
        '''
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
                    with open (fname,'w') as f:
                        f.write(fit.text)
                    f_df=pd.read_csv(fname, skiprows=3,delim_whitespace=True, names=['Q','Ie','Ib','E'])
                    if abs(f_df.iloc[22,2]-f_df.iloc[22,1])>abs(f_df.iloc[22,3]-f_df.iloc[22,1]):
                        f_df.rename(columns={'Q':'Q','Ie':'Ie','Ib':'E','E':'Ib'},inplace=True)
                    f_df['logIe']=np.log(f_df['Ie'])
                    f_df['logIb']=np.log(f_df['Ib'])
                    f_df['r']=f_df['Ie']-f_df['Ib']
                    if f_df['E'].isnull().values.any():
                        f_df['rsigma']=0
                    else:
                        f_df['rsigma']=f_df['r']/f_df['E']
                    print (key)
                    print (f_df['rsigma'])
                    #f_df['rsigma']=f_df['r']/f_df['E']
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
        '''
        '''
        r2=df['r2a'].sum()/df['r2b'].sum()
        return round(r2,2)
    
    def get_total_fits(self):
        '''
        get number of fits
        '''
        data_dic=self.get_data_from_SASBDB()
        num=0
        for key,val in data_dic.items():
            num += len(val['fits'])
        return num

    def get_fit_image(self):
        '''
        get fit image from fit, deprecated
        '''
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
                    filename = os.path.abspath(os.path.join(os.getcwd(), 'static/images/',self.ID+key+str(i)+'fit.png'))
                    with open (filename,'wb') as f:
                        f.write(fit.content)

