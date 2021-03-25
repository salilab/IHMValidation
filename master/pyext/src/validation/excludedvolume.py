###################################
# Script : 
# 1) Contains class for excluded 
# volume calculation
#
# ganesans - Salilab - UCSF
# ganesans@salilab.org
###################################
from validation import get_input_information
import ihm
import multiprocessing as mp
import pandas as pd
import numpy as np
import math,os
from multiprocessing import Process, Queue, Pool, Manager


class get_excluded_volume(get_input_information):
    def __init__(self,mmcif_file):
        super().__init__(mmcif_file)
        self.ID=str(get_input_information.get_id(self))
        self.nos=get_input_information.get_number_of_models(self)
        
    def get_all_spheres (self,filetemp=None):
        """get information on all spheres for each model"""
        if filetemp is None:
            Model_object=[b for i in self.system.state_groups for j in i for a in j for b in a]
            model_dict={i+1:j._spheres for i,j in enumerate(Model_object)}
        else:
            system, = ihm.reader.read(filetemp,
                                  model_class=ihm.model.Model)
            Model_object=[b for i in system.state_groups for j in i for a in j for b in a]
            model_dict={i+1:j._spheres for i,j in enumerate(Model_object)}
            #print (model_dict)
        return model_dict
    
    def get_nCr(self,n,r):
        """get all combinations"""
        f=math.factorial
        return (f(n)/(f(r)*f(n-r)))

    def get_violation_percentage(self,models_spheres_df:pd.DataFrame,viols:dict)->float:
        """get information on all spheres for each model"""
        number_of_violations=sum(list(viols.values()))
        number_of_combinations=self.get_nCr(models_spheres_df.shape[1],2)
        return (1-number_of_violations/number_of_combinations)*100

    def get_violation_normalized(self,models_spheres_df:pd.DataFrame,viols:dict)->float:
        """ """
        number_of_violations=sum(list(viols.values()))
        normalization_constant=models_spheres_df.shape[1]*math.log(models_spheres_df.shape[1],10)
        return (1-number_of_violations/normalization_constant)*100

    def get_xyzr(self,spheres:pd.DataFrame)->pd.DataFrame:
        """ get X,Y, Z coords from sphere objects"""
        model_spheres={i+1:[j.x,j.y,j.z,j.radius] for i,j in enumerate(spheres)}
        model_spheres_df=pd.DataFrame(model_spheres, index=['X','Y','Z','R'])
        return model_spheres_df

    def get_xyzr_complete(self,model_ID,spheres:list)->pd.DataFrame:
        """ get X,Y,Z,R, chain and model ID from sphere objects"""

        model_spheres={i+1:[j.x,j.y,j.z,j.radius,j.asym_unit._id,model_ID] for i,j in enumerate(spheres)}
        model_spheres_df=pd.DataFrame(model_spheres, index=['X','Y','Z','R','Chain_ID','Model_ID'])
        return model_spheres_df

    def get_violation_dict(self,model_spheres_df:pd.DataFrame)->dict:
        """ get violation from model_sphere df"""
        viols={}
        for indx,col in model_spheres_df.iteritems():
            if indx< model_spheres_df.shape[1]:
                sphere_R=model_spheres_df.iloc[-1,indx:]
                remaining=model_spheres_df.iloc[:-1,indx:]
                subt_alone=remaining.sub(col[:-1], axis=0)
                final_df=np.square(subt_alone)
                final_df.loc['sqrt']=np.sqrt(final_df.sum(axis=0))
                final_df.loc['R_tot']=sphere_R.add(col[[-1]].tolist()[0]).to_list()
                final_df.loc['distances']=final_df.loc['sqrt']-final_df.loc['R_tot']
                final_df.loc['violations']=final_df.loc['distances'].apply(lambda x: 1 if x <0 else 0)
                viols[indx]=final_df.loc['violations'].sum(axis=0)
        return viols 
    
    def get_exc_vol_for_models(self,model_dict:dict)->dict:
        excluded_volume={'Models':[],'Excluded Volume Satisfaction':[], 'Number of violations':[]}
        for indx, model in model_dict.items():
            excluded_volume['Models'].append(indx)
            df=self.get_xyzr(model)
            excluded_volume['Excluded Volume Satisfaction'].append(round(self.get_violation_percentage(df,self.get_violation_dict(df)),2))
            excluded_volume['Number of violations'].append(sum(list(self.get_violation_dict(df).values())))
        f_exv=open(os.path.join(os.getcwd(),'static/results/',self.ID+'exv.txt'),'w+')
        return excluded_volume

    def get_exc_vol_for_models_normalized(self,model_dict:dict)->dict:
        excluded_volume={'Models':[],'Excluded Volume Satisfaction':[]}
        for indx, model in model_dict.items():
            excluded_volume['Models'].append(indx)
            df=self.get_xyzr(model)
            satisfaction=self.get_violation_percentage(df,self.get_violation_dict(df))
            excluded_volume['Excluded Volume Satisfaction'].append(round(satisfaction,2))
        return excluded_volume

    def get_exc_vol_given_sphere_parallel(self,sphere_list:list)->(float,int):
        """ """
        df=self.get_xyzr(sphere_list)
        violation_dict=self.get_violation_dict(df)
        satisfaction=round(self.get_violation_percentage(df,violation_dict),2)
        violations=sum(list(violation_dict.values()))
        return (satisfaction,violations)

    def run_exc_vol_parallel(self,model_dict:dict)->dict:
        """ get exc vol info in parallel """
        list_of_sphere_list=list(model_dict.values())
        pool=mp.Pool(processes=len(list(model_dict.keys())))
        complete_list=pool.map(self.get_exc_vol_given_sphere_parallel,list(model_dict.values()))
        excluded_volume={'Models':list(model_dict.keys()),'Excluded Volume Satisfaction':[i[0] for i in complete_list], 'Number of violations':[i[1] for i in complete_list]}
        f_exv=open(os.path.join(os.getcwd(),'../static/results/',self.ID+'exv.txt'),'w+')
        return excluded_volume

    def exv_readable_format(self,exv:dict)->str:
        fin_string=''
        for i in exv['Models']:
            fin_string+'Model-'+str(i)+': '+'Number of violations-' + str(exv['Number of violations'])
        return fin_string


