###################################
# Script : 
# 1) Contains Base class to 
# generate specific information 
# from mmcif files
# 2) Contains infromation for 
#  entry composition section
#
# ganesans - Salilab - UCSF
# ganesans@salilab.org
###################################

import ihm
import ihm.reader
import pandas as pd
import glob
import os
import numpy as np
import re
from collections import defaultdict
#########################
#Useful functions
#########################
def dict_to_list(d):
    L = []
    if bool(d):
        for k,v in d.items():
            if isinstance(v, list):
                L.append([k,v])
            else:
                L.append([k,[v]])
    return L

def dict_to_JSlist(d):
    L = []
    if bool(d):
        if len(list(d.keys()))>0:
            L.append(list(d.keys()))
            target=list(d.values())
            for i in range(len(target[0])):
                ltt=[]
                for j in target:
                    ltt.append(str(j[i]))
                L.append(ltt)
    return L


def format_RB_text(tex):
    val=''
    for a in tex:
        for b in a:
            if b==a[-1] and a==tex[-1]:
                val+=str(b)+'. '
            elif b==a[-1] and a!= tex[-1]:
                val+=str(b)+', '
            else:
                val+=str(b)+':'
    if val=='':
        val='-'
    return val

def format_flex_text(tex):
    val=''
    for a in tex:
        for b in a:
            if b==a[-1] and a==tex[-1]:
                val+=str(b)+'. '
            else:
                val+=str(b)+', '

    if val=='':
        val='-'

    return val


def format_tupple(tex):
    return str(tex[0])+'-'+str(tex[1])


def dict_to_JSlist_rows(d1,d2):
    L=[]
    L.append(['Chain ID','Rigid bodies','Non-rigid segments'])
    for i,j in d1.items():
        L.append([i,format_RB_text(j),format_flex_text(d2[i])])
    return L

def islistempty(inlist):
    if isinstance (inlist,list):
        return all(map(islistempty,inlist))
    return False

def cat_list_string(listn):
    result=' '
    for i in range(len(listn)):
        if i==0:
            result += str(listn[i])
        else:
            result += ','
            result += str(listn[i])
    return result

def get_key_from_val(dict1,val1):
    return dict1.keys()[dict1.values().index(val1)]

def get_val_from_key(dict1,key1):
    return dict1[key1]

def get_name(name):
    #if str(name) in ['?','',1,'.']:
    #    return 'None Listed'
    #else:
    return str(name)

def get_copy(name):
    if str(name)=='?':
        copy='None listed'
    elif '.' in name:
        copy=(name.split('.')[1]).split('@')[0]
    else:
        copy=1
    return copy

#########################
#Get information from IHM reader
#########################

class get_input_information(object):
    def __init__(self,mmcif_file):
        self.mmcif_file = mmcif_file
        self.datasets = {}
        self.entities = {}
        self.model = ihm.model.Model
        self.system, = ihm.reader.read(open(self.mmcif_file),
                                  model_class=self.model)
        
    def get_databases(self):
        dbs=self.system.orphan_datasets
        return dbs

    def get_id(self):
        if self.system.id is 'model':
            id=self.get_id_from_entry()
        else:
            id=self.system.id.split('_')[0]+self.system.id.split('_')[1]
        return id

    def get_id_from_entry(self):
        sf=open(self.mmcif_file,'r')
        for i,ln in enumerate(sf.readlines()):
            line =ln.strip().split(' ')
            if '_entry' in line[0]:
                entry_init=line[-1]
                entry=entry_init.split('_')[0]+entry_init.split('_')[1]
        return entry

    def get_title(self):
        cit=self.system.citations
        tit=cit[0].title
        return tit

    def get_authors(self):
        cit=self.system.citations
        aut=cit[0].authors
        for i in range(0,len(aut)):
            if i==0:
                authors=str(aut[i])
            else:
                authors+=';'+str(aut[i])
        return authors

    def get_struc_title(self):
        strc=self.system.title
        if strc is None:
            e=self.system.entities
            mol_name=e.description
        else:
            mol_name=strc
        return mol_name
    
    def check_sphere(self):
        spheres=[len(b._spheres) for i in self.system.state_groups for j in i for a in j for b in a]
        if 0 not in spheres:
            return 1
        else:
            return 0 

    def get_assembly_ID_of_models(self):
        assembly_id=[b.assembly._id for i in self.system.state_groups for j in i for a in j for b in a]
        return assembly_id

    def get_representation_ID_of_models(self):
        representation_id=[b.representation._id for i in self.system.state_groups for j in i for a in j for b in a]
        return representation_id

    def get_model_names(self):
        model_name1=[a.name for i in self.system.state_groups for j in i for a in j]
        model_name2=[b.name for i in self.system.state_groups for j in i for a in j for b in a]
        if len(model_name1)==len(model_name2):
            model_name =[str(t[0])+'/'+str(t[1]) for t in zip(model_name1,model_name2)]
        else:
            model_name=model_name2
        return model_name
    
    def get_model_assem_dict(self):
        model_id=map(int,[b._id for i in self.system.state_groups for j in i for a in j for b in a])
        assembly_id=map(int,self.get_assembly_ID_of_models())
        model_assembly=dict(zip(model_id,assembly_id))
        return model_assembly

    def get_model_rep_dict(self):
        model_id=map(int,[b._id for i in self.system.state_groups for j in i for a in j for b in a])
        rep_id=map(int,self.get_representation_ID_of_models())
        model_rep=dict(zip(model_id,rep_id))
        return model_rep

    def get_number_of_models(self):
        models=[b._id for i in self.system.state_groups for j in i for a in j for b in a]
        return len(models)

    def get_residues(self,asym):
        if asym.seq_id_range[0] is not None:
            residues=asym.seq_id_range[1]-asym.seq_id_range[0]+1
        elif asym.seq_id_range[0] is None:
            residues='None listed'
        return residues
    
    def get_composition(self):
        entry_comp={'Model ID':[],'Subunit number':[],'Subunit ID':[],
                    'Subunit name':[],'Chain ID':[],
                    'Total residues':[]}
        for i,j in self.get_model_assem_dict().items():
            for m in self.system.orphan_assemblies:
                if int(m._id)==int(j):
                    count=0
                    for n in m:
                        try:
                            count += 1
                            entry_comp['Model ID'].append(i)
                            entry_comp['Subunit number'].append(count)
                            entry_comp['Subunit ID'].append(n.entity._id)
                            entry_comp['Subunit name'].append(str(n.entity.description))
                            entry_comp['Chain ID'].append(n._id)
                            entry_comp['Total residues'].append(self.get_residues(n))
                        except AttributeError:
                            break
        return entry_comp

    def get_protocol_number(self):
        return len(self.system.orphan_protocols)

    def get_sampling(self):
        sampling_comp={'Step number':[], 'Protocol ID':[],'Method name':[],'Method type':[],'Number of computed models':[],'Multi state modeling':[],
                        'Multi scale modeling':[]}
        for prot in self.system.orphan_protocols:
            #print (prot,prot._id)
            for step in prot.steps:
                #print (step._id)
                sampling_comp['Step number'].append(step._id)
                sampling_comp['Multi state modeling'].append(str(step.multi_state))
                sampling_comp['Multi scale modeling'].append(str(step.multi_scale))                
                sampling_comp['Protocol ID'].append(self.system.orphan_protocols.index(prot)+1)
                sampling_comp['Method name'].append(step.method)
                sampling_comp['Method type'].append(step.name)
                sampling_comp['Number of computed models'].append(step.num_models_end)
        return sampling_comp

    def get_representation(self):
        representation_comp={'Chain ID':[],'Subunit name':[],'Rigid bodies':[],
                    'Non-rigid regions':[]}
        for i in self.system.orphan_representations:
            print (["%s:%d-%d" % ((x.asym_unit._id,) + x.asym_unit.seq_id_range) for x in i if not x.rigid])

    def get_RB_flex_dict(self):
        RB=self.get_empty_chain_dict();RB_nos=[];all_nos=[];flex=self.get_empty_chain_dict()
        for i in self.system.orphan_representations:
            for j in i:
                all_nos.append(j.asym_unit.seq_id_range)
                if j.rigid and j.starting_model:
                    RB_nos.append(j.asym_unit.seq_id_range)
                    RB[j.starting_model.asym_unit._id].append([format_tupple(j.asym_unit.seq_id_range),
                        get_val_from_key(self.get_dataset_dict(),j.starting_model.dataset._id)])
                elif j.rigid and not j.starting_model:
                    RB_nos.append(j.asym_unit.seq_id_range)
                    RB[j.asym_unit._id].append([format_tupple(j.asym_unit.seq_id_range),'None'])
                else:
                    flex[j.asym_unit._id].append([format_tupple(j.asym_unit.seq_id_range)])
        return RB,flex,len(RB_nos),len(all_nos)      


    def get_RB_dict(self):
        RB=self.get_empty_chain_dict();RB_nos=[];all_nos=[]
        for i in self.system.orphan_representations:
            for j in i:
                all_nos.append(j.asym_unit.seq_id_range)
                if j.rigid: #starting_model:
                    RB_nos.append(j.asym_unit.seq_id_range)
                    RB[j.starting_model.asym_unit._id].append([get_val_from_key(self.get_chain_subunit_dict(),j.starting_model.asym_unit._id),j.asym_unit.seq_id_range,
                        get_val_from_key(self.get_dataset_dict(),j.starting_model.dataset._id)])
        return RB,len(RB_nos),len(all_nos)      

    def get_number_of_assemblies(self):
        return (len(self.system.orphan_assemblies))

    def get_number_of_entities (self):
        return (len(self.system.entities))

    def get_number_of_chains(self):
        used=[];assembly={}
        lists= self.system.orphan_assemblies
        for i,k in enumerate(self.system.orphan_assemblies):
            chain=[]
            for l in k:
                chain.append(l._id)
            assembly[i]=chain
            unique=[used.append(x) for x in chain if x not in used]
        number_of_chains=[len(i) for i in assembly.values()]
        return number_of_chains

    def get_all_asym(self):
        parents=[(a._id,a.details,a.entity.description,a.entity._id,i) for i,a in enumerate(self.system.asym_units)]
        return parents

    def get_empty_chain_dict(self):
        empty_chain_dict={}
        for i,a in enumerate(self.system.asym_units):
            empty_chain_dict[a._id]=[]
        return empty_chain_dict

    def get_chain_subunit_dict(self):
        chain_subunit_dict={}
        for i,a in enumerate(self.system.asym_units):
            chain_subunit_dict[a._id]=a.details.split('.')[0]
        return chain_subunit_dict

    def get_residues_subunit_dict(self):
        residues_subunit_dict={}
        for i,a in enumerate(self.system.asym_units):
            residues_subunit_dict[a._id]=self.get_residues(a)
        return residues_subunit_dict

    def get_dataset_xl_info(self,id):
        restraints=self.get_restraints()
        print (restraints)
        print (id)
        lt='Linker name and number of cross-links: %s' % (restraints['Restraint info'][restraints['Dataset ID'].index(id)])
        return (lt)

    def get_dataset_dict(self):
        dataset_dict={}
        lists=self.system.orphan_datasets
        if len(lists)>0:
            for i in lists:
                try:
                    loc=i.location.db_name
                except AttributeError as error:
                    loc=str('Not listed')
                try:
                    acc=i.location.access_code
                except AttributeError as error:
                    acc=str('None')
                dataset_dict[i._id]=str(i.data_type)+'/'+str(acc)
        return dataset_dict 

    def get_software_length(self):
        lists=self.system.software
        if lists is None:
            return 0
        else:
            return len(lists)

    def get_software_comp (self):
        software_comp={'ID':[],'Software name':[],'Software version':[],'Software classification':[],'Software location':[]}
        lists=self.system.software
        if len(lists)>0:
            for i in lists:
                software_comp['ID'].append(i._id)
                software_comp['Software name'].append(i.name)
                software_comp['Software location'].append(i.location)
                if str(i.version) is '?':
                    vers='None'
                else:
                    vers=str(i.version)
                software_comp['Software version'].append(vers)
                software_comp['Software classification'].append(i.classification)
            final_software_composition=software_comp
        else:
            final_software_composition={}
        return final_software_composition

    def get_dataset_length(self):
        lists=self.system.orphan_datasets
        if lists is None:
            return 0
        else:
            return len(lists)

    def check_ensembles(self):
        a=self.system.ensembles
        return len(a)

    def get_ensembles(self):
        s=self.system.ensembles
        if len(s)>0:
            ensemble_comp={'Ensemble number':[],'Ensemble name':[],'Model ID':[],'Number of models':[],
                        'Clustering method': [], 'Clustering feature': [], 'Cluster precision':[]}
            for i in s:
                ensemble_comp['Ensemble number'].append(str(i._id))
                ensemble_comp['Ensemble name'].append(str(i.name))
                ensemble_comp['Model ID'].append(str(i.model_group._id))
                ensemble_comp['Number of models'].append(str(i.num_models))
                ensemble_comp['Clustering method'].append(str(i.clustering_method))
                ensemble_comp['Clustering feature'].append(str(i.clustering_feature))
                ensemble_comp['Cluster precision'].append(str(i.precision))
            return ensemble_comp
        else:
            return None

    def get_dataset_comp (self):
        dataset_comp={'ID':[],'Dataset type':[],'Database name':[],'Data access code':[]}
        lists=self.system.orphan_datasets
        if len(lists)>0:
            for i in lists:
                try:
                    loc=i.location.db_name
                except AttributeError as error:
                    loc=str('Not listed')
                try:
                    acc=i.location.access_code
                except AttributeError as error:
                    acc=str('None')
                dataset_comp['ID'].append(i._id)
                #if i.data_type=='unspecified' and 'None' not in i.details:
                #    dataset_comp['Dataset type'].append(i.details)
                #else:
                dataset_comp['Dataset type'].append(i.data_type)
                dataset_comp['Database name'].append(str(loc))
                dataset_comp['Data access code'].append(acc)
                #print (dataset_comp)
        return dataset_comp

    def dataset_id_type_dic(self):
        dataset_dic={}
        if len(self.system.orphan_datasets)>0:
            for i in self.system.orphan_datasets:
                if i.data_type =='unspecified':
                    dataset_dic[str(i._id)]=str(i.details)
                else:
                    dataset_dic[str(i._id)]=str(i.data_type)
        return dataset_dic
        
    def get_restraints(self):
        r=self.system.restraints
        restraints_comp={'ID':[],'Dataset ID':[],'Restraint type':[],'Restraint info':[]}
        for j,i in enumerate(r):
            restraints_comp['ID'].append(j+1)
            restraints_comp['Restraint type'].append(str(i.__class__.__name__))
            restraints_comp['Dataset ID'].append(str(i.dataset._id))
            if 'CrossLink' in str(i.__class__.__name__):
                restraints_comp['Restraint info'].append(str(i.linker.auth_name)+', '+str(len(i.experimental_cross_links))+' cross-links')
            if 'EM3D' in str(i.__class__.__name__):
                restraints_comp['Restraint info'].append(str(i.fitting_method)+ ', '+str(i.number_of_gaussians))
            if 'EM2D' in str(i.__class__.__name__):
                restraints_comp['Restraint info'].append('Number of micrographs: '+str(i.number_raw_micrographs)+','+' Image resolution: '+str(i.image_resolution))
            if 'SAS' in str(i.__class__.__name__):
                restraints_comp['Restraint info'].append('Assembly name: '+str(i.assembly.name)+' Fitting method: '+str(i.fitting_method)+ ' Multi-state: '+str(i.multi_state))
            if 'UpperBound' in str(i.__class__.__name__):

                restraints_comp['Restraint info'].append('Distance: '+str(i.distance))
            if 'Mutagenesis' in str(i.__class__.__name__):
                restraints_comp['Restraint info'].append('Details: '+str(i.details))
            if 'DerivedDistance' in str(i.__class__.__name__):
                dic=self.dataset_id_type_dic()
                ID=str(i.dataset._id)
                #restraints_comp['Restraint info'].append(dic[ID])
                if 'UpperBound' in str(i.distance.__class__.__name__):
                    restraints_comp['Restraint info'].append(('Upper Bound Distance: '+str(i.distance.distance)))
                else:
                    restraints_comp['Restraint info'].append(['restraint type '+ str(i.distance.__class__.__name__),dic[ID]])
        return (restraints_comp) 
    
    def get_dataset_details(self):
        dataset_comp={'ID':[],'Dataset type':[],'Database name':[],'Details':[]}
        lists=self.system.orphan_datasets
        if len(lists)>0:
            for i in lists:
                try:
                    loc=i.location.db_name
                except AttributeError as error:
                    loc=str('')
                try:
                    acc=i.location.access_code
                except AttributeError as error:
                    acc=str('Not listed')
                dataset_comp['ID'].append(i._id)
                if i.data_type=='unspecified' and i.details is not None:
                    dataset_comp['Dataset type'].append(i.details)
                else:
                    dataset_comp['Dataset type'].append(i.data_type)
                dataset_comp['Database name'].append(str(loc))
                if 'ComparativeModel' in str(i.__class__.__name__):
                    acc1='template PDB ID: '+ acc
                    dataset_comp['Details'].append(acc1)
                elif 'PDB' in str(i.__class__.__name__):
                    acc1='PDB ID: '+ acc
                    dataset_comp['Details'].append(acc1)
                elif 'CX' in str(i.__class__.__name__):
                    acc1=self.get_dataset_xl_info(i._id)
                    dataset_comp['Details'].append(acc1)
                elif 'EM' in str(i.__class__.__name__):
                    acc1='EMDB ID: '+acc
                    dataset_comp['Details'].append(acc1)
                else:
                    dataset_comp['Details'].append(acc)
                    
        return dataset_comp
    
    def get_atomic_coverage(self):
        for i in self.system.orphan_representations:
            if self.check_sphere()==1:
                flex=sum([int(x.asym_unit.seq_id_range[1])-int(x.asym_unit.seq_id_range[0])+1 for x in i if not x.rigid])
                rigid=sum([int(x.asym_unit.seq_id_range[1])-int(x.asym_unit.seq_id_range[0])+1 for x in i if x.rigid])
                percentage=str(round(rigid/(rigid+flex)*100))+'%'
            else:
                percentage='100%'
        return percentage

    def check_for_sas(self):
        dataset=self.get_dataset_comp()
        data_type=dataset['Dataset type']
        print (data_type)
        if 'SAS' in str(data_type):
            return True
        else:
            return False

    def check_for_cx(self,dataset):
        dataset=self.get_dataset_comp()
        data_type=dataset['Dataset type']
        print (data_type)
        if 'CX' in str(data_type):
            return True
        else:
            return False

    def check_for_em(self,dataset):
        dataset=self.get_dataset_comp()
        data_type=dataset['Dataset type']
        print (data_type)
        if 'EM' in str(data_type):
            return True
        else:
            return False

    def get_resolution(self):
        pass


