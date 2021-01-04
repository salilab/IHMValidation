###################################
# Script : 
# 1) Contains class to 
# generate molprobity assessments
#
# ganesans - Salilab - UCSF
# ganesans@salilab.org
###################################

import pickle,os
from subprocess import run, call, PIPE
from validation import get_input_information,utility
import ihm
import ihm.reader


class get_molprobity_information(get_input_information):
    def __init__(self,mmcif_file):
        super().__init__(mmcif_file)
        self.ID=str(get_input_information.get_id(self))
        self.nos=get_input_information.get_number_of_models(self)

    def check_for_molprobity(self,filetemp=None):
        """ Check the biso and occupancy columns for mmcif files"""
        if filetemp is not None:
            model=ihm.model.Model
            system, = ihm.reader.read(filetemp,
                                  model_class=model)
            models=[b for i in system.state_groups for j in i for a in j for b in a]
        else:
            """check if file is in the right format for molprobity analysis """
            models=[b for i in self.system.state_groups for j in i for a in j for b in a]
        if models[0]._atoms[0].biso is not None: #and  models[0]._atoms[0].occupancy is not None:
            print ("File in the appropriate format for molprobity")
            return True
            # uncomment the lines below to run MP analysis
            #call([r"/home/ganesans/PDB-dev-test/MolProbity-master/build/bin/molprobity.molprobity",self.mmcif_file,"prefix="+out])
            #print ("printing result", result_mp.stdout)
        else:
            #print ("bisoval",models[0]._atoms[0].biso)
            #print ("occ",models[0]._atoms[0].occupancy)
            print ("File is not in the appropriate format for molprobity")  
            return False

    def run_ramalyze(self,d):
        """run ramalyze to get outliers """
        f=str(self.ID)+'_temp_rama.txt'
        f1=open(f,'w+')
        with f1 as outfile:
            run([r"/home/ganesans/PDB-dev-test/MolProbity-master/build/bin/molprobity.ramalyze",self.mmcif_file],stdout=outfile)
        with open(f,'r+') as inf:
            line=[ln.strip() for ln in inf.readlines()]
        d['rama']=line
        filename = open(os.path.join(os.getcwd(), 'static/results/',self.ID+'_temp_rama.txt'))
        with open(filename,'wb') as fp:
            pickle.dump(d['rama'],fp)

    def run_molprobity(self,d):
        """run molprobity"""
        f=str(self.ID)+'_temp_mp.txt'
        f1=open(f,'w+')
        with f1 as outfile:
            run([r"/home/ganesans/PDB-dev-test/MolProbity-master/build/bin/molprobity.molprobity",self.mmcif_file],stdout=outfile)
        with open(f,'r+') as inf:
            line=[ln.strip() for ln in inf.readlines()]
        d['molprobity']=line
        filename = open(os.path.join(os.getcwd(), 'static/results/',self.ID+'_temp_mp.txt'))
        with open(filename,'wb') as fp:
            pickle.dump(d['molprobity'],fp)

    def run_clashscore(self,d):
        """run clashscore to get information on steric clashes"""
        f=str(self.ID)+'_temp_clash.txt'
        f1=open(f,'w+')
        with f1 as outfile:
            run([r"/home/ganesans/PDB-dev-test/MolProbity-master/build/bin/molprobity.clashscore",self.mmcif_file],stdout=outfile)
        with open(f,'r+') as inf:
            line=[ln.strip() for ln in inf.readlines()]
        d['clash']=line
        filename = open(os.path.join(os.getcwd(), 'static/results/',self.ID+'_temp_clash.txt'))
        with open(filename,'wb') as fp:
            pickle.dump(d['clash'],fp)
        
    def run_rotalyze(self,d):
        """run rotalyZe to get rotameric outliers"""
        f=str(self.ID)+'_temp_rota.txt'
        f1=open(f,'w+')
        with f1 as outfile:
            run([r"/home/ganesans/PDB-dev-test/MolProbity-master/build/bin/molprobity.rotalyze",self.mmcif_file],stdout=outfile)
        with open(f,'r+') as inf:
            line=[ln.strip() for ln in inf.readlines()]
        d['rota']=line
        filename = open(os.path.join(os.getcwd(), 'Output/static/results/',self.ID+'_temp_rota.txt'))
        with open(filename,'wb') as fp:
            pickle.dump(d['rota'],fp)

    def write_all_lines(self,file_handle):
        with open(file_handle.name, 'r+') as inf:
            line=[ln.strip() for ln in inf.readlines()] 
        return line

    def process_rama(self,line):
        """ reading and processing molprobity output from rama outliers.
        Outputs information specific to models """
        line_new=line[1:-3];count=1
        models={i:[] for i in range(1,self.nos+1)}
        cutoff=len(line_new)/self.nos
        for i,k in enumerate(line_new):
            if self.nos>1 and i <len(line_new)-1:
                if i<count*cutoff:
                #if (abs(int(line_new[i+1].strip().split()[1])-int(line_new[i].strip().split()[1]))<10 or (line_new[i+1].strip().split()[0] != line_new[i].strip().split()[0])) and i<count*cutoff :
                    #print (i,count,cutoff, int(line_new[i+1].strip().split()[1]),int(line_new[i].strip().split()[1]),abs(int(line_new[i+1].strip().split()[1])-int(line_new[i].strip().split()[1])))
                    models[count].append(k)
                    #print (line_new[i+1].strip().split()[0])
                else:
                    count=count+1
                    #print (i,count,int(line_new[i+1].strip().split()[1]),int(line_new[i].strip().split()[1]),abs(int(line_new[i+1].strip().split()[1])-int(line_new[i].strip().split()[1])))
                    models[count].append(k)
            else:
                models[count].append(k) 
        return models

    def process_molprobity(self,line):
        """ process molprobity files to extract relevant information """
        bond_index=angle_index=None
        for i,j in enumerate(line):
            if 'Bond outliers' in j:
                bond_index=i
            if 'Angle outliers' in j:
                angle_index=i
            if 'Molprobity validation' in j:
                end=i
        if angle_index is None and bond_index is not None:
            bond_outliers=line[bond_index+2:end-1]
            return (bond_outliers,[])
        elif bond_index is None and angle_index is not None:
            angle_outliers=line[angle_index+2:end-1]
            return ([],angle_outliers)
        elif bond_index is None and angle_index is None:
            return ([],[])
        else:
            bond_outliers=line[bond_index+2:angle_index-1]
            angle_outliers=line[angle_index+2:end-1]
            return (bond_outliers,angle_outliers)


    def process_clash(self,line):
        """ process clash files to extract relevant information """
        count=[i for i,j in enumerate(line) if 'Bad Clashes' in j]
        if self.nos>1:
            vals=[j.split(' ')[5] for i,j in enumerate(line) if 'Bad Clashes' in j]
            clashes={j[-1]:[] for j in vals}
        else:
            clashes={'Model 1':[]}
        count.append(len(line)-self.nos)
        for i in range(0,len(clashes.keys())):
            l=[j for k,j in enumerate(line) if k>int(count[i]) and k<int(count[i+1])]
            clashes[list(clashes.keys())[i]].append(l)
        clashes_ordered=dict(sorted(clashes.items()))
        return clashes_ordered

    def process_rota(self,line):
        """ process rota files to extract relevant information """
        line_new=line[1:-1];count=1
        models={i:[] for i in range(1,self.nos+1)}
        cutoff=len(line_new)/self.nos
        for i,k in enumerate(line_new):
            if self.nos>1 and i <len(line_new)-1:
                if i<count*cutoff:
                    models[count].append(k)
                else:
                    count=count+1
                    models[count].append(k)
            else:
                models[count].append(k)
        return models

    def rama_summary_table(self,models):
        """ write out summary table from rama, clash and other tables"""
        f_rama=open(os.path.join(os.getcwd(),'static/results/',self.ID+'_rama_summary.txt'),'w+')
        dict1={'Model ID':[],'Analyzed':[],'Favored':[],'Allowed':[],'Outliers':[]}
        for i,j in models.items():
            dict1['Model ID'].append(i)
            F=[];A=[];U=[]
            for k in j:
                if k.strip().split()[-2]=='or':
                    a=':'.join(k.strip().split()[-3:])
                else:
                    a=k.strip().split()[-1]
                if a.split(':')[4] =='Favored':
                    F.append(a.split(':')[4])
                elif a.split(':')[4]=='Allowed':
                    A.append(a.split(':')[4])
                else:
                    U.append(a.split(':')[4])
            dict1['Analyzed'].append(len(F)+len(A)+len(U))
            dict1['Favored'].append(len(F))
            dict1['Allowed'].append(len(A)) 
            dict1['Outliers'].append(len(U))
        print (dict1['Model ID'], file=f_rama)
        print (dict1['Analyzed'],file=f_rama)
        print (dict1['Favored'], file=f_rama)
        print (dict1['Allowed'],file=f_rama)
        print (dict1['Outliers'],file=f_rama)
        return dict1

    def clash_summary_table(self,line):
        """ format clash data to print to file """
        f_clash=open(os.path.join(os.getcwd(),'static/results/',self.ID+'_clash_summary.txt'),'w+')
        clashes=self.process_clash(line)
        if self.nos>1:
            clashscore_list=line[len(line)-self.nos:]
        else:
            clashscore_list=['Model 1 ' + (line[len(line)-self.nos:])[0]]
        dict1={'Model ID':[],'Clash score':[],'Number of clashes':[]}
        for i in clashscore_list:
            dict1['Model ID'].append(str(i.split(' ')[0].title()+' '+i.split(' ')[1]))
            dict1['Clash score'].append(i.split(' ')[-1])
        for j in list(clashes.values()):
            dict1['Number of clashes'].append(len(j[0]))
        clash_total=(sum(dict1['Number of clashes']))
        print (dict1['Model ID'], file=f_clash)
        print (dict1['Clash score'],file=f_clash)
        print (dict1['Number of clashes'], file=f_clash)
        return dict1,clash_total

    def rota_summary_table(self,models):
        """ format rota data to print to file """
        dict1={'Model ID':[],'Analyzed':[],'Favored':[],'Allowed':[],'Outliers':[]}
        f_rota=open(os.path.join(os.getcwd(),'static/results/',self.ID+'_rota_summary.txt'),'w+')        
        for i,j in models.items():
            dict1['Model ID'].append(i)
            F=[];A=[];U=[]
            for k in j:
                if k.strip().split()[-1].split(':')[-2] =='Favored':
                    F.append(k.strip().split()[-1].split(':')[-2])
                elif k.strip().split()[-1].split(':')[-2]=='Allowed':
                    A.append(k.strip().split()[-1].split(':')[-2])
                else:
                    U.append(k.strip().split()[-1].split(':')[-2])
            dict1['Analyzed'].append(len(F)+len(A)+len(U))
            dict1['Favored'].append(len(F))
            dict1['Allowed'].append(len(A))
            dict1['Outliers'].append(len(U))
        print (dict1['Model ID'], file=f_rota)
        print (dict1['Analyzed'],file=f_rota)
        print (dict1['Favored'], file=f_rota)
        print (dict1['Allowed'],file=f_rota)
        print (dict1['Outliers'],file=f_rota)
        return dict1

    def rama_detailed_table(self,models):
        """ format rama information to print to file"""
        dict1={'Model ID':[],'Chain and res ID':[],'Residue type':[]}
        f_rama_D=open(os.path.join(os.getcwd(),'static/results/',self.ID+'_rama_detail.txt'),'w+')  
        for i,j in models.items():
            for k in j:
                if k.strip().split()[-2]=='or':
                    a=':'.join(k.strip().split()[-3:])
                else:
                    a=k.strip().split()[-1]
                if a.split(':')[4] =='OUTLIER':
                    dict1['Model ID'].append(i)
                    dict1['Residue type'].append(a.split(':')[0])
                    if len(k.strip().split()[0])>2:
                        dict1['Chain and res ID'].append(k.strip().split()[0])
                    else:
                        dict1['Chain and res ID'].append(':'.join(k.strip().split()[:2]))
        print (dict1['Model ID'], file=f_rama_D)
        print (dict1['Chain and res ID'],file=f_rama_D)
        print (dict1['Residue type'], file=f_rama_D)
        return dict1

    def molprobity_detailed_table_bonds(self,bonds):
        """process molprobity bonded information and format to table"""
        f_mp_D=open(os.path.join(os.getcwd(),'static/results/',self.ID+'_molprobity_bond.txt'),'w+')
        if len(bonds)>0:
            bondtype_diff={};bondtype_dist={};bondtype_ex={};
            dict1={'Bond type':[],'Observed distance (&#8491)':[],'Ideal distance (&#8491)':[],'Number of outliers':[]}
            for i,j in enumerate(bonds):
                if len(j.replace(',','').replace(':','').split()[0])>2:
                    bondtype_dist[j.replace(',','').replace(':','').split()[4]]=float(j.replace(',','').replace(':','').split()[6])
                    bondtype_diff[j.replace(',','').replace(':','').split()[4]]=float(j.replace(',','').replace(':','').split()[10])
                    if j.replace(',','').replace(':','').split()[4] not in list(bondtype_ex.keys()):
                        bondtype_ex[j.replace(',','').replace(':','').split()[4]]=[]
                    else:
                        bondtype_ex[j.replace(',','').replace(':','').split()[4]].append(float(j.replace(',','').replace(':','').split()[6]))
                else:
                    bondtype_dist[j.replace(',','').replace(':','').split()[5]]=float(j.replace(',','').replace(':','').split()[7])
                    bondtype_diff[j.replace(',','').replace(':','').split()[5]]=float(j.replace(',','').replace(':','').split()[11])
                    if j.replace(',','').replace(':','').split()[5] not in list(bondtype_ex.keys()):
                        bondtype_ex[j.replace(',','').replace(':','').split()[5]]=[]
                    else:
                        bondtype_ex[j.replace(',','').replace(':','').split()[5]].append(float(j.replace(',','').replace(':','').split()[7]))
            for m,n in bondtype_dist.items():
                dict1['Bond type'].append(m)
                dict1['Observed distance (&#8491)'].append(n)
                #ideal=sum(abs(bondtype_diff[i]))/len(bondtype_diff[i])
                dict1['Ideal distance (&#8491)'].append(round(n+bondtype_diff[m],2))
                dict1['Number of outliers'].append(len(bondtype_ex[m]))
            print (dict1['Bond type'], file=f_mp_D)
            print (dict1['Observed distance (&#8491)'],file=f_mp_D)
            print (dict1['Ideal distance (&#8491)'], file=f_mp_D)
            print (dict1['Number of outliers'], file=f_mp_D)
            return dict1
        else:
            return {}

    def molprobity_detailed_table_angles(self,angles):
        """process molprobity angles information and format to table"""
        f_mp_A=open(os.path.join(os.getcwd(),'static/results/',self.ID+'_molprobity_angles.txt'),'w+')
        if len(angles)>0:
            angle_diff={};angle_dist={};angle_ex={};
            dict1={'Angle type':[],'Observed angle (&#176)':[],'Ideal angle (&#176)':[],'Number of outliers':[]}
            for i,j in enumerate(angles):
                if len(j.replace(',','').replace(':','').split()[0])>2:
                    angle_dist[j.replace(',','').replace(':','').split()[4]]=float(j.replace(',','').replace(':','').split()[6])
                    angle_diff[j.replace(',','').replace(':','').split()[4]]=float(j.replace(',','').replace(':','').split()[10])
                    if j.replace(',','').replace(':','').split()[4] not in list(angle_ex.keys()):
                        angle_ex[j.replace(',','').replace(':','').split()[4]]=[]
                    else:
                        angle_ex[j.replace(',','').replace(':','').split()[4]].append(float(j.replace(',','').replace(':','').split()[6]))
                else:
                    angle_dist[j.replace(',','').replace(':','').split()[5]]=float(j.replace(',','').replace(':','').split()[7])
                    angle_diff[j.replace(',','').replace(':','').split()[5]]=float(j.replace(',','').replace(':','').split()[11])
                    if j.replace(',','').replace(':','').split()[5] not in list(angle_ex.keys()):
                        angle_ex[j.replace(',','').replace(':','').split()[5]]=[]
                    else:
                        angle_ex[j.replace(',','').replace(':','').split()[5]].append(float(j.replace(',','').replace(':','').split()[7]))
            for i,j in angle_dist.items():
                dict1['Angle type'].append(i)
                dict1['Observed angle (&#176)'].append(j)
                #ideal=sum(abs(bondtype_diff[i]))/len(bondtype_diff[i])
                dict1['Ideal angle (&#176)'].append(round(j+angle_diff[i],2))
                dict1['Number of outliers'].append(len(angle_ex[i]))
            print (dict1['Angle type'], file=f_mp_A)
            print (dict1['Observed angle (&#176)'],file=f_mp_A)
            print (dict1['Ideal angle (&#176)'], file=f_mp_A)
            print (dict1['Number of outliers'], file=f_mp_A)           
            return dict1
        else:
            return {}

    def clash_detailed_table(self,line):
        """process molprobity clash information and format to table"""
        f_clash_D=open(os.path.join(os.getcwd(),'static/results/',self.ID+'_clash_detailed.txt'),'w+')
        dict1={'Model ID':[],'Atom-1':[],'Atom-2':[],'Clash overlap (&#8491)':[]}
        clashes=self.process_clash(line)
        #clashes=self.cleanup_clashes(clashes,nos)
        for i,j in clashes.items():
            for k in j[0]:
                k1=[i for i in k.split(' ') if i not in '']
                dict1['Model ID'].append(i.title()[-1])
                if len(k1)<9 and len(k1[0])>2:
                    dict1['Atom-1'].append(':'.join(k1[0:3]))
                else:
                    dict1['Atom-1'].append(':'.join(k1[0:4]))
                if len(k1)<9 and len(k1[3])>4:
                    dict1['Atom-2'].append(':'.join(k1[3:-1]))
                else: 
                    dict1['Atom-2'].append(':'.join(k1[4:-1]))
                dict1['Clash overlap (&#8491)'].append(k1[-1].replace(':',''))

        print (dict1['Model ID'], file=f_clash_D)
        print (dict1['Atom-1'],file=f_clash_D)
        print (dict1['Atom-2'], file=f_clash_D)
        print (dict1['Clash overlap (&#8491)'], file=f_clash_D)
        return dict1 

    def rota_detailed_table(self,models):
        """process molprobity rotamers information and format to table"""
        f_rota_D=open(os.path.join(os.getcwd(),'static/results/',self.ID+'_rota_detailed.txt'),'w+')
        dict1={'Model ID':[],'Chain and res ID':[],'Residue type':[]}
        for i,j in models.items():
            for k in j:
                if k.strip().split()[-1].split(':')[-2] =='OUTLIER':
                    dict1['Model ID'].append(i)
                    dict1['Residue type'].append(k.strip().split()[-1].split(':')[0])
                    if len(k.strip().split()[0])>2:
                        dict1['Chain and res ID'].append(k.strip().split()[0])
                    else:
                        dict1['Chain and res ID'].append(':'.join(k.strip().split()[:2]))
        print (dict1['Model ID'], file=f_rota_D)
        print (dict1['Chain and res ID'],file=f_rota_D)
        print (dict1['Residue type'], file=f_rota_D)
        return dict1

    def get_data_for_quality_at_glance(self,line):
        """format mean information of models for quality at glance plots, read from temp_mp file"""
        count=[i for i,j in enumerate(line) if 'Summary' in j]
        info=[(i,line[i+1],line[i+2],line[i+3]) for i,j in enumerate(line) if 'Summary' in j]
        print (info)
        for i in range(count[0]+2,len(line)):
            m=line[i].strip().split('=')
            if 'Ramachandran outliers' in m[0]:
                rama=float(m[1].replace(' ','').replace('%',''))
            if 'Rotamer outliers' in m[0]:
                sidechain=float(m[1].replace(' ','').replace('%',''))
            if 'Clashscore' in m[0]:
                clashscore=float(m[1].replace(' ',''))
        return clashscore,rama,sidechain

