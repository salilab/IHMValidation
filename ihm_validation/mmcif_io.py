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
import os
import re
from collections import defaultdict
from itertools import chain
import utility

import logging
from typing import Final
from enum import Enum
from pathlib import Path

#########################
# Setup operational mode
#########################

logging.basicConfig(level=logging.INFO)


class IHMVAvailableModes(Enum):
    PRODUCTION = 0
    DEVELOPMENT = 1


def get_operational_mode() -> IHMVAvailableModes:
    """ Check environment variables and set the operational mode"""
    ihmv_env_name = "IHMV_MODE"  # Name of the environment variable

    try:
        ihmv_mode = IHMVAvailableModes[os.environ.get(ihmv_env_name)]  # Check the environment variable
        logging.info(f"Picked up environment variable: {ihmv_env_name}={ihmv_mode.name}")
    except KeyError:
        ihmv_mode = IHMVAvailableModes.PRODUCTION  # Default mode

    return ihmv_mode


IHMV_MODE: Final = get_operational_mode()  # Set constant for the operational mode
logging.info(f"Current operational mode is: {IHMV_MODE.name}")

# Setup default values for variables
__max_num_models = 100000  # Hopefully this value is large enough

# Alter variables for the DEVELOPMENT mode
if IHMV_MODE == IHMVAvailableModes.DEVELOPMENT:
    __max_num_models = 20  # Cap number of structures for development purposes

# Setup final values for constants
MAX_NUM_MODELS: Final = __max_num_models  # Set constant for maximum number of models in a file


#########################
# Get information from IHM reader
#########################

class GetInputInformation(object):
    def __init__(self, mmcif_file):
        self.mmcif_file = mmcif_file
        self.encoding = None
        self.datasets = {}
        self.entities = {}
        self.model = ihm.model.Model
        encoding = 'utf8'
        try:
            with open(self.mmcif_file, encoding=encoding) as fh:
                self.system, = ihm.reader.read(fh, model_class=self.model)
        except UnicodeDecodeError:
            encoding = 'ascii'
            with open(self.mmcif_file, encoding=encoding, errors='ignore') as fh:
                self.system, = ihm.reader.read(fh, model_class=self.model)

        self.encoding = encoding

    def get_databases(self):
        """ get all datasets from the mmcif file"""
        dbs = self.system.orphan_datasets
        return dbs

    def get_id(self):
        """ get id from model name, eg: PDBDEV_00XX will be PDBDEV00XX"""
        # if self.system.id == 'model':
        #    id = self.get_id_from_entry()
        # else:
        #    id = self.system.id.split('_')[0] + self.system.id.split('_')[1]
        return self.get_id_from_entry()

    def get_id_from_entry(self) -> str:
        """ get id name from entry for cif files
            deprecated """
        sf = open(self.mmcif_file, 'r', encoding='latin1')
        for ind, ln in enumerate(sf.readlines()):
            line = ln.strip().split(' ')
            if '_entry.id' in line[0]:
                entry_init = line[-1]
                try:
                    entry = entry_init.split('_')[0] + \
                        entry_init.split('_')[1]
                except:
                    entry = entry_init.strip()
        return entry

    def get_primary_citation_info(self) -> tuple:
        '''get title and authors for the primary citation'''
        title, authors = None, None
        for citation in self.system.citations:
            if citation.is_primary:
                try:
                    title = citation.title
                except AttributeError:
                    title = 'Title not available/Citation not provided'

                try:
                    authors =  '; '.join(citation.authors)
                except AttributeError:
                    authors = 'Authors are not available/Citation not provided'

        return (title, authors)


    def get_authors(self) -> str:
        """get authors of the structure; fallback to authors of primary citation """
        output = None
        if len(self.system.authors) > 0:
            output = '; '.join(self.system.authors)
        elif len(self.system.citations) > 0:
            for citation in self.system.citations:
                if citation.is_primary and len(citation.authors) > 0:
                    output =  '; '.join(citation.authors)

        if output is None:
            output = 'Authors are not available'

        return output

    def get_struc_title(self) -> str:
        """get name of molecule"""
        strc = self.system.title
        if strc is None:
            entities = self.system.entities
            mol_name = entities[0].description
        else:
            mol_name = strc
        return mol_name

    def check_sphere(self) -> int:
        """check resolution of structure,
        returns 0 if its atomic and 1 if the model is multires"""
        spheres = [len(b._spheres) for i in self.system.state_groups
                   for j in i for a in j for b in a]
        if 0 not in spheres:
            return 1
        else:
            return 0

    def get_assembly_ID_of_models(self) -> list:
        """Assembly info i.e. model assemblies in the file """
        assembly_id = [
            b.assembly._id for i in self.system.state_groups
            for j in i for a in j for b in a]
        return assembly_id

    def get_representation_ID_of_models(self) -> list:
        """Number of representations in model """
        representation_id = [
            b.representation._id for i in self.system.state_groups
            for j in i for a in j for b in a]
        return representation_id

    def get_model_names(self) -> list:
        """ Names of models"""
        model_name1 = [
            a.name for i in self.system.state_groups for j in i for a in j]
        model_name2 = [
            b.name for i in self.system.state_groups
            for j in i for a in j for b in a]
        if len(model_name1) == len(model_name2):
            model_name = [str(t[0])+'/'+str(t[1])
                          for t in zip(model_name1, model_name2)]
        else:
            model_name = model_name2
        return model_name

    def get_model_assem_dict(self) -> dict:
        """Map models to assemblies """
        model_id = [int(b._id) for i in self.system.state_groups
                    for j in i for a in j for b in a]
        assembly_id = [int(asmb) for asmb in self.get_assembly_ID_of_models()]
        model_assembly = dict(zip(model_id, assembly_id))
        return model_assembly

    def get_model_rep_dict(self) -> dict:
        """Map models to representations
            useful especially for multi-state systems"""
        model_id = [int(b._id) for i in self.system.state_groups
                    for j in i for a in j for b in a]
        rep_id = [int(repr) for repr in self.get_representation_ID_of_models()]
        model_rep = dict(zip(model_id, rep_id))
        return model_rep

    def get_number_of_models(self) -> int:
        """ Get total number of models """
        models = [
            b._id for i in self.system.state_groups
            for j in i for a in j for b in a]
        return len(models)

    def get_residues(self, asym):
        """Get residues per chain """
        if asym.seq_id_range[0] is not None:
            residues = asym.seq_id_range[1]-asym.seq_id_range[0]+1
        elif asym.seq_id_range[0] is None:
            residues = utility.NA
        return residues

    def get_composition(self) -> dict:
        """Get composition dictionary"""
        entry_comp = {'Model ID': [], 'Subunit number': [], 'Subunit ID': [],
                      'Subunit name': [], 'Chain ID': [], 'Chain ID [auth]': [],
                      'Total residues': []}
        for i, j in self.get_model_assem_dict().items():
            for m in self.system.orphan_assemblies:
                if int(m._id) == int(j):
                    count = 0
                    for n in m:
                        try:
                            count += 1
                            entry_comp['Model ID'].append(i)
                            entry_comp['Subunit number'].append(count)
                            entry_comp['Subunit ID'].append(n.entity._id)
                            entry_comp['Subunit name'].append(
                                str(n.entity.description))
                            if isinstance(n, ihm.AsymUnit):
                                aid = n.id
                                sid = n.strand_id
                            elif isinstance(n, ihm.AsymUnitRange):
                                aid = n.asym.id
                                sid = n.asym.strand_id
                            else:
                                raise(ValueError('Unexpected entity type. Only AsymUnit and AsymUnitRange are allowed'))
                            entry_comp['Chain ID'].append(aid)
                            entry_comp['Chain ID [auth]'].append(sid)
                            entry_comp['Total residues'].append(
                                self.get_residues(n))
                        except AttributeError:
                            break

        return entry_comp

    def get_protocol_number(self) -> int:
        """ number of protocols/methods used to create model"""
        return len(self.system.orphan_protocols)

    def get_sampling(self) -> dict:
        """ sampling composition/details """
        sampling_comp = {'Step number': [], 'Protocol ID': [],
                         'Method name': [], 'Method type': [],
                         'Method description': [],
                         'Number of computed models': [],
                         'Multi state modeling': [],
                         'Multi scale modeling': []}
        for prot in self.system.orphan_protocols:
            for step in prot.steps:
                sampling_comp['Step number'].append(step._id)
                sampling_comp['Multi state modeling'].append(
                    str(step.multi_state))
                sampling_comp['Multi scale modeling'].append(
                    str(step.multi_scale))
                sampling_comp['Protocol ID'].append(
                    self.system.orphan_protocols.index(prot)+1)
                # cit = self.system.citations[0].pmid
                # link = 'https://pubmed.ncbi.nlm.nih.gov/'+str(cit)+'/'
                # if step.name:
                #     method_link = '<a href='+link+'>'+str(step.method)+'</a>'
                # else:
                method_link = str(step.method)

                sampling_comp['Method description'].append(step.description)
                sampling_comp['Method name'].append(method_link)
                sampling_comp['Method type'].append(step.name)
                sampling_comp['Number of computed models'].append(
                    step.num_models_end)
        return sampling_comp

    def get_representation(self):
        """ get details on number of model composition based on
        types of representation listed """
        # representation_comp = {'Chain ID': [],
        # 'Subunit name': [], 'Rigid bodies': [],
        #  'Non-rigid regions': []}
        for rep in self.system.orphan_representations:
            # print (rep,rep[0].rigid,rep[0].
            # asym_unit.seq_id_range,rep[0].asym_unit._id)
            # print(["%s:%d-%d" % ((x.asym_unit._id,) + x.asym_unit.seq_id_range)
            #        for x in rep if not x.rigid])
            pass

    def get_RB_flex_dict(self) -> (dict, dict, int, int):
        """ get RB and flexible segments from model information"""
        RB = self.get_empty_chain_dict()
        RB_nos = []
        all_nos = []
        flex = self.get_empty_chain_dict()
        for rep in self.system.orphan_representations:
            for el in rep:
                all_nos.append(el.asym_unit.seq_id_range)
                if el.rigid and el.starting_model:
                    RB_nos.append(el.asym_unit.seq_id_range)
                    RB[el.starting_model.asym_unit._id].append(
                      [utility.format_tuple(el.asym_unit.seq_id_range)]
                    #,
                    #   utility.get_val_from_key(self.get_dataset_dict(),
                    #                            el.starting_model.dataset._id)]
                    )
                elif el.rigid and not el.starting_model:
                    RB_nos.append(el.asym_unit.seq_id_range)
                    RB[el.asym_unit._id].append(
                        [utility.format_tuple(el.asym_unit.seq_id_range),
                         'None'])
                else:
                    flex[el.asym_unit._id].append(
                        [utility.format_tuple(el.asym_unit.seq_id_range)])
        return RB, flex, len(RB_nos), len(all_nos)

    def get_number_of_assemblies(self) -> int:
        return (len(self.system.orphan_assemblies))

    def get_number_of_entities(self) -> int:
        return (len(self.system.entities))

    def get_number_of_chains(self) -> int:
        """get number of chains per protein per assembly """
        assembly = defaultdict()
        for ind, ass in enumerate(self.system.orphan_assemblies):
            chain = [el._id for el in ass]
            assembly[ind] = chain
        number_of_chains = [len(i) for i in assembly.values()]
        return number_of_chains

    def get_all_asym(self) -> list:
        """ get all asym units"""
        parents = [(a._id, a.details, a.entity.description, a.entity._id, i)
                   for i, a in enumerate(self.system.asym_units)]
        return parents

    def get_empty_chain_dict(self) -> dict:
        empty_chain_dict = defaultdict()
        for ind, el in enumerate(self.system.asym_units):
            empty_chain_dict[el._id] = []
        return empty_chain_dict

    def get_chain_subunit_dict(self) -> dict:
        """ Get chains of subunits"""
        chain_subunit_dict = defaultdict()
        for ind, el in enumerate(self.system.asym_units):
            chain_subunit_dict[el._id] = el.details.split('.')[0]
        return chain_subunit_dict

    def get_residues_subunit_dict(self) -> dict:
        """Get residues of chains in subunits"""
        residues_subunit_dict = defaultdict()
        for el in self.system.asym_units:
            residues_subunit_dict[el._id] = self.get_residues(el)
        return residues_subunit_dict

    def get_software_length(self) -> int:
        lists = self.system.software
        if lists is None:
            return 0
        else:
            return len(lists)

    def get_software_comp(self) -> dict:
        """get software composition to write out as a table"""
        software_comp = {'ID': [], 'Software name': [], 'Software version': [
        ], 'Software classification': [], 'Software location': []}
        lists = self.system.software
        self.read_all_references()
        if len(lists) > 0:
            for software in lists:
                software_comp['ID'].append(software._id)

                # Get the name of software from mmCIF
                ref_name = software.name.lower()

                # Check if URL for sotware was defined in the mmCIF and parsed
                # and parsed by ihm
                if software.location is None:
                    try:
                        self.ref_link[ref_name]
                    except KeyError:
                        logging.warning(
                            # LookupError(
                            f'The URL for {software.name} is missing from '
                            'from both mmCIF and references.csv')
                            # )
                    else:
                        logging.debug(
                            f'Filling the url for {software.name} from '
                            'references.csv'
                            )

                        software.location = self.ref_link[ref_name]

                if software.location is None:
                    ref_tot = f'{software.name}'
                    ref_loc = utility.NA
                else:

                    ref_tot = f'<a href="{software.location}">' \
                                f'{software.name}</a>'

                    ref_loc = f'<a href="{software.location}">' \
                                f'{software.location}</a>'

                software_comp['Software name'].append(ref_tot)
                software_comp['Software location'].append(ref_loc)
                if str(software.version) == ihm.unknown:
                    vers = utility.NA
                elif str(software.version) == '?':
                    vers = utility.NA
                elif software.version is None:
                    vers = utility.NA
                else:
                    vers = str(software.version)
                software_comp['Software version'].append(vers)
                software_comp['Software classification'].append(
                    software.classification)
            final_software_composition = software_comp
        else:
            final_software_composition = {}
        return final_software_composition


    def read_all_references(self) -> None:
        self.ref_link = dict()
        self.ref_cit = dict()

        template_path = Path(Path(__file__).parent.parent.resolve(), 'templates')
        reference_filename = str(Path(template_path, 'references.csv'))

        with open(reference_filename, 'r') as f:
            allref = [_.strip().split('|') for _ in f.readlines()]

        for line in allref:
            self.ref_link[line[0].lower().rstrip()] = line[1].rstrip().lstrip()
            self.ref_cit[line[0].lower()] = line[2]

    def check_ensembles(self) -> int:
        """check if ensembles exist"""
        return len(self.system.ensembles)

    def get_ensembles(self):
        """details on ensembles, if it exists"""
        if len(self.system.ensembles) > 0:
            ensemble_comp = {'Ensemble number': [],
                             'Ensemble name': [],
                             'Model ID': [],
                             'Number of models': [],
                             'Clustering method': [],
                             'Clustering feature': [],
                             'Cluster precision': []}
            for ensm in self.system.ensembles:
                ensemble_comp['Ensemble number'].append(str(ensm._id))
                ensemble_comp['Ensemble name'].append(str(ensm.name))
                try:
                    ensemble_comp['Model ID'].append(str(ensm.model_group._id))
                except AttributeError:
                    ensemble_comp['Model ID'].append(utility.NA)

                ensemble_comp['Number of models'].append(str(ensm.num_models))
                ensemble_comp['Clustering method'].append(
                    str(ensm.clustering_method))
                ensemble_comp['Clustering feature'].append(
                    str(ensm.clustering_feature))
                ensemble_comp['Cluster precision'].append(str(ensm.precision))
            return ensemble_comp
        else:
            return None

    def get_dataset_xl_info(self, id: int) -> str:
        """Get dataset XL info given dataset ID"""
        restraints = self.get_restraints()
        return 'Linker name and number of cross-links: %s' % (restraints
                                                              ['Restraint info']
                                                              [restraints['Dataset ID']
                                                               .index(id)])

    def get_dataset_dict(self):
        """get dataset dictionary """
        dataset_dict = defaultdict()
        lists = self.system.orphan_datasets
        if len(lists) > 0:
            for _ in lists:
                try:
                    acc = _.location.access_code
                except AttributeError:
                    acc = str('None')
                dataset_dict[_._id] = str(_.data_type)+'/'+str(acc)
        return dataset_dict

    def get_dataset_length(self) -> int:
        lists = self.system.orphan_datasets
        if lists is None:
            return 0
        else:
            return len(lists)

    def get_dataset_comp(self) -> dict:
        """detailed dataset composition"""
        dataset_comp = {'ID': [], 'Dataset type': [],
                        'Database name': [], 'Data access code': []}
        lists = self.system.orphan_datasets
        if len(lists) > 0:
            for _ in lists:
                try:
                    loc = _.location.db_name
                except AttributeError:
                    if isinstance(_.location, ihm.location.InputFileLocation) \
                            or isinstance(_.location, ihm.location.FileLocation):
                        loc = 'File'
                    else:
                        loc = utility.NA
                try:
                    acc = _.location.access_code
                except AttributeError:
                    if isinstance(_.location, ihm.location.InputFileLocation) \
                            or isinstance(_.location, ihm.location.FileLocation):
                        try:
                            acc = _.location.repo.doi
                        except AttributeError:
                            acc = utility.NA
                    else:
                        acc = utility.NA
                dataset_comp['ID'].append(_._id)
                # if i.data_type=='unspecified' and 'None' not in i.details:
                #    dataset_comp['Dataset type'].append(i.details)
                # else:
                dataset_comp['Dataset type'].append(_.data_type)
                dataset_comp['Database name'].append(str(loc))
                dataset_comp['Data access code'].append(acc)
        return dataset_comp

    def dataset_id_type_dic(self) -> dict:
        """dataset id and data items"""
        dataset_dic = defaultdict()
        if len(self.system.orphan_datasets) > 0:
            for i in self.system.orphan_datasets:
                if i.data_type == 'other':
                    dataset_dic[str(i._id)] = str(i.details)
                else:
                    try:
                        dataset_dic[str(i._id)] = str(i.data_type)
                    except TypeError:
                        dataset_dic[str(i._id)] = 'None'
        return dataset_dic

    def get_restraints(self) -> dict:
        """ get restraints table from cif file"""
        r = self.system.restraints
        restraints_comp = {'ID': [], 'Dataset ID': [],
                           'Restraint type': [], 'Restraint info': []}
        for j, i in enumerate(r):
            restraints_comp['ID'].append(j+1)
            restraints_comp['Restraint type'].append(str(i.__class__.__name__))
            try:
                restraints_comp['Dataset ID'].append(str(i.dataset._id))
            except AttributeError:
                restraints_comp['Dataset ID'].append('None')
            if isinstance(i, ihm.restraint.CrossLinkRestraint):
                restraints_comp['Restraint info'].append(
                    str(i.linker.auth_name) + ', ' +
                    str(len(i.experimental_cross_links)) + ' cross-links')
            elif isinstance(i, ihm.restraint.EM3DRestraint):
                restraints_comp['Restraint info'].append(
                    str(i.fitting_method)) # + ', '+str(i.number_of_gaussians) + ' components')

            elif isinstance(i, ihm.restraint.PredictedContactRestraint):
                # Temporary fix for Entry 161;
                # Should be unified with DerivedRestraint
                if isinstance(i.distance, ihm.restraint.LowerUpperBoundDistanceRestraint):
                    restraints_comp['Restraint info'].append(
                        ('Lower Upper Bound Distance: '+str(i.distance.distance_lower_limit)+'-' +
                         str(i.distance.distance_upper_limit)))
                else:
                    restraints_comp['Restraint info'].append('Distance: '+str(i.distance.distance))
#                                                         + ' between residues ' +
#                                                         str(i.resatom1.seq_id)
#                                                         + ' and ' + str(i.resatom2.seq_id))

            elif isinstance(i, ihm.restraint.EM2DRestraint):
                restraints_comp['Restraint info'].append('Number of micrographs: '
                                                         + str(i.number_raw_micrographs)
                                                         + ',' + ' Image resolution: '
                                                         + str(i.image_resolution))
            elif isinstance(i, ihm.restraint.SASRestraint):
                restraints_comp['Restraint info'].append('Assembly name: '+str(
                    i.assembly.name)+' Fitting method: ' +
                    str(i.fitting_method) + ' Multi-state: ' + str(i.multi_state))
            elif isinstance(i, ihm.restraint.UpperBoundDistanceRestraint):
                restraints_comp['Restraint info'].append(
                    'Distance: '+str(i.distance))
            elif 'Mutagenesis' in str(i.__class__.__name__):
                restraints_comp['Restraint info'].append(
                    'Details: '+str(i.details))
            elif isinstance(i, ihm.restraint.DerivedDistanceRestraint):
                dic = self.dataset_id_type_dic()
                try:
                    ID = str(i.dataset._id)
                except AttributeError:
                    ID = 'None'
                # restraints_comp['Restraint info'].append(dic[ID])
                if isinstance(i.distance, ihm.restraint.UpperBoundDistanceRestraint):
                    restraints_comp['Restraint info'].append(
                        ('Upper Bound Distance: '+str(i.distance.distance)))
                elif isinstance(i.distance, ihm.restraint.LowerUpperBoundDistanceRestraint):
                    restraints_comp['Restraint info'].append(
                        ('Lower Upper Bound Distance: '+str(i.distance.distance_lower_limit)+'-' +
                         str(i.distance.distance_upper_limit)))
                else:
                    restraints_comp['Restraint info'].append(
                        'restraint type ' + str(i.distance.__class__.__name__) + str(dic[ID]))
                '''
                if 'UpperBound' in str(i.distance.__class__.__name__):
                    print (i.distance,i.distance.distance_lower_limit)
                    restraints_comp['Restraint info'].append(
                        ('Upper Bound Distance: '+str(i.distance.distance)))
                '''
        return restraints_comp

    def get_dataset_details(self) -> dict:
        """get information on dataset and databases"""
        dataset_comp = {'ID': [], 'Dataset type': [],
                        'Database name': [], 'Details': []}
        lists = self.system.orphan_datasets
        if len(lists) > 0:
            for i in lists:
                try:
                    loc = i.location.db_name
                except AttributeError:
                    loc = str('')
                try:
                    acc = i.location.access_code
                except AttributeError:
                    acc = utility.NA

                dataset_comp['ID'].append(i._id)
                if i.data_type == 'unspecified' and i.details is not None:
                    dataset_comp['Dataset type'].append(i.details)
                else:
                    dataset_comp['Dataset type'].append(i.data_type)
                dataset_comp['Database name'].append(str(loc))
                if 'ComparativeModel' in str(i.__class__.__name__):
                    acc1 = 'template PDB ID: ' + acc
                    dataset_comp['Details'].append(acc1)
                elif 'PDB' in str(i.__class__.__name__):
                    acc1 = 'PDB ID: ' + acc
                    dataset_comp['Details'].append(acc1)
                elif 'CX' in str(i.__class__.__name__):
                    acc1 = self.get_dataset_xl_info(i._id)
                    dataset_comp['Details'].append(acc1)
                elif 'EMDB' in str(i.__class__.__name__):
                    acc1 = 'EMDB ID: ' + acc
                    dataset_comp['Details'].append(acc1)
                else:
                    if acc is not utility.NA:
                        dataset_comp['Details'].append(f'{loc}: {acc}')
                    elif isinstance(i.location, ihm.location.FileLocation):
                        doi = utility.NA
                        try:
                            doi = i.location.repo.doi
                        except AttributeError:
                            doi = utility.NA

                        acc1 = f'File: {doi}'
                        dataset_comp['Details'].append(acc1)

                    else:
                        dataset_comp['Details'].append(acc)


        return dataset_comp

    def get_atomic_coverage(self) -> str:
        """Measure amount of atomic residues"""
        for _ in self.system.orphan_representations:
            if self.check_sphere() == 1:
                flex = sum([(x.asym_unit.seq_id_range[1]) -
                            (x.asym_unit.seq_id_range[0])+1 for x in _ if not x.rigid])
                rigid = sum([(x.asym_unit.seq_id_range[1]) -
                             (x.asym_unit.seq_id_range[0])+1 for x in _ if x.rigid])

                if rigid > 0 or flex > 0:
                    percentage = str(round(rigid/(rigid+flex)*100))+'%'
                else:
                    percentage = '0%'
            else:
                percentage = '100%'
        return percentage

    def check_for_sas(self, dataset: dict) -> bool:
        """check if sas is in the dataset"""
        dataset = self.get_dataset_comp()
        data_type = dataset['Dataset type']
        database = dataset['Database name']
        return 'SAS' in str(data_type) and 'SAS' in str(database)

    def check_for_cx(self, dataset: dict) -> bool:
        """check if CX-XL is in the dataset"""
        dataset = self.get_dataset_comp()
        data_type = dataset['Dataset type']
        if 'CX' in str(data_type):
            return True
        else:
            return False

    def check_for_em(self, dataset: dict) -> bool:
        """check if em is in the dataset"""
        dataset = self.get_dataset_comp()
        data_type = dataset['Dataset type']
        if 'EM' in str(data_type):
            return True
        else:
            return False

    def mmcif_get_lists(self, filetemp=None) -> (list, dict, dict, list):
        """function to help re-write mmcif file for molprobity
        this function reads the atom_site dictionary terms and returns a list"""
        if filetemp is None:
            file = open(self.mmcif_file, 'r', encoding='latin1')
        else:
            file = filetemp
            filetemp.seek(0)
        all_lines = []
        for i, j in enumerate(file.readlines()):
            all_lines.append(j.strip().split())
        atom_site = {}
        atoms = {}
        before_atom_site = []
        after_atom = []
        for i, j in enumerate(all_lines):
            if len(j) > 0 and '_atom_site.' in j[0]:
                if len(before_atom_site) == 0:
                    before_atom_site = all_lines[:i+1]
                atom_site[i] = j[0]
            elif ('_atom_site.B_iso_or_equiv' not in list(atom_site.values())) and len(list(atom_site.values())) > 0:
                atom_site[i] = '_atom_site.B_iso_or_equiv'
            elif ('_atom_site.occupancy' not in list(atom_site.values())) and len(list(atom_site.values())) > 0:
                atom_site[i] = '_atom_site.occupancy'

            elif ('_atom_site.label_seq_id' not in list(atom_site.values())) and len(list(atom_site.values())) > 0:
                atom_site[i] = '_atom_site.label_seq_id'

        total_list = list(atom_site.values())
        index_biso = total_list.index('_atom_site.B_iso_or_equiv')
        index_occu = total_list.index('_atom_site.occupancy')
        index_label_seq = total_list.index('_atom_site.label_seq_id')
        for i, j in enumerate(all_lines):
            if len(j) > 0 and ('ATOM' in j[0] or 'HETATM' in j[0]) and (i > list(atom_site.keys())[-1]):
                if len(j) <= index_occu:
                    j.extend(['1'])
                elif j[index_occu] == '.':
                    j[index_occu] = '0.67'
                if len(j) <= index_biso:
                    j.extend(['1'])
                elif j[index_biso] == '.':
                    j[index_biso] = '0.00'
                if len(j) <= index_label_seq:
                    j.extend(['1'])
                elif j[index_label_seq] == '.':
                    j[index_label_seq] = str(i)
                atoms[i] = j
            elif len(j) > 0 and (i > list(atom_site.keys())[-1]):
                if len(after_atom) == 0:
                    after_atom = all_lines[i:]
        return before_atom_site, atom_site, atoms, after_atom

    def delete_extra_loops(self, some_text=list()) -> list:
        """function to help re-write mmcif file for molprobity
        this cleans up extra loops in the cif file"""
        new_text = []
        line = 0
        total_lines = len(some_text)
        while line < total_lines:
            if some_text[line] and 'loop_' in some_text[line][0] and '_' not in some_text[line+1][0]:
                skip = line+2
                line = min(total_lines-1, skip)
            new_text.append(some_text[line])
            line += 1
        return new_text

    def remove_flr(self, some_text=list()) -> list:
        """function to help re-write mmcif file for molprobity
        this deletes all flr related text/key-value pairs"""
        new_text = []
        line = 0
        total_lines = len(some_text)
        while line < total_lines:
            if some_text[line] and '_flr' in some_text[line][0]:
                while '#' not in some_text[line][0]:
                    line += 1
            new_text.append(some_text[line])
            line += 1
        return new_text

    def get_model_id_column(self, atom_site: dict()) -> int:
        for colnum, val in enumerate(atom_site.values()):
            if val.lower() in ('_atom_site.ihm_model_id',
                               '_atom_site.pdbx_pdb_model_num'):
                return colnum
        return None

    def get_representation_details(self) -> dict:
        """Extract details about representation (atomic/coarse-grained)"""
        reprs = {'atomic': False, 'coarse-grained': False, 'coarse-grain_levels': []}
        for rep in chain.from_iterable(self.system.orphan_representations):
            if rep.granularity == 'by-atom':
                reprs['atomic'] = True
            elif rep.granularity == 'by-residue':
                reprs['coarse-grained'] = True
                reprs['coarse-grain_levels'].append(1)
            elif rep.granularity == 'by-feature':
                reprs['coarse-grained'] = True
                reprs['coarse-grain_levels'].append(rep.count)
            else:
                raise ValueError('Wrong representation granularity')

        if len(reprs['coarse-grain_levels']) > 1:
            levels = sorted(set(reprs['coarse-grain_levels']))
            reprs['coarse-grain_levels'] = levels

        return reprs

    def get_auth_label_map(self, system=None) -> dict:
        """get map between auth_seq_id and label_seq_id"""

        logging.info('Building auth <-> label map')
        if system is None:
            system = self.system

        emap = {}
        for st in system.state_groups:
            for s in st:
                for mg in s:
                    for m in mg:
                        for a in m.get_atoms():
                            label_asym_id = a.asym_unit.id
                            auth_asym_id = a.asym_unit.strand_id

                            label_seq_id = str(a.seq_id)
                            auth_seq_id = str(a.asym_unit.residue(a.seq_id).auth_seq_id)

                            key = (auth_asym_id, auth_seq_id)
                            val = (label_asym_id, label_seq_id)

                            if key not in emap:
                                emap[key] = val
                            else:
                                if emap[key] != val:
                                    logging.warning(f"Conflicting residue ids. auth, label1, label2 {key}, {emap[key]}, {val}. Keeping first {key} -> {emap[key]}")

        return emap
