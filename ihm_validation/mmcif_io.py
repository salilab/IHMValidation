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

import logging
import os, sys
import re
from pathlib import Path
from enum import Enum
from typing import Final
from collections import defaultdict
from itertools import chain

import ihm
import utility

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
    nocache = False
    def __init__(self, mmcif_file, cache='.', nocache=False):
        self.mmcif_file = mmcif_file
        self.datasets = {}
        self.entities = {}

        if not Path(cache).is_dir():
            try:
                os.makedirs(cache)
                logging.info(f"Created cache directory {cache}")
            except OSError:
                logging.error(f"Couldn't create cache directory {cache}")
                sys.exit(1)

        self.cache = cache
        self.nocache = nocache

        self.system, self.encoding = utility.parse_ihm_cif(mmcif_file)
        self.stem = Path(self.mmcif_file).stem
        self.ID = self.get_id()
        self.ID_f = self.get_file_id()

    def get_databases(self):
        """ get all datasets from the mmcif file"""
        dbs = self.system.orphan_datasets
        return dbs

    def get_id(self):
        """Return _entry.id; Requires compliant CIF file"""
        ids = self.get_ranked_id_list()

        if len(ids) == 0:
            raise(ValueError('Missing system ID'))

        id_type, entry_id = ids[0]

        return entry_id

    def get_file_id(self):
        """Return _entry.id; Requires compliant CIF file"""
        ids = self.get_ranked_id_list()

        if len(ids) == 0:
            raise(ValueError('Missing system ID'))

        id_type, entry_id = ids[0]

        if id_type == 'PDB ID':
            # PDB filenames have to be lowercase
            entry_id = entry_id.lower()
        elif id_type == 'PDB-Dev ID':
            # PDB-Dev filenames have to be uppercase
            entry_id = entry_id.upper()
        else:
            # Use entry ID as is
            pass

        return entry_id

    def get_pdb_id(self) -> str:
        """Check database2 table for PDB ID"""
        entry_id = None
        if len(self.system.databases) > 0:
            for db in self.system.databases:
                if db.id == 'PDB':
                    entry_id =  db.code.upper()

        return entry_id

    def get_pdb_dev_id(self) -> str:
        """Check database2 table for PDB ID"""
        entry_id = None
        if len(self.system.databases) > 0:
            for db in self.system.databases:
                if db.id == 'PDB-Dev':
                    entry_id =  db.code.upper()

        return entry_id

    def get_ranked_id_list(self) -> list:
        """Get sorted list of multiple ids"""
        ids = []
        entry_id = self.system.id

        # Sort primary/secondary IDs
        # If we have database2 table
        if len(self.system.databases) > 0:
            pdb_id = self.get_pdb_id()
            pdb_dev_id = self.get_pdb_dev_id()
            # if PDB is a primary
            if pdb_id is not None and entry_id == pdb_id:
                k = ('PDB ID', pdb_id)
                ids.append(k)

                if pdb_dev_id is not None:
                    k = ('PDB-Dev ID', pdb_dev_id)
                    ids.append(k)

            # if PDB-Dev is a primary
            elif pdb_dev_id is not None and entry_id == pdb_dev_id:
                k = ('PDB-Dev ID', pdb_dev_id)
                ids.append(k)

                if pdb_id is not None:
                    k = ('PDB ID', pdb_id)
                    ids.append(k)
            # Else entity is a primary
            else:
                k = ('Entry ID', entry_id)
                ids.append(k)

                if pdb_id is not None:
                    k = ('PDB ID', pdb_id)
                    ids.append(k)

                if pdb_dev_id is not None:
                    k = ('PDB-Dev ID', pdb_dev_id)
                    ids.append(k)

        # Compatibility mode for old ID scheme
        else:
            if re.match('PDBDEV_', entry_id):
                k = ('PDB-Dev ID', entry_id)
                ids.append(k)
            else:
                k = ('Entry ID', entry_id)
                ids.append(k)

        return ids


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
                if el.rigid:
                    RB_nos.append(el.asym_unit.seq_id_range)
                    RB[el.asym_unit._id].append(
                      [utility.format_tuple(el.asym_unit.seq_id_range)]
                    )
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

                try:
                    p = float(ensm.precision)
                except (ValueError, TypeError):
                    p = None
                ensemble_comp['Cluster precision'].append(p)
            return ensemble_comp
        else:
            return None

    def get_dataset_xl_info(self, id: int) -> str:
        """Get dataset XL info given dataset ID"""
        restraints = self.get_restraints()
        return 'Linker name and number of crosslinks: %s' % (restraints
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
                    acc = utility.NA
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
                if isinstance(_.location, ihm.location.DatabaseLocation):
                    try:
                        loc = _.location.db_name
                    except AttributeError:
                        loc = utility.NA

                    try:
                        acc = _.location.access_code
                    except AttributeError:
                        acc = utility.NA

                    if isinstance(_.location, ihm.location.PDBDevLocation) and acc != utility.NA:
                        acc = f"<a href=https://pdb-dev.wwpdb.org/entry.html?{acc}>{acc}</a>"

                    if isinstance(_.location, ihm.location.PDBLocation) and acc != utility.NA:
                        acc = f"<a href=https://www.rcsb.org/structure/{acc}>{acc}</a>"

                    if isinstance(_.location, ihm.location.ModelArchiveLocation) and acc != utility.NA:
                        acc = f"<a href=https://doi.org/10.5452/{acc}>{acc}</a>"

                    if isinstance(_.location, ihm.location.AlphaFoldDBLocation) and acc != utility.NA:
                        acc = f"<a href=https://alphafold.ebi.ac.uk/entry/{acc}>{acc}</a>"

                    if isinstance(_.location, ihm.location.EMPIARLocation) and acc != utility.NA:
                        acc = f"<a href=https://www.ebi.ac.uk/empiar/{acc}>{acc}</a>"

                    if isinstance(_.location, ihm.location.EMDBLocation) and acc != utility.NA:
                        acc = f"<a href=https://emdb-empiar.org/{acc}>{acc}</a>"

                    if isinstance(_.location, ihm.location.SASBDBLocation) and acc != utility.NA:
                        acc = f"<a href=https://www.sasbdb.org/data/{acc}>{acc}</a>"

                    if isinstance(_.location, ihm.location.PRIDELocation) and acc != utility.NA:
                        acc = f"<a href=https://www.ebi.ac.uk/pride/archive/projects/{acc}>{acc}</a>"

                    if isinstance(_.location, ihm.location.MassIVELocation) and acc != utility.NA:
                        acc = f"<a href=https://www.omicsdi.org/dataset/massive/{acc}>{acc}</a>"

                    if isinstance(_.location, ihm.location.ProteomeXchangeLocation) and acc != utility.NA:
                        acc = f"<a href=https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID={acc}>{acc}</a>"

                elif isinstance(_.location, ihm.location.FileLocation):
                    try:
                        loc = _.location.repo.reference_provider
                        if loc is None:
                            loc = utility.NA
                        doi = _.location.repo.doi
                        acc = f"<a href=https://doi.org/{doi}>{doi}</a>"
                    except AttributeError as e:
                        loc = 'File'
                        acc = utility.NA
                    except TypeError as e:
                        loc = utility.NA
                        acc = utility.NA
                        logging.warning('Missing repository data')
                        logging.warning(e)
                    except Exception as e:
                        logging.error(f'Unexepcted dataset error')
                        logging.error(e)

                else:
                    loc = utility.NA
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
                        dataset_dic[str(i._id)] = utility.NA
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
                restraints_comp['Dataset ID'].append(utility.NA)
            if isinstance(i, ihm.restraint.CrossLinkRestraint):
                restraints_comp['Restraint info'].append(
                    str(i.linker.auth_name) + ', ' +
                    str(len(i.experimental_cross_links)) + ' crosslinks')
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
                    ID = utility.NA
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
                elif i.data_type == 'CX-MS data':
                    dataset_comp['Dataset type'].append('Crosslinking-MS')
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
        # Martini-like coarse-granining (multiple beads per residue)
        # is not supported by the dictionary
        reprs = []
        for rep in self.system.orphan_representations:
            reprs_ = {'atomic': False, 'coarse-grained': False, 'coarse-grain_levels': []}
            for rep_ in rep:
                if rep_.granularity == 'by-atom':
                    reprs_['atomic'] = True
                elif rep_.granularity == 'by-residue':
                    reprs_['coarse-grained'] = True
                    reprs_['coarse-grain_levels'].append(1)
                elif rep_.granularity == 'by-feature':
                    reprs_['coarse-grained'] = True
                else:
                    raise ValueError('Wrong representation granularity')
            if reprs_['coarse-grained']:
                for stg in self.system.state_groups:
                    for st in stg:
                        for mg in st:
                            for m in mg:
                                if m.representation == rep:
                                    for s in m.get_spheres():
                                        s_size = utility.get_bead_size(s)
                                        reprs_['coarse-grain_levels'].append(s_size)


            if len(reprs_['coarse-grain_levels']) > 1:
                levels = sorted(set(reprs_['coarse-grain_levels']))
                reprs_['coarse-grain_levels'] = levels
            reprs.append(reprs_)

        return reprs

    def get_representation_scale(self, rep, chid=None) -> dict:
        """Extract details about representation (atomic/coarse-grained)"""
        # Martini-like coarse-granining (multiple beads per residue)
        # is not supported by the dictionary
        reprs_ = {'atomic': False, 'coarse-grained': False, 'coarse-grain_levels': []}
        for rep_ in rep:
            if chid:
                x = rep_.asym_unit

                if isinstance(x, ihm.AsymUnitRange):
                    x = x.asym

                aid = x.id
                sid = x.strand_id

                chid_ = (aid, sid)

                if chid_ != chid:
                    continue

            if rep_.granularity == 'by-atom':
                reprs_['atomic'] = True
            elif rep_.granularity == 'by-residue':
                reprs_['coarse-grained'] = True
                reprs_['coarse-grain_levels'].append(1)
            elif rep_.granularity == 'by-feature':
                reprs_['coarse-grained'] = True
            else:
                raise ValueError('Wrong representation granularity')

        if reprs_['coarse-grained']:
            for stg in self.system.state_groups:
                for st in stg:
                    for mg in st:
                        for m in mg:
                            if m.representation == rep:
                                for s in m.get_spheres():
                                    if chid:
                                        x = s.asym_unit

                                        if isinstance(x, ihm.AsymUnitRange):
                                            x = x.asym

                                        aid = x.id
                                        sid = x.strand_id

                                        chid_ = (aid, sid)

                                        if chid_ != chid:
                                            continue

                                    s_size = utility.get_bead_size(s)
                                    reprs_['coarse-grain_levels'].append(s_size)


        if len(reprs_['coarse-grain_levels']) > 1:
            levels = sorted(set(reprs_['coarse-grain_levels']))
            reprs_['coarse-grain_levels'] = levels

        return reprs_

    def pretty_print_scale(scale, reprs_: dict) -> list:
        """Pretty print information about representation scales"""

        out = ''
        if ((reprs_['atomic'] and reprs_['coarse-grained']) or
            (reprs_['coarse-grained'] and len(reprs_['coarse-grain_levels']) > 1)):
            out += 'Multiscale: '

        if reprs_['atomic']:
            out += 'Atomic'

        if reprs_['coarse-grained']:
            if out != '':
                if out[-1] != ' ':
                    out += '; '

            out += 'Coarse-grained: '
            levels = reprs_['coarse-grain_levels']
            if len(levels) == 1:
                out += f'{levels[0]:d}'
            else:
                min_level = min(levels)
                max_level = max(levels)
                out += f'{min_level:d} - {max_level:d}'
            out += ' residue(s) per bead'

        return out

    def get_representation_info(self):
        output = []

        for i, rep in enumerate(self.system.orphan_representations, 1):
            rep_ = {'ID': i, 'Models':[], 'Chains': {}}

            models_ = []

            for stg in self.system.state_groups:
                for st in stg:
                    for mg in st:
                        for m in mg:
                            if m.representation == rep:
                                models_.append(int(m._id))

            models__ = []

            for x in utility.ranges(models_):
                if x[0] == x[1]:
                    models__.append(x[0])
                else:
                    models__.append(utility.format_tuple(x))

            rep_['Models'] = models__

            chain_ids_raw_ = []
            for s in rep:
                x = s.asym_unit

                if isinstance(x, ihm.AsymUnitRange):
                    x = x.asym

                aid = x.id
                sid = x.strand_id

                chid = (aid, sid)

                chain_ids_raw_.append(chid)

            chain_ids_ = sorted(set(chain_ids_raw_), key=lambda x: f'{len(x[0])}{x[0]}')
            chains_ = {x: {
                # 'Subunit ID': None,
                'Entity ID': None,
                'Molecule name': None,
                'Chains': None,
                'Total residues': None,
                'Rigid segments': [],
                'Flexible segments': [],
                'Identical subunits': [],
                'Model coverage': None,
                'Starting model segments': [],
                'Starting model coverage': None,
                } for x in chain_ids_}

            for k in chains_.keys():
                aid, sid = k
                if aid == sid:
                    chid = aid
                else:
                    chid = f'{aid} [{sid}]'
                chains_[k]['Chains'] = chid

            for s in rep:
                x = s.asym_unit

                if isinstance(x, ihm.AsymUnitRange):
                    x = x.asym

                aid = x.id
                sid = x.strand_id

                chid = (aid, sid)

                if x.entity.is_polymeric():
                    tr_ = x.seq_id_range
                    tr = int(tr_[1] - tr_[0] + 1)
                else:
                    tr = 'Non-polymeric'

                attributes = {
                    # 'Subunit ID':
                    'Entity ID': s.asym_unit.entity._id,
                    'Molecule name': s.asym_unit.entity.description,
                    'Total residues': tr,
                }

                for k, v in attributes.items():
                    if chains_[chid][k] is None:
                        chains_[chid][k] = v
                    else:
                        assert chains_[chid][k] == v

                if x.entity.is_polymeric():
                    seq_range_ = s.asym_unit.seq_id_range

                    if s.rigid:
                        chains_[chid]['Rigid segments'].append(seq_range_)
                    else:
                        chains_[chid]['Flexible segments'].append(seq_range_)

                    if s.starting_model:
                        chains_[chid]['Starting model segments'].append(seq_range_)

            subkeys = ['Rigid segments', 'Flexible segments', 'Starting model segments']

            for ci, ki in enumerate(chain_ids_):
                for cj, kj in enumerate(chain_ids_[ci + 1:]):
                    if chains_[ki]['Entity ID'] == chains_[kj]['Entity ID']:
                        sub_ki = {k: chains_[ki][k] for k in subkeys}
                        sub_kj = {k: chains_[kj][k] for k in subkeys}
                        # print(sub_ki, sub_kj)
                        if sub_ki == sub_kj:
                            chains_[ki]['Identical subunits'].append(kj)

            duplicates = []

            for k in chain_ids_:
                if k in duplicates:
                    try:
                        del chains_[k]
                    except KeyError:
                        pass
                else:
                    chains_dup_ = sorted(set(chains_[k].pop('Identical subunits')), key=lambda x: f'{len(x[0])}{x[0]}')
                    chains_dup__ = [chains_[x]['Chains'] for x in chains_dup_]
                    chains_[k]['Chains'] = [chains_[k]['Chains']] + chains_dup__

                    duplicates.extend(chains_dup_)
                    duplicates = list(set(duplicates))

            for k, chain_ in chains_.items():
                tr = chain_['Total residues']

                modeled = []

                for sgn in ['Rigid segments', 'Flexible segments']:
                    segs = chain_[sgn]
                    segs_ = []
                    modeled_ = []

                    for s_ in sorted(segs):
                        if s_[0] == s_[1]:
                            modeled_.append(s_[0])
                            segs_.append(f'{s_[0]:d}')
                        else:
                            modeled_.extend(list(range(s_[0], s_[1] + 1)))
                            segs_.append(utility.format_tuple(s_))

                    chains_[k][sgn] = segs_
                    modeled.extend(modeled_)

                modeled_ = []
                for s_ in chains_[k].pop('Starting model segments'):
                    if s_[0] == s_[1]:
                        modeled_.append(s_[0])
                    else:
                        modeled_.extend(list(range(s_[0], s_[1] + 1)))

                modeled_count = len(set(sorted(modeled)))

                if isinstance(tr, (int, float)):
                    coverage = modeled_count / tr * 100.0
                    starting_model_coverage = len(set(sorted(modeled_))) / modeled_count * 100
                else:
                    coverage = utility.NA
                    starting_model_coverage = utility.NA
                chains_[k]['Model coverage'] = coverage
                chains_[k]['Starting model coverage'] = starting_model_coverage

                scale = self.get_representation_scale(rep, k)
                pretty_scale = self.pretty_print_scale(scale)
                chains_[k]['Scale'] = pretty_scale

            rep_['Chains'] = chains_

            output.append(rep_)

        return output

    @property
    def num_models(self):
        return sum(1 for group, model in self.system._all_models())

    @property
    def cg(self):
        return utility.is_cg(self.system)

    @property
    def atomic(self):
        return utility.is_atomic(self.system)

    def _has_dataset_type(self, dataset_type):
        flag = False
        for i, d in enumerate(self.system.orphan_datasets):
            if isinstance(d, dataset_type):
                flag = True
                break

        return flag

    @property
    def has_crosslinking_ms_dataset(self):
        return self._has_dataset_type(ihm.dataset.CXMSDataset)

    @property
    def has_sas_dataset(self):
        return self._has_dataset_type(ihm.dataset.SASDataset)

    @property
    def deposition_date(self):
        """Return initial deposition date"""
        date = utility.NA
        try:
            date = self.system._database_status['recvd_initial_deposition_date']
        except KeyError:
            logging.error('Deposition date is missing')

        return date
