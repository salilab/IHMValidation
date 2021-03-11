import os,sys
import unittest  
from io import StringIO, BytesIO
sys.path.insert(0, "../master/pyext/src/")
from validation import get_input_information,utility
from validation.sas import sas_validation
import warnings,tempfile


def ignore_warnings(test_func):
	def do_test(self, *args, **kwargs):
		with warnings.catch_warnings():
			warnings.simplefilter("ignore", ResourceWarning)
			test_func(self, *args, **kwargs)
	return do_test

class Testing(unittest.TestCase):
	def __init__(self, *args, **kwargs):
		super(Testing, self).__init__(*args, **kwargs)
		#self.mmcif_test_file='test.cif'
		#self.IO=cx_validation(self.mmcif_test_file)

	def test_get_SASBDB_code(self):
		with tempfile.TemporaryDirectory() as tempdir:
			tmpfilepath = os.path.join(tempdir, 'someFileInTmpDir.txt')

			with open(tmpfilepath, 'w') as tmpfile:
				tmpfile.write("""

data_PDBDEV_test
_entry.id PDBDEV_test
_struct.title 'Structure irrelevant'
loop_
_ihm_dataset_list.id
_ihm_dataset_list.data_type
_ihm_dataset_list.database_hosted
_ihm_dataset_list.details
1 'SAS data' YES .
loop_
_ihm_dataset_related_db_reference.id
_ihm_dataset_related_db_reference.dataset_list_id
_ihm_dataset_related_db_reference.db_name
_ihm_dataset_related_db_reference.accession_code
_ihm_dataset_related_db_reference.version
_ihm_dataset_related_db_reference.details
1 1 SASBDB SASDC29 . .
					""")
			temp=sas_validation(tmpfilepath)
			self.assertEqual(['SASDC29'],temp.get_SASBDB_code())
			self.assertEqual(['SASDC29'],temp.clean_SASBDB_code())

			data_dic=temp.get_data_from_SASBDB()
			self.assertEqual('SASDC29',data_dic['SASDC29']['code'])
			self.assertEqual(55.07,data_dic['SASDC29']['fits'][0]['models'][0]['model_mw'])

			temp.get_sascif_file()
			self.assertRaises(TypeError,lambda:temp.get_all_sascif())
			all_lines=temp.get_all_sascif(sasbdb='SASDC29')
			self.assertEqual(6249,len(all_lines))

			intensities=temp.get_intensities()
			self.assertEqual(['Q', 'I', 'E'],list(intensities['SASDC29'].columns))

			mod_intensities=temp.modify_intensity()
			cols=['Q', 'I', 'E', 'err_x', 'err_y', 'logI', 'logQ', 'logX', 'Ky', 'Kx', 'Px', 'Py']
			self.assertEqual(cols,list(mod_intensities['SASDC29'].columns))
			self.assertEqual(2e-05,round(mod_intensities['SASDC29']['Py'][0],5))

	def test_get_SASBDB_code_negative(self):
		with tempfile.TemporaryDirectory() as tempdir:
			tmpfilepath = os.path.join(tempdir, 'someFileInTmpDir.txt')
			with open(tmpfilepath, 'w') as tmpfile:
				tmpfile.write("""

data_PDBDEV_test
_entry.id PDBDEV_test
_struct.title 'Structure irrelevant'
loop_
_ihm_dataset_list.id
_ihm_dataset_list.data_type
_ihm_dataset_list.database_hosted
_ihm_dataset_list.details
1 'SAS data' YES .
loop_
_ihm_dataset_related_db_reference.id
_ihm_dataset_related_db_reference.dataset_list_id
_ihm_dataset_related_db_reference.db_name
_ihm_dataset_related_db_reference.accession_code
_ihm_dataset_related_db_reference.version
_ihm_dataset_related_db_reference.details
1 1 . . . .
					""")
			temp=sas_validation(tmpfilepath)
			self.assertEqual([None],temp.get_SASBDB_code())


	def test_get_rg_data(self):
		with tempfile.TemporaryDirectory() as tempdir:
			tmpfilepath = os.path.join(tempdir, 'someFileInTmpDir.txt')
			with open(tmpfilepath, 'w') as tmpfile:
				tmpfile.write("""

data_PDBDEV_test
_entry.id PDBDEV_test
_struct.title 'Structure irrelevant'
loop_
_ihm_dataset_list.id
_ihm_dataset_list.data_type
_ihm_dataset_list.database_hosted
_ihm_dataset_list.details
1 'SAS data' YES .
loop_
_ihm_dataset_related_db_reference.id
_ihm_dataset_related_db_reference.dataset_list_id
_ihm_dataset_related_db_reference.db_name
_ihm_dataset_related_db_reference.accession_code
_ihm_dataset_related_db_reference.version
_ihm_dataset_related_db_reference.details
1 1 SASBDB SASDC29 . .
					""")
			temp=sas_validation(tmpfilepath)
			self.assertEqual({'SASDC29': [3.01, 2.93]},temp.get_rg_for_plot())
			self.assertEqual({'SASDC29': (3.01, 0.369)},temp.get_rg_and_io())
			self.assertEqual(['2.93 nm'],temp.get_rg_table_many()['Rg'])
			self.assertEqual([25.13],temp.get_fits_for_plot()['SASDC29'])
			self.assertEqual(['R', 'P', 'E'],list(temp.get_pofr()['SASDC29'].columns))
			#self.assertEqual(0,int(temp.get_pofr()['SASDC29']['E'][0]))


if __name__ == '__main__':
    unittest.main(warnings='ignore')



