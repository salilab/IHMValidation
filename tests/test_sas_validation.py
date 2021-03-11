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
			print (temp.get_SASBDB_code())

if __name__ == '__main__':
    unittest.main(warnings='ignore')



