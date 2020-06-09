import os,sys
import unittest

sys.path.append("../../")
from FlaskApp import app
  
 
class BasicTests(unittest.TestCase):
 
    ############################
    #### setup and teardown ####
    ############################
 
    # executed prior to each test
    def setUp(self):
        app.config['TESTING'] = True
        app.config['DEBUG'] = False
        self.app = app.test_client()

    # executed after each test
    def tearDown(self):
        pass
 
###############
#### tests ####
###############
 
    def test_main_page(self):
        response = self.app.get('/', follow_redirects=True)
        self.assertEqual(response.status_code, 200)
 
    def test_home_page(self):
        response = self.app.get('/home/', follow_redirects=True)
        self.assertEqual(response.status_code, 200)

    def test_intro_page(self):
        response = self.app.get('/introduction.html', follow_redirects=True)
        self.assertEqual(response.status_code, 200)

    def test_homeintro_page(self):
        response = self.app.get('/home/introduction.html', follow_redirects=True)
        self.assertEqual(response.status_code, 200)

    def test_definitions_page(self):
        response = self.app.get('/definitions.html', follow_redirects=True)
        self.assertEqual(response.status_code, 200)

    def test_homedefinitions_page(self):
        response = self.app.get('/home/definitions.html', follow_redirects=True)
        self.assertEqual(response.status_code, 200)

    def test_main_all(self):
        for i in [1,2,4,9]:
            name='PDBDEV_0000000'+str(i)
            response = self.app.get('/'+name+'/main.html', follow_redirects=True)
            self.assertEqual(response.status_code, 200)
        for i in [17,27]:   
            name='PDBDEV_000000'+str(i)
            response = self.app.get('/'+name+'/main.html', follow_redirects=True)
            self.assertEqual(response.status_code, 200)

    def test_model_composition_all(self):
        for i in [1,2,4,9]:
            name='PDBDEV_0000000'+str(i)
            response = self.app.get('/'+name+'/model_composition.html', follow_redirects=True)
            self.assertEqual(response.status_code, 200)
        for i in [17,27]:   
            name='PDBDEV_000000'+str(i)
            response = self.app.get('/'+name+'/model_composition.html', follow_redirects=True)
            self.assertEqual(response.status_code, 200)

    def test_data_quality_all(self):
        for i in [1,2,4,9]:
            name='PDBDEV_0000000'+str(i)
            response = self.app.get('/'+name+'/data_quality.html', follow_redirects=True)
            self.assertEqual(response.status_code, 200)
        for i in [17,27]:   
            name='PDBDEV_000000'+str(i)
            response = self.app.get('/'+name+'/data_quality.html', follow_redirects=True)
            self.assertEqual(response.status_code, 200)

    def test_model_quality_all(self):
        for i in [1,2,4,9]:
            name='PDBDEV_0000000'+str(i)
            response = self.app.get('/'+name+'/model_quality.html', follow_redirects=True)
            self.assertEqual(response.status_code, 200)
        for i in [17,27]:   
            name='PDBDEV_000000'+str(i)
            response = self.app.get('/'+name+'/model_quality.html', follow_redirects=True)
            self.assertEqual(response.status_code, 200)

    def test_formodeling_all(self):
        for i in [1,2,4,9]:
            name='PDBDEV_0000000'+str(i)
            response = self.app.get('/'+name+'/formodeling.html', follow_redirects=True)
            self.assertEqual(response.status_code, 200)
        for i in [17,27]:   
            name='PDBDEV_000000'+str(i)
            response = self.app.get('/'+name+'/formodeling.html', follow_redirects=True)
            self.assertEqual(response.status_code, 200)

    def test_notformodeling_all(self):
        for i in [1,2,4,9]:
            name='PDBDEV_0000000'+str(i)
            response = self.app.get('/'+name+'/notformodeling.html', follow_redirects=True)
            self.assertEqual(response.status_code, 200)
        for i in [17,27]:   
            name='PDBDEV_000000'+str(i)
            response = self.app.get('/'+name+'/notformodeling.html', follow_redirects=True)
            self.assertEqual(response.status_code, 200)

    def test_uncertainty_all(self):
        for i in [1,2,4,9]:
            name='PDBDEV_0000000'+str(i)
            response = self.app.get('/'+name+'/main.html', follow_redirects=True)
            self.assertEqual(response.status_code, 200)
        for i in [17,27]:   
            name='PDBDEV_000000'+str(i)
            response = self.app.get('/'+name+'/main.html', follow_redirects=True)
            self.assertEqual(response.status_code, 200)

    def test_download_all(self):
        for i in [1,2,4,9]:
            name='PDBDEV_0000000'+str(i)
            response = self.app.get('/'+name+'/download', follow_redirects=True)
            self.assertEqual(response.status_code, 200)
            response.close()
        for i in [17,27]:   
            name='PDBDEV_000000'+str(i)
            response = self.app.get('/'+name+'/download', follow_redirects=True)
            self.assertEqual(response.status_code, 200)
            response.close()

    def test_download_table(self):
        for i in [1,2,4,9]:
            name='PDBDEV_0000000'+str(i)
            response = self.app.get('/'+name+'/downloadTable', follow_redirects=True)
            self.assertEqual(response.status_code, 200)
            response.close()
        for i in [17,27]:   
            name='PDBDEV_000000'+str(i)
            response = self.app.get('/'+name+'/downloadTable', follow_redirects=True)
            self.assertEqual(response.status_code, 200)
            response.close()


if __name__ == "__main__":
    unittest.main()
