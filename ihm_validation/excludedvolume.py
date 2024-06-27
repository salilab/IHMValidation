###################################
# Script :
# 1) Contains class for excluded
# volume calculation
#
# ganesans - Salilab - UCSF
# ganesans@salilab.org
###################################
from pathlib import Path
from mmcif_io import GetInputInformation
import ihm
import multiprocessing as mp
import pandas as pd
import numpy as np
from scipy.spatial import KDTree
import math
import pickle
import logging

class GetExcludedVolume(GetInputInformation):
    ID = None

    def __init__(self, mmcif_file, cache):
        super().__init__(mmcif_file)
        self.ID = str(self.get_id())
        self.cache = cache

    def get_all_spheres(self, filetemp=None):
        """get information on all spheres for each model"""
        if filetemp is None:
            model_object = [
                b for i in self.system.state_groups for j in i for a in j for b in a]
            model_dict = {i+1: j._spheres for i, j in enumerate(model_object)}
        else:
            system, = ihm.reader.read(filetemp,
                                      model_class=ihm.model.Model)
            model_object = [
                b for i in system.state_groups for j in i for a in j for b in a]
            model_dict = {i+1: j._spheres for i, j in enumerate(model_object)}
        return model_dict

    def get_nCr(self, n, r):
        """get all combinations"""
        f = math.factorial
        return int(f(n)/(f(r)*f(n-r)))

    def get_violation_percentage(self, models_spheres_df: pd.DataFrame, viols: dict) -> float:
        """get information on all spheres for each model"""
        number_of_violations = sum(viols.values())
        number_of_combinations = self.get_nCr(models_spheres_df.shape[1], 2)
        return (1-number_of_violations/number_of_combinations)*100

    def get_violation_normalized(self, models_spheres_df: pd.DataFrame, viols: dict) -> float:
        """
        normalize violations, not currently used.
        """
        number_of_violations = sum(list(viols.values()))
        normalization_constant = models_spheres_df.shape[1]*math.log(
            models_spheres_df.shape[1], 10)
        return (1-number_of_violations/normalization_constant)*100

    def get_xyzr(self, spheres: pd.DataFrame) -> pd.DataFrame:
        """ get X,Y, Z coords from sphere objects"""
        # model_spheres={i+1:[j.x,j.y,j.z,j.radius] for i,j in enumerate(spheres)}
        # model_spheres_df=pd.DataFrame(model_spheres, index=['X','Y','Z','R'])
        model_spheres_df = pd.DataFrame.from_records([(j.x, j.y, j.z, j.radius)
                                                      for i, j in enumerate(spheres)], columns=['X', 'Y', 'Z', 'R'])
        model_spheres_df.index += 1
        return model_spheres_df.T

    def get_xyzr_complete(self, model_ID, spheres: list) -> pd.DataFrame:
        """ get X,Y,Z,R, chain and model ID from sphere objects"""

        # model_spheres={i+1:[j.x,j.y,j.z,j.radius,j.asym_unit._id,model_ID] for i,j in enumerate(spheres)}
        # model_spheres_df=pd.DataFrame(model_spheres, index=['X','Y','Z','R','Chain_ID','Model_ID'])
        model_spheres_df = pd.DataFrame.from_records([(j.x, j.y, j.z, j.radius, j.asym_unit._id, model_ID)
                                                      for i, j in enumerate(spheres)], columns=['X', 'Y', 'Z', 'R',
                                                                                                'Chain_ID', 'Model_ID'])
        model_spheres_df.index += 1
        return model_spheres_df

    def get_violation_dict(self, model_spheres_df: pd.DataFrame) -> dict:
        """ get violation from model_sphere df"""
        viols = {}

        # Get coordinates
        xyz = model_spheres_df.T[['X', 'Y', 'Z']].to_numpy()
        # Get radii
        radii = model_spheres_df.T[['R']].to_numpy()
        # Get maximum radius
        maxr = np.max(radii)
        # Build tree
        t = KDTree(xyz)

        # The enumeration is done to preveserve
        # compatibility as it's a drop-in replacement
        for indx, i in enumerate(range(len(xyz) - 1), 1):
            viols_ = 0
            # Get neighours in R1+R2 radius, where
            # R1 is particle's i radius
            # and R2 is the maxium radius
            # Thus it's a greedy search
            nb = t.query_ball_point(xyz[i], radii[i] + maxr)

            # Check each neighbour
            for j in nb:
                # Only check pairs in a triangle
                if j > i:
                    # np.linalg.norm is slow, but convenient
                    d = np.linalg.norm(xyz[i] - xyz[j])
                    if d < (radii[i] + radii[j]):
                        viols_ += 1

            viols[indx] = viols_

        return viols

    def get_exc_vol_given_sphere_parallel(self, sphere_list: list) -> (float, int):
        """
        get violations from cart coords
        """
        total = self.get_nCr(len(sphere_list), 2)
        df = self.get_xyzr(sphere_list)
        violation_dict = self.get_violation_dict(df)
        satisfaction = round(
            self.get_violation_percentage(df, violation_dict), 2)
        violations = sum(violation_dict.values())
        return (total, violations, satisfaction)

    def run_exc_vol_parallel(self, model_dict: dict) -> dict:
        """ get exc vol info in parallel """
        pool = mp.Pool(processes=min(4, len(model_dict.keys())))
        complete_list = pool.map(
            self.get_exc_vol_given_sphere_parallel, list(model_dict.values()))
        excluded_volume = {'Models': list(model_dict.keys()),
                           'Analysed': [i[0] for i in complete_list],
                           'Number of violations': [i[1] for i in complete_list],
                           'Excluded Volume Satisfaction (%)': [i[2] for i in complete_list],
                           }

        return excluded_volume

    def get_excluded_volume(self):
        cache_fn = Path(self.cache, f'{self.ID}_exv.pickle')
        data = None

        # Check if we already requested the data
        if Path(cache_fn).is_file():
            logging.info(f'Found {self.ID} in cache: {cache_fn}')
            with open(cache_fn, 'rb') as f:
                data = pickle.load(f)

        elif not Path(cache_fn).is_file():
            spheres = self.get_all_spheres()
            model_dict = self.get_all_spheres()
            data = self.run_exc_vol_parallel(model_dict)

            with open(cache_fn, 'wb') as f:
                pickle.dump(data, f)

        return data
