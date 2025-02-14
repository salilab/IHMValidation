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
import mendeleev

class GetExcludedVolume(GetInputInformation):
    def __init__(self, mmcif_file: str, cache: str='.', nocache: bool=False):
        super().__init__(mmcif_file, cache=cache, nocache=nocache)
        self.cache = cache

    def get_all_spheres(self, fname: str | None = None) -> dict:
        """ Get spheres and atoms for each model """

        if fname is not None:
            system, encoding = utility.parse_ihm_cif(fname)
        else:
            system = self.system

        model_dict = {}

        for st in system.state_groups:
            for s in st:
                for mg in s:
                    for m in mg:
                        model_dict[m._id] = list(m.get_spheres()) + list(m.get_atoms())

        return model_dict

    def get_nCr(self, n: int, r: int) -> int:
        """get all combinations"""
        f = math.factorial
        return int(f(n)/(f(r)*f(n-r)))

    def get_violation_percentage(self, viols: dict, n: int) -> float:
        """get information on all spheres for each model"""
        number_of_violations = sum(viols.values())
        number_of_combinations = self.get_nCr(n, 2)
        return (1 - number_of_violations / number_of_combinations ) * 100.

    @staticmethod
    def get_radii(atoms: list) -> dict:
        """ Get VdW radii for all atom types """
        elems = []
        radii = {}

        for atom in atoms:
            if isinstance(atom, ihm.model.Atom):
                elem_ = atom.type_symbol
                elems.append(elem_)

        elems = set(elems)

        for elem in elems:
            if not isinstance(elem, str):
                logging.warning(f"Skipping missing element")
                continue

            # Try to guess elem:
            elem_ = elem

            if len(elem) == 1:
                elem_ = elem.upper()
            else:
                elem_ = elem[0].upper() + elem[1].lower()

            try:
                vdwr = getattr(mendeleev, elem_).vdw_radius
                # Convert picometers to angstroms
                vdwr /= 100.0
                radii[elem] = vdwr
            except AttributeError as e:
                logging.warning(f"Skipping unknown element {elem}")

        return radii


    def get_xyzr(self, spheres: list) -> pd.DataFrame:
        """ get X,Y, Z coords from sphere objects"""
        raw = []

        vdw = self.get_radii(spheres)

        for sphere in spheres:
            if isinstance(sphere, ihm.model.Sphere):
                data_ = (sphere.x, sphere.y, sphere.z, sphere.radius)
                raw.append(data_)
            elif isinstance(sphere, ihm.model.Atom):
                try:
                    radius = vdw[sphere.type_symbol]
                except KeyError as e:
                    # We've already put a message in the log
                    # so just silently skip
                    continue

                data_ = (sphere.x, sphere.y, sphere.z, radius)
                raw.append(data_)

        model_spheres_df = pd.DataFrame(raw, columns=['X', 'Y', 'Z', 'R'])

        return model_spheres_df

    def get_violation_dict(self, model_spheres_df: pd.DataFrame) -> dict:
        """ get violation from model_sphere df"""
        viols = {}

        # Get coordinates
        xyz = model_spheres_df[['X', 'Y', 'Z']].to_numpy()
        # Get radii
        radii = model_spheres_df[['R']].to_numpy()
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

    def _get_exc_vol(self, sphere_list: list) -> (float, int):
        """
        get violations from cart coords
        """
        total = self.get_nCr(len(sphere_list), 2)
        df = self.get_xyzr(sphere_list)
        violation_dict = self.get_violation_dict(df)
        violations = sum(violation_dict.values())
        satisfaction = round(
            self.get_violation_percentage(violation_dict, len(df)), 2)
        return (total, violations, satisfaction)

    def _run_exc_vol(self, model_dict: dict) -> dict:
        """ get exc vol info in parallel """
        complete_list = [
            self._get_exc_vol(x) for x in model_dict.values()]

        excluded_volume = {'Model ID': list(model_dict.keys()),
                           'Analysed': [i[0] for i in complete_list],
                           'Number of violations': [i[1] for i in complete_list],
                           'Excluded Volume Satisfaction (%)': [i[2] for i in complete_list],
                           }

        return excluded_volume

    def get_excluded_volume(self):
        cache_fn = Path(self.cache, f'{self.stem}.exv.pkl')
        data = None

        # Check if we already requested the data
        if Path(cache_fn).is_file() and not self.nocache:
            logging.info(f'Found {self.stem} in cache: {cache_fn}')
            with open(cache_fn, 'rb') as f:
                data = pickle.load(f)

        elif not Path(cache_fn).is_file() or self.nocache:
            logging.info("Excluded volume is being calculated...")
            model_dict = self.get_all_spheres()
            data = self._run_exc_vol(model_dict)

            with open(cache_fn, 'wb') as f:
                pickle.dump(data, f)

        return data
