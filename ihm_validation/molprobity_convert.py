#!/bin/env python3
# -*- coding: utf-8 -*-
#
# molprobity_convert.py - Convert MolProbity output from json to a pickled object
#
# Copyright (C) 2025 Arthur Zalevsky <aozalevsky@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
Convert MolProbity output from json to a pickled object
"""

import argparse
import pickle
import json
from pathlib import Path

def get_atom_tuple(atom):
    chain = atom.chain_id.strip()
    resid = atom.resid.strip()
    resname = atom.resname.strip()
    name = atom.name.strip()

    r_ = (chain, resid, resname, name)

    return r_

def rota_rama_result_to_tuple(r):
    chid = r.chain_id.strip()
    resid = r.resid.strip()
    resname = r.resname.strip()
    score = r.score

    r_ = (chid, resid, resname, score)

    return r_


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input',
                    required=True)

    parser.add_argument('-o', '--output',
                    required=True)

    args = parser.parse_args()

    with open(args.input, 'rb') as f:
        data = pickle.load(f)

    clashscore = data.clashes.clashscore

    clashes_list = []
    for r in data.clashes.results:
        a1, a2 = r.atoms_info
        a1_ = get_atom_tuple(a1)
        a2_ = get_atom_tuple(a2)
        overlap = abs(r.overlap)

        r_ = (a1_, a2_, overlap)
        clashes_list.append(r_)

    bonds_total = data.restraints.bonds.n_total
    bonds_outliers = data.restraints.bonds.n_outliers

    bonds_outliers_list = []
    for r in data.restraints.bonds.results:
        a1, a2 = r.atoms_info
        a1_ = get_atom_tuple(a1)
        a2_ = get_atom_tuple(a2)
        ideal = r.target
        observed = r.model
        score = r.score

        r_ = (a1_, a2_, observed, ideal, score)
        bonds_outliers_list.append(r_)

    angles_total = data.restraints.angles.n_total
    angles_outliers = data.restraints.angles.n_outliers

    angles_outliers_list = []
    for r in data.restraints.angles.results:
        a1, a2, a3 = r.atoms_info
        a1_ = get_atom_tuple(a1)
        a2_ = get_atom_tuple(a2)
        a3_ = get_atom_tuple(a3)
        ideal = r.target
        observed = r.model
        score = r.score

        r_ = (a1_, a2_, a3_, observed, ideal, score)
        angles_outliers_list.append(r_)

    rama_total = data.ramalyze.n_total
    rama_favored = data.ramalyze.n_favored
    rama_allowed = data.ramalyze.n_allowed
    rama_outliers = data.ramalyze.n_outliers

    rama_list = []

    for r in data.ramalyze.results:
        if r.is_outlier():
            r_ = rota_rama_result_to_tuple(r)
            rama_list.append(r_)

    rota_total = data.rotalyze.n_total
    rota_favored = data.rotalyze.n_favored
    rota_allowed = data.rotalyze.n_allowed
    rota_outliers = data.rotalyze.n_outliers

    rota_list = []

    for r in data.rotalyze.results:
        if r.is_outlier():
            r_ = rota_rama_result_to_tuple(r)
            rota_list.append(r_)

    out = {
        'clash': {
            'clashscore': clashscore,
            'clashes_list': clashes_list,
        },

        'bonds': {
            'total': bonds_total,
            'outliers': bonds_outliers,
            'outliers_list': bonds_outliers_list,
        },

        'angles': {
            'total': angles_total,
            'outliers': angles_outliers,
            'outliers_list': angles_outliers_list,
        },

        'rama': {
            'total': rama_total,
            'favored': rama_favored,
            'allowed': rama_allowed,
            'outliers': rama_outliers,
            'outliers_list': rama_list,
        },

        'rota': {
            'total': rota_total,
            'favored': rota_favored,
            'allowed': rota_allowed,
            'outliers': rota_outliers,
            'outliers_list': rota_list,
        }
    }

    with open(args.output, 'wb') as f:
        pickle.dump(out, f)
