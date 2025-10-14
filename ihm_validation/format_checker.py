#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# format_checker.py - Check residue and atom names in IHMCIF file
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
Check residue and atom names in IHMCIF file
"""

import sys
import logging
import argparse as ag
import ihm, ihm.reader, ihm.util.make_mmcif

# Non-standard histidine names (protonation states)
HISTIDINES = frozenset(('HIP', 'HID', 'HIE'))

def parse_ihm_cif(fname, encoding='utf8') -> tuple:
    try:
        with open(fname, encoding=encoding) as fh:
            system, = ihm.reader.read(fh)
    except UnicodeDecodeError:
        encoding = 'ascii'
        with open(fname, encoding=encoding, errors='ignore') as fh:
            system, = ihm.reader.read(fh)

    return(system, encoding)

def check_entities_histidines(system: ihm.System, histidines=HISTIDINES):
    """Find any non-standard histidine chemical components"""
    out = []
    his = ihm.LPeptideAlphabet()['H']
    for e in system.entities:
        for c in e.sequence:
            if c.id in histidines:
                out.append(c.id)

    if len(out) > 0:
        raise(ValueError(f"Non-canonical histidine variant found: {', '.join(set(out))}"))

def check_models(system: ihm.System):
    """Find any non-standard histidine chemical components"""
    for state_group in system.state_groups:
        for state in state_group:
            for model_group in state:
                for model in model_group:
                    ihm.util.make_mmcif._check_atom_names(model, check_all=True)

def check_all_exception(system: ihm.System):
    """Perform all checks. Throw an exception if a check fails."""
    # Disable atom check until python-ihm fixes
    # checks = [check_entities_histidines, check_models]
    checks = [check_entities_histidines]

    for check in checks:
        check(system)

def check_all_log(system: ihm.System) -> int:
    """Perform all checks. Throw a message in the log if a check fails and return a non-zero exit code"""
    # Disable atom check until python-ihm fixes
    # checks = [check_entities_histidines, check_models]
    checks = [check_entities_histidines]
    exit_code = 0
    for check in checks:
        try:
            check(system)
        except ValueError as e:
            logging.error(e)
            exit_code = 127

    return exit_code

def check_file_exception(fname: str):
    """Parse a file, do all checks, throw an exception if a check fails."""
    system, encoding = parse_ihm_cif(fname)
    check_all_exception(system)

def check_file_log(fname: str) -> int:
    """Parse a file, do all checks, throw a log message if a check fails and return a non-zero exit code"""
    system, encoding = parse_ihm_cif(fname)
    exit_code = check_all_log(system)
    return exit_code

if __name__ == "__main__":
    parser = ag.ArgumentParser(description="Check residue and atom names in IHMCIF file")
    parser.add_argument("-i", "--input_file", help="Path to the input file")

    args = parser.parse_args()
    check_file_exception(args.input_file)
