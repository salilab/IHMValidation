#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# format_checker.py - Check residue and atom names in IHMCIF file
#
# Copyright (C) 2025 Arthur Zalevsky <aozalevsky@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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
