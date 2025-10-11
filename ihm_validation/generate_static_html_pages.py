#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# generate_static_html_pages.py - A script for generation of validation_help.html and about_validation.html
#
# Copyright (C) 2019-2025 Arthur Zalevsky, Sai Ganesan, Benjamin M. Webb, Brinda Vallat
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
A script for generation of validation_help.html and about_validation.html
"""

import logging
import argparse
from pathlib import Path
import jinja2
import utility
from report import REPORT_VERSION

#############################################################################################################################
# Jinja scripts
#############################################################################################################################

def createdirs(dirNames: dict):
    for name in list(dirNames.values()):
        if Path(name).is_dir():
            logging.info(f"Directory {name} already exists")
        else:
            Path(name).mkdir(parents=True)
            logging.info(f"Directory {name} created ")

def write_html(template_dict: dict, template_list: list, dirName: str, overwrite=False):
    """Generate HTML pages from templates"""

    for template_file in template_list:
        template = templateEnv.get_template(template_file)
        outputText = template.render(template_dict)

        fn = Path(dirName, template_file)

        if fn.is_file() and not overwrite:
            logging.info(f'Output directory {output_path} exists. '
                         'Use --force to overwright')
        else:
            with open(fn, "w") as fh:
                fh.write(outputText)

############################################################################################################################
# Run script
#################################################

if __name__ == "__main__":
    ###################################################################################################################
    # Parser
    #####################################################################

    parser = argparse.ArgumentParser()
    parser.add_argument('-v', dest='verbose', action='store_true',
                        help="Verbose output")
    parser.add_argument('--about-validation', dest='about_validation', action='store_true',
                        help="Generate about validation page")
    parser.add_argument('--validation-help', dest='validation_help', action='store_true',
                        help="Generate validation help page")
    parser.add_argument('--output-root', type=str, default=str(Path(Path(__file__).parent.resolve(), 'Validation')),
                        help="Path to a directory where the output will be written")
    # parser.add_argument('--html-mode', type=str, default='pdb-ihm',
    #                     choices=['local', 'pdb-ihm'],
    #                     help="HTML mode affects paths to various statis resources")
    # parser.add_argument('--html-resources',
    #                     type=str,
    #                     default=str(Path(Path(__file__).parent.parent.resolve(), 'static')),
    #                     help="Path to static HTML resources")
    parser.add_argument('--force', action='store_true', default=False,
                        help="Overwright output files")

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING)

    #############################################################################################################################
    # Input for Jinja
    ####################################################################################

    templates = [
    #    "about_validation.html",
    #    "validation_help.html",
    ]

    if args.about_validation:
        templates.append('about_validation.html')

    if args.validation_help:
        templates.append('validation_help.html')

    output_path = args.output_root

    dirNames = {
        'root': str(output_path),
    }

    # This is a temporary hack for ../templates
    template_path = Path(Path(__file__).parent.parent.resolve(), 'templates')
    templateLoader = jinja2.FileSystemLoader(searchpath=template_path)
    templateEnv = jinja2.Environment(loader=templateLoader)
    Template_Dict = {}
    Template_Dict['version'] = REPORT_VERSION
    Template_Dict['html_mode'] = 'pdb-ihm'

    logging.info("Clean up and create output directories")

    createdirs(dirNames)

    logging.info("Generating pages")
    write_html(Template_Dict, templates, dirNames['root'], args.force)
