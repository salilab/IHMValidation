###################################
# Script to run validation and write
# to HTML and PDF
# ganesans - Salilab - UCSF
# ganesans@salilab.org
###################################
import os
import datetime
import decouple
import json
import argparse
from multiprocessing import Manager
import pdfkit
import jinja2
import pytz
import sys
sys.path.insert(0, "../master/pyext/src/")
from validation.Report import WriteReport
from validation import utility

# from validation.WKhtmlToPdf import  wkhtmltopdf
# import utility

####################################################################################################################
# Add input arguments for supp table
#####################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-p', type=str, default='No',
                    help="Physical principles used in modeling yes/no?")
parser.add_argument('-f', default='PDBDEV_00000001.cif',
                    help="Input mmcif file")
parser.add_argument(
    '-ls', type=list, default=['No location specified'], help="add location of your scripts")
parser.add_argument(
    '-ld', type=list, default=['No location specified'], help="add location of your analysis files")
parser.add_argument(
    '-m', type=list, default=['Method details unspecified'], help="add information on your method")
parser.add_argument('-models', type=str, default='1',
                    help="number of models in an ensemble, if you have multiple ensembles, add comma-separated string")
parser.add_argument('-c', type=str, default='Distance threshold-based clustering used if ensembles are deposited',
                    help="The type of clustering algorithm used to analyze the ensemble")
parser.add_argument('-mp', type=str, default='10 &#8491 (average RMSF of the solution ensemble with respect to the centroid structure)',
                    help="add model precision. Model precision is defined as average RMSF of the solution ensemble with respect to the centroid structure")
parser.add_argument('-sv', type=list, default=['Information related to sampling validation has not been provided'],
                    help="add model precision. Model precision is defined as average RMSF of the solution ensemble with respect to the centroid structure")
parser.add_argument('-v1', type=list, default=['Fit of model to information used to compute it has not been determined'],
                    help="Add information on satisfaction of input data/restraints")
parser.add_argument('-v2', type=list, default=['Fit of model to information not used to compute it has not been determined'],
                    help="Add information on satisfaction of data not used for modeling")
parser.add_argument('-dv', type=list, default=[
                    'Quality of input data has not be assessed'], help="Add information on quality of input data")
parser.add_argument('-res', type=list, default=['Rigid bodies: 1 residue per bead.',
                                                'Flexible regions: N/A'], help="Add information on model quality (molprobity or excluded volume)")

args = parser.parse_args()
if args.p.upper() == 'YES':
    physics = 'Excluded volume and Sequence connectivity.'
else:
    physics = 'Information about physical principles was not provided'
#############################################################################################################################
# Input for Jinja
####################################################################################
config = pdfkit.configuration(wkhtmltopdf=decouple.config('wkhtmltopdf'))
options = {
    'page-size': 'Letter',
    'margin-top': '0.5in',
    'margin-right': '0.5in',
    'margin-bottom': '0.5in',
    'margin-left': '0.5in',
    'enable-javascript': None,
    'javascript-delay': '50000',
    'header-left': '[page] of [topage]',
    'footer-center': 'Full RCSB IM Structure Validation Report',
    'footer-line': '',
    'header-line': '',
    'footer-spacing': '5',
    'header-spacing': '5',
    "enable-local-file-access": "",
}

options_supp = {
    'page-size': 'A4',
    'margin-top': '0.75in',
    'margin-right': '0.75in',
    'margin-bottom': '0.75in',
    'margin-left': '0.75in',
    'enable-javascript': None,
    'javascript-delay': '500',
    'header-left': '[page] of [topage]',
    'footer-center': 'RCSB IM Methods Table',
    'footer-line': '',
    'header-line': '',
    'footer-spacing': '5',
    'header-spacing': '5'
}

template_flask = ["main.html",
                  "data_quality.html",
                  "model_quality.html",
                  "model_composition.html",
                  "formodeling.html"]


d = datetime.datetime.now()
t = pytz.timezone("America/Los_Angeles")
d1 = t.localize(d)
timestamp = d1.strftime("%B %d, %Y --  %I:%M %p")

dirNames = {
    'images': '../static/images',
    'pdf': '../static/pdf',
    'json': '../static/json',
            'supp': '../static/supplementary',
            'template': '../static/htmls'
}

templateLoader = jinja2.FileSystemLoader(searchpath="../templates/")
templateEnv = jinja2.Environment(loader=templateLoader)
template_pdf = "template_pdf.html"
template_file_supp = "supplementary_template.html"


Template_Dict = {}
Template_Dict['date'] = timestamp
#############################################################################################################################
# Jinja scripts
#############################################################################################################################


def createdirs(dirNames: dict):
    for name in list(dirNames.values()):
        try:
            os.mkdir(name)
            print("Directory ", name,  " Created ")
        except FileExistsError:
            print("Directory ", name,  " already exists")


def write_html(mmcif_file, Template_Dict, template_list, dirName):
    for template_file in template_list:
        template = templateEnv.get_template(template_file)
        outputText = template.render(Template_Dict)
        with open(os.path.join(os.path.join(dirName, template_file)), "w") as fh:
            fh.write(outputText)


def write_pdf(mmcif_file: str, Template_Dict: dict, template_file: str, dirName: str, dirName_Output: str):
    template = templateEnv.get_template(template_file)
    outputText = template.render(Template_Dict)
    with open(os.path.join(os.path.join(dirName, utility.get_output_file_temp_html(mmcif_file))), "w") as fh:
        fh.write(outputText)
    pdfkit.from_file(os.path.join(os.path.join(dirName, utility.get_output_file_temp_html(mmcif_file))),
                     os.path.join(os.path.join(dirName_Output,
                                               utility.get_output_file_pdf(mmcif_file))),
                     options=options)
    # os.remove(os.path.join(os.path.join(dirName,utility.get_output_file_temp_html(mmcif_file))))


def write_supplementary_table(mmcif_file, Template_Dict, template_file, dirName, dirName_supp):
    template = templateEnv.get_template(template_file)
    outputText = template.render(Template_Dict)
    with open(os.path.join(os.path.join(dirName, utility.get_supp_file_html(mmcif_file))), "w") as fh:
        fh.write(outputText)
    pdfkit.from_file(os.path.join(os.path.join(dirName, utility.get_supp_file_html(mmcif_file))),
                     os.path.join(os.path.join(
                         dirName_supp, utility.get_supp_file_pdf(mmcif_file))),
                     options=options_supp)


def write_json(mmcif_file, Template_Dict, dirName, dirName_Output):
    j = json.dumps([{'Category': k, 'Itemized_List': v}
                    for k, v in Template_Dict.items()], indent=4)
    with open(os.path.join(os.path.join(dirName, utility.get_output_file_json(mmcif_file))), "w") as fh:
        fh.write(j)
    fh.close()


def convert_html_to_pdf(template_file, pdf_name, dirName, dirName_Output):
    pdfkit.from_file(os.path.join(os.path.join(dirName, template_file)),
                     os.path.join(os.path.join(dirName_Output, pdf_name)),
                     options=options_supp)


############################################################################################################################
# Run script
#################################################

if __name__ == "__main__":
    utility.clean_all()
    createdirs(dirNames)
    manager = Manager()  # create only 1 mgr
    d = manager.dict()  # create only 1 dict
    report = WriteReport(args.f)
    template_dict = report.run_entry_composition(Template_Dict)
    template_dict, clashscore, rama, sidechain, exv_data = report.run_model_quality(
        template_dict)
    template_dict, sas_data, sas_fit = report.run_sas_validation(template_dict)
    cx_fit, template_dict = report.run_cx_validation(template_dict)
    report.run_quality_glance(
        clashscore, rama, sidechain, exv_data, sas_data, sas_fit, cx_fit)
    report.run_sas_validation_plots(template_dict)
    write_pdf(args.f, template_dict, template_pdf,
              dirNames['pdf'], dirNames['pdf'])
    template_dict = report.run_supplementary_table(template_dict,
                                                   location=args.ls,
                                                   physics=physics,
                                                   method_details=args.m,
                                                   sampling_validation=args.sv,
                                                   validation_input=args.v1,
                                                   cross_validation=args.v2,
                                                   Data_quality=args.dv,
                                                   clustering=args.c,
                                                   resolution=args.res)
    write_supplementary_table(
        args.f, template_dict, template_file_supp, dirNames['supp'], dirNames['supp'])
    write_json(args.f, template_dict, dirNames['json'], dirNames['json'])
    write_html(args.f, template_dict, template_flask, dirNames['template'])
    utility.clean_all()
