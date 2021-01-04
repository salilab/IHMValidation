from flask import json,Flask, flash, render_template, request, url_for, redirect,send_file
import jinja2
import os,sys
import pandas as pd
pd.set_option('mode.chained_assignment', None)
import numpy as np
import json


app = Flask(__name__)
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0

def get_dictionary(file_name):
	with open(file_name) as handle:
		dictdump = json.loads(handle.read())
		data=dict()
		for item in dictdump:
			for key,val in item.items():
				data[key]=val
	return data


@app.route("/")
@app.route("/home/")
@app.route("/introduction.html")
@app.route("/home/introduction.html")
def home():
	return render_template("introduction.html")

@app.route('/test')
def test():
	#data={'id':1, 'name':'Josh'}
	data=get_dictionary(file_name='static/json/ValidationReport_PDBDEV_00000001.cif.json')
	#print (data)
	return render_template("index.html", data = data)

@app.route("/<id>/main.html")
def main(id):
	name='static/json/ValidationReport_'+id+''+'.cif.json'
	data=get_dictionary(file_name=name)
	return render_template("main_json.html", data = data)

@app.route("/<id>/model_composition.html")
def model_composition(id):
	name='static/json/ValidationReport_'+id+''+'.cif.json'
	data=get_dictionary(file_name=name)
	return render_template('model_composition_json.html',data=data)

@app.route("/<id>/data_quality.html")
def data_quality(id):
	name='static/json/ValidationReport_'+id+''+'.cif.json'
	data=get_dictionary(file_name=name)
	return render_template('data_quality_json.html',data=data)

@app.route("/<id>/model_quality.html")
def model_quality(id):
	name='static/json/ValidationReport_'+id+''+'.cif.json'
	data=get_dictionary(file_name=name)
	return render_template('model_quality_json.html',data=data)

@app.route("/<id>/formodeling.html")
def formodeling(id):
	name='static/json/ValidationReport_'+id+''+'.cif.json'
	data=get_dictionary(file_name=name)
	return render_template('formodeling_json.html',data=data,enumerate=enumerate)

@app.route("/notformodeling.html")
def notformodeling():
	return render_template('notformodeling.html')

@app.route("/uncertainty.html")
def uncertainty():
	return render_template('uncertainty.html')

@app.route('/<id>/download')
def downloadFile (id):
	print (id)
	path='static/pdf/'+id+'.cif.pdf'
	return send_file(path, as_attachment=True)

@app.route('/<id>/downloadTable')
def downloadTable (id):
	print (id)
	path='static/pdf/Supplementary_'+id+'.cif.pdf'
	return send_file(path, as_attachment=True)



if __name__ == "__main__":
	app.run(debug=True)
	#d=get_dictionary(file_name='static/json/ValidationReport_PDBDEV_00000001.cif.json')
	#print (d)


