from flask import json,Flask, flash, render_template, request, url_for, redirect
import jinja2
import os,sys
import pandas as pd
pd.set_option('mode.chained_assignment', None)
import numpy as np


app = Flask(__name__)
@app.context_processor
def override_url_for():
    return dict(url_for=dated_url_for)

def dated_url_for(endpoint, **values):
    if endpoint == 'static':
        filename = values.get('filename', None)
        if filename:
            file_path = os.path.join(app.root_path,
                                 endpoint, filename)
            values['q'] = int(os.stat(file_path).st_mtime)
    return url_for(endpoint, **values)

@app.route("/")
def home():
    return render_template("main.html")

@app.route("/main.html")
def main():
    return render_template("main.html")

@app.route("/model_composition.html")
def model_composition():
    return render_template("model_composition.html")

@app.route("/data_quality.html")
def data_quality():
    return render_template("data_quality.html")

@app.route("/model_quality.html")
def model_quality():
    return render_template("model_quality.html")

@app.route("/formodeling.html")
def for_modeling():
    return render_template("formodeling.html")

@app.route("/notmodeling.html")
def not_modeling():
    return render_template("notmodeling.html")

@app.route("/uncertainty.html")
def uncertainty():
    return render_template("uncertainty.html")


if __name__ == "__main__":
    app.run(debug=True)


