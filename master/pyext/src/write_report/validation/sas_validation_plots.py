import pandas as pd
import glob
import sys,os,math
import numpy as np
import pandas as pd
import validation
import ihm
import ihm.reader
import re,pickle,requests,json
from bokeh.io import output_file, show, curdoc, export_png, export_svgs
from bokeh.models import ColumnDataSource, LinearAxis, Legend
from bokeh.palettes import GnBu3, RdBu,OrRd3,Blues,YlOrBr, Spectral6, Set1
from bokeh.plotting import figure, output_file, save
from bokeh.models.widgets import Tabs, Panel
from bokeh.layouts import row,column
import multiprocessing as mp

class sas_validation_plots(validation.sas_validation.sas_validation):
    def __init__(self,mmcif_file):
        super().__init__(mmcif_file)
        self.ID=str(validation.get_input_information.get_id(self))
        self.df=validation.sas_validation.sas_validation.get_intensities(self)
 
    def  plot_intensities(self):
        output_file(self.ID+"intensities.html",mode="inline")
        source = ColumnDataSource(self.df)
        p = figure(plot_height=350, plot_width=350, title="Log I(s) vs s with errorbars")
        p.circle(x='Q',y='logI',source=source, color='blue',fill_alpha=0.1,size=5)
        p.multi_line('err_x','err_y',source=source, color='gray',line_width=0.5)
        p.xaxis.major_label_text_font_size="14pt"
        p.yaxis.major_label_text_font_size="14pt"
        p.title.text_font_size='12pt'
        p.title.align="center"
        p.title.vertical_align='top'
        p.xaxis.axis_label = "s nm\u207B\u00B9"
        p.xaxis.axis_label_text_font_size='14pt'
        p.yaxis.axis_label = 'Log I(s)'
        p.yaxis.axis_label_text_font_size='14pt'
        dirname=os.path.dirname(os.path.abspath(__file__))
        filename = os.path.abspath(os.path.join(dirname, '../static/images/'))
        print (filename)
        save(p,filename=filename+'/'+self.ID+"intensities.html")
        p.output_backend="svg"
        export_svgs(p,filename=filename+'/'+self.ID+"intensities.svg")
        export_png(p,filename=filename+'/'+self.ID+"intensities.png")

    def plot_intensities_log(self):
        output_file(self.ID+"intensities_log.html",mode="inline")
        source = ColumnDataSource(self.df)
        p = figure(plot_height=350, plot_width=350, title="Log I(s) vs Log s with errorbars")
        p.circle(x='logQ',y='logI',source=source,color='blue',fill_alpha=0.1,size=5)
        p.multi_line('logX','err_y',source=source, color='gray',line_width=0.5)
        p.xaxis.major_label_text_font_size="14pt"
        p.yaxis.major_label_text_font_size="14pt"
        p.title.text_font_size='12pt'
        p.title.align="center"
        p.title.vertical_align='top' 
        p.xaxis.axis_label = 'Log s nm\u207B\u00B9'
        p.xaxis.axis_label_text_font_size='14pt'
        p.yaxis.axis_label = 'Log I(s)'
        p.yaxis.axis_label_text_font_size='14pt'
        dirname=os.path.dirname(os.path.abspath(__file__))
        filename = os.path.abspath(os.path.join(dirname, '../static/images/'))
        save(p,filename=filename+'/'+self.ID+"intensities_log.html")
        p.output_backend="svg"
        export_svgs(p,filename=filename+'/'+self.ID+"intensities_log.svg")
        export_png(p,filename=filename+'/'+self.ID+"intensities_log.png")

    def plot_kratky(self):
        output_file(self.ID+"Kratky.html",mode="inline")
        source = ColumnDataSource(self.df)
        p = figure(plot_height=350, plot_width=350, title="Kratky plot")
        p.circle(x='Q',y='Ky',source=source,color='blue',fill_alpha=0.1,size=5)
        p.xaxis.major_label_text_font_size="14pt"
        p.yaxis.major_label_text_font_size="14pt"
        p.title.text_font_size='12pt'
        p.title.align="center"
        p.title.vertical_align='top'
        p.xaxis.axis_label = 'Log s nm\u207B\u00B9'
        p.xaxis.axis_label_text_font_size='14pt'
        p.yaxis.axis_label = 's\u00B2 I(s)'
        p.yaxis.axis_label_text_font_size='14pt'
        dirname=os.path.dirname(os.path.abspath(__file__))
        filename = os.path.abspath(os.path.join(dirname, '../static/images/'))
        save(p,filename=filename+'/'+self.ID+"Kratky.html")
        p.output_backend="svg"
        export_svgs(p,filename=filename+'/'+self.ID+"Kratky.svg")
        export_png(p,filename=filename+'/'+self.ID+"Kratky.png")
    
    def plot_porod_debye(self):
        output_file(self.ID+"porod.html",mode="inline")
        source = ColumnDataSource(self.df)
        p = figure(plot_height=350, plot_width=350, title="Porod-Debye plot")
        p.circle(x='Px',y='Py',source=source,color='blue',fill_alpha=0.1,size=5)
        p.xaxis.major_label_text_font_size="14pt"
        p.yaxis.major_label_text_font_size="14pt"
        p.title.text_font_size='12pt'
        p.title.align="center"
        p.title.vertical_align='top'
        p.xaxis.axis_label = 's \u2074'
        p.xaxis.axis_label_text_font_size='14pt'
        p.yaxis.axis_label = 's\u2074 I(s)'
        p.yaxis.axis_label_text_font_size='14pt'
        dirname=os.path.dirname(os.path.abspath(__file__))
        filename = os.path.abspath(os.path.join(dirname, '../static/images/'))
        save(p,filename=filename+'/'+self.ID+"porod.html")
        p.output_backend="svg"
        export_svgs(p,filename=filename+'/'+self.ID+"porod.svg")
        export_png(p,filename=filename+'/'+self.ID+"porod.png")

