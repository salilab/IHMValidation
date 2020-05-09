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
from bokeh.models import Span,ColumnDataSource, LinearAxis, Legend
from bokeh.palettes import GnBu3, RdBu,OrRd3,Blues,YlOrBr, Spectral6, Set1
from bokeh.plotting import figure, output_file, save
from bokeh.models.widgets import Tabs, Panel
from bokeh.layouts import row,column
import multiprocessing as mp

class sas_validation_plots(validation.sas_validation.sas_validation):
    def __init__(self,mmcif_file):
        super().__init__(mmcif_file)
        self.ID=str(validation.get_input_information.get_id(self))
        self.df_dict=validation.sas_validation.sas_validation.modify_intensity(self)
        self.pdf_dict=validation.sas_validation.sas_validation.get_pddf(self)
        self.fdf_dict=validation.sas_validation.sas_validation.get_fit_data(self)
        self.score,self.gdf=validation.sas_validation.sas_validation.get_Guinier_data(self)

    def plot_intensities(self,sasbdb,df):
        output_file(self.ID+sasbdb+"intensities.html",mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=500, plot_width=500, title="Log I(q) vs q with errorbars ("+sasbdb+")")
        p.circle(x='Q',y='logI',source=source, color='blue',fill_alpha=0.3,size=5)
        p.multi_line('err_x','err_y',source=source, color='gray',line_width=0.5)
        p.xaxis.major_label_text_font_size="14pt"
        p.yaxis.major_label_text_font_size="14pt"
        p.title.text_font_size='12pt'
        p.title.align="center"
        p.title.vertical_align='top'
        p.xaxis.axis_label = "q nm\u207B\u00B9"
        p.xaxis.axis_label_text_font_size='14pt'
        p.yaxis.axis_label = 'Log I(q)'
        p.yaxis.axis_label_text_font_size='14pt'
        dirname=os.path.dirname(os.path.abspath(__file__))
        filename = os.path.abspath(os.path.join(dirname, '../static/images/'))
        print (filename)
        save(p,filename=filename+'/'+self.ID+sasbdb+"intensities.html")
        p.output_backend="svg"
        export_svgs(p,filename=filename+'/'+self.ID+sasbdb+"intensities.svg")
        export_png(p,filename=filename+'/'+self.ID+sasbdb+"intensities.png")

    def plot_multiple(self):
        for key,val in self.df_dict.items():
            self.plot_intensities(key,val)
            self.plot_intensities_log(key,val)
            self.plot_kratky(key,val)
            self.plot_porod_debye(key,val)
 
    def plot_Guinier(self):
        for key,val in self.gdf.items():
            self.Guinier_plot_fit(key,val,self.score[key])
            self.Guinier_plot_residuals(key,val)

    def plot_pf(self):
        for key,val in self.pdf_dict.items():
            self.plot_pddf(key,val)

    def plot_fits(self):
        for key,val in self.fdf_dict.items():
            for key_m,val_m in val.items():
                if val_m[1].empty==False:
                    self.plot_fit(key,key_m,val_m[0],val_m[1])
                    self.plot_fit_residuals(key,key_m,val_m[1])


    def plot_intensities_log(self,sasbdb,df):
        output_file(self.ID+sasbdb+"intensities_log.html",mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=500, plot_width=500, title="Log I(q) vs Log q with errorbars ("+sasbdb+")")
        p.circle(x='logQ',y='logI',source=source,color='blue',fill_alpha=0.3,size=5)
        p.multi_line('logX','err_y',source=source, color='gray',line_width=0.5)
        p.xaxis.major_label_text_font_size="14pt"
        p.yaxis.major_label_text_font_size="14pt"
        p.title.text_font_size='12pt'
        p.title.align="center"
        p.title.vertical_align='top' 
        p.xaxis.axis_label = 'Log q nm\u207B\u00B9'
        p.xaxis.axis_label_text_font_size='14pt'
        p.yaxis.axis_label = 'Log I(q)'
        p.yaxis.axis_label_text_font_size='14pt'
        dirname=os.path.dirname(os.path.abspath(__file__))
        filename = os.path.abspath(os.path.join(dirname, '../static/images/'))
        save(p,filename=filename+'/'+self.ID+sasbdb+"intensities_log.html")
        p.output_backend="svg"
        export_svgs(p,filename=filename+'/'+self.ID+sasbdb+"intensities_log.svg")
        export_png(p,filename=filename+'/'+self.ID+sasbdb+"intensities_log.png")

    def plot_kratky(self,sasbdb,df):
        output_file(self.ID+sasbdb+"Kratky.html",mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=500, plot_width=500, title="Kratky plot ("+sasbdb+")")
        p.circle(x='Q',y='Ky',source=source,color='blue',fill_alpha=0.3,size=5)
        #vline = Span(location=0.1732, dimension='height', line_color='red', line_width=3)
        #hline = Span(location=0.1104, dimension='width', line_color='green', line_width=3)
        #p.renderers.extend([vline, hline])
        p.xaxis.major_label_text_font_size="14pt"
        p.yaxis.major_label_text_font_size="14pt"
        p.title.text_font_size='12pt'
        p.title.align="center"
        p.title.vertical_align='top'
        p.xaxis.axis_label = 'Log q nm\u207B\u00B9'
        p.xaxis.axis_label_text_font_size='14pt'
        p.yaxis.axis_label = 'q\u00B2 I(q)'
        p.yaxis.axis_label_text_font_size='14pt'
        dirname=os.path.dirname(os.path.abspath(__file__))
        filename = os.path.abspath(os.path.join(dirname, '../static/images/'))
        save(p,filename=filename+'/'+self.ID+sasbdb+"Kratky.html")
        p.output_backend="svg"
        export_svgs(p,filename=filename+'/'+self.ID+sasbdb+"Kratky.svg")
        export_png(p,filename=filename+'/'+self.ID+sasbdb+"Kratky.png")
    
    def plot_porod_debye(self,sasbdb,df):
        output_file(self.ID+sasbdb+"porod.html",mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=500, plot_width=500, title="Porod-Debye plot ("+sasbdb+")")
        p.circle(x='Px',y='Py',source=source,color='blue',fill_alpha=0.3,size=5)
        p.xaxis.major_label_text_font_size="14pt"
        p.yaxis.major_label_text_font_size="14pt"
        p.title.text_font_size='12pt'
        p.title.align="center"
        p.title.vertical_align='top'
        p.xaxis.axis_label = 'q \u2074'
        p.xaxis.axis_label_text_font_size='14pt'
        p.yaxis.axis_label = 'q\u2074 I(q)'
        p.yaxis.axis_label_text_font_size='14pt'
        dirname=os.path.dirname(os.path.abspath(__file__))
        filename = os.path.abspath(os.path.join(dirname, '../static/images/'))
        save(p,filename=filename+'/'+self.ID+sasbdb+"porod.html")
        p.output_backend="svg"
        export_svgs(p,filename=filename+'/'+self.ID+sasbdb+"porod.svg")
        export_png(p,filename=filename+'/'+self.ID+sasbdb+"porod.png")

    def plot_pddf(self,sasbdb,df):
        output_file(self.ID+sasbdb+"pddf.html",mode="inline") 
        source = ColumnDataSource(df)
        p = figure(plot_height=500, plot_width=500, title="Pair distance distribution function ("+sasbdb+")")
        p.circle(x='R',y='P',source=source,color='blue',fill_alpha=0.3,size=5)
        p.multi_line('err_x','err_y',source=source, color='gray',line_width=1.5)
        p.xaxis.major_label_text_font_size="14pt"
        p.yaxis.major_label_text_font_size="14pt"
        p.title.text_font_size='12pt'
        p.title.align="center"
        p.title.vertical_align='top'
        p.xaxis.axis_label = "r nm"
        p.xaxis.axis_label_text_font_size='14pt'
        p.yaxis.axis_label = 'P(r)'
        p.yaxis.axis_label_text_font_size='14pt'
        dirname=os.path.dirname(os.path.abspath(__file__))
        filename = os.path.abspath(os.path.join(dirname, '../static/images/'))
        print (filename)
        save(p,filename=filename+'/'+self.ID+sasbdb+"pddf.html")
        p.output_backend="svg"
        export_svgs(p,filename=filename+'/'+self.ID+sasbdb+"pddf.svg")
        export_png(p,filename=filename+'/'+self.ID+sasbdb+"pddf.png")

    def Guinier_plot_fit(self,sasbdb,df,score):
        output_file(self.ID+sasbdb+"guinier.html",mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=500, plot_width=500, title="Guinier plot for "+sasbdb+" (R\u00B2="+str(score)+")")
        legend1='Observed';legend2="Predicted"
        p.circle(x='Q2A',y='logI',source=source, color='blue',line_width=1,fill_alpha=0.3,size=5,legend=legend1)
        p.line(x='Q2A',y='y_pred',source=source, color='red',line_width=3,legend=legend2)
        p.legend.orientation = "vertical"
        p.legend.location = "top_right"
        p.xaxis.major_label_text_font_size="14pt"
        p.yaxis.major_label_text_font_size="14pt"
        p.title.text_font_size='12pt'
        p.title.align="center"
        p.title.vertical_align='top'
        p.xaxis.axis_label = "q\u00B2 nm\u00B2" #\u212B\u207B\u00B2"
        p.xaxis.axis_label_text_font_size='14pt'
        p.yaxis.axis_label = 'Log I(q)'
        p.yaxis.axis_label_text_font_size='14pt'
        dirname=os.path.dirname(os.path.abspath(__file__))
        filename = os.path.abspath(os.path.join(dirname, '../static/images/'))
        save(p,filename=filename+'/'+self.ID+sasbdb+"guinier.html")
        p.output_backend="svg"
        export_svgs(p,filename=filename+'/'+self.ID+sasbdb+"guinier.svg")
        export_png(p,filename=filename+'/'+self.ID+sasbdb+"guinier.png")

    def Guinier_plot_residuals(self,sasbdb,df):
        output_file(self.ID+sasbdb+"guinier_residuals.html",mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=500, plot_width=500, title="Residuals for Guinier plot fit ("+sasbdb+")")
        p.circle(x='Q2A',y='res',source=source, color='blue',fill_alpha=0.3,size=5)
        hline = Span(location=0, dimension='width', line_color='red', line_width=3)
        p.renderers.extend([hline])
        p.xaxis.major_label_text_font_size="14pt"
        p.yaxis.major_label_text_font_size="14pt"
        p.title.text_font_size='12pt'
        p.title.align="center"
        p.title.vertical_align='top'
        p.xaxis.axis_label = "q\u00B2 nm\u00B2"#\u212B\u207B\u00B2"
        p.xaxis.axis_label_text_font_size='14pt'
        p.yaxis.axis_label = 'R'
        p.yaxis.axis_label_text_font_size='14pt'
        dirname=os.path.dirname(os.path.abspath(__file__))
        filename = os.path.abspath(os.path.join(dirname, '../static/images/'))
        save(p,filename=filename+'/'+self.ID+sasbdb+"guinier_residuals.html")
        p.output_backend="svg"
        export_svgs(p,filename=filename+'/'+self.ID+sasbdb+"guinier_residuals.svg")
        export_png(p,filename=filename+'/'+self.ID+sasbdb+"guinier_residuals.png")

    def plot_fit(self,sasbdb,fit,score,df):
        output_file(self.ID+sasbdb+str(fit)+"fit1.html",mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=500, plot_width=500, title="Model fit for "+sasbdb)        
        legend1='Observed';legend2="Predicted"
        p.circle(x='Q',y='logIe',source=source, color='blue',line_width=1,fill_alpha=0.3,size=3,legend=legend1)
        p.line(x='Q',y='logIb',source=source, color='red',line_width=3,legend=legend2)        
        #p.circle(x='Q',y='logIb',source=source, color='red',line_width=1,fill_alpha=0.1,size=3,legend=legend2)
        p.legend.orientation = "vertical"
        p.legend.location = "top_right"
        p.xaxis.major_label_text_font_size="14pt"
        p.yaxis.major_label_text_font_size="14pt"
        p.title.text_font_size='12pt'
        p.title.align="center"
        p.title.vertical_align='top'
        p.xaxis.axis_label = "q \u212B\u207B\u207B\u00B9"
        p.xaxis.axis_label_text_font_size='14pt'
        p.yaxis.axis_label = 'Log I(q)'
        p.yaxis.axis_label_text_font_size='14pt'
        dirname=os.path.dirname(os.path.abspath(__file__))
        filename = os.path.abspath(os.path.join(dirname, '../static/images/'))
        save(p,filename=filename+'/'+self.ID+sasbdb+str(fit)+"fit1.html")
        p.output_backend="svg"
        export_svgs(p,filename=filename+'/'+self.ID+sasbdb+str(fit)+"fit1.svg")
        export_png(p,filename=filename+'/'+self.ID+sasbdb+str(fit)+"fit1.png")

    def plot_fit_residuals(self,sasbdb,fit,df):
        output_file(self.ID+sasbdb+str(fit)+"residuals.html",mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=500, plot_width=500, title="Residuals for model fit ("+sasbdb+")")
        p.circle(x='Q',y='r',source=source, color='blue',fill_alpha=0.3,size=5)
        hline = Span(location=0, dimension='width', line_color='red', line_width=3)
        p.renderers.extend([hline])
        p.xaxis.major_label_text_font_size="14pt"
        p.yaxis.major_label_text_font_size="14pt"
        p.title.text_font_size='12pt'
        p.title.align="center"
        p.title.vertical_align='top'
        p.xaxis.axis_label = "q \u212B\u207B\u207B\u00B9"
        p.xaxis.axis_label_text_font_size='14pt'
        p.yaxis.axis_label = 'R'
        p.yaxis.axis_label_text_font_size='14pt'
        dirname=os.path.dirname(os.path.abspath(__file__))
        filename = os.path.abspath(os.path.join(dirname, '../static/images/'))
        save(p,filename=filename+'/'+self.ID+sasbdb+str(fit)+"residuals.html")
        p.output_backend="svg"
        export_svgs(p,filename=filename+'/'+self.ID+sasbdb+str(fit)+"residuals.svg")
        export_png(p,filename=filename+'/'+self.ID+sasbdb+str(fit)+"residuals.png")


