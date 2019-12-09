import glob
import sys,os,re
import ihm
import ihm.reader
import validation
from bokeh.io import output_file, show, curdoc, export_png, export_svgs
from bokeh.models import ColumnDataSource, LinearAxis, Legend
from bokeh.palettes import GnBu3, Blues, OrRd3, Spectral6, Set1
from bokeh.plotting import figure, output_file, save
from bokeh.models.widgets import Tabs, Panel
from bokeh.layouts import row,column

class plots_mp(validation.get_molprobity_information.get_molprobity_information):
    def __init__(self,mmcif):
        super().__init__(mmcif)
        self.ID=str(validation.get_input_information.get_id(self))

    def plot_quality_at_glance_mp(self,clashscore,rama,sidechain):
        output_file(self.ID+"quality_at_glance.html",mode="inline")
        counts=[clashscore,rama,sidechain]
        Scores = ['Clashscore', 'Ramachandran Outliers', 'Sidechain Outliers']
        legends=[str(counts[0]),str(counts[1])+'%',str(counts[2])+'%']
        source = ColumnDataSource(data=dict(Scores=Scores, counts=counts, legends=legends, color=Blues[6]))
        p = figure(y_range=Scores, x_range=(0,max(counts)+1), plot_height=250, plot_width=700, title="Model Quality: Molprobity Analysis")
        p.hbar(y='Scores',right='counts', height=0.5, color='color', legend="legends", source=source,alpha=0.8)
        p.legend.orientation = "vertical"
        p.legend.location = "top_right"
        p.ygrid.grid_line_color = None
        p.xaxis.major_label_text_font_size="14pt"
        p.yaxis.major_label_text_font_size="14pt"
        p.title.text_font_size='14pt'
        p.title.align="center"
        p.title.vertical_align='top'
        Scores=[''];counts = [''];legends=['']
        source = ColumnDataSource(data=dict(Scores=Scores, counts=counts, legends=legends))
        p1 = figure(y_range=Scores, x_range=(0,1), plot_height=250, plot_width=700)
        p1.legend.orientation = "vertical"
        p1.legend.location = "top_right"
        p1.ygrid.grid_line_color = None
        p1.xaxis.major_label_text_font_size="14pt"
        p1.yaxis.major_label_text_font_size="14pt"
        p1.title.text_font_size='14pt'
        p1.title.align="center"
        first = Panel(child=row(p), title='Model quality')
        second = Panel(child=row(p1), title='Data quality')
        third = Panel(child=row(p1), title='Fit to input data')
        fourth = Panel(child=row(p1), title='Fit to data not used for modeling')
        fifth=Panel(child=row(p1), title='Uncertainty in models')
        tabs = Tabs(tabs=[first,second,third,fourth,fifth])
        curdoc().add_root(tabs)
        dirname=os.path.dirname(os.path.abspath(__file__))
        filename = os.path.abspath(os.path.join(dirname, '../static/images//'))
        print (filename)
        save(tabs,filename=filename+'/'+self.ID+"quality_at_glance.html")
        p.output_backend="svg"
        export_svgs(p,filename=filename+'/'+self.ID+"quality_at_glance.svg")
        export_png(p,filename=filename+'/'+self.ID+"quality_at_glance.png")


class plots_exv(validation.get_excluded_volume.get_excluded_volume):
    def __init__(self,mmcif):
        super().__init__(mmcif)
        self.ID=str(validation.get_input_information.get_id(self))


    def plot_quality_at_glance_ev(self,exv_data):
        output_file(self.ID+"quality_at_glance.html",mode="inline")
        model=exv_data['Models']
        satisfaction=exv_data['Number of violations']
        print (satisfaction)
        counts=[float(i) for i in satisfaction]
        Scores = ['Model '+ str(i+1) for i,j in enumerate (model)]
        legends=['Model '+ str(i+1)+ ': ' +str(int(j))  for i,j in enumerate(counts)]
        n=3 if len(model)< 3 else len(model)
        source = ColumnDataSource(data=dict(Scores=Scores, counts=counts, legends=legends, color=Blues[n]))
        if len(model)>1:
            upper=max(counts)+ (max(counts)-min(counts))
            lower=min(counts)-(max(counts)-min(counts))
        else:
            upper=max(counts)+ 20
            lower=min(counts)-20
        p = figure(y_range=Scores, x_range=(lower,upper), plot_height=250, plot_width=700, title='Model quality: Excluded Volume Analysis')
        p.hbar(y='Scores',right='counts', height=0.5, color='color', legend="legends", source=source,alpha=0.8)
        p.xaxis.axis_label = 'Number of violations'
        p.xaxis.axis_label_text_font_size='14pt'
        p.legend.orientation = "vertical"
        p.legend.location = "top_right"
        p.ygrid.grid_line_color = None
        p.xaxis.major_label_text_font_size="14pt"
        p.yaxis.major_label_text_font_size="14pt"
        p.title.text_font_size='14pt'
        p.title.align="center"
        p.title.vertical_align='top'
        Scores=[''];counts = [''];legends=['']
        source = ColumnDataSource(data=dict(Scores=Scores, counts=counts, legends=legends))
        p1 = figure(y_range=Scores, x_range=(0,1), plot_height=250, plot_width=700)
        p1.legend.orientation = "vertical"
        p1.legend.location = "top_right"
        p1.ygrid.grid_line_color = None
        p1.xaxis.major_label_text_font_size="14pt"
        p1.yaxis.major_label_text_font_size="14pt"
        p1.title.text_font_size='14pt'
        p1.title.align="center"
        first = Panel(child=row(p), title='Model quality')
        second = Panel(child=row(p1), title='Data quality')
        third = Panel(child=row(p1), title='Fit to input data')
        fourth = Panel(child=row(p1), title='Fit to data not used for modeling')
        fifth=Panel(child=row(p1), title='Uncertainty in models')
        tabs = Tabs(tabs=[first,second,third,fourth,fifth])
        curdoc().add_root(tabs)
        dirname=os.path.dirname(os.path.abspath(__file__))
        filename = os.path.abspath(os.path.join(dirname, '../static/images/'))
        print (filename)
        save(tabs,filename=filename+'/'+self.ID+"quality_at_glance.html")
        p.output_backend="svg"
        export_svgs(p,filename=filename+'/'+self.ID+"quality_at_glance.svg")
        export_png(p,filename=filename+'/'+self.ID+"quality_at_glance.png")

