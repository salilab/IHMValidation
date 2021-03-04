###################################
# Script : 
# 1) Contains class for plots that
# combines infromation from multiple
# datasets 
#
# ganesans - Salilab - UCSF
# ganesans@salilab.org
###################################
import glob
import sys,os,re
import ihm
import ihm.reader
import validation
import bokeh
from bokeh.io import output_file, show, curdoc, export_png, export_svgs
from bokeh.models import ColumnDataSource, LinearAxis, Legend
from bokeh.palettes import GnBu3, Blues, OrRd3, Spectral, Set1,RdYlBu,Paired,Viridis256,viridis
from bokeh.plotting import figure, output_file, save
from bokeh.models.widgets import Tabs, Panel
from bokeh.layouts import row,column,gridplot

from bokeh.core.validation import silence
from bokeh.core.validation.warnings import MISSING_RENDERERS,EMPTY_LAYOUT
silence(MISSING_RENDERERS, True)
silence(EMPTY_LAYOUT, True)

class plots(validation.get_input_information):
	def __init__(self,mmcif):
		super().__init__(mmcif)
		self.ID=str(validation.get_input_information.get_id(self))
		self.dirname=os.path.dirname(os.path.abspath(__file__))
		self.filename = os.path.join('../static/images//')

	def plot_quality_at_glance(self,clashscore: list,rama: list,sidechain: list, \
		exv_data: dict,sas_data: dict,sas_fit: dict,cx_fit:dict)->bokeh.plotting.figure:
				
		# create tabs list to add all the panel figures (model quality, data quality.. etc)
		output_file(self.ID+"quality_at_glance.html",mode="inline")
		tabsI=[]

		#MODEL QUALITY
		## check for molprobity or excluded volume data
		if clashscore or rama or sidechain:
			#if molprobity data, plot that
			counts=[clashscore,rama,sidechain]
			Scores = ['Clashscore', 'Ramachandran Outliers', 'Sidechain Outliers']
			legends=[str(counts[0]),str(counts[1])+'%',str(counts[2])+'%']
			source = ColumnDataSource(data=dict(Scores=Scores, counts=counts, legends=legends, color=viridis(len(Scores))))
			p = figure(y_range=Scores, x_range=(0,max(counts)+1), plot_height=250, plot_width=700, title="Model Quality: Molprobity Analysis")
			p.hbar(y='Scores',right='counts', height=0.75, color='color', legend_group='legends', source=source,alpha=0.8)
		elif exv_data:
			# if excluded vol data, plot that
			model=exv_data['Models']
			satisfaction=exv_data['Number of violations']
			counts=[float(i) for i in satisfaction]
			Scores = ['Model '+ str(i+1) for i,j in enumerate (model)]
			legends=['Model '+ str(i+1)+ ': ' +str(int(j))  for i,j in enumerate(counts)]
			n=3 if len(model)< 3 else len(model)
			source = ColumnDataSource(data=dict(Scores=Scores, counts=counts, legends=legends, color=viridis(n)))
			if len(model)>1:
				upper=max(counts)+ (max(counts)-min(counts))
				lower=min(counts)-(max(counts)-min(counts))
			else:
				upper=max(counts)+ 20
				lower=min(counts)-20
			p = figure(y_range=Scores, x_range=(lower,upper), plot_height=250, plot_width=700, title='Model quality: Excluded Volume Analysis',sizing_mode='stretch_both')
			p.hbar(y='Scores',right='counts', height=0.5, color='color', legend_group='legends', source=source,alpha=0.8)
			p.xaxis.axis_label = 'Number of violations'
		else:
			Scores=[''];counts = [''];legends=['']
			source = ColumnDataSource(data=dict(Scores=Scores, counts=counts, legends=legends))
			p = figure(y_range=Scores, x_range=(0,1), plot_height=250, plot_width=700)
		# if neither: ie. atomic but unable to calculate molprobity, return empty plot
		Scores=[''];counts = [''];legends=['']
		source = ColumnDataSource(data=dict(Scores=Scores, counts=counts, legends=legends))
		p1 = figure(y_range=Scores, x_range=(0,1), plot_height=250, plot_width=700)
		p.legend.orientation=p1.legend.orientation = "vertical"
		p.legend.location=p1.legend.location = "top_right"
		p.ygrid.grid_line_color=p1.ygrid.grid_line_color = None
		p.xaxis.axis_label_text_font_size=p1.xaxis.axis_label_text_font_size="14pt"
		p.yaxis.axis_label_text_font_size=p1.yaxis.axis_label_text_font_size="14pt"
		p.title.text_font_size=p1.title.text_font_size='14pt'
		p.title.align=p1.title.align="center"
		p.output_backend=p1.output_backend="svg"
		p.title.vertical_align='top'
		# make panel figures
		# first panel is model quality
		first = Panel(child=row(p), title='Model quality')
		tabsI.append(first)
		export_svgs(p,filename=self.filename+'/'+self.ID+"quality_at_glance_MQ.svg")
		export_png(p,filename=self.filename+'/'+self.ID+"quality_at_glance_MQ.png")

		#DATA QUALITY
		#check for sas data, if exists, plot
		# this section will be updated with more data assessments, as and when it is complete
		if len(sas_data.keys())>0:
			Rgl={0:'p(r)',1:'Guinier'}
			Scores=['Rg from '+ Rgl[m]+ ' ('+i+')' for i,j in sas_data.items() for m,n in enumerate(j) ]
			counts=[float(n)for i,j in sas_data.items() for m,n in enumerate(j) ];
			legends=[str(i)+' nm' for i in counts]
			source = ColumnDataSource(data=dict(Scores=Scores, counts=counts, legends=legends, color=viridis(len(legends))))
			pd = figure(y_range=Scores, x_range=(0,max(counts)+1), plot_height=250, plot_width=700, title="Data Quality: Rg Analysis",sizing_mode='stretch_both')
			pd.hbar(y='Scores',right='counts', height=0.5, color='color', legend_group="legends", source=source,alpha=0.8)
			pd.legend.orientation = "vertical"
			pd.legend.location = "top_right"
			pd.sizing_mode = 'stretch_both'
			pd.legend.label_text_font_size = "8px"
			pd.ygrid.grid_line_color = None
			pd.xaxis.major_label_text_font_size="14pt"
			pd.yaxis.major_label_text_font_size="14pt"
			pd.title.text_font_size='14pt'
			pd.title.align="center"
			second = Panel(child=row(pd), title='Data quality: SAS') 
			tabsI.append(second)
			pd.output_backend="svg"
			p = gridplot([[pd]], sizing_mode='stretch_both',plot_height=450, plot_width=800)
			export_svgs(p,filename=self.filename+'/'+self.ID+"quality_at_glance_DQ.svg")
			export_png(p,filename=self.filename+'/'+self.ID+"quality_at_glance_DQ.png")
		else:
			second = Panel(child=row(p1), title='Data quality')
			tabsI.append(second)

		#FIT TO DATA QUALITY
		#check for sas data, if exists, plot
		# this section will be updated with more data assessments, as and when it is complete
		if len(sas_fit.keys())>0:
			Scores=['\u03C7\u00b2 value fit '+ str(int(m+1)) +' ('+i+')' for i,j in sas_fit.items() for m,n in enumerate(j)]
			counts=[ float(n) for i,j in sas_fit.items() for m,n in enumerate(j)]
			legends=[str(i) for i in counts]
			source = ColumnDataSource(data=dict(Scores=Scores, counts=counts, legends=legends, color=viridis(len(legends))))
			pf = figure(y_range=Scores, x_range=(0,max(counts)+1), plot_height=250, plot_width=700, title="Fit to SAS input: \u03C7\u00b2 Fit",sizing_mode='stretch_both')
			pf.hbar(y='Scores',right='counts', height=0.25*len(list(sas_fit.values())[0]), color='color', legend_group="legends", source=source,alpha=0.8)
			pf.legend.orientation = "vertical"
			pf.legend.location = "top_right"
			pf.ygrid.grid_line_color = None
			pf.xaxis.major_label_text_font_size="14pt"
			pf.yaxis.major_label_text_font_size="14pt"
			pf.title.text_font_size='14pt'
			pf.legend.label_text_font_size = "8px"
			pf.title.align="center"
			pf.output_backend="svg"
			third = Panel(child=row(pf), title='Fit to data: SAS')
			tabsI.append(third)
			pf.output_backend="svg"
			export_svgs(pf,filename=self.filename+'/'+self.ID+"quality_at_glance_FQ.svg")
			export_png(pf,filename=self.filename+'/'+self.ID+"quality_at_glance_FQ.png")

		#check for XL_MS data, if exists, plot
		if len(cx_fit.keys())>0:
			Scores=['Model '+str(i) for i,j in cx_fit.items()]
			counts=[round(float(j),2) for i,j in cx_fit.items()]
			#legends=[str(i) for i in counts]
			legends=['Model '+ str(i+1)+ ': ' +str(j)+'%'  for i,j in enumerate(counts)]
			source = ColumnDataSource(data=dict(Scores=Scores, counts=counts, legends=legends, color=viridis(len(legends))))
			pf1 = figure(y_range=Scores, x_range=(0,max(counts)+1), plot_height=250, plot_width=700, title="Fit to XL-MS input",sizing_mode='stretch_both')
			pf1.hbar(y='Scores',right='counts', height=0.15*len(list(cx_fit.values())), color='color', legend_group="legends", source=source,alpha=0.8)
			pf1.legend.orientation = "vertical"
			pf1.legend.location = "top_right"
			pf1.ygrid.grid_line_color = None
			pf1.xaxis.major_label_text_font_size="14pt"
			pf1.yaxis.major_label_text_font_size="14pt"
			pf1.title.text_font_size='14pt'
			pf1.title.align="center"
			pf1.output_backend="svg"
			third1 = Panel(child=row(pf1), title='Fit to data: XL-MS')
			tabsI.append(third1)
			pf1.output_backend="svg"
			export_svgs(pf1,filename=self.filename+'/'+self.ID+"quality_at_glance_FQ1.svg")
			export_png(pf1,filename=self.filename+'/'+self.ID+"quality_at_glance_FQ1.png")

		#check for XL_MS or SAS data, if exists, plot third panel as FIT TO INPUT
		if (len(cx_fit.keys())<1) or (len(sas_fit.keys())<1) :
			third = Panel(child=row(p1), title='Fit to input data')
			tabsI.append(third)

		# remaining panels are empty for now, will be updated
		fourth = Panel(child=row(p1), title='Fit to data not used for modeling')
		tabsI.append(fourth)
		fifth=Panel(child=row(p1), title='Uncertainty in models')
		tabsI.append(fifth)
		tabs = Tabs(tabs=tabsI)
		curdoc().add_root(tabs)
		save(tabs,filename=self.filename+'/'+self.ID+"quality_at_glance.html")


