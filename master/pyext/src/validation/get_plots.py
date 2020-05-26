import glob
import sys,os,re
import ihm
import ihm.reader
import validation
from bokeh.io import output_file, show, curdoc, export_png, export_svgs
from bokeh.models import ColumnDataSource, LinearAxis, Legend
from bokeh.palettes import GnBu3, Blues, OrRd3, Spectral, Set1,RdYlBu,Paired
from bokeh.plotting import figure, output_file, save
from bokeh.models.widgets import Tabs, Panel
from bokeh.layouts import row,column

class plots(validation.get_input_information):
	def __init__(self,mmcif):
		super().__init__(mmcif)
		self.ID=str(validation.get_input_information.get_id(self))

	def plot_quality_at_glance(self,clashscore,rama,sidechain,exv_data,sas_data,sas_fit):
		output_file(self.ID+"quality_at_glance.html",mode="inline")
		print (clashscore,rama,sidechain,exv_data,sas_data,sas_fit)
		tabsI=[]
		if clashscore or rama or sidechain:
			counts=[clashscore,rama,sidechain]
			Scores = ['Clashscore', 'Ramachandran Outliers', 'Sidechain Outliers']
			legends=[str(counts[0]),str(counts[1])+'%',str(counts[2])+'%']
			source = ColumnDataSource(data=dict(Scores=Scores, counts=counts, legends=legends, color=RdYlBu[9]))
			p = figure(y_range=Scores, x_range=(0,max(counts)+1), plot_height=250, plot_width=700, title="Model Quality: Molprobity Analysis")
			p.hbar(y='Scores',right='counts', height=0.75, color='color', legend="legends", source=source,alpha=0.8)
		elif exv_data:
			model=exv_data['Models']
			satisfaction=exv_data['Number of violations']
			print (satisfaction)
			counts=[float(i) for i in satisfaction]
			Scores = ['Model '+ str(i+1) for i,j in enumerate (model)]
			legends=['Model '+ str(i+1)+ ': ' +str(int(j))  for i,j in enumerate(counts)]
			n=3 if len(model)< 3 else len(model)
			source = ColumnDataSource(data=dict(Scores=Scores, counts=counts, legends=legends, color=RdYlBu[n]))
			if len(model)>1:
				upper=max(counts)+ (max(counts)-min(counts))
				lower=min(counts)-(max(counts)-min(counts))
			else:
				upper=max(counts)+ 20
				lower=min(counts)-20
			p = figure(y_range=Scores, x_range=(lower,upper), plot_height=250, plot_width=700, title='Model quality: Excluded Volume Analysis')
			p.hbar(y='Scores',right='counts', height=0.5, color='color', legend="legends", source=source,alpha=0.8)
			p.xaxis.axis_label = 'Number of violations'
		
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

		p.xaxis.axis_label_text_font_size='14pt'
		p.legend.orientation = "vertical"
		p.legend.location = "top_right"
		p.ygrid.grid_line_color = None
		p.xaxis.major_label_text_font_size="14pt"
		p.yaxis.major_label_text_font_size="14pt"
		p.title.text_font_size='14pt'
		p.title.align="center"
		p.title.vertical_align='top'
		first = Panel(child=row(p), title='Model quality')
		tabsI.append(first)
		if len(sas_data.keys())>0:
			Rgl={0:'p(r)',1:'Guinier'}
			Scores=['Rg from '+ Rgl[m]+ ' ('+i+')' for i,j in sas_data.items() for m,n in enumerate(j) ]
			counts=[float(n)for i,j in sas_data.items() for m,n in enumerate(j) ];
			legends=[str(i)+' nm' for i in counts]
			print (Scores)
			print (counts)
			print (legends)
			source = ColumnDataSource(data=dict(Scores=Scores, counts=counts, legends=legends, color=RdYlBu[9]))
			pd = figure(y_range=Scores, x_range=(0,max(counts)+1), plot_height=250, plot_width=700, title="Data Quality: Rg Analysis")
			pd.hbar(y='Scores',right='counts', height=0.5, color='color', legend="legends", source=source,alpha=0.8)
			pd.legend.orientation = "vertical"
			pd.legend.location = "top_right"
			pd.ygrid.grid_line_color = None
			pd.xaxis.major_label_text_font_size="14pt"
			pd.yaxis.major_label_text_font_size="14pt"
			pd.title.text_font_size='14pt'
			pd.title.align="center"
			second = Panel(child=row(pd), title='Data quality: SAS') 
			tabsI.append(second)
		else:
			second = Panel(child=row(p1), title='Data quality')
			tabsI.append(second)
		if len(sas_fit.keys())>0:
			Scores=['\u03C7\u00b2 value fit '+ str(int(m+1)) +' ('+i+')' for i,j in sas_fit.items() for m,n in enumerate(j)]
			counts=[ float(n) for i,j in sas_fit.items() for m,n in enumerate(j)]
			legends=[str(i) for i in counts]
			print (Scores)
			print (counts)
			print (legends)
			source = ColumnDataSource(data=dict(Scores=Scores, counts=counts, legends=legends, color=RdYlBu[9]))
			pf = figure(y_range=Scores, x_range=(0,max(counts)+1), plot_height=250, plot_width=700, title="Fit to SAS input: \u03C7\u00b2 Fit")
			pf.hbar(y='Scores',right='counts', height=0.25*len(list(sas_fit.values())[0]), color='color', legend="legends", source=source,alpha=0.8)
			pf.legend.orientation = "vertical"
			pf.legend.location = "top_right"
			pf.ygrid.grid_line_color = None
			pf.xaxis.major_label_text_font_size="14pt"
			pf.yaxis.major_label_text_font_size="14pt"
			pf.title.text_font_size='14pt'
			pf.title.align="center"
			third = Panel(child=row(pf), title='Fit to data: SAS')
			tabsI.append(third)
		else:
			third = Panel(child=row(p1), title='Fit to input data')
			tabsI.append(third)
		fourth = Panel(child=row(p1), title='Fit to data not used for modeling')
		tabsI.append(fourth)
		fifth=Panel(child=row(p1), title='Uncertainty in models')
		tabsI.append(fifth)
		tabs = Tabs(tabs=tabsI)
		curdoc().add_root(tabs)
		dirname=os.path.dirname(os.path.abspath(__file__))
		filename = os.path.abspath(os.path.join(dirname, '../static/images//'))
		print (filename)
		save(tabs,filename=filename+'/'+self.ID+"quality_at_glance.html")
		p.output_backend="svg"
		export_svgs(p,filename=filename+'/'+self.ID+"quality_at_glance.svg")
		export_png(p,filename=filename+'/'+self.ID+"quality_at_glance.png")

