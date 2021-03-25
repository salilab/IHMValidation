###################################
# Script : 
# 1) Contains class to generate XL-MS
# plots 
# 2) Inherits from CX class
#
# ganesans - Salilab - UCSF
# ganesans@salilab.org
###################################
import pandas as pd
import glob
import sys,os,math,itertools
import numpy as np
import pandas as pd
from validation import sas, get_input_information,cx
from bokeh.io import output_file, show, curdoc, export_png, export_svgs
from bokeh.models import Span,ColumnDataSource, LinearAxis, Legend, FactorRange
from bokeh.palettes import GnBu3, RdBu,OrRd3,Blues,YlOrBr, Spectral6, Set1
from bokeh.plotting import figure, output_file, save
from bokeh.models.widgets import Tabs, Panel
from bokeh.layouts import row,column,gridplot


class CxValidationPlots(cx.CxValidation):
	def __init__(self,mmcif_file):
		super().__init__(mmcif_file)
		self.ID=str(GetInputInformation.get_id(self))
		self.xl_df=cx.CxValidation.get_xl_data(self)
		self.model_df=cx.CxValidation.get_df_for_models(self)        
		self.filename = os.path.join('../static/images//')


	def plot_linker_dist_I(self,df,intra=1,key='Intra'):
		'''
		plot distance distribution per linker
		based on inter and intra links
		'''
		for i in df['Linker'].unique():
			df_c=df[df['Linker']==i]
			if i=='DSS':
				loc=30
			elif i=='EDC':
				loc=20
			else:
				loc=30
			output_file(self.ID+i+"linker.html",mode="inline")
			measured=df_c[df_c['Intra']==intra]['dist']
			hist, edges = np.histogram(measured, density=False, bins=50)
			#hist_l, edges_l = np.histogram(measured, density=False, bins=25)
			p = figure(title=key+'-molecular distances/Linker '+i, 
				plot_height=400, plot_width=400)
			p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
						fill_color="navy", line_color="white", alpha=0.3)
			#p.line(edges, hist, line_color="navy", line_width=4, alpha=0.7, legend_label=key+"/"+i)
			vline = Span(location=loc, dimension='height', line_color='red', line_width=3,line_dash='dashed')
			p.renderers.extend([vline])
			p.xaxis.major_label_text_font_size="14pt"
			p.yaxis.major_label_text_font_size="14pt"
			p.title.text_font_size='12pt'
			p.title.align="center"
			p.title.vertical_align='top'
			p.xaxis.axis_label = 'Distance \u212B'
			p.xaxis.axis_label_text_font_size='14pt'
			p.yaxis.axis_label = 'Number of cross-links'
			p.yaxis.axis_label_text_font_size='14pt'
			p.output_backend="svg"
			save(p,filename=self.filename+'/'+self.ID+i+key+"linker.html")
			export_svgs(p,filename=self.filename+'/'+self.ID+i+key+"linker.svg")

	def plot_linker_dist_S(self,df,struc=1,key='Structured'):
		'''
		plot distance distribution per linker
		based on structured/unstrcutured/between struc&unstruc
		'''
		for i in df['Linker'].unique():
			df_c=df[df['Linker']==i]
			if i=='DSS':
				loc=30
			elif i=='EDC':
				loc=20
			else:
				loc=30
			output_file(self.ID+i+"linker.html",mode="inline")
			measured=df_c[df_c['Structured']==struc]['dist']
			hist, edges = np.histogram(measured, density=False, bins=50)
			#hist_l, edges_l = np.histogram(measured, density=False, bins=25)
			p = figure(title=key+ ' regions/Linker '+i, 
				plot_height=350, plot_width=350)
			p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
						fill_color="navy", line_color="white", alpha=0.3)
			#p.line(edges, hist, line_color="navy", line_width=4, alpha=0.7, legend_label=key+"/"+i)
			vline = Span(location=loc, dimension='height', line_color='red', line_width=3,line_dash='dashed')
			p.renderers.extend([vline])
			p.xaxis.major_label_text_font_size="14pt"
			p.yaxis.major_label_text_font_size="14pt"
			p.title.text_font_size='12pt'
			p.title.align="center"
			p.title.vertical_align='top'
			p.xaxis.axis_label = 'Distance \u212B'
			p.xaxis.axis_label_text_font_size='14pt'
			p.yaxis.axis_label = 'Number of cross-links'
			p.yaxis.axis_label_text_font_size='14pt'
			p.output_backend="svg"
			save(p,filename=self.filename+'/'+self.ID+i+key+"linker.html")
			export_svgs(p,filename=self.filename+'/'+self.ID+i+key+"linker.svg")

	def plot_intra_summary_deprecated(self):
		'''
		plot summary of intra/inter xl-ms satisfaction
		'''
		for model_id,df in self.model_df.items():
			factors=self.get_factors(df)
			regions = ['Satisfied', 'Violated']
			source = ColumnDataSource(data=dict(
	    				x=factors,
	    				Satisfied=self.get_satisfied(df,factors),
	    				Violated=self.get_violated(df,factors), ))
			fig_id = figure(x_range=FactorRange(*factors), plot_height=400,plot_width=500,
	           			title='CX-MS Satisfaction for model:'+str(model_id))

			fig_id.vbar_stack(regions, x='x', width=0.9, alpha=0.5, color=["blue", "red"], source=source,
	             				legend_label=regions)
			fig_id.xaxis.major_label_text_font_size="14pt"
			fig_id.yaxis.major_label_text_font_size="14pt"
			fig_id.yaxis.axis_label_text_font_size='14pt'
			fig_id.title.text_font_size='12pt'
			fig_id.title.align="center"
			fig_id.title.vertical_align='top'
			fig_id.yaxis.axis_label = 'Number of cross-links'
			fig_id.y_range.start = 0
			fig_id.y_range.end = df.shape[0]
			fig_id.x_range.range_padding = 0.1
			fig_id.xaxis.major_label_orientation = 1
			fig_id.xgrid.grid_line_color = None
			fig_id.legend.location = "top_center"
			fig_id.legend.orientation = "horizontal"
			fig_id.output_backend="svg"
			save(fig_id,filename=self.filename+'/'+self.ID+str(model_id)+"IS.html")
			export_svgs(fig_id,filename=self.filename+'/'+self.ID+str(model_id)+"IS.svg")


	def plot_distributions(self):
		'''
		plot inter and intra distance distributions
		'''
		for model_id,df in self.model_df.items():
			self.plot_linker_dist_I(df,intra=1,key='Intra')
			self.plot_linker_dist_I(df,intra=0,key='Inter')
			self.plot_linker_dist_S(df,struc=1,key='Structured')
			self.plot_linker_dist_S(df,struc=0,key='Unstructured')
			self.plot_linker_dist_S(df,struc=2,key='Intermediate')

	def get_factors(self,df):
		'''
		get grouped inter/intra factors for stacked bar plot
		'''
		link=df['Linker'].unique()
		xl=['Inter','Intra']
		factors=list(itertools.product(link, xl))
		return factors

	def get_factors_struc(self,df):
		'''
		get grouped struc factors for stacked bar plot
		'''
		link=df['Linker'].unique()
		xl=['Structured','Unstructured','Intermediate']
		factors=list(itertools.product(link, xl))
		return factors

	def get_satisfied(self,df,factors):
		'''
		get satisfied list for stacked bar plot;inter/intra info
		'''
		Satisfied=[]
		for i in factors:
			df_1=df[df['Linker']==i[0]]
			df_2=df_1[df_1[i[1]]==1]
			Satisfied.append(df_2[df_2['Satisfied']==1].shape[0])
		return Satisfied

	def get_satisfied_struc(self,df,factors):
		'''
		get satisfied list for stacked bar plot;struc info
		'''
		Satisfied=[]
		struc_dict={'Structured':1,'Unstructured':0,'Intermediate':2}
		for i in factors:
			df_1=df[df['Linker']==i[0]]
			df_2=df_1[df_1['Structured']==struc_dict[i[1]]]
			Satisfied.append(df_2[df_2['Satisfied']==1].shape[0])
		return Satisfied

	def get_violated(self,df,factors):
		'''
		get violated list for stacked bar plot;inter/intra info
		'''
		Violated=[]
		for i in factors:
			df_1=df[df['Linker']==i[0]]
			df_2=df_1[df_1[i[1]]==1]
			Violated.append(df_2[df_2['Satisfied']==0].shape[0])
		return Violated

	def get_violated_struc(self,df,factors):
		'''
		get violated list for stacked bar plot;struc info
		'''
		struc_dict={'Structured':1,'Unstructured':0,'Intermediate':2}
		Violated=[]
		for i in factors:
			df_1=df[df['Linker']==i[0]]
			df_2=df_1[df_1['Structured']==struc_dict[i[1]]]
			Violated.append(df_2[df_2['Satisfied']==0].shape[0])
		return Violated

	def make_gridplot_intra(self):
		'''
		make gridplot;inter/intra info
		'''
		grid=[]
		for model_id,df in self.model_df.items():
			grid.append(self.plot_intra_summary(df,model_id))
		gridP=gridplot(grid, ncols=len(grid))
		save(gridP,filename=self.filename+'/'+self.ID+"IS.html")
		export_png(gridP,filename=self.filename+'/'+self.ID+"IS.png")


	def plot_intra_summary(self,df,model_id):
		'''
		plot summary stats for inter/intra data
		'''
		factors=self.get_factors(df)
		regions = ['Satisfied', 'Violated']
		source = ColumnDataSource(data=dict(
    				x=factors,
    				Satisfied=self.get_satisfied(df,factors),
    				Violated=self.get_violated(df,factors), ))
		fig_id = figure(x_range=FactorRange(*factors), plot_height=300,plot_width=350,
           			title='Model:'+str(model_id))

		fig_id.vbar_stack(regions, x='x', width=0.9, alpha=0.5, color=["blue", "red"], source=source,
             				legend_label=regions)
		fig_id.xaxis.major_label_text_font_size="12pt"
		fig_id.yaxis.major_label_text_font_size="12pt"
		fig_id.yaxis.axis_label_text_font_size='12pt'
		fig_id.title.text_font_size='12pt'
		fig_id.title.align="center"
		fig_id.title.vertical_align='top'
		fig_id.yaxis.axis_label = 'Number of cross-links'
		#fig_id.y_range.start = 0
		#fig_id.y_range.end = df.shape[0]
		fig_id.x_range.range_padding = 0.1
		fig_id.xaxis.major_label_orientation = 1
		fig_id.xgrid.grid_line_color = None
		fig_id.legend.location = "top_center"
		fig_id.legend.orientation = "horizontal"
		return fig_id

	def make_gridplot_struc(self):
		'''
		plot grid plot for struc info
		'''
		grid=[]
		for model_id,df in self.model_df.items():
			grid.append(self.plot_struc_summary(df,model_id))
		gridP=gridplot(grid, ncols=len(grid))
		save(gridP,filename=self.filename+'/'+self.ID+"SS.html")
		export_png(gridP,filename=self.filename+'/'+self.ID+"SS.png")
		save(gridP,filename=self.filename_add+'/'+self.ID+"SS.html")
		export_png(gridP,filename=self.filename_add+'/'+self.ID+"SS.png")


	def plot_struc_summary(self,df,model_id):
		'''
		plot summary stats: struc/unstruc/intermediate
		'''
		factors=self.get_factors_struc(df)
		regions = ['Satisfied', 'Violated']
		source = ColumnDataSource(data=dict(
    				x=factors,
    				Satisfied=self.get_satisfied_struc(df,factors),
    				Violated=self.get_violated_struc(df,factors), ))
		fig_id = figure(x_range=FactorRange(*factors), plot_height=350, plot_width=400,
           			title='Model:'+str(model_id))

		fig_id.vbar_stack(regions, x='x', width=0.9, alpha=0.5, color=["blue", "red"], source=source,
             				legend_label=regions)
		fig_id.xaxis.major_label_text_font_size="12pt"
		fig_id.yaxis.major_label_text_font_size="12pt"
		fig_id.yaxis.axis_label_text_font_size='12pt'
		fig_id.title.text_font_size='12pt'
		fig_id.title.align="center"
		fig_id.title.vertical_align='top'
		fig_id.yaxis.axis_label = 'Number of cross-links'
		#fig_id.y_range.start = 0
		#fig_id.y_range.end = df.shape[0]
		fig_id.x_range.range_padding = 0.1
		fig_id.xaxis.major_label_orientation = 1
		fig_id.xgrid.grid_line_color = None
		fig_id.legend.location = "top_center"
		fig_id.legend.orientation = "horizontal"
		return fig_id

