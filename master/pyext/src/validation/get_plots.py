###################################
# Script :
# 1) Contains class for plots that
# combines infromation from multiple
# datasets
#
# ganesans - Salilab - UCSF
# ganesans@salilab.org
###################################
import os
import validation
import bokeh
import numpy as np
from bokeh.io import output_file, curdoc, export_png, export_svgs, show
from bokeh.models import ColumnDataSource, Legend, LegendItem, FactorRange
from bokeh.palettes import viridis, Reds256, linear_palette
from bokeh.plotting import figure, save
from bokeh.models.widgets import Tabs, Panel
from bokeh.layouts import row
from bokeh.core.validation import silence
from bokeh.core.validation.warnings import MISSING_RENDERERS, EMPTY_LAYOUT
from bokeh.transform import factor_cmap
silence(MISSING_RENDERERS, True)
silence(EMPTY_LAYOUT, True)


class Plots(validation.GetInputInformation):
    def __init__(self, mmcif, imageDirName):
        super().__init__(mmcif)
        self.ID = str(validation.GetInputInformation.get_id(self))
        self.dirname = os.path.dirname(os.path.abspath(__file__))
        self.imageDirName = imageDirName
        self.filename = os.path.join(self.imageDirName)

    def plot_quality_at_glance(self, molprobity_data: dict, exv_data: dict,
                               sas_data: dict, sas_fit: dict, cx_fit: dict) -> bokeh.plotting.figure:

        # create tabs list to add all the panel figures (model quality, data quality.. etc)
        output_file(self.ID+"quality_at_glance.html", mode="inline")

        # MODEL QUALITY
        # check for molprobity or excluded volume data
        if molprobity_data:
            # if molprobity data, plot that
            Models = molprobity_data['Names']
            Scores = ['Clashscore', 'Ramachandran outliers',
                      'Sidechain outliers']
            data = {'models': Models,
                    'Clashscore': molprobity_data['Clashscore'],
                    'Ramachandran outliers': molprobity_data['Ramachandran outliers'],
                    'Sidechain outliers': molprobity_data['Sidechain outliers']}

            y = [(model, score) for model in Models for score in Scores]
            counts = sum(zip(data['Clashscore'], data['Ramachandran outliers'],
                         data['Sidechain outliers']), ())  # like an hstack
            source = ColumnDataSource(data=dict(y=y, counts=counts))

            if len(molprobity_data['Names']) > 7:
                p = figure(y_range=FactorRange(*y), plot_height=1200,
                           plot_width=700, title="Model Quality: Molprobity Analysis")
            else:

                p = figure(y_range=FactorRange(*y), plot_height=450,
                           plot_width=800, title="Model Quality: Molprobity Analysis")

            p.hbar(y='y', right='counts', width=0.9, source=source, line_color="white",
                   fill_color=factor_cmap('y', palette=viridis(len(Scores)), factors=Scores, start=1, end=2))
            p.xaxis.major_label_text_font_size = "12pt"
            p.yaxis.major_label_text_font_size = "12pt"
            p.xaxis.axis_label = 'Outliers'
            p.xaxis.axis_label_text_font_style = 'italic'
            p.left[0].group_text_font_size = '14px'
            p.left[0].group_label_orientation = 'horizontal'
            p.title.vertical_align = 'top'
            p.title.align = "center"

        elif exv_data:
            model = exv_data['Models']
            satisfaction = exv_data['Number of violations']

            try:
                counts = [float(i) for i in satisfaction]
            except (ValueError):
                return
            violations = exv_data['Excluded Volume Satisfaction']
            Scores = ['Model ' + str(i+1) for i, j in enumerate(model)]
            legends = ['Model ' + str(i+1) + ': ' + str(int(j)) +
                       '('+str(violations[i])+' %)' for i, j in enumerate(counts)]
            n = 3 if len(model) < 3 else len(model)
            source = ColumnDataSource(
                data=dict(Scores=Scores, counts=counts, legends=legends, color=viridis(n)))
            if len(model) > 1:
                upper = max(counts) + (max(counts)-min(counts))
                lower = min(counts)-(max(counts)-min(counts))
            else:
                upper = max(counts) + 20
                lower = min(counts)-20
            if upper == lower:
                upper = lower+20

            p = figure(y_range=Scores, x_range=(lower, upper), plot_height=450,
                       plot_width=800, title='Model Quality: Excluded Volume Analysis')
            r = p.hbar(y='Scores', right='counts', color='color', height=0.5,
                       source=source, alpha=0.8, line_color='black')
            p.xaxis.axis_label = 'Number of violations'
            legend = Legend(items=[LegendItem(label=legends[i], renderers=[
                            r], index=i) for i in range(len(legends))], location='center',
                            label_text_font_size='12px', orientation='vertical')
            p.add_layout(legend, 'right')
            p.xaxis.major_label_text_font_size = "12pt"
            p.yaxis.major_label_text_font_size = "12pt"
            p.title.vertical_align = 'top'
            p.title.align = "center"

        else:
            Scores = ['']
            counts = ['']
            legends = ['']
            source = ColumnDataSource(
                data=dict(Scores=Scores, counts=counts, legends=legends))
            p = figure(y_range=Scores, x_range=(0, 1),
                       plot_height=300, plot_width=800)
        # if neither: ie. atomic but unable to calculate molprobity, return empty plot
        Scores = ['']
        counts = ['']
        legends = ['']
        source = ColumnDataSource(
            data=dict(Scores=Scores, counts=counts, legends=legends))
        p1 = figure(y_range=Scores, x_range=(0, 1),
                    plot_height=300, plot_width=700, title='Under Development')
        p.ygrid.grid_line_color = p1.ygrid.grid_line_color = None
        p.xaxis.axis_label_text_font_size = p1.xaxis.axis_label_text_font_size = "14pt"
        p.yaxis.axis_label_text_font_size = p1.yaxis.axis_label_text_font_size = "14pt"
        p.title.text_font_size = p1.title.text_font_size = '14pt'
        p.title.align = p1.title.align = "center"
        p.output_backend = p1.output_backend = "svg"
        p.title.vertical_align = 'top'
        # make panel figures
        # first panel is model quality
        export_svgs(p, filename=self.filename+'/' +
                    self.ID+"quality_at_glance_MQ.svg")
        export_png(p, filename=self.filename+'/' +
                   self.ID+"quality_at_glance_MQ.png")
        save(p, filename=self.filename+'/'+self.ID+"quality_at_glance_MQ.html")

        # DATA QUALITY
        # check for sas data, if exists, plot
        # this section will be updated with more data assessments, as and when it is complete
        if len(sas_data.keys()) > 0:
            Rgl = {0: 'P(r)', 1: 'Guinier'}
            Scores = [Rgl[m] + ' ('+i+')' for i, j in sas_data.items()
                      for m, n in enumerate(j)]
            counts = [float(n)for i, j in sas_data.items()
                      for m, n in enumerate(j)]
            legends = [str(i)+' nm' for i in counts]
            source = ColumnDataSource(data=dict(
                Scores=Scores, counts=counts, legends=legends, color=viridis(len(legends))))
            pd = figure(y_range=Scores, x_range=(0, max(
                counts)+1), plot_height=450, plot_width=800, title="Data Quality for SAS: Rg Analysis",)
            rd = pd.hbar(y='Scores', right='counts', color='color', height=0.5,
                         source=source, alpha=0.8, line_color='black')
            pd.ygrid.grid_line_color = None
            pd.xaxis.axis_label = 'Distance (nm)'
            pd.xaxis.major_label_text_font_size = "12pt"
            pd.yaxis.major_label_text_font_size = "12pt"
            pd.title.text_font_size = '14pt'
            legend = Legend(items=[LegendItem(label=legends[i], renderers=[
                            rd], index=i) for i in range(len(legends))], location='center',
                            orientation='vertical', label_text_font_size="12px")
            pd.add_layout(legend, 'right')
            pd.legend.label_text_font_size = "12px"
            pd.xaxis.axis_label_text_font_style = 'italic'
            pd.yaxis.axis_label_text_font_style = 'italic'
            pd.xaxis.axis_label_text_font_size = "14pt"
            pd.yaxis.major_label_text_font_size = "14pt"
            pd.title.vertical_align = 'top'
            pd.title.align = "center"
            pd.output_backend = "svg"
            export_svgs(pd, filename=self.filename+'/' +
                        self.ID+"quality_at_glance_DQ.svg")
            export_png(pd, filename=self.filename+'/' +
                       self.ID+"quality_at_glance_DQ.png")
            save(pd, filename=self.filename+'/' +
                 self.ID+"quality_at_glance_DQ.html")

        # FIT TO DATA QUALITY
        # check for sas data, if exists, plot
        # this section will be updated with more data assessments, as and when it is complete
        if len(sas_fit.keys()) > 0:
            Scores = [' \u03C7\u00b2 Fit ' +
                      str(int(m+1)) + ' ('+i+')' for i, j in sas_fit.items() for m, n in enumerate(j)]
            counts = [float(n) for i, j in sas_fit.items()
                      for m, n in enumerate(j)]
            legends = [str(i) for i in counts]
            source = ColumnDataSource(data=dict(
                Scores=Scores, counts=counts, legends=legends, color=viridis(len(legends))))
            pf = figure(y_range=Scores, x_range=(0, max(counts)+1), plot_height=450,
                        plot_width=800, title="Fit to SAS Data:  \u03C7\u00b2 Fit")
            rf = pf.hbar(y='Scores', right='counts', color='color', height=0.5,
                         source=source, alpha=0.8, line_color='black')
            pf.ygrid.grid_line_color = None
            pf.title.text_font_size = '14pt'
            pf.xaxis.axis_label = 'Fit value'
            pf.xaxis.major_label_text_font_size = "12pt"
            pf.yaxis.major_label_text_font_size = "12pt"

            legend = Legend(items=[LegendItem(label=legends[i], renderers=[
                            rf], index=i) for i in range(len(legends))], location="center",
                            orientation='vertical', label_text_font_size="12px")
            pf.add_layout(legend, 'right')
            pf.title.vertical_align = 'top'
            pf.title.align = "center"

            pf.output_backend = "svg"
            pf.legend.label_text_font_size = "12px"
            pf.xaxis.axis_label_text_font_style = 'italic'
            pf.yaxis.axis_label_text_font_style = 'italic'
            pf.xaxis.axis_label_text_font_size = "14pt"
            pf.yaxis.major_label_text_font_size = "14pt"

            pf.title.vertical_align = 'top'
            pf.title.align = "center"
            pf.output_backend = "svg"

            export_svgs(pf, filename=self.filename+'/' +
                        self.ID+"quality_at_glance_FQ.svg")
            export_png(pf, filename=self.filename+'/' +
                       self.ID+"quality_at_glance_FQ.png")
            save(pf, filename=self.filename+'/' +
                 self.ID+"quality_at_glance_FQ.html")

        # check for XL_MS data, if exists, plot
        if len(cx_fit.keys()) > 0:
            Scores = ['Model '+str(i) for i, j in cx_fit.items()]
            counts = [round(float(j), 2) for i, j in cx_fit.items()]
            # legends=[str(i) for i in counts]
            legends = ['Model ' + str(i+1) + ': ' +
                       str(j)+'%' for i, j in enumerate(counts)]
            source = ColumnDataSource(data=dict(
                Scores=Scores, counts=counts, legends=legends, color=viridis(len(legends))))
            pf1 = figure(y_range=Scores, x_range=(0, max(counts)+1),
                         plot_height=450, plot_width=800, title="Fit to XL-MS Input")
            rf1 = pf1.hbar(y='Scores', right='counts', color='color',
                           source=source, alpha=0.8, line_color='black')
            pf1.ygrid.grid_line_color = None
            pf1.xaxis.major_label_text_font_size = "12pt"
            pf1.yaxis.major_label_text_font_size = "12pt"
            pf1.title.text_font_size = '12pt'
            pf1.title.align = "center"
            pf1.title.vertical_align = 'top'

            legend = Legend(items=[LegendItem(label=legends[i], renderers=[
                            rf1], index=i) for i in range(len(legends))], location='center',
                            orientation='vertical', label_text_font_size="12px")
            pf1.add_layout(legend, 'right')
            pf1.output_backend = "svg"
            export_svgs(pf1, filename=self.filename+'/' +
                        self.ID+"quality_at_glance_FQ1.svg")
            export_png(pf1, filename=self.filename+'/' +
                       self.ID+"quality_at_glance_FQ1.png")
            save(pf1, filename=self.filename+'/' +
                 self.ID+"quality_at_glance_FQ1.html")
