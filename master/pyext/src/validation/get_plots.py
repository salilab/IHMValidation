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
from bokeh.io import output_file, curdoc, export_png, export_svgs
from bokeh.models import ColumnDataSource, Legend, LegendItem
from bokeh.palettes import viridis
from bokeh.plotting import figure, save
from bokeh.models.widgets import Tabs, Panel
from bokeh.layouts import row
from bokeh.core.validation import silence
from bokeh.core.validation.warnings import MISSING_RENDERERS, EMPTY_LAYOUT
silence(MISSING_RENDERERS, True)
silence(EMPTY_LAYOUT, True)


class Plots(validation.GetInputInformation):
    def __init__(self, mmcif, imageDirName):
        super().__init__(mmcif)
        self.ID = str(validation.GetInputInformation.get_id(self))
        self.dirname = os.path.dirname(os.path.abspath(__file__))
        self.imageDirName = imageDirName
        self.filename = os.path.join(self.imageDirName)

    def plot_quality_at_glance(self, clashscore: list, rama: list, sidechain: list,
                               exv_data: dict, sas_data: dict, sas_fit: dict, cx_fit: dict) -> bokeh.plotting.figure:

        # create tabs list to add all the panel figures (model quality, data quality.. etc)
        output_file(self.ID+"quality_at_glance.html", mode="inline")
        tabsI = []

        # MODEL QUALITY
        # check for molprobity or excluded volume data
        if clashscore or rama or sidechain:
            # if molprobity data, plot that
            counts = [clashscore, rama, sidechain]
            Scores = ['Clashscore', 'Ramachandran Outliers',
                      'Sidechain Outliers']
            legends = [str(counts[0]), str(counts[1])+'%', str(counts[2])+'%']
            source = ColumnDataSource(data=dict(
                Scores=Scores, counts=counts, legends=legends, color=viridis(len(Scores))))
            p = figure(y_range=Scores, x_range=(0, max(counts)+1), plot_height=300,
                       plot_width=800, title="Model Quality: Molprobity Analysis")
            r = p.hbar(y='Scores', right='counts', height=0.5,
                       color='color', source=source, alpha=0.8, line_color='black')
            legend = Legend(items=[LegendItem(label=legends[i], renderers=[
                            r], index=i) for i in range(len(legends))])
            p.add_layout(legend, 'below')
            p.legend.orientation = "horizontal"
            p.legend.location = 'bottom_center'
            p.legend.label_text_font_size = "8px"
            p.xaxis.major_label_text_font_size = "14pt"
            p.yaxis.major_label_text_font_size = "14pt"

        elif exv_data:
            # if excluded vol data, plot that
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

            p = figure(y_range=Scores, x_range=(lower, upper), plot_height=300,
                       plot_width=800, title='Model quality: Excluded Volume Analysis')
            r = p.hbar(y='Scores', right='counts', color='color',
                       source=source, alpha=0.8, line_color='black')
            p.xaxis.axis_label = 'Number of violations'
            legend = Legend(items=[LegendItem(label=legends[i], renderers=[
                            r], index=i) for i in range(len(legends))])
            p.add_layout(legend, 'below')
            p.legend.orientation = "horizontal"
            p.legend.location = 'bottom_center'
            p.legend.label_text_font_size = "8px"
            p.xaxis.major_label_text_font_size = "14pt"
            p.yaxis.major_label_text_font_size = "14pt"

        else:
            Scores = ['']
            counts = ['']
            legends = ['']
            source = ColumnDataSource(
                data=dict(Scores=Scores, counts=counts, legends=legends))
            p = figure(y_range=Scores, x_range=(0, 1),
                       plot_height=300, plot_width=700)
        # if neither: ie. atomic but unable to calculate molprobity, return empty plot
        Scores = ['']
        counts = ['']
        legends = ['']
        source = ColumnDataSource(
            data=dict(Scores=Scores, counts=counts, legends=legends))
        p1 = figure(y_range=Scores, x_range=(0, 1),
                    plot_height=300, plot_width=700,title='Under Development')
        p.ygrid.grid_line_color = p1.ygrid.grid_line_color = None
        p.xaxis.axis_label_text_font_size = p1.xaxis.axis_label_text_font_size = "14pt"
        p.yaxis.axis_label_text_font_size = p1.yaxis.axis_label_text_font_size = "14pt"
        p.title.text_font_size = p1.title.text_font_size = '14pt'
        p.title.align = p1.title.align = "center"
        p.output_backend = p1.output_backend = "svg"
        p.title.vertical_align = 'top'
        # make panel figures
        # first panel is model quality
        first = Panel(child=row(p), title='Model quality')
        tabsI.append(first)
        export_svgs(p, filename=self.filename+'/' +
                    self.ID+"quality_at_glance_MQ.svg")
        export_png(p, filename=self.filename+'/' +
                   self.ID+"quality_at_glance_MQ.png")

        # DATA QUALITY
        # check for sas data, if exists, plot
        # this section will be updated with more data assessments, as and when it is complete
        if len(sas_data.keys()) > 0:
            Rgl = {0: 'p(r)', 1: 'Guinier'}
            Scores = [Rgl[m] + ' ('+i+')' for i, j in sas_data.items()
                      for m, n in enumerate(j)]
            counts = [float(n)for i, j in sas_data.items()
                      for m, n in enumerate(j)]
            legends = [str(i)+' nm' for i in counts]
            source = ColumnDataSource(data=dict(
                Scores=Scores, counts=counts, legends=legends, color=viridis(len(legends))))
            pd = figure(y_range=Scores, x_range=(0, max(
                counts)+1), plot_height=300, plot_width=800, title="Data Quality: Rg Analysis",)
            rd = pd.hbar(y='Scores', right='counts', color='color',
                         source=source, alpha=0.8, line_color='black')
            pd.ygrid.grid_line_color = None
            pd.xaxis.major_label_text_font_size = "14pt"
            pd.yaxis.major_label_text_font_size = "14pt"
            pd.title.text_font_size = '14pt'
            pd.title.align = "center"
            legend = Legend(items=[LegendItem(label=legends[i], renderers=[
                            rd], index=i) for i in range(len(legends))])
            pd.add_layout(legend, 'below')
            pd.legend.orientation = "horizontal"
            pd.legend.location = 'bottom_center'
            pd.legend.label_text_font_size = "8px"
            pd.xaxis.axis_label_text_font_style = 'italic'
            pd.yaxis.axis_label_text_font_style = 'italic'
            second = Panel(child=row(pd), title='Data quality: SAS')
            tabsI.append(second)
            pd.output_backend = "svg"
            export_svgs(pd, filename=self.filename+'/' +
                        self.ID+"quality_at_glance_DQ.svg")
            export_png(pd, filename=self.filename+'/' +
                       self.ID+"quality_at_glance_DQ.png")
        else:
            second = Panel(child=row(p1), title='Data quality')
            tabsI.append(second)

        # FIT TO DATA QUALITY
        # check for sas data, if exists, plot
        # this section will be updated with more data assessments, as and when it is complete
        if len(sas_fit.keys()) > 0:
            Scores = ['\uab53 \u00b2 value fit ' +
                      str(int(m+1)) + ' ('+i+')' for i, j in sas_fit.items() for m, n in enumerate(j)]
            counts = [float(n) for i, j in sas_fit.items()
                      for m, n in enumerate(j)]
            legends = [str(i) for i in counts]
            source = ColumnDataSource(data=dict(
                Scores=Scores, counts=counts, legends=legends, color=viridis(len(legends))))
            pf = figure(y_range=Scores, x_range=(0, max(counts)+1), plot_height=300,
                        plot_width=800, title="Fit to SAS input: \uab53 \u00b2 Fit")
            rf = pf.hbar(y='Scores', right='counts', color='color',
                         source=source, alpha=0.8, line_color='black')
            pf.ygrid.grid_line_color = None
            pf.xaxis.major_label_text_font_size = "14pt"
            pf.yaxis.major_label_text_font_size = "14pt"
            pf.title.text_font_size = '14pt'
            pf.title.align = "center"
            legend = Legend(items=[LegendItem(label=legends[i], renderers=[
                            rf], index=i) for i in range(len(legends))])
            pf.add_layout(legend, 'below')
            pf.legend.orientation = "horizontal"
            pf.legend.location = 'bottom_center'
            pf.legend.label_text_font_size = "8px"
            third = Panel(child=row(pf), title='Fit to data: SAS')
            tabsI.append(third)
            pf.output_backend = "svg"
            export_svgs(pf, filename=self.filename+'/' +
                        self.ID+"quality_at_glance_FQ.svg")
            export_png(pf, filename=self.filename+'/' +
                       self.ID+"quality_at_glance_FQ.png")

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
                         plot_height=300, plot_width=800, title="Fit to XL-MS input")
            rf1 = pf1.hbar(y='Scores', right='counts', color='color',
                           source=source, alpha=0.8, line_color='black')
            pf1.ygrid.grid_line_color = None
            pf1.xaxis.major_label_text_font_size = "14pt"
            pf1.yaxis.major_label_text_font_size = "14pt"
            pf1.title.text_font_size = '14pt'
            pf1.title.align = "center"
            legend = Legend(items=[LegendItem(label=legends[i], renderers=[
                            rf1], index=i) for i in range(len(legends))])
            pf1.add_layout(legend, 'below')
            pf1.legend.orientation = "horizontal"
            pf1.legend.label_text_font_size = "8px"
            pf1.legend.location = 'bottom_center'
            third1 = Panel(child=row(pf1), title='Fit to data: XL-MS')
            tabsI.append(third1)
            pf1.output_backend = "svg"
            export_svgs(pf1, filename=self.filename+'/' +
                        self.ID+"quality_at_glance_FQ1.svg")
            export_png(pf1, filename=self.filename+'/' +
                       self.ID+"quality_at_glance_FQ1.png")

        # check for XL_MS or SAS data, if exists, plot third panel as FIT TO INPUT
        if (len(cx_fit.keys()) < 1) or (len(sas_fit.keys()) < 1):
            third = Panel(child=row(p1), title='Fit to input data')
            tabsI.append(third)

        # remaining panels are empty for now, will be updated
        # fourth = Panel(
        #   child=row(p1), title='Fit to data not used for modeling')
        # tabsI.append(fourth)
        fifth = Panel(child=row(p1), title='Uncertainty in models')
        tabsI.append(fifth)
        tabs = Tabs(tabs=tabsI)
        curdoc().add_root(tabs)
        save(tabs, filename=self.filename+'/'+self.ID+"quality_at_glance.html")
