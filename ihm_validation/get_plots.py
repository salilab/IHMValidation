# -*- coding: utf-8 -*-
#
# get_plots.py.py - Generate overview plots
#
# Copyright (C) 2019-2025 Arthur Zalevsky, Sai Ganesan, Benjamin M. Webb, Brinda Vallat
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
Generate overview plots
"""

import os
from pathlib import Path
import utility
from mmcif_io import GetInputInformation
import bokeh
import numpy as np
from bokeh.io import output_file, curdoc, export_svg, show
from bokeh.models import (ColumnDataSource, Legend, LegendItem, FactorRange,
                          Div, BasicTickFormatter)
from bokeh.palettes import linear_palette, Greys256, Blues256, Oranges256, Greens256
from bokeh.plotting import figure, save
from bokeh.layouts import row
from bokeh.core.validation import silence
from bokeh.core.validation.warnings import MISSING_RENDERERS, EMPTY_LAYOUT
from bokeh.transform import factor_cmap
from bokeh.layouts import gridplot, column
silence(MISSING_RENDERERS, True)
silence(EMPTY_LAYOUT, True)

import json
from bokeh.embed import json_item

MAXPLOTS = 256

# Total canvas width of the "quality at a glance" plots.
#
# These figures are sized by their total width, so everything bokeh puts around
# the data area - axis labels, the categorical group labels and the legend -
# eats into the plot itself. Bokeh 3 reserves more room for that chrome than
# bokeh 2 did, which left the plots cramped; the crosslinking-MS data quality
# plot lost the most because its category labels are long (its data area went
# from 192 to 119 px). Keep this comfortably below the printable width of the
# letter-sized PDF report (540 pt, i.e. about 910 px); note that the model
# quality plots are wrapped in a column layout that adds another 30 px.
GLANCE_PLOT_WIDTH = 840

# Height of the data area per horizontal bar. Set as a frame (data area) height
# instead of being folded into the total height, so the bars keep their
# thickness regardless of how much vertical space the title and axes take.
# These plots are meant to be tight, so this is deliberately only a little more
# than the height of the 14pt category labels next to them (about 19 px).
GLANCE_BAR_HEIGHT = 21

class Plots(GetInputInformation):
    def __init__(self, *args, imageDirName, **kwargs):
        super().__init__(*args, **kwargs)

        self.imageDirName = imageDirName

    def plot_quality_at_glance(self,
                               molprobity_data: dict=None, exv_data: dict=None,
                               sas_data_quality: dict=None, sas_fit: dict=None,
                               cx_data_quality: dict=None, cx_fit: dict=None,
                               em_data_quality: dict=None, em_fit: dict=None
                               ) -> dict:
        """
        Generate the "quality at a glance" summary plots.

        Produces up to three horizontal bar charts - model quality (MQ), data
        quality (DQ) and fit to data used for modeling (FQ) - saved as SVG and
        JSON, and returns a dict telling the report templates which of them
        exist.

        All of these figures are sized by their total width
        (:data:`GLANCE_PLOT_WIDTH`), so the axis labels, categorical group
        labels and legend are subtracted from the data area rather than added
        around it. Their data area is therefore only as wide as that width
        minus whatever bokeh decides the labels need; the bar thickness, on the
        other hand, is pinned by :data:`GLANCE_BAR_HEIGHT`. See the comments on
        both constants before changing either.
        """

        glance_plots = {
            'MQ': False,
            'DQ': False,
            'FQ': False,
        }

        # MODEL QUALITY
        mq_plots = []

        # check for MolProbity or excluded volume data
        if molprobity_data is not None:
            # if MolProbity data, plot that
            # every model has clashscore, rama outliers, and rota outliers
            Models = list(molprobity_data.keys())
            Scores = [
                'Sidechain outliers',
                'Ramachandran outliers',
                'Clashscore',
            ]
            data = {'models': Models}

            for s in Scores:
                data[s] = [x[s] for x in molprobity_data.values()]

            y = [(f"Model {model}", score) for model in Models for score in Scores]
            counts = sum(
                zip(
                    data['Sidechain outliers'],
                    data['Ramachandran outliers'],
                    data['Clashscore'],
                ),
                ()
            )
            source = ColumnDataSource(data=dict(y=y, counts=counts))

            # if there are more than 7 models, we will increase the size of the plots
            # this is important, else the plots look ugly

            plots = []

            # get data ranges
            lower, upper = utility.calc_optimal_range(counts)

            # create plot
            for i, name_ in enumerate(Models):

                p = figure(
                    y_range=FactorRange(*y[i * 3: (i + 1) * 3]),
                    # Force left limit at zero
                    x_range=(0, upper),
                    frame_height=3 * GLANCE_BAR_HEIGHT,
                    width=GLANCE_PLOT_WIDTH
                )

                p_ = p.hbar(y=source.data['y'][i * 3: (i + 1) * 3],
                       right=source.data['counts'][i * 3: (i + 1) * 3],
                       line_color="white",
                            fill_color=factor_cmap('y', palette=linear_palette(Greys256, len(Scores) + 2)[1:-1],
                                              factors=Scores,
                                              start=1, end=2)
                       )

                # set labels and fonts
                p.xaxis.axis_label_text_font_size = "14pt"
                p.yaxis.axis_label_text_font_size = "14pt"
                p.xaxis.major_label_text_font_size = "14pt"
                p.yaxis.major_label_text_font_size = "14pt"
                p.xaxis.axis_label_text_font_style = 'normal'
                p.yaxis.axis_label_text_font_style = 'normal'
                p.xaxis.axis_label = 'Outliers'
                p.left[0].group_text_font_size = '14pt'
                p.left[0].group_text_color = p.left[0].major_label_text_color
                p.left[0].group_text_font_style = p.left[0].major_label_text_font_style
                p.left[0].group_text_font_size = p.left[0].major_label_text_font_size
                p.left[0].group_label_orientation = 'horizontal'
                p.title.vertical_align = 'top'
                p.title.align = "center"
                p.output_backend = "svg"
                p.border_fill_color = None
                p.background_fill_color = None

                fname = Path(self.imageDirName, f"{self.ID_f}_{name_}_quality_at_glance_MQ_mp.svg")
                export_svg(p, filename=fname)
                utility.strip_bokeh_invisible_outline(fname)

                plots.append(p)

            grid = gridplot(plots, ncols=1,
                            merge_tools=True,
                            toolbar_location='right')
            # grid.children[1].css_classes = ['scrollable']
            # grid.children[1].sizing_mode = 'fixed'
            # grid.children[1].height = 450
            # grid.children[1].width = 800

            title = Div(text="<p>Model Quality: Molprobity Analysis</p>",
                        styles={"font-size": "1.5em", "font-weight": "bold",
                               "text-align": "center", "width": '100%'}, width=GLANCE_PLOT_WIDTH
                        )

            fullplot = column(title, grid)
            mq_plots.append(fullplot)

        # If excluded volume data is available
        if exv_data:
            numplots = min(len(exv_data['Model ID']), MAXPLOTS)

            Models = exv_data['Model ID'][:numplots]
            # analysed = exv_data['Analysed'][:numplots]
            # violations = exv_data['Number of violations'][:numplots]
            # satisfaction = exv_data['Number of violations'][:numplots]
            satisfaction = exv_data['Excluded Volume Satisfaction (%)'][:numplots]
            # make sure data is plot-able
            try:
                counts = [float(i) for i in satisfaction]
            except (ValueError):
                return

            Scores = [f'Model {m}' for m in Models]
            legends = [f'{s:.2f} %' for s in satisfaction]

            # set the size of the axis
            # n = 3 if len(model) < 3 else len(model)
            n = len(counts)
            color = linear_palette(Greys256, 5)[1] # match with darkest molprobity color
            colors = [color for x in range(n)]
            source = ColumnDataSource(
                data=dict(Scores=Scores, counts=counts, legends=legends, color=colors))

            #  build plots
            plots = []

            # get ranges
            lower, upper = 0, 102

            for i, name_ in enumerate(Models):
                p = figure(y_range=FactorRange(factors=source.data['Scores'][i: i + 1]), x_range=(lower, upper), frame_height=GLANCE_BAR_HEIGHT,
                           width=GLANCE_PLOT_WIDTH)  # , title='Model Quality: Excluded Volume Analysis')
                # p.xaxis.formatter = BasicTickFormatter(use_scientific=True, power_limit_high=3)
                p.xaxis.ticker.desired_num_ticks = 3

                p_ = p.hbar(y=source.data['Scores'][i:i + 1], right=source.data['counts'][i: i + 1], color=source.data['color'][i:i + 1], height=1.0,
                           line_color='white')

                p.xaxis.axis_label = 'Satisfaction rate [%]'
                legend = Legend(items=[LegendItem(label=legends[i:i + 1][j], renderers=[
                    p_], index=j) for j in range(len(legends[i:i + 1]))], location='center',
                    label_text_font_size='14pt', orientation='vertical')
                p.add_layout(legend, 'right')
                p.legend.border_line_width = 0
                p.xaxis.major_label_text_font_size = "14pt"
                p.yaxis.major_label_text_font_size = "14pt"
                p.xaxis.axis_label_text_font_style = 'normal'
                p.yaxis.axis_label_text_font_style = 'normal'
                p.title.vertical_align = 'top'
                p.title.align = "center"
                p.output_backend = "svg"
                p.border_fill_color = None
                p.background_fill_color = None
                p.min_border_top = 20

                fname = Path(self.imageDirName, f"{self.ID_f}_{name_}_quality_at_glance_MQ_exv.svg")
                export_svg(p, filename=fname)
                utility.strip_bokeh_invisible_outline(fname)

                plots.append(p)

            grid = gridplot(plots, ncols=1,
                            merge_tools=True,
                            toolbar_location='right')
            # grid.children[1].css_classes = ['scrollable']
            # grid.children[1].sizing_mode = 'fixed'
            # grid.children[1].height = 450
            # grid.children[1].width = 800

            title = Div(text='<p>Model Quality: Excluded Volume Analysis</p>',
                        styles={"font-size": "1.5em", "font-weight": "bold",
                               "text-align": "center", "width": '100%'}, width=GLANCE_PLOT_WIDTH
                        )

            fullplot = column(title, grid)
            mq_plots.append(fullplot)

        # if neither exc vol nor molp data exists, we create a blank plot
        # pdb-ihm visuals keep changing, so this plot might or might not make sense
        # we are keeping it, just in case the visuals change again
        # else:
        #    Scores = ['']
        #    counts = ['']
        #    legends = ['']
        #    source = ColumnDataSource(
        #        data=dict(Scores=Scores, counts=counts, legends=legends))
        #    p = figure(y_range=Scores, x_range=(0, 1),
        #               height=300, width=800)
        #
        #    p.ygrid.grid_line_color = None
        #    p.xaxis.axis_label_text_font_size = "14pt"
        #    p.yaxis.axis_label_text_font_size = "14pt"
        #    p.title.text_font_size = '14pt'
        #    p.title.align = "center"
        #
        #    p.output_backend = "svg"
        #    p.border_fill_color = None
        #    p.background_fill_color = None
        #    p.title.vertical_align = 'top'
        #    fullplot = p
        #
        #    mq_plots.append(fullplot)
        # make panel figures
        # first panel is model quality
#        export_svg(fullplot, filename=self.filename+'/' +
#                    self.ID_f + "_quality_at_glance_MQ.svg")
#        save(fullplot, filename=self.filename+'/' +
#             self.ID_f +"_quality_at_glance_MQ.html")

        if len(mq_plots) > 0:
            pd = gridplot(mq_plots, ncols=1,
                merge_tools=True,
                     toolbar_location="above",
                     # sizing_mode='stretch_width'
                     )
            # pd = column(*dq_plots)

            plots = self.save_plots(pd, 'quality_at_glance_MQ')
            glance_plots['MQ'] = True

        # DATA QUALITY
        # check for sas data, if exists, plot
        # this section will be updated with more data assessments, as and when it is complete
        dq_plots = []

        if len(sas_data_quality.keys()) > 0:
            # Don't forget to update palette if adding new metric
            Rgl = {0: 'P(r)', 1: 'Guinier'}
            Scores = [Rgl[m] + ' ('+i+')' for i, j in sas_data_quality.items()
                      for m, n in enumerate(j)]
            counts = [float(n)for i, j in sas_data_quality.items()
                      for m, n in enumerate(j)]
            legends = [str(i)+' nm' for i in counts]
            # Update palette if we need more colors
            colors_ = linear_palette(Blues256, len(Rgl) + 2)[1:-1]
            colors = [colors_[m] for i, j in sas_data_quality.items()
                      for m, n in enumerate(j)]
            source = ColumnDataSource(data=dict(
                Scores=Scores, counts=counts, legends=legends, color=colors))
            pd = figure(y_range=Scores, x_range=(0, max(
                counts)+1), frame_height=len(counts) * GLANCE_BAR_HEIGHT, width=GLANCE_PLOT_WIDTH, title="Data Quality for SAS: Rg Analysis",)
            rd = pd.hbar(y='Scores', right='counts', color='color', height=1.0,
                         source=source, line_color='white')
            pd.ygrid.grid_line_color = None
            pd.xaxis.axis_label = 'Distance [nm]'
            pd.title.text_font_size = '14pt'
            legend = Legend(items=[LegendItem(label=legends[i], renderers=[
                            rd], index=i) for i in range(len(legends))], location='center',
                            orientation='vertical', label_text_font_size="14pt")
            pd.add_layout(legend, 'right')
            pd.legend.items.reverse()
            pd.legend.border_line_width = 0
            pd.legend.label_text_font_size = "14pt"
            pd.xaxis.axis_label_text_font_size = "14pt"
            pd.yaxis.axis_label_text_font_size = "14pt"
            pd.xaxis.major_label_text_font_size = "14pt"
            pd.yaxis.major_label_text_font_size = "14pt"
            pd.xaxis.axis_label_text_font_style = 'normal'
            pd.yaxis.axis_label_text_font_style = 'normal'
            pd.title.vertical_align = 'top'
            pd.title.align = "center"
            pd.output_backend = "svg"
            pd.border_fill_color = None
            pd.background_fill_color = None

            dq_plots.append(pd)

        # If crosslinking-MS data is available
        if cx_data_quality is not None and len(cx_data_quality) > 0:
            Models = [data["pride_id"] for data in cx_data_quality]
            # Don't forget to update palette if adding new metric
            Scores = ['Total', 'Mapped to matching entities', 'Matched']
            legends = []
            for data in cx_data_quality:
                legends.append(f'{data["stats"]["ms"]["total"]}')
                legends.append(f'{data["stats"]["ms"]["mapped_entities"]} ({data["stats"]["ms"]["mapped_entities_pct"]:.2f}%)')
                legends.append(f'{data["stats"]["ms"]["matched"]} ({data["stats"]["ms"]["matched_pct"]:.2f}%)')

            data = {'models': Models,
                    'Total': [data['stats']['ms']['total'] for data in cx_data_quality],
                    'Mapped to matching entities': [data['stats']['ms']['mapped_entities'] for data in cx_data_quality],
                    'Matched': [data['stats']['ms']['matched'] for data in cx_data_quality]
                    }
            y = [(model, score) for model in Models for score in Scores]
            counts = sum(zip(data['Total'], data['Mapped to matching entities'],
                         data['Matched']), ())
            source = ColumnDataSource(data=dict(y=y, counts=counts, legends=legends))

            # if there are more than 7 models, we will increase the size of the plots
            # this is important, else the plots look ugly

            plots = []

            # get data ranges
            lower, upper = utility.calc_optimal_range(counts)

            # create plot
            for i, name_ in enumerate(Models):

                title_txt = "Crosslinking-MS Data Quality"
                title = Div(text=f"<p>{title_txt}</p>",
                            styles={"font-size": "1.5em", "font-weight": "bold",
                                   "text-align": "center"}
                            )


                p = figure(
                    y_range=FactorRange(*y[i * 3: (i + 1) * 3]),
                    # Force left limit at zero
                    x_range=(lower, upper),
                    frame_height=3 * GLANCE_BAR_HEIGHT,
                    width=GLANCE_PLOT_WIDTH,
                    title=title_txt
                )
                p.xaxis.ticker.desired_num_ticks = 3

                rd = p.hbar(y=source.data['y'][i * 3: (i + 1) * 3],
                       right=source.data['counts'][i * 3: (i + 1) * 3],
                       line_color="white",
                        fill_color=factor_cmap('y', palette=linear_palette(Oranges256, len(Scores) + 2)[1:-1],
                                              factors=Scores,
                                              start=1, end=2)
                       )
                legends_ = source.data['legends'][i * 3: (i + 1) * 3]
                legend = Legend(items=[LegendItem(label=legends_[j], renderers=[
                            rd], index=j) for j in range(len(legends_))], location='center',
                            orientation='vertical', label_text_font_size="14pt")
                legend.items = legend.items[::-1]
                p.add_layout(legend, 'right')
                p.legend.border_line_width = 0
                # set labels and fonts
                p.xaxis.axis_label_text_font_size = "14pt"
                p.yaxis.axis_label_text_font_size = "14pt"
                p.xaxis.major_label_text_font_size = "14pt"
                p.yaxis.major_label_text_font_size = "14pt"
                p.xaxis.axis_label_text_font_style = 'normal'
                p.yaxis.axis_label_text_font_style = 'normal'
                p.yaxis.major_label_text_align='right'
                p.xaxis.axis_label = 'Residue pairs'
                p.left[0].group_text_font_size = '14px'
                p.left[0].group_label_orientation = 'horizontal'
                p.title.text_font_size = '14pt'
                p.title.vertical_align = 'top'
                p.title.align = "center"
                p.output_backend = "svg"
                p.border_fill_color = None
                p.background_fill_color = None
                plots.append(p)

            grid = gridplot(plots, ncols=1,
                            merge_tools=True,
                            toolbar_location=None,
                            )
#            grid.children[1].css_classes = ['scrollable']
#            grid.children[1].sizing_mode = 'fixed'
#            grid.children[1].height = 450
#            grid.children[1].width = 800


            # fullplot = column(title, grid)

            dq_plots.append(grid)

        if em_data_quality is not None and len(em_data_quality) > 0:
            Scores = []
            counts = []
            for dataset in em_data_quality:
                emdbid = dataset['emdbid']
                try:
                    s = float(dataset['data_stats']['resolution'])
                except (ValueError, TypeError):
                    continue
                Scores.append(emdbid)
                counts.append(s)

            if len(counts) > 0:
                legends = [f'{i:.2f} Å' for i in counts]
                source = ColumnDataSource(data=dict(
                    Scores=Scores, counts=counts, legends=legends, color=linear_palette(Greens256, len(legends) + 2)[1:-1]))
                pf = figure(y_range=Scores, x_range=(0, 80), frame_height=len(counts) * GLANCE_BAR_HEIGHT,
                            width=GLANCE_PLOT_WIDTH, title="3DEM resolution")
                rf = pf.hbar(y='Scores', right='counts', color='color', height=1.0,
                             source=source, line_color='white')
                pf.ygrid.grid_line_color = None
                pf.title.text_font_size = '14pt'
                pf.xaxis.axis_label = 'Resolution [Å]'
                legend = Legend(items=[LegendItem(label=legends[i], renderers=[
                                rf], index=i) for i in range(len(legends))], location="center",
                                orientation='vertical', label_text_font_size="14pt")
                pf.add_layout(legend, 'right')
                pf.legend.items.reverse()
                pf.legend.border_line_width = 0
                pf.title.vertical_align = 'top'
                pf.title.align = "center"
                pf.output_backend = "svg"
                pf.border_fill_color = None
                pf.background_fill_color = None
                pf.legend.label_text_font_size = "14pt"
                pf.xaxis.axis_label_text_font_size = "14pt"
                pf.yaxis.axis_label_text_font_size = "14pt"
                pf.xaxis.major_label_text_font_size = "14pt"
                pf.yaxis.major_label_text_font_size = "14pt"
                pf.xaxis.axis_label_text_font_style = 'normal'
                pf.yaxis.axis_label_text_font_style = 'normal'
                dq_plots.append(pf)

        if len(dq_plots) > 0:
            pd = gridplot(dq_plots, ncols=1,
                merge_tools=True,
                     toolbar_location="above",
                     # sizing_mode='stretch_width'
                     )
            # pd = column(*dq_plots)

            plots = self.save_plots(pd, 'quality_at_glance_DQ')
            glance_plots['DQ'] = True

        # FIT TO DATA QUALITY
        # check for sas data, if exists, plot
        # this section will be updated with more data assessments, as and when it is complete
        fq_plots = []
        if len(sas_fit.keys()) > 0:
            Scores = [' \u03C7\u00b2 Fit ' +
                      str(int(m+1)) + ' ('+i+')' for i, j in sas_fit.items() for m, n in enumerate(j)]
            counts = [float(n) for i, j in sas_fit.items()
                      for m, n in enumerate(j)]
            legends = [str(i) for i in counts]
            source = ColumnDataSource(data=dict(
                Scores=Scores, counts=counts, legends=legends, color=linear_palette(Blues256, len(legends) + 2)[1:-1]))
            pf = figure(y_range=Scores, x_range=(0, max(counts)+1), frame_height=len(counts) * GLANCE_BAR_HEIGHT,
                        width=GLANCE_PLOT_WIDTH, title="Fit to SAS Data:  \u03C7\u00b2 Fit")
            rf = pf.hbar(y='Scores', right='counts', color='color', height=1.0,
                         source=source, line_color='white')
            pf.ygrid.grid_line_color = None
            pf.xaxis.axis_label = 'Fit value'
            legend = Legend(items=[LegendItem(label=legends[i], renderers=[
                            rf], index=i) for i in range(len(legends))], location="center",
                            orientation='vertical', label_text_font_size="14pt")
            pf.add_layout(legend, 'right')
            pf.legend.items.reverse()
            pf.legend.border_line_width = 0
            pf.title.vertical_align = 'top'
            pf.title.align = "center"
            pf.output_backend = "svg"
            pf.border_fill_color = None
            pf.background_fill_color = None
            pf.title.text_font_size = '14pt'
            pf.legend.label_text_font_size = "14pt"
            pf.xaxis.axis_label_text_font_size = "14pt"
            pf.yaxis.axis_label_text_font_size = "14pt"
            pf.xaxis.major_label_text_font_size = "14pt"
            pf.yaxis.major_label_text_font_size = "14pt"
            pf.xaxis.axis_label_text_font_style = 'normal'
            pf.yaxis.axis_label_text_font_style = 'normal'
            fq_plots.append(pf)

        if cx_fit is not None and len(cx_fit) > 0:
            Scores = []
            counts = []
            i = 0
            for sg, sgv in cx_fit.items():
                for st, stv in sgv.items():
                    for mg, mgv in stv.items():
                        i += 1
                        try:
                            s = float(mgv['cx_stats']['All']['Satisfied'])
                        except ValueError:
                            continue
                        Scores.append(f'Model group/Ensemble {i}')
                        counts.append(s)

            if len(counts) > 0:
                legends = [f'{i} %' for i in counts]
                # Select dark orange;
                # identical in all plots because they're separated
                source = ColumnDataSource(data=dict(
                    Scores=Scores, counts=counts, legends=legends, color=linear_palette(Oranges256, len(legends) + 2)[1:-1]))
                pf = figure(y_range=Scores, x_range=(0, 102), frame_height=len(counts) * GLANCE_BAR_HEIGHT,
                            width=GLANCE_PLOT_WIDTH, title="Crosslink satisfaction")
                rf = pf.hbar(y='Scores', right='counts', color='color', height=1.0,
                             source=source, line_color='white')
                pf.ygrid.grid_line_color = None
                pf.xaxis.axis_label = 'Satisfaction rate [%]'
                legend = Legend(items=[LegendItem(label=legends[i], renderers=[
                                rf], index=i) for i in range(len(legends))], location="center",
                                orientation='vertical', label_text_font_size="14pt")
                pf.add_layout(legend, 'right')
                pf.legend.items.reverse()
                pf.legend.border_line_width = 0
                pf.title.vertical_align = 'top'
                pf.title.align = "center"
                pf.output_backend = "svg"
                pf.border_fill_color = None
                pf.background_fill_color = None
                pf.title.text_font_size = '14pt'
                pf.legend.label_text_font_size = "14pt"
                pf.xaxis.axis_label_text_font_size = "14pt"
                pf.yaxis.axis_label_text_font_size = "14pt"
                pf.xaxis.major_label_text_font_size = "14pt"
                pf.yaxis.major_label_text_font_size = "14pt"
                pf.xaxis.axis_label_text_font_style = 'normal'
                pf.yaxis.axis_label_text_font_style = 'normal'
                fq_plots.append(pf)

        if em_fit is not None and len(em_fit) > 0:
            Scores = []
            counts = []
            for dataset in em_fit:
                emdbid = dataset['emdbid']
                for mid, data_ in dataset['fit_stats'].items():
                    try:
                        s = float(data_['q_score']['average'])
                    except ValueError:
                        continue
                    Scores.append(f'Model {mid}/{emdbid}')
                    counts.append(s)

            if len(counts) > 0:
                legends = [f'{i:.3f}' for i in counts]
                source = ColumnDataSource(data=dict(
                    Scores=Scores, counts=counts, legends=legends, color=linear_palette(Greens256, len(legends) + 2)[1:-1]))
                pf = figure(y_range=Scores, x_range=(-1, 1), frame_height=len(counts) * GLANCE_BAR_HEIGHT,
                            width=GLANCE_PLOT_WIDTH, title="Q-score")
                rf = pf.hbar(y='Scores', right='counts', color='color', height=1.0,
                             source=source, line_color='white')
                pf.ygrid.grid_line_color = None
                pf.xaxis.axis_label = 'Q-score'
                legend = Legend(items=[LegendItem(label=legends[i], renderers=[
                                rf], index=i) for i in range(len(legends))], location="center",
                                orientation='vertical', label_text_font_size="14pt")
                pf.add_layout(legend, 'right')
                pf.legend.items.reverse()
                pf.legend.border_line_width = 0
                pf.title.vertical_align = 'top'
                pf.title.align = "center"
                pf.output_backend = "svg"
                pf.border_fill_color = None
                pf.background_fill_color = None
                pf.title.text_font_size = '14pt'
                pf.legend.label_text_font_size = "14pt"
                pf.xaxis.axis_label_text_font_size = "14pt"
                pf.yaxis.axis_label_text_font_size = "14pt"
                pf.xaxis.major_label_text_font_size = "14pt"
                pf.yaxis.major_label_text_font_size = "14pt"
                pf.xaxis.axis_label_text_font_style = 'normal'
                pf.yaxis.axis_label_text_font_style = 'normal'
                fq_plots.append(pf)

        if len(fq_plots) > 0:
            pd = gridplot(fq_plots, ncols=1,
                     toolbar_location="above",
                     # sizing_mode='stretch_width'
                     )

            plots = self.save_plots(pd, 'quality_at_glance_FQ')
            glance_plots['FQ'] = True

        return glance_plots

    def save_plots(self, p, plot_name: str) -> dict:
        """Save html and svg plots"""
        fname_html = Path(self.imageDirName, f"{self.ID_f}_{plot_name}.html")
        fname_svg = fname_html.with_suffix('.svg')
        fname_json = fname_html.with_suffix('.json')
        # output_file(filename=fname_html, mode='inline')
        # save(p)
        export_svg(p, filename=fname_svg)
        utility.strip_bokeh_invisible_outline(fname_svg)

        with open(fname_json, 'w') as f:
            json.dump(json_item(p), f)

        return {'html': fname_html, 'svg': fname_svg, 'json': fname_json}
