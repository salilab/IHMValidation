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
from pathlib import Path
import utility
from mmcif_io import GetInputInformation
import bokeh
import numpy as np
from bokeh.io import output_file, curdoc, export_svg, show
from bokeh.models import (ColumnDataSource, Legend, LegendItem, FactorRange,
                          Div, BasicTickFormatter)
from bokeh.palettes import viridis, Reds256, linear_palette
from bokeh.plotting import figure, save
from bokeh.models.widgets import Tabs, Panel
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

class Plots(GetInputInformation):
    def __init__(self, mmcif, imageDirName, driver):
        super().__init__(mmcif)
        # self.dirname = os.path.dirname(os.path.abspath(__file__))
        self.imageDirName = imageDirName
        self.filename = imageDirName
        self.driver=driver

    def plot_quality_at_glance(self, molprobity_data: dict, exv_data: dict,
                               sas_data: dict, sas_fit: dict, cx_data_quality: dict, cx_fit: dict) -> dict:

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
                    plot_height=120,
                    plot_width=700
                )

                p_ = p.hbar(y=source.data['y'][i * 3: (i + 1) * 3],
                       right=source.data['counts'][i * 3: (i + 1) * 3],
                       line_color="white",
                       fill_color=factor_cmap('y', palette=viridis(len(Scores)),
                                              factors=Scores,
                                              start=1, end=2)
                       )

                # set labels and fonts
                p.xaxis.major_label_text_font_size = "12pt"
                p.yaxis.major_label_text_font_size = "12pt"
                p.xaxis.axis_label = 'Outliers'
                p.xaxis.axis_label_text_font_style = 'italic'
                p.left[0].group_text_font_size = '14px'
                p.left[0].group_text_color = p.left[0].major_label_text_color
                p.left[0].group_text_font_style = p.left[0].major_label_text_font_style
                p.left[0].group_text_font_size = p.left[0].major_label_text_font_size
                p.left[0].group_label_orientation = 'horizontal'
                p.title.vertical_align = 'top'
                p.title.align = "center"
                p.output_backend = "svg"

                fname = Path(self.imageDirName, f"{self.ID_f}_{name_}_quality_at_glance_MQ_mp.svg")
                export_svg(p, filename=fname, webdriver=self.driver)

                plots.append(p)

            grid = gridplot(plots, ncols=1,
                            merge_tools=True,
                            toolbar_location='right')
            grid.children[1].css_classes = ['scrollable']
            grid.children[1].sizing_mode = 'fixed'
            grid.children[1].height = 450
            grid.children[1].width = 800

            title = Div(text="<p>Model Quality: Molprobity Analysis</p>",
                        style={"font-size": "1.5em", "font-weight": "bold",
                               "text-align": "center", "width": '100%'}, width=800
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
            legends = [f'{s:.2f}%' for s in satisfaction]

            # set the size of the axis
            # n = 3 if len(model) < 3 else len(model)
            n = len(counts)
            source = ColumnDataSource(
                data=dict(Scores=Scores, counts=counts, legends=legends, color=viridis(n)))

            #  build plots
            plots = []

            # get ranges
            lower, upper = 0, 102

            for i, name_ in enumerate(Models):
                p = figure(y_range=FactorRange(factors=source.data['Scores'][i: i + 1]), x_range=(lower, upper), plot_height=90,
                           plot_width=700)  # , title='Model Quality: Excluded Volume Analysis')
                # p.xaxis.formatter = BasicTickFormatter(use_scientific=True, power_limit_high=3)
                p.xaxis.ticker.desired_num_ticks = 3

                p_ = p.hbar(y=source.data['Scores'][i:i + 1], right=source.data['counts'][i: i + 1], color=source.data['color'][i:i + 1], height=1.0,
                           alpha=0.8, line_color='white')

                p.xaxis.axis_label = 'Satisfaction rate, %'
                legend = Legend(items=[LegendItem(label=legends[i:i + 1][j], renderers=[
                    p_], index=j) for j in range(len(legends[i:i + 1]))], location='center',
                    label_text_font_size='12px', orientation='vertical')
                p.add_layout(legend, 'right')
                p.legend.border_line_width = 0
                p.xaxis.major_label_text_font_size = "12pt"
                p.yaxis.major_label_text_font_size = "12pt"
                p.title.vertical_align = 'top'
                p.title.align = "center"
                p.output_backend = "svg"
                p.min_border_top = 20

                fname = Path(self.imageDirName, f"{self.ID_f}_{name_}_quality_at_glance_MQ_exv.svg")
                export_svg(p, filename=fname, webdriver=self.driver)

                plots.append(p)

            grid = gridplot(plots, ncols=1,
                            merge_tools=True,
                            toolbar_location='right')
            grid.children[1].css_classes = ['scrollable']
            grid.children[1].sizing_mode = 'fixed'
            grid.children[1].height = 450
            grid.children[1].width = 800

            title = Div(text='<p>Model Quality: Excluded Volume Analysis</p>',
                        style={"font-size": "1.5em", "font-weight": "bold",
                               "text-align": "center", "width": '100%'}, width=800
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
        #               plot_height=300, plot_width=800)
        #
        #    p.ygrid.grid_line_color = None
        #    p.xaxis.axis_label_text_font_size = "14pt"
        #    p.yaxis.axis_label_text_font_size = "14pt"
        #    p.title.text_font_size = '14pt'
        #    p.title.align = "center"
        #
        #    p.output_backend = "svg"
        #    p.title.vertical_align = 'top'
        #    fullplot = p
        #
        #    mq_plots.append(fullplot)
        # make panel figures
        # first panel is model quality
#        export_svg(fullplot, filename=self.filename+'/' +
#                    self.ID_f + "_quality_at_glance_MQ.svg", webdriver=self.driver)
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
                counts)+1), plot_height=90 + len(counts) * 20, plot_width=800, title="Data Quality for SAS: Rg Analysis",)
            rd = pd.hbar(y='Scores', right='counts', color='color', height=1.0,
                         source=source, alpha=0.8, line_color='white')
            pd.ygrid.grid_line_color = None
            pd.xaxis.axis_label = 'Distance (nm)'
            pd.xaxis.major_label_text_font_size = "12pt"
            pd.yaxis.major_label_text_font_size = "12pt"
            pd.title.text_font_size = '14pt'
            legend = Legend(items=[LegendItem(label=legends[i], renderers=[
                            rd], index=i) for i in range(len(legends))], location='center',
                            orientation='vertical', label_text_font_size="12px")
            pd.add_layout(legend, 'right')
            pd.legend.items.reverse()
            pd.legend.border_line_width = 0
            pd.legend.label_text_font_size = "12px"
            pd.xaxis.axis_label_text_font_style = 'italic'
            pd.yaxis.axis_label_text_font_style = 'italic'
            pd.xaxis.axis_label_text_font_size = "14pt"
            pd.yaxis.major_label_text_font_size = "14pt"
            pd.title.vertical_align = 'top'
            pd.title.align = "center"
            pd.output_backend = "svg"

            dq_plots.append(pd)

        # If crosslinking-MS data is available
        if cx_data_quality is not None and len(cx_data_quality) > 0:
            Models = [data["pride_id"] for data in cx_data_quality]
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
                            style={"font-size": "1.5em", "font-weight": "bold",
                                   "text-align": "center"}
                            )


                p = figure(
                    y_range=FactorRange(*y[i * 3: (i + 1) * 3]),
                    # Force left limit at zero
                    x_range=(lower, upper),
                    plot_height=95 + 3 * 20,
                    plot_width=700,
                    title=title_txt
                )
                p.xaxis.ticker.desired_num_ticks = 3

                rd = p.hbar(y=source.data['y'][i * 3: (i + 1) * 3],
                       right=source.data['counts'][i * 3: (i + 1) * 3],
                       line_color="white",
                       fill_color=factor_cmap('y', palette=viridis(len(Scores)),
                                              factors=Scores,
                                              start=1, end=2)
                       )
                legends_ = source.data['legends'][i * 3: (i + 1) * 3]
                legend = Legend(items=[LegendItem(label=legends_[j], renderers=[
                            rd], index=j) for j in range(len(legends_))], location='center',
                            orientation='vertical', label_text_font_size="12px")
                legend.items = legend.items[::-1]
                p.add_layout(legend, 'right')
                p.legend.border_line_width = 0
                # set labels and fonts
                p.xaxis.major_label_text_font_size = "12pt"
                p.yaxis.major_label_text_font_size = "12pt"
                p.yaxis.major_label_text_align='right'
                p.xaxis.axis_label = 'Residue pairs'
                p.xaxis.axis_label_text_font_style = 'italic'
                p.left[0].group_text_font_size = '14px'
                p.left[0].group_label_orientation = 'horizontal'
                p.title.text_font_size = '14pt'
                p.title.vertical_align = 'top'
                p.title.align = "center"
                p.output_backend = "svg"
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
                Scores=Scores, counts=counts, legends=legends, color=viridis(len(legends))))
            pf = figure(y_range=Scores, x_range=(0, max(counts)+1), plot_height=100 + len(counts) * 20,
                        plot_width=800, title="Fit to SAS Data:  \u03C7\u00b2 Fit")
            rf = pf.hbar(y='Scores', right='counts', color='color', height=1.0,
                         source=source, alpha=0.8, line_color='white')
            pf.ygrid.grid_line_color = None
            pf.title.text_font_size = '14pt'
            pf.xaxis.axis_label = 'Fit value'
            pf.xaxis.major_label_text_font_size = "12pt"
            pf.yaxis.major_label_text_font_size = "12pt"
            legend = Legend(items=[LegendItem(label=legends[i], renderers=[
                            rf], index=i) for i in range(len(legends))], location="center",
                            orientation='vertical', label_text_font_size="12px")
            pf.add_layout(legend, 'right')
            pf.legend.items.reverse()
            pf.legend.border_line_width = 0
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
                source = ColumnDataSource(data=dict(
                    Scores=Scores, counts=counts, legends=legends, color=viridis(len(legends))))
                pf = figure(y_range=Scores, x_range=(0, 102), plot_height=95 + len(counts) * 20,
                            plot_width=800, title="Crosslink satisfaction")
                rf = pf.hbar(y='Scores', right='counts', color='color', height=1.0,
                             source=source, alpha=0.8, line_color='white')
                pf.ygrid.grid_line_color = None
                pf.title.text_font_size = '14pt'
                pf.xaxis.axis_label = 'Satisfaction rate, %'
                pf.xaxis.major_label_text_font_size = "12pt"
                pf.yaxis.major_label_text_font_size = "12pt"
                legend = Legend(items=[LegendItem(label=legends[i], renderers=[
                                rf], index=i) for i in range(len(legends))], location="center",
                                orientation='vertical', label_text_font_size="12px")
                pf.add_layout(legend, 'right')
                pf.legend.items.reverse()
                pf.legend.border_line_width = 0
                pf.title.vertical_align = 'top'
                pf.title.align = "center"
                pf.output_backend = "svg"
                pf.legend.label_text_font_size = "12px"
                pf.xaxis.axis_label_text_font_style = 'italic'
                pf.yaxis.axis_label_text_font_style = 'italic'
                pf.xaxis.axis_label_text_font_size = "14pt"
                pf.yaxis.major_label_text_font_size = "14pt"
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
        export_svg(p, filename=fname_svg, webdriver=self.driver)

        with open(fname_json, 'w') as f:
            json.dump(json_item(p), f)

        return {'html': fname_html, 'svg': fname_svg, 'json': fname_json}
