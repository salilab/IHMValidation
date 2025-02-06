###################################
# Script :
# 1) Contains class to generate SAS
# plots
# 2) Inherits from SAS class
#
# ganesans - Salilab - UCSF
# ganesans@salilab.org
###################################
import pandas as pd
from pathlib import Path
from mmcif_io import GetInputInformation
import sas
from bokeh.io import output_file, export_svg
from bokeh.models import Span, ColumnDataSource
from bokeh.plotting import figure, save
from bokeh.layouts import column, gridplot
from bokeh.embed import json_item
import utility
import json

class SasValidationPlots(sas.SasValidation):
    def __init__(self, mmcif_file, imageDirName, driver, db='.'):
        super().__init__(mmcif_file, db=db)
        self.df_dict = self.intensities
        self.pdf_dict = sas.SasValidation.get_pddf(self)
        self.fdf_dict = sas.SasValidation.get_fit_data(self)
        self.pdf_ext_dict = sas.SasValidation.get_pofr_ext(self)
        self.pdf_dict_err = sas.SasValidation.get_pofr_errors(self)
        self.score, self.gdf = sas.SasValidation.get_Guinier_data(self)
        self.imageDirName = imageDirName
        self.driver = driver


    @staticmethod
    def set_plot_style(p):
        p.xaxis.major_label_text_font_size = "14pt"
        p.yaxis.major_label_text_font_size = "14pt"
        p.xaxis.axis_label_text_font_size = '14pt'
        p.yaxis.axis_label_text_font_size = '14pt'
        p.xaxis.axis_label_text_font_style = 'italic'
        p.yaxis.axis_label_text_font_style = 'italic'

        if p.title is not None:
            p.title.text_font_size = '12pt'
            p.title.align = "center"
            p.title.vertical_align = 'top'

        p.output_backend = "svg"

        return p

    def plot_intensities(self, sasbdb: str, df: pd.DataFrame):
        '''
        plot intensities with errors
        '''
        source = ColumnDataSource(df)
        p = figure(plot_height=350, plot_width=350,
                   title=f"Log I(q) vs q with for {sasbdb}",
                   toolbar_location="above",
                   toolbar_sticky=False)
        p.circle(x='Q', y='logI', source=source,
                 # color='blue',
                 fill_alpha=0.3,
                 size=5)
        p.multi_line('err_x', 'err_y', source=source,
                     color='gray', line_width=0.5)

        p.xaxis.axis_label = "q [nm\u207B\u00B9]"
        p.yaxis.axis_label = 'Log I(q) [a.u]'

        p = self.set_plot_style(p)

        self.save_plots(p, f"{sasbdb}_intensities")

    def plot_intensities_log(self, sasbdb: str, df: pd.DataFrame):
        '''
        plot intensities on a log scale with errors
        '''
        source = ColumnDataSource(df)
        p = figure(plot_height=350, plot_width=350,
                   title=f"Log I(q) vs Log q for {sasbdb}",
                   toolbar_location="above",
                   toolbar_sticky=False)
        p.circle(x='logQ', y='logI', source=source,
                 # color='blue',
                 fill_alpha=0.3,
                 size=5)
        p.multi_line('logX', 'err_y', source=source,
                     color='gray', line_width=0.5)

        p.xaxis.axis_label = 'Log q [nm\u207B\u00B9]'
        p.yaxis.axis_label = 'Log I(q) [a.u]'

        p = self.set_plot_style(p)

        self.save_plots(p, f"{sasbdb}_intensities_log")

    def plot_kratky(self, sasbdb: str, df: pd.DataFrame):
        '''
        plot dimensionless kratky
        '''
        source = ColumnDataSource(df)
        p = figure(plot_height=350, plot_width=350,
                   title=f"Dimensionless Kratky plot {sasbdb}",
                   toolbar_location="above",
                   toolbar_sticky=False)
        p.circle(x='Kx', y='Ky', source=source,
                 # color='blue',
                 fill_alpha=0.3,
                 size=5)
        # vline = Span(location=0.1732, dimension='height', line_color="crimson", line_width=3)
        # hline = Span(location=0.1104, dimension='width', line_color='green', line_width=3)
        # p.renderers.extend([vline, hline])
        p.xaxis.axis_label = 'qRg'
        p.yaxis.axis_label = 'q\u00B2 Rg\u00B2 I(q)/I(0)'

        p = self.set_plot_style(p)

        self.save_plots(p, f"{sasbdb}_kratky")

    def plot_porod_debye(self, sasbdb: str, df: pd.DataFrame):
        '''
        Porod-Debye plot for flexibility
        '''
        source = ColumnDataSource(df)
        p = figure(plot_height=350, plot_width=350,
                   title=f"Porod-Debye plot {sasbdb}",
                   toolbar_location="above",
                   toolbar_sticky=False)
        p.circle(x='Px', y='Py', source=source,
                 #color='blue',
                 fill_alpha=0.3,
                 size=5)

        p.xaxis.axis_label = 'q \u2074'
        p.yaxis.axis_label = 'q\u2074 I(q)'

        p = self.set_plot_style(p)

        self.save_plots(p, f"{sasbdb}_porod")

    def plot_pddf(self, sasbdb: str, df: pd.DataFrame, Rg=None, Dmax=None):
        '''
        p(r) plot, deprecated function
        '''
        source = ColumnDataSource(df)
        ymax = max(df['P']) * 1.1

        p = figure(plot_height=350, plot_width=350,
                   y_range=(0, ymax),
                   title=f"P(r) {sasbdb}",
                   toolbar_location="above",
                   toolbar_sticky=False)
        p.circle(x='R', y='P', source=source,
                 # color='blue',
                 fill_alpha=0.3,
                 size=5)
        p.multi_line('err_x', 'err_y', source=source,
                     color='gray', line_width=1.5)

        if Rg is not None or Rg != utility.NA:
            p.line([Rg, Rg], [0, ymax], legend_label='Rg',
                          line_color="crimson", line_width=3)

        if Dmax is not None or Dmax != utility.NA:
            p.line([Dmax, Dmax], [0, ymax], legend_label=r"Dmax",
                          line_color="grey", line_width=3)

        p.xaxis.axis_label = "r [nm]"
        p.yaxis.axis_label = 'P(r)'

        p.legend.orientation = "vertical"
        p.legend.location = "top_right"
        p.legend.border_line_width = 0

        p = self.set_plot_style(p)

        self.save_plots(p, f"{sasbdb}_pddf")

    def plot_pddf_int_rwt(self, sasbdb: str,
                      df_int: pd.DataFrame,
                      df_pofr: pd.DataFrame,
                      df_error: pd.DataFrame,):
        '''
        p(r) with fit
        '''
        source1 = ColumnDataSource(df_int)
        source2 = ColumnDataSource(df_pofr)
        p1 = figure(plot_height=350, plot_width=350,
                   title=f"P(r) extrapolated fit for {sasbdb}")
        legend1 = 'Experimental data'
        legend2 = "Extrapolated fit"

        p1.circle(x='Q', y='logI', source=source1,
                  #color='blue',
                  line_width=1, fill_alpha=0.3, size=3, legend_label=legend1)
        p1.line(x='Q', y='logI', source=source2, color="crimson",
               line_width=3, legend_label=legend2)
        # p.circle(x='Q',y='logIb',source=source, color="crimson",line_width=1,fill_alpha=0.1,size=3,legend_label=legend2)

        p1.xaxis.axis_label = "q [nm\u207B\u00B9]"
        p1.yaxis.axis_label = 'Log I(q) [a.u]'

        p1.legend.orientation = "vertical"
        p1.legend.location = "top_right"
        p1.legend.border_line_width = 0

        p1 = self.set_plot_style(p1)
        p1.xaxis.visible = False

        source3 = ColumnDataSource(df_error)
        p2 = figure(plot_height=150, plot_width=350,
                   x_range=p1.x_range,
                   # title="Error weighted residuals for P(r) fit ("+sasbdb+")"
                   )
        p2.circle(x='Q', y='WR', source=source3,
                  # color='blue',
                  fill_alpha=0.3, size=5)
        hline = Span(location=0, dimension='width',
                     line_color="crimson", line_width=3)
        p2.renderers.extend([hline])

        p2.xaxis.axis_label = 'q [nm\u207B\u00B9]'
        p2.yaxis.axis_label = 'R/\u03C3'

        p2 = self.set_plot_style(p2)
        p2.yaxis.ticker = [-3, 0, 3]

        p = gridplot([p1, p2], ncols=1,
                     toolbar_location="above",
                     # sizing_mode='stretch_width'
                     )

        p.width = 350
        p.height = 550

        self.save_plots(p, f"{sasbdb}_pddf_int")

    def Guinier_plot_fit_rwt(self, sasbdb: str, df: pd.DataFrame, score: float):
        '''
        Gunier plot with fit
        '''
        source = ColumnDataSource(df)
        p1 = figure(plot_height=350, plot_width=350,
                   title=f"Guinier plot for {sasbdb}")
        legend1 = 'Experimental data'
        legend2 = f"Linear fit (R\u00B2={score})"
        p1.circle(x='Q2', y='logI', source=source,
                  # color='blue',
                 line_width=1, fill_alpha=0.3, size=5, legend_label=legend1)
        p1.line(x='Q2', y='y_pred', source=source,
               color="crimson", line_width=3, legend_label=legend2)

        p1.xaxis.axis_label = "q [nm\u207B\u00B2]"  # \u212B\u207B\u00B2"
        p1.yaxis.axis_label = 'Log I(q)'

        p1.legend.orientation = "vertical"
        p1.legend.location = "top_right"
        p1.legend.border_line_width = 0

        p1 = self.set_plot_style(p1)
        p1.yaxis.ticker.desired_num_ticks = 3
        p1.xaxis.ticker.desired_num_ticks = 3
        p1.xaxis.visible = False

        p2 = figure(plot_height=150, plot_width=350,
                    x_range=p1.x_range
                   # title="Residuals for Guinier plot fit ("+sasbdb+")"
                    )
        p2.circle(x='Q2', y='res', source=source,
                  # color='blue', fill_alpha=0.3, size=5,
                  )
        hline = Span(location=0, dimension='width',
                     line_color="crimson", line_width=3)
        p2.renderers.extend([hline])

        p2.xaxis.axis_label = "q [nm\u207B\u00B2]"  # \u212B\u207B\u00B2"
        p2.yaxis.axis_label = 'R'

        p2 = self.set_plot_style(p2)

        p2.xaxis.ticker.desired_num_ticks = 3
        p2.yaxis.ticker.desired_num_ticks = 3

        p = gridplot([p1, p2], ncols=1,
                     toolbar_location="above",
                     # sizing_mode='stretch_width'
                     )

        p.width = 350
        p.height = 550

        self.save_plots(p, f"{sasbdb}_guinier")

    def plot_multiple(self):
        for sasbdb, df in self.df_dict.items():
            self.plot_intensities(sasbdb, df)
            self.plot_intensities_log(sasbdb, df)
            self.plot_kratky(sasbdb, df)
            self.plot_porod_debye(sasbdb, df)
            self.plot_pddf_int_rwt(sasbdb, df,
                                   self.pdf_ext_dict[sasbdb],
                                   self.pdf_dict_err[sasbdb])

    def plot_Guinier(self):
        for sasbdb, df in self.gdf.items():
            self.Guinier_plot_fit_rwt(sasbdb, df, self.score[sasbdb])

    def plot_pf(self):
        for sasbdb, df in self.pdf_dict.items():
            sascif = self.sascif_dicts[sasbdb]
            main = f'{sasbdb}_MAIN'
            data = sascif[main]['_sas_result']

            try:
                Dmax = float(data["D_max"])
            except ValueError:
                Dmax = None
            try:
                Rg = float(data["Rg_from_PR"])
            except ValueError:
                Rg = None

            self.plot_pddf(sasbdb, df,
                           Rg=Rg, Dmax=Dmax)

    def plot_fits(self):
        for sasbdb, df in self.fdf_dict.items():
            for sasdb_m, df_m in df.items():
                if not df_m[1].empty:
                    self.plot_fit_rwt(sasbdb, sasdb_m, df_m[0], df_m[1])

    def plot_fit_rwt(self, sasbdb: str, fit: int, score: float, df: pd.DataFrame):
        '''
        plot chi-squared fit
        '''
        source = ColumnDataSource(df)
        p1 = figure(plot_height=350, plot_width=350,
                   title=f"Model fit for {sasbdb}, model {fit + 1}")
        legend1 = 'Experimental data'
        legend2 = "Model fit"
        p1.circle(x='Q', y='logIe', source=source,
                  # color='blue',
                  line_width=1, fill_alpha=0.3, size=3, legend_label=legend1)
        p1.line(x='Q', y='logIb', source=source, color="crimson",
                line_width=3, legend_label=legend2)
        # p.circle(x='Q',y='logIb',source=source, color="crimson",line_width=1,fill_alpha=0.1,size=3,legend_label=legend2)

        p1.yaxis.axis_label = 'Log I(q)'

        p1.legend.orientation = "vertical"
        p1.legend.location = "top_right"
        p1.legend.border_line_width = 0

        p1 = self.set_plot_style(p1)

        p1.xaxis.visible = False
        p1.yaxis.ticker.desired_num_ticks = 3
        p1.xaxis.ticker.desired_num_ticks = 3

        p2 = figure(plot_height=150, plot_width=350,
                   title=None, x_range=p1.x_range)
        p2.circle(x='Q', y='rsigma', source=source,
                  # color='blue',
                  fill_alpha=0.3, size=3)

        p2.xaxis.axis_label = "q [\u212B\u207B\u00B9]"
        p2.yaxis.axis_label = 'R/\u03C3'

        hline = Span(location=0, dimension='width',
                     line_color="crimson", line_width=3)
        p2.renderers.append(hline)
        hline = Span(location=-3, dimension='width',
                     line_color='grey', line_width=2)
        p2.renderers.append(hline)
        hline = Span(location=3, dimension='width',
                     line_color='grey', line_width=2)
        p2.renderers.append(hline)

        p2 = self.set_plot_style(p2)

        p2.yaxis.ticker = [-3, 0, 3]
        p2.xaxis.ticker.desired_num_ticks = 3

        p = gridplot([p1, p2], ncols=1,
                     toolbar_location="above")

        self.save_plots(p, f"{sasbdb}_{fit:d}_fit1")

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
