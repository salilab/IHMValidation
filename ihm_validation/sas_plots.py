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
import os
from mmcif_io import GetInputInformation
import sas
from bokeh.io import output_file, export_svgs, export_png, export_svg
from bokeh.models import Span, ColumnDataSource
from bokeh.plotting import figure, save
from bokeh.layouts import column, gridplot

class SasValidationPlots(sas.SasValidation):
    def __init__(self, mmcif_file, imageDirName, driver):
        super().__init__(mmcif_file)
        self.ID = str(GetInputInformation.get_id(self))
        # self.intensities = self.get_intensities()
        # self.intensities = self.modify_intensity()

        # self.df_dict = sas.SasValidation.get_intensities(self)
        # self.df_dict = sas.SasValidation.modify_intensity(self)
        self.df_dict = self.intensities
        self.pdf_dict = sas.SasValidation.get_pddf(self)
        self.fdf_dict = sas.SasValidation.get_fit_data(self)
        self.pdf_ext_dict = sas.SasValidation.get_pofr_ext(self)
        self.pdf_dict_err = sas.SasValidation.get_pofr_errors(self)
        self.score, self.gdf = sas.SasValidation.get_Guinier_data(self)
        self.imageDirName = imageDirName
        self.filename = os.path.join(self.imageDirName)
        self.driver = driver

    def plot_intensities(self, sasbdb: str, df: pd.DataFrame):
        '''
        plot intensities with errors
        '''
        output_file(self.ID+sasbdb+"intensities.html", mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=350, plot_width=350,
                   title=f"Log I(q) vs q with error bars\n{sasbdb}",
                   toolbar_location="above",
                   toolbar_sticky=False)
        p.circle(x='Q', y='logI', source=source,
                 color='blue', fill_alpha=0.3, size=5)
        p.multi_line('err_x', 'err_y', source=source,
                     color='gray', line_width=0.5)
        p.xaxis.major_label_text_font_size = "14pt"
        p.yaxis.major_label_text_font_size = "14pt"
        p.title.text_font_size = '12pt'
        p.title.align = "center"
        p.title.vertical_align = 'top'
        p.xaxis.axis_label = "q [nm\u207B\u00B9]"
        p.xaxis.axis_label_text_font_size = '14pt'
        p.yaxis.axis_label = 'Log I(q) [a.u]'
        p.yaxis.axis_label_text_font_size = '14pt'
        p.xaxis.axis_label_text_font_style = 'italic'
        p.yaxis.axis_label_text_font_style = 'italic'

        save(p, filename=self.filename+'/'+self.ID+sasbdb+"intensities.html")
        p.output_backend = "svg"
        export_svg(p, height=350, width=350, filename=self.filename +
                    '/'+self.ID+sasbdb+"intensities.svg", webdriver=self.driver)
        export_png(p, height=350, width=350, filename=self.filename +
                   '/'+self.ID+sasbdb+"intensities.png", webdriver=self.driver)

    def plot_intensities_log(self, sasbdb: str, df: pd.DataFrame):
        '''
        plot intensities on a log scale with errors
        '''
        output_file(self.ID+sasbdb+"intensities_log.html", mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=350, plot_width=350,
                   title=f"Log I(q) vs Log q with error bars\n{sasbdb}",
                   toolbar_location="above",
                   toolbar_sticky=False)
        p.circle(x='logQ', y='logI', source=source,
                 color='blue', fill_alpha=0.3, size=5)
        p.multi_line('logX', 'err_y', source=source,
                     color='gray', line_width=0.5)
        p.xaxis.major_label_text_font_size = "14pt"
        p.yaxis.major_label_text_font_size = "14pt"
        p.title.text_font_size = '12pt'
        p.title.align = "center"
        p.title.vertical_align = 'top'
        p.xaxis.axis_label = 'Log q [nm\u207B\u00B9]'
        p.xaxis.axis_label_text_font_size = '14pt'
        p.yaxis.axis_label = 'Log I(q) [a.u]'
        p.yaxis.axis_label_text_font_size = '14pt'
        p.xaxis.axis_label_text_font_style = 'italic'
        p.yaxis.axis_label_text_font_style = 'italic'

        save(p, filename=self.filename+'/' +
             self.ID+sasbdb+"intensities_log.html")
        p.output_backend = "svg"
        export_svg(p, height=350, width=350, filename=self.filename +
                    '/'+self.ID+sasbdb+"intensities_log.svg", webdriver=self.driver)
        export_png(p, height=350, width=350, filename=self.filename +
                   '/'+self.ID+sasbdb+"intensities_log.png", webdriver=self.driver)

    def plot_kratky(self, sasbdb: str, df: pd.DataFrame):
        '''
        plot dimensionless kratky
        '''
        output_file(self.ID+sasbdb+"Kratky.html", mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=350, plot_width=350,
                   title=f"Dimensionless Kratky plot\n{sasbdb}",
                   toolbar_location="above",
                   toolbar_sticky=False)
        p.circle(x='Kx', y='Ky', source=source,
                 color='blue', fill_alpha=0.3, size=5)
        # vline = Span(location=0.1732, dimension='height', line_color='red', line_width=3)
        # hline = Span(location=0.1104, dimension='width', line_color='green', line_width=3)
        # p.renderers.extend([vline, hline])
        p.xaxis.major_label_text_font_size = "14pt"
        p.yaxis.major_label_text_font_size = "14pt"
        p.title.text_font_size = '12pt'
        p.title.align = "center"
        p.title.vertical_align = 'top'
        p.xaxis.axis_label = 'qRg'
        p.xaxis.axis_label_text_font_size = '14pt'
        p.yaxis.axis_label = 'q\u00B2 Rg\u00B2 I(q)/I(0)'
        p.yaxis.axis_label_text_font_size = '14pt'
        p.xaxis.axis_label_text_font_style = 'italic'
        p.yaxis.axis_label_text_font_style = 'italic'
        save(p, filename=self.filename+'/'+self.ID+sasbdb+"Kratky.html")
        p.output_backend = "svg"
        export_svg(p, filename=self.filename+'/'+self.ID+sasbdb+"Kratky.svg", webdriver=self.driver)
        export_png(p,  height=350, width=350,
                   filename=self.filename+'/'+self.ID+sasbdb+"Kratky.png", webdriver=self.driver)

    def plot_porod_debye(self, sasbdb: str, df: pd.DataFrame):
        '''
        porod debye plot for flexibility
        '''
        output_file(self.ID+sasbdb+"porod.html", mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=350, plot_width=350,
                   title=f"Porod-Debye plot\n{sasbdb}",
                   toolbar_location="above",
                   toolbar_sticky=False)
        p.circle(x='Px', y='Py', source=source,
                 color='blue', fill_alpha=0.3, size=5)
        p.xaxis.major_label_text_font_size = "14pt"
        p.yaxis.major_label_text_font_size = "14pt"
        p.title.text_font_size = '12pt'
        p.title.align = "center"
        p.title.vertical_align = 'top'
        p.xaxis.axis_label = 'q \u2074'
        p.xaxis.axis_label_text_font_size = '14pt'
        p.yaxis.axis_label = 'q\u2074 I(q)'
        p.yaxis.axis_label_text_font_size = '14pt'
        p.xaxis.axis_label_text_font_style = 'italic'
        p.yaxis.axis_label_text_font_style = 'italic'
        p.output_backend = "svg"
        save(p, filename=self.filename+'/'+self.ID+sasbdb+"porod.html")
        export_svg(p, filename=self.filename+'/'+self.ID+sasbdb+"porod.svg", webdriver=self.driver)
        export_png(p,  height=350, width=350,
                   filename=self.filename+'/'+self.ID+sasbdb+"porod.png", webdriver=self.driver)

    def plot_pddf(self, sasbdb: str, df: pd.DataFrame, Rg=None, Dmax=None):
        '''
        p(r) plot, deprecated function
        '''
        output_file(self.ID+sasbdb+"pddf.html", mode="inline")
        source = ColumnDataSource(df)
        ymax = max(df['P']) * 1.1

        p = figure(plot_height=350, plot_width=350,
                   y_range=(0, ymax),
                   title=f"Pair distance distribution function\n{sasbdb}",
                   toolbar_location="above",
                   toolbar_sticky=False)
        p.circle(x='R', y='P', source=source,
                 color='blue', fill_alpha=0.3, size=5)
        p.multi_line('err_x', 'err_y', source=source,
                     color='gray', line_width=1.5)


        if Rg is not None or Rg != 'N/A':
            p.line([Rg, Rg], [0, ymax], legend_label='Rg',
                          line_color="red", line_width=3)

        if Dmax is not None or Dmax != 'N/A':
            p.line([Dmax, Dmax], [0, ymax], legend_label='Dmax',
                          line_color="grey", line_width=3)

        p.legend.orientation = "vertical"
        p.legend.location = "top_right"

        p.xaxis.major_label_text_font_size = "14pt"
        p.yaxis.major_label_text_font_size = "14pt"
        p.title.text_font_size = '12pt'
        p.title.align = "center"
        p.title.vertical_align = 'top'
        p.xaxis.axis_label = "r [nm]"
        p.xaxis.axis_label_text_font_size = '14pt'
        p.yaxis.axis_label = 'P(r)'
        p.yaxis.axis_label_text_font_size = '14pt'
        p.output_backend = "svg"
        p.xaxis.axis_label_text_font_style = 'italic'
        p.yaxis.axis_label_text_font_style = 'italic'
        save(p, filename=self.filename+'/'+self.ID+sasbdb+"pddf.html")
        export_svg(p, height=350, width=350, filename=self.filename+'/'+self.ID+sasbdb+"pddf.svg", webdriver=self.driver)
        export_png(p,  height=350, width=350,
                   filename=self.filename+'/'+self.ID+sasbdb+"pddf.png", webdriver=self.driver)

    def plot_pddf_int_rwt(self, sasbdb: str,
                      df_int: pd.DataFrame,
                      df_pofr: pd.DataFrame,
                      df_error: pd.DataFrame,):
        '''
        p(r) with fit
        '''
        output_file(self.ID+sasbdb+"pddf_int.html", mode="inline")
        source1 = ColumnDataSource(df_int)
        source2 = ColumnDataSource(df_pofr)
        p1 = figure(plot_height=350, plot_width=350,
                   title="P(r) extrapolated fit for "+sasbdb)
        legend1 = 'Experimental data'
        legend2 = "Linear fit"
        p1.circle(x='Q', y='logI', source=source1, color='blue',
                 line_width=1, fill_alpha=0.3, size=3, legend_label=legend1)
        p1.line(x='Q', y='logI', source=source2, color='red',
               line_width=3, legend_label=legend2)
        # p.circle(x='Q',y='logIb',source=source, color='red',line_width=1,fill_alpha=0.1,size=3,legend_label=legend2)
        p1.legend.orientation = "vertical"
        p1.legend.location = "top_right"
        p1.xaxis.major_label_text_font_size = "14pt"
        p1.yaxis.major_label_text_font_size = "14pt"
        p1.title.text_font_size = '12pt'
        p1.title.align = "center"
        p1.title.vertical_align = 'top'
        p1.xaxis.axis_label = "q [nm\u207B\u00B9]"
        p1.xaxis.axis_label_text_font_size = '14pt'
        p1.yaxis.axis_label = 'Log I(q) [a.u]'
        p1.yaxis.axis_label_text_font_size = '14pt'
        p1.xaxis.axis_label_text_font_style = 'italic'
        p1.yaxis.axis_label_text_font_style = 'italic'
        p1.xaxis.visible = False

        source3 = ColumnDataSource(df_error)
        p2 = figure(plot_height=150, plot_width=350,
                   x_range=p1.x_range,
                   # title="Error weighted residuals for P(r) fit ("+sasbdb+")"
                   )
        p2.circle(x='Q', y='WR', source=source3,
                 color='blue', fill_alpha=0.3, size=5)
        hline = Span(location=0, dimension='width',
                     line_color='red', line_width=3)
        p2.renderers.extend([hline])
        p2.xaxis.major_label_text_font_size = "14pt"
        p2.yaxis.major_label_text_font_size = "14pt"
        p2.title.text_font_size = '12pt'
        p2.title.align = "center"
        p2.title.vertical_align = 'top'
        p2.xaxis.axis_label = 'q [nm\u207B\u00B9]'
        p2.xaxis.axis_label_text_font_size = '14pt'
        p2.yaxis.axis_label = 'R/\u03C3'
        p2.yaxis.ticker = [-3, 0, 3]
        p2.yaxis.axis_label_text_font_size = '14pt'
        p2.xaxis.axis_label_text_font_style = 'italic'
        p2.yaxis.axis_label_text_font_style = 'italic'

        p = gridplot([p1, p2], ncols=1,
                     toolbar_location="above",
                     sizing_mode='stretch_width')

        p.width = 350
        p.height = 550

        save(p, filename=self.filename+'/'+self.ID+sasbdb+"pddf_int.html")
        # p.output_backend = "svg"
        export_svg(p, filename=self.filename+'/' +
                    self.ID+sasbdb+"pddf_int.svg", webdriver=self.driver)
        export_png(p, filename=self.filename+'/' +
                   self.ID+sasbdb+"pddf_int.png", webdriver=self.driver)

    def Guinier_plot_fit_rwt(self, sasbdb: str, df: pd.DataFrame, score: float):
        '''
        Gunier plot with fit
        '''
        output_file(self.ID+sasbdb+"guinier.html", mode="inline")
        source = ColumnDataSource(df)
        p1 = figure(plot_height=350, plot_width=350,
                   title=f"Guinier plot for {sasbdb}")
        legend1 = 'Experimental data'
        legend2 = f"Linear fit (R\u00B2={score})"
        p1.circle(x='Q2', y='logI', source=source, color='blue',
                 line_width=1, fill_alpha=0.3, size=5, legend_label=legend1)
        p1.line(x='Q2', y='y_pred', source=source,
               color='red', line_width=3, legend_label=legend2)
        p1.legend.orientation = "vertical"
        p1.legend.location = "top_right"
        p1.xaxis.major_label_text_font_size = "14pt"
        p1.yaxis.major_label_text_font_size = "14pt"
        p1.title.text_font_size = '12pt'
        p1.title.align = "center"
        p1.title.vertical_align = 'top'
        p1.xaxis.axis_label = "q [nm\u207B\u00B2]"  # \u212B\u207B\u00B2"
        p1.xaxis.axis_label_text_font_size = '14pt'
        p1.yaxis.axis_label = 'Log I(q)'
        p1.yaxis.axis_label_text_font_size = '14pt'
        p1.xaxis.axis_label_text_font_style = 'italic'
        p1.yaxis.axis_label_text_font_style = 'italic'
        p1.xaxis.visible = False

        p2 = figure(plot_height=150, plot_width=350,
                    x_range=p1.x_range
                   # title="Residuals for Guinier plot fit ("+sasbdb+")"
                    )
        p2.circle(x='Q2', y='res', source=source,
                 color='blue', fill_alpha=0.3, size=5,
                  )
        hline = Span(location=0, dimension='width',
                     line_color='red', line_width=3)
        p2.renderers.extend([hline])
        p2.xaxis.major_label_text_font_size = "14pt"
        p2.yaxis.major_label_text_font_size = "14pt"
        p2.title.text_font_size = '12pt'
        p2.title.align = "center"
        p2.title.vertical_align = 'top'
        p2.xaxis.axis_label = "q [nm\u207B\u00B2]"  # \u212B\u207B\u00B2"
        p2.xaxis.axis_label_text_font_size = '14pt'
        p2.yaxis.axis_label = 'R'
        p2.yaxis.axis_label_text_font_size = '14pt'
        p2.xaxis.axis_label_text_font_style = 'italic'
        p2.yaxis.axis_label_text_font_style = 'italic'
        p2.yaxis.ticker.desired_num_ticks = 3

        p = gridplot([p1, p2], ncols=1,
                     toolbar_location="above",
                     sizing_mode='stretch_width')

        p.width = 350
        p.height = 550


        save(p, filename=self.filename+'/'+self.ID+sasbdb+"guinier.html")
        # p.output_backend = "svg"
        export_svg(p, filename=self.filename+'/'+self.ID+sasbdb+"guinier.svg", webdriver=self.driver)
        export_png(p, filename=self.filename +
                   '/'+self.ID+sasbdb+"guinier.png", webdriver=self.driver)

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
        output_file(self.ID+sasbdb+str(fit)+"fit1.html", mode="inline")
        source = ColumnDataSource(df)
        p1 = figure(plot_height=350, plot_width=350,
                   title=f"Model fit for {sasbdb}, model {fit + 1}")
        legend1 = 'Experimental data'
        legend2 = "Linear fit"
        p1.circle(x='Q', y='logIe', source=source, color='blue',
                 line_width=1, fill_alpha=0.3, size=3, legend_label=legend1)
        p1.line(x='Q', y='logIb', source=source, color='red',
               line_width=3, legend_label=legend2)
        # p.circle(x='Q',y='logIb',source=source, color='red',line_width=1,fill_alpha=0.1,size=3,legend_label=legend2)
        p1.legend.orientation = "vertical"
        p1.legend.location = "top_right"
        p1.xaxis.visible = False

        p1.yaxis.major_label_text_font_size = "14pt"
        p1.title.text_font_size = '12pt'
        p1.title.align = "center"
        p1.title.vertical_align = 'top'
        p1.yaxis.axis_label = 'Log I(q)'
        p1.yaxis.axis_label_text_font_size = '14pt'
        p1.yaxis.axis_label_text_font_style = 'italic'

        p2 = figure(plot_height=150, plot_width=350,
                   title=None, x_range=p1.x_range)
        p2.circle(x='Q', y='rsigma', source=source,
                 color='blue', fill_alpha=0.3, size=3)

        p2.xaxis.axis_label = "q [\u212B\u207B\u00B9]"
        p2.xaxis.axis_label_text_font_size = '14pt'
        p2.xaxis.axis_label_text_font_style = 'italic'
        p2.xaxis.major_label_text_font_size = "14pt"

        p2.yaxis.axis_label = 'R/\u03C3'
        p2.yaxis.ticker = [-3, 0, 3]
        p2.yaxis.axis_label_text_font_size = '14pt'
        p2.yaxis.axis_label_text_font_style = 'italic'

        hline = Span(location=0, dimension='width',
                     line_color='red', line_width=3)
        p2.renderers.append(hline)
        hline = Span(location=-3, dimension='width',
                     line_color='grey', line_width=2)
        p2.renderers.append(hline)
        hline = Span(location=3, dimension='width',
                     line_color='grey', line_width=2)
        p2.renderers.append(hline)

        p = gridplot([p1, p2], ncols=1,
                     toolbar_location="above")

        save(p, filename=self.filename+'/'+self.ID+sasbdb+str(fit)+"fit1.html")
        # p.output_backend = "svg"
        export_svg(p, height=500, width=500, filename=self.filename+'/' +
                    self.ID+sasbdb+str(fit)+"fit1.svg", webdriver=self.driver)
        export_png(p,  height=500, width=500, filename=self.filename+'/' +
                   self.ID+sasbdb+str(fit)+"fit1.png", webdriver=self.driver)
