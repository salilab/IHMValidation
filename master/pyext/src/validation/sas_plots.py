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
from validation import sas, GetInputInformation
from bokeh.io import output_file, export_svgs, export_png
from bokeh.models import Span, ColumnDataSource
from bokeh.plotting import figure, save


class SasValidationPlots(sas.SasValidation):
    def __init__(self, mmcif_file):
        super().__init__(mmcif_file)
        self.ID = str(GetInputInformation.get_id(self))
        self.df_dict = sas.SasValidation.modify_intensity(self)
        self.pdf_dict = sas.SasValidation.get_pddf(self)
        self.fdf_dict = sas.SasValidation.get_fit_data(self)
        self.pdf_ext_dict = sas.SasValidation.get_pofr_ext(self)
        self.pdf_dict_err = sas.SasValidation.get_pofr_errors(self)
        self.score, self.gdf = sas.SasValidation.get_Guinier_data(self)
        self.filename = os.path.join('../static/images//')

    def plot_intensities(self, sasbdb: str, df: pd.DataFrame):
        '''
        plot intensities with errors
        '''
        output_file(self.ID+sasbdb+"intensities.html", mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=500, plot_width=500,
                   title="Log I(q) vs q with error bars ("+sasbdb+")")
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
        save(p, filename=self.filename+'/'+self.ID+sasbdb+"intensities.html")
        p.output_backend = "svg"
        export_svgs(p, height=500, width=500, filename=self.filename +
                    '/'+self.ID+sasbdb+"intensities.svg")
        export_png(p, height=500, width=500, filename=self.filename +
                   '/'+self.ID+sasbdb+"intensities.png")

    def plot_intensities_log(self, sasbdb: str, df: pd.DataFrame):
        '''
        plot intensities on a log scale with errors
        '''
        output_file(self.ID+sasbdb+"intensities_log.html", mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=500, plot_width=500,
                   title="Log I(q) vs Log q with error bars ("+sasbdb+")")
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
        save(p, filename=self.filename+'/' +
             self.ID+sasbdb+"intensities_log.html")
        p.output_backend = "svg"
        export_svgs(p, height=500, width=500, filename=self.filename +
                    '/'+self.ID+sasbdb+"intensities_log.svg")
        export_png(p, height=500, width=500, filename=self.filename +
                   '/'+self.ID+sasbdb+"intensities_log.png")

    def plot_kratky_dep(self, sasbdb: str, df: pd.DataFrame):
        '''
        plot kratky plot, deprecated function
        '''
        output_file(self.ID+sasbdb+"Kratky_dep.html", mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=500, plot_width=500,
                   title="Kratky plot ("+sasbdb+")")
        p.circle(x='Q', y='Ky', source=source,
                 color='blue', fill_alpha=0.3, size=5)
        p.xaxis.major_label_text_font_size = "14pt"
        p.yaxis.major_label_text_font_size = "14pt"
        p.title.text_font_size = '12pt'
        p.title.align = "center"
        p.title.vertical_align = 'top'
        p.xaxis.axis_label = 'Log q [nm\u207B\u00B9]'
        p.xaxis.axis_label_text_font_size = '14pt'
        p.yaxis.axis_label = 'q\u00B2 I(q)'
        p.yaxis.axis_label_text_font_size = '14pt'
        save(p, filename=self.filename+'/'+self.ID+sasbdb+"Kratky_dep.html")
        p.output_backend = "svg"
        export_svgs(p, filename=self.filename+'/' +
                    self.ID+sasbdb+"Kratky_dep.svg")
        export_png(p,  height=500, width=500, filename=self.filename+'/' +
                   self.ID+sasbdb+"Kratky_dep.png")

    def plot_kratky(self, sasbdb: str, df: pd.DataFrame):
        '''
        plot dimensionless kratky
        '''
        output_file(self.ID+sasbdb+"Kratky.html", mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=500, plot_width=500,
                   title="Dimensionless Kratky plot ("+sasbdb+")")
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
        save(p, filename=self.filename+'/'+self.ID+sasbdb+"Kratky.html")
        p.output_backend = "svg"
        export_svgs(p, filename=self.filename+'/'+self.ID+sasbdb+"Kratky.svg")
        export_png(p,  height=500, width=500,
                   filename=self.filename+'/'+self.ID+sasbdb+"Kratky.png")

    def plot_porod_debye(self, sasbdb: str, df: pd.DataFrame):
        '''
        porod debye plot for flexibility
        '''
        output_file(self.ID+sasbdb+"porod.html", mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=500, plot_width=500,
                   title="Porod-Debye plot ("+sasbdb+")")
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
        p.output_backend = "svg"
        save(p, filename=self.filename+'/'+self.ID+sasbdb+"porod.html")
        export_svgs(p, filename=self.filename+'/'+self.ID+sasbdb+"porod.svg")
        export_png(p,  height=500, width=500,
                   filename=self.filename+'/'+self.ID+sasbdb+"porod.png")

    def plot_pddf(self, sasbdb: str, df: pd.DataFrame):
        '''
        p(r) plot, deprecated function
        '''
        output_file(self.ID+sasbdb+"pddf.html", mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=500, plot_width=500,
                   title="Pair distance distribution function ("+sasbdb+")")
        p.circle(x='R', y='P', source=source,
                 color='blue', fill_alpha=0.3, size=5)
        p.multi_line('err_x', 'err_y', source=source,
                     color='gray', line_width=1.5)
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
        save(p, filename=self.filename+'/'+self.ID+sasbdb+"pddf.html")
        export_svgs(p, filename=self.filename+'/'+self.ID+sasbdb+"pddf.svg")
        export_png(p,  height=500, width=500,
                   filename=self.filename+'/'+self.ID+sasbdb+"pddf.png")

    def plot_pddf_residuals(self, sasbdb: str, df: pd.DataFrame):
        '''
        p(r) residuals
        '''
        output_file(self.ID+sasbdb+"pddf_residuals.html", mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=500, plot_width=500,
                   title="Residuals for P(r) fit ("+sasbdb+")")
        p.circle(x='Q', y='R', source=source,
                 color='blue', fill_alpha=0.3, size=5)
        hline = Span(location=0, dimension='width',
                     line_color='red', line_width=3)
        p.renderers.extend([hline])
        p.xaxis.major_label_text_font_size = "14pt"
        p.yaxis.major_label_text_font_size = "14pt"
        p.title.text_font_size = '12pt'
        p.title.align = "center"
        p.title.vertical_align = 'top'
        p.xaxis.axis_label = "q\u00B2 [nm \u00B2]"  # \u212B\u207B\u00B2"
        p.xaxis.axis_label_text_font_size = '14pt'
        p.yaxis.axis_label = 'R'
        p.yaxis.axis_label_text_font_size = '14pt'
        save(p, filename=self.filename+'/'+self.ID+sasbdb+"pddf_residuals.html")
        p.output_backend = "svg"
        export_svgs(p, filename=self.filename+'/' +
                    self.ID+sasbdb+"pddf_residuals.svg")
        export_png(p,  height=500, width=500, filename=self.filename+'/' +
                   self.ID+sasbdb+"pddf_residuals.png")

    def plot_pddf_residuals_wt(self, sasbdb: str, df: pd.DataFrame):
        '''
        p(r) error weighted residuals
        '''
        output_file(self.ID+sasbdb+"pddf_residuals_wt.html", mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=500, plot_width=500,
                   title="Error weighted residuals for P(r) fit ("+sasbdb+")")
        p.circle(x='Q', y='WR', source=source,
                 color='blue', fill_alpha=0.3, size=5)
        hline = Span(location=0, dimension='width',
                     line_color='red', line_width=3)
        p.renderers.extend([hline])
        p.xaxis.major_label_text_font_size = "14pt"
        p.yaxis.major_label_text_font_size = "14pt"
        p.title.text_font_size = '12pt'
        p.title.align = "center"
        p.title.vertical_align = 'top'
        p.xaxis.axis_label = "q\u00B2 [nm \u00B2]"  # \u212B\u207B\u00B2"
        p.xaxis.axis_label_text_font_size = '14pt'
        p.yaxis.axis_label = 'R/\u03C3'
        p.yaxis.axis_label_text_font_size = '14pt'
        save(p, filename=self.filename+'/' +
             self.ID+sasbdb+"pddf_residuals_wt.html")
        p.output_backend = "svg"
        export_svgs(p, filename=self.filename+'/' +
                    self.ID+sasbdb+"pddf_residuals_wt.svg")
        export_png(p, height=500, width=500, filename=self.filename+'/' +
                   self.ID+sasbdb+"pddf_residuals_wt.png")

    def plot_pddf_int(self, sasbdb: str, df_int: pd.DataFrame, df_pofr: pd.DataFrame):
        '''
        p(r) with fit
        '''
        output_file(self.ID+sasbdb+"pddf_int.html", mode="inline")
        source1 = ColumnDataSource(df_int)
        source2 = ColumnDataSource(df_pofr)
        p = figure(plot_height=500, plot_width=500,
                   title="P(r) extrapolated fit for "+sasbdb)
        legend1 = 'Experimental data'
        legend2 = "Linear fit"
        p.circle(x='Q', y='logI', source=source1, color='blue',
                 line_width=1, fill_alpha=0.3, size=3, legend_label=legend1)
        p.line(x='Q', y='logI', source=source2, color='red',
               line_width=3, legend_label=legend2)
        # p.circle(x='Q',y='logIb',source=source, color='red',line_width=1,fill_alpha=0.1,size=3,legend_label=legend2)
        p.legend.orientation = "vertical"
        p.legend.location = "top_right"
        p.xaxis.major_label_text_font_size = "14pt"
        p.yaxis.major_label_text_font_size = "14pt"
        p.title.text_font_size = '12pt'
        p.title.align = "center"
        p.title.vertical_align = 'top'
        p.xaxis.axis_label = "q [\u212B\u207B\u207B\u00B9]"
        p.xaxis.axis_label_text_font_size = '14pt'
        p.yaxis.axis_label = 'Log I(q) [a.u]'
        p.yaxis.axis_label_text_font_size = '14pt'
        save(p, filename=self.filename+'/'+self.ID+sasbdb+"pddf_int.html")
        p.output_backend = "svg"
        export_svgs(p, filename=self.filename+'/' +
                    self.ID+sasbdb+"pddf_int.svg")
        export_png(p,  height=500, width=500, filename=self.filename+'/' +
                   self.ID+sasbdb+"pddf_int.png")

    def Guinier_plot_fit(self, sasbdb: str, df: pd.DataFrame, score: float):
        '''
        Gunier plot with fit
        '''
        output_file(self.ID+sasbdb+"guinier.html", mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=500, plot_width=500,
                   title="Guinier plot for "+sasbdb+" (R\u00B2="+str(score)+")")
        legend1 = 'Experimental data'
        legend2 = "Linear fit"
        p.circle(x='Q2A', y='logI', source=source, color='blue',
                 line_width=1, fill_alpha=0.3, size=5, legend_label=legend1)
        p.line(x='Q2A', y='y_pred', source=source,
               color='red', line_width=3, legend_label=legend2)
        p.legend.orientation = "vertical"
        p.legend.location = "top_right"
        p.xaxis.major_label_text_font_size = "14pt"
        p.yaxis.major_label_text_font_size = "14pt"
        p.title.text_font_size = '12pt'
        p.title.align = "center"
        p.title.vertical_align = 'top'
        p.xaxis.axis_label = "q\u00B2 [nm \u00B2]"  # \u212B\u207B\u00B2"
        p.xaxis.axis_label_text_font_size = '14pt'
        p.yaxis.axis_label = 'Log I(q)'
        p.yaxis.axis_label_text_font_size = '14pt'
        save(p, filename=self.filename+'/'+self.ID+sasbdb+"guinier.html")
        p.output_backend = "svg"
        export_svgs(p, filename=self.filename+'/'+self.ID+sasbdb+"guinier.svg")
        export_png(p, height=500, width=500, filename=self.filename +
                   '/'+self.ID+sasbdb+"guinier.png")

    def Guinier_plot_residuals(self, sasbdb: str, df: pd.DataFrame):
        '''
        Guinier plot residuals
        '''
        output_file(self.ID+sasbdb+"guinier_residuals.html", mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=500, plot_width=500,
                   title="Residuals for Guinier plot fit ("+sasbdb+")")
        p.circle(x='Q2A', y='res', source=source,
                 color='blue', fill_alpha=0.3, size=5)
        hline = Span(location=0, dimension='width',
                     line_color='red', line_width=3)
        p.renderers.extend([hline])
        p.xaxis.major_label_text_font_size = "14pt"
        p.yaxis.major_label_text_font_size = "14pt"
        p.title.text_font_size = '12pt'
        p.title.align = "center"
        p.title.vertical_align = 'top'
        p.xaxis.axis_label = "q\u00B2 [nm \u00B2]"  # \u212B\u207B\u00B2"
        p.xaxis.axis_label_text_font_size = '14pt'
        p.yaxis.axis_label = 'R'
        p.yaxis.axis_label_text_font_size = '14pt'
        save(p, filename=self.filename+'/' +
             self.ID+sasbdb+"guinier_residuals.html")
        p.output_backend = "svg"
        export_svgs(p, filename=self.filename+'/' +
                    self.ID+sasbdb+"guinier_residuals.svg")
        export_png(p, height=500, width=500, filename=self.filename+'/' +
                   self.ID+sasbdb+"guinier_residuals.png")

    def plot_fit(self, sasbdb: str, fit: int, score: float, df: pd.DataFrame):
        '''
        plot chi-squared fit
        '''
        output_file(self.ID+sasbdb+str(fit)+"fit1.html", mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=500, plot_width=500,
                   title="Model fit for "+sasbdb)
        legend1 = 'Experimental data'
        legend2 = "Linear fit"
        p.circle(x='Q', y='logIe', source=source, color='blue',
                 line_width=1, fill_alpha=0.3, size=3, legend_label=legend1)
        p.line(x='Q', y='logIb', source=source, color='red',
               line_width=3, legend_label=legend2)
        # p.circle(x='Q',y='logIb',source=source, color='red',line_width=1,fill_alpha=0.1,size=3,legend_label=legend2)
        p.legend.orientation = "vertical"
        p.legend.location = "top_right"
        p.xaxis.major_label_text_font_size = "14pt"
        p.yaxis.major_label_text_font_size = "14pt"
        p.title.text_font_size = '12pt'
        p.title.align = "center"
        p.title.vertical_align = 'top'
        p.xaxis.axis_label = "q [\u212B\u207B\u207B\u00B9]"
        p.xaxis.axis_label_text_font_size = '14pt'
        p.yaxis.axis_label = 'Log I(q)'
        p.yaxis.axis_label_text_font_size = '14pt'
        save(p, filename=self.filename+'/'+self.ID+sasbdb+str(fit)+"fit1.html")
        p.output_backend = "svg"
        export_svgs(p, filename=self.filename+'/' +
                    self.ID+sasbdb+str(fit)+"fit1.svg")
        export_png(p,  height=500, width=500, filename=self.filename+'/' +
                   self.ID+sasbdb+str(fit)+"fit1.png")

    def plot_fit_residuals(self, sasbdb: str, fit: int, df: pd.DataFrame):
        '''
        plot residuals for each fit
        '''
        output_file(self.ID+sasbdb+str(fit)+"residuals.html", mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=500, plot_width=500,
                   title="Residuals for model fit ("+sasbdb+")")
        p.circle(x='Q', y='r', source=source,
                 color='blue', fill_alpha=0.3, size=5)
        hline = Span(location=0, dimension='width',
                     line_color='red', line_width=3)
        p.renderers.extend([hline])
        p.xaxis.major_label_text_font_size = "14pt"
        p.yaxis.major_label_text_font_size = "14pt"
        p.title.text_font_size = '12pt'
        p.title.align = "center"
        p.title.vertical_align = 'top'
        p.xaxis.axis_label = "q [\u212B\u207B\u207B\u00B9]"
        p.xaxis.axis_label_text_font_size = '14pt'
        p.yaxis.axis_label = 'R'
        p.yaxis.axis_label_text_font_size = '14pt'
        save(p, filename=self.filename+'/'+self.ID +
             sasbdb+str(fit)+"residuals.html")
        p.output_backend = "svg"
        export_svgs(p, filename=self.filename+'/' +
                    self.ID+sasbdb+str(fit)+"residuals.svg")
        export_png(p, height=500, width=500, filename=self.filename+'/' +
                   self.ID+sasbdb+str(fit)+"residuals.png")

    def plot_fit_residuals_wt(self, sasbdb: str, fit: int, df: pd.DataFrame):
        '''
        plot error weighted residuals for each fit
        '''
        output_file(self.ID+sasbdb+str(fit)+"residuals_wt.html", mode="inline")
        source = ColumnDataSource(df)
        p = figure(plot_height=500, plot_width=500,
                   title="Error-weighted residuals for model fit ("+sasbdb+")")
        p.circle(x='Q', y='rsigma', source=source,
                 color='blue', fill_alpha=0.3, size=5)
        hline = Span(location=0, dimension='width',
                     line_color='red', line_width=3)
        p.renderers.extend([hline])
        p.xaxis.major_label_text_font_size = "14pt"
        p.yaxis.major_label_text_font_size = "14pt"
        p.title.text_font_size = '12pt'
        p.title.align = "center"
        p.title.vertical_align = 'top'
        p.xaxis.axis_label = "q [\u212B\u207B\u207B\u00B9]"
        p.xaxis.axis_label_text_font_size = '14pt'
        p.yaxis.axis_label = 'R/\u03C3'
        p.yaxis.axis_label_text_font_size = '14pt'
        save(p, filename=self.filename+'/'+self.ID +
             sasbdb+str(fit)+"residuals_wt.html")
        p.output_backend = "svg"
        export_svgs(p, filename=self.filename+'/'+self.ID +
                    sasbdb+str(fit)+"residuals_wt.svg")
        export_png(p,  height=500, width=500, filename=self.filename+'/'+self.ID +
                   sasbdb+str(fit)+"residuals_wt.png")

    def plot_multiple(self):
        for sasbdb, df in self.df_dict.items():
            self.plot_intensities(sasbdb, df)
            self.plot_intensities_log(sasbdb, df)
            self.plot_kratky(sasbdb, df)
            self.plot_porod_debye(sasbdb, df)
            self.plot_pddf_int(sasbdb, df, self.pdf_ext_dict[sasbdb])

    def plot_Guinier(self):
        for sasbdb, df in self.gdf.items():
            self.Guinier_plot_fit(sasbdb, df, self.score[sasbdb])
            self.Guinier_plot_residuals(sasbdb, df)

    def plot_pf(self):
        for sasbdb, df in self.pdf_dict.items():
            self.plot_pddf(sasbdb, df)
            self.plot_pddf_residuals(sasbdb, self.pdf_dict_err[sasbdb])
            self.plot_pddf_residuals_wt(sasbdb, self.pdf_dict_err[sasbdb])

    def plot_fits(self):
        for sasbdb, df in self.fdf_dict.items():
            for sasdb_m, df_m in df.items():
                if not df_m[1].empty:
                    self.plot_fit(sasbdb, sasdb_m, df_m[0], df_m[1])
                    self.plot_fit_residuals(sasbdb, sasdb_m, df_m[1])
                    self.plot_fit_residuals_wt(sasbdb, sasdb_m, df_m[1])
