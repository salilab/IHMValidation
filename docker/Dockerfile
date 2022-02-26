FROM centos:7
MAINTAINER Ben Webb <ben@salilab.org>
COPY ATSAS-3.0.5-2.el7.x86_64.rpm /
RUN yum -y localinstall /ATSAS-3.0.5-2.el7.x86_64.rpm && yum -y install libgfortran5 git && rm -f /ATSAS*.rpm
RUN curl -L https://github.com/conda-forge/miniforge/releases/download/4.11.0-0/Miniforge3-4.11.0-0-Linux-x86_64.sh -o miniforge.sh && echo "4c24b38969ac413efa3a984290296f72578340d06004f2c7ba5efcbf828ec86b  miniforge.sh" sha256sum -c --status && bash ./miniforge.sh -b -p /root/miniforge && rm -f miniforge.sh && /root/miniforge/bin/conda init -q bash
RUN source /root/.bashrc && eval "$(conda shell.bash hook)" && conda install -y cctbx-base pandas flask ihm scikit-learn selenium bokeh requests wkhtmltopdf firefox geckodriver && pip install python-decouple pdfkit
RUN git clone -b development https://github.com/salilab/IHMValidation.git
RUN source /root/.bashrc && eval "$(conda shell.bash hook)" && yum -y install which && echo "ATSAS=\"`which datcmp`\"" > IHMValidation/.env && echo "Molprobity_ramalyze=\"`which molprobity.ramalyze`\"" >> IHMValidation/.env && echo "Molprobity_molprobity=\"`which molprobity.molprobity`\"" >> IHMValidation/.env && echo "Molprobity_clashscore=\"`which molprobity.clashscore`\"" >> IHMValidation/.env && echo "Molprobity_rotalyze=\"`which molprobity.rotalyze`\"" >> IHMValidation/.env && echo "wkhtmltopdf=\"`which wkhtmltopdf`\"" >> IHMValidation/.env
