FROM rocker/shiny:4.1.2

LABEL maintainer=TODO
LABEL description=TODO
LABEL version=TODO
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get -y install xvfb libgtk2.0-0 libgconf-2-4 curl gpg apt-utils dialog libnss3-dev libgtk-3-dev
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
RUN conda install -c plotly plotly-orca

RUN apt-get install -y libxss1 libasound-dev

#devtools::install_github("rstudio/fontawesome"); \
#devtools::install_version("htmltools", version="0.5.2"); \
#install.packages("https://cran.r-project.org/src/contrib/shiny_1.7.1.tar.gz", repos=NULL, type="source"); \
RUN Rscript -e 'install.packages("shiny"); \
                install.packages("shinycssloaders"); \
                install.packages("DT"); \
                install.packages("VennDiagram"); \
                install.packages("UpSetR"); \
                install.packages("ggplot2"); \
                install.packages("svglite"); \
                install.packages("plotly"); \
                install.packages("pheatmap"); \
                install.packages("gtools"); \
                install.packages("openxlsx"); \
                install.packages("vcfR"); \
                install.packages("tidytext"); \
                install.packages("dplyr"); \
                install.packages("heatmaply"); \
                install.packages("https://cran.r-project.org/src/contrib/Archive/upsetjs/upsetjs_1.6.0.tar.gz", repos=NULL, type="source"); \
                install.packages("BiocManager"); \
                BiocManager::install(c("ComplexHeatmap", "InteractiveComplexHeatmap")); \
                '

# Delete default sample app
RUN rm -rf /srv/shiny-server/*
# Add the virus app
COPY env/.Renviron /home/shiny/
COPY env/.Renviron /srv/shiny-server/
COPY app.R /srv/shiny-server/
COPY www /srv/shiny-server/www/

RUN chmod 777 /home/shiny/.Renviron
RUN chmod 777 /srv/shiny-server/
RUN chmod 777 /srv/shiny-server/.Renviron
RUN chmod 777 -R /root/miniconda3/bin

RUN ls -alh /srv/shiny-server
RUN ls -alh /root/miniconda3
RUN ls -alh /root/miniconda3/bin
RUN which orca
