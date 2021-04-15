FROM rocker/shiny:4.0.3

LABEL maintainer=TODO
LABEL description=TODO
LABEL version=TODO

RUN Rscript -e 'install.packages(“https://cran.r-project.org/src/contrib/Archive/upsetjs/upsetjs_1.6.0.tar.gz", repos=NULL, type=”source”); \
                install.packages("shiny"); \
                install.packages("DT"); \
                install.packages("VennDiagram"); \
                install.packages("UpSetR"); \
                install.packages("ggplot2"); \
                install.packages("plotly"); \
                install.packages("pheatmap"); \
                install.packages("gtools"); \
                install.packages("openxlsx"); \
                install.packages("vcfR"); \
                install.packages("tidytext"); \
                install.packages("dplyr"); \
                install.packages("heatmaply"); \
                '

# Delete default sample app
RUN rm -rf /srv/shiny-server/*
# Add the virus app
COPY env/Renviron /home/shiny/.Renviron
COPY app.R /srv/shiny-server/