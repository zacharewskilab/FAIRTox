# Install R version 3.6.1
FROM r-base:3.6.1

# Install Ubuntu packages
RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    #libcurl4-gnutls-dev \
    #libcairo2-dev/unstable \
    libxt-dev \
    libssl-dev \
        sqlite3 \
        libsqlite3-dev \
    libv8-dev \
	xml-core \
	libcurl4-openssl-dev -o APT::Immediate-Configure=0 \
	libxml2-dev

# Download and install ShinyServer (latest version)
RUN wget --no-verbose https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb
	#wget --save-cookies cookies.txt 'https://docs.google.com/uc?export=download&id='1FCqYkAFmQ46QqqsgnniJ3d582ACN0Vm3 -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1/p' > confirm.txt && \

# Install R packages that are required
RUN R -e "install.packages('rvest')"
RUN R -e "install.packages('BiocManager'); BiocManager::install(c('graph', 'fgsea', 'mixOmics')); install.packages(c('shiny', 'plotly', 'ggm', 'shinycssloaders', 'ggplot2', 'dplyr', 'reshape2', 'RSQLite', 'markdown', 'plyr', 'tibble', 'UpSetR', 'stringr', 'data.table', 'ggnewscale', 'V8', 'tidyverse', 'rlist', 'readxl', 'devtools', 'factoextra', 'randomcoloR', 'Seurat'), repos='http://cran.rstudio.com/', quiet = TRUE)"
#RUN R -e "devtools::install_version('shinyjs', version = '1.1')"
RUN R -e "install.packages('reticulate')"
RUN R -e "reticulate::install_miniconda(path = reticulate::miniconda_path(), update = TRUE, force = FALSE)"
RUN R -e "reticulate::py_install(packages = 'umap-learn')"
RUN R -e "install.packages('ggridges')"
RUN R -e "install.packages('shinyjs')"
RUN R -e "install.packages('pheatmap')"

# Copy configuration files into the Docker image
COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY /app /srv/shiny-server/
COPY README.md /srv/

# Make the ShinyApp available at port 80
EXPOSE 80

# Copy further configuration files into the Docker image
COPY shiny-server.sh /usr/bin/shiny-server.sh

CMD ["/usr/bin/shiny-server.sh"]

#Acknowledgements: https://www.bjoern-hartmann.de/post/learn-how-to-dockerize-a-shinyapp-in-7-steps/
