# Zacharewski Lab Expression DataViewer (ZED) ShinyApp
Version 4.0

Initiatives that promote Findable, Accessible, Interoperable, and Reusable (FAIR) data principles 
are becoming a key criteria of public research funding agencies. FAIR principles enable new knowledge 
to be derived from existing data with minimal additional investment in data generation. The Michigan 
State University (MSU) Superfund Research Center (SRC) has produced vast amounts of data characterizing 
the mechanisms and impact of exposure to aryl hydrocarbon receptor (AHR) ligands. To promote their 
FAIRness, we have developed FAIRTox, an open-source web-based data exploration, visualization, and 
analysis application for toxicogenomic datasets. FAIRTox is built using a R Shiny framework, chosen 
due to its wide use as a tool for omics data analysis, and the availability of other R packages for 
enrichment and multidimensional analyses such as principal component analysis (PCA). Unprocessed 
and analyzed datasets are stored in a SQL relational database enabling querying through a user-friendly 
interface. Metadata filters allow users to visualize and compare gene expression responses to various 
experimental factors such as zeitgeber time, dose, and duration of exposure to various environmental 
contaminants. Enrichment and advanced analysis features also facilitate the integration of datasets 
furthering the development of novel hypotheses. In addition, the implementation of FAIRTox makes publicly 
available MSU SRC toxicogenomic data more accessible to researchers without transcriptomic expertise. 
FAIRTox is implemented as a Docker container ensuring portability and reproducibility for groups looking 
to execute their own local version. Ultimately, the goal of FAIRTox is to improve data sharing, reuse, 
and reproducibility through an intuitive interface that serves bioinformaticians and novices alike. 

## Getting Started
1) Pull ZED to a fresh repository.

2) Download data files from GDrive link below and place in /app folder within repo
https://drive.google.com/drive/folders/11RZwfEpDIJs8gYE1neV0loOpz041_tRX?usp=sharing

### Prerequisites
The following software is required and is available for download following 
the corresponding links.

R - [Version 3.6.0](https://www.r-project.org/)
```
Packages:
- shiny #App framework
- RSQLite #SQL database import
- ggplot2 #Figure generation
- plotly #Interactive figures
- shinyjs #Tableau implementation
- shinycssloaders #Loading animations
- Seurat #For single cell data
- ggm #Enables powersets
- reshape2 #Array/vector operations
- plyr #Array/vector operations
- tibble #Array operations
- dplyr #Advanced array operations
- UpSetR #UpSet plots
- ggnewscale #tzheatmap requirement
- stringr #Advanced string ops
- fgsea #GSEA analysis
- data.table #Advanced dataframe ops
- readxl #Read .xl files
- tidyverse #Data organization and manipulation
- factoextra #Dimensionality reduction for PCA
- matrixStats #Required for PCA
- FactoMineR #Required for PCA
- randomcoloR #Large randomized color palettes
- mixOmics #Required for PLS
- reticulate #Required for single cell plots
- ggridges #Single cell ridge plots
```
See Dockerfile for specific versions and installation commands.

## Running Locally
3a) Install correct R version and required packages (listed above)

4a) Open in RStudio and run

## Deployment in Docker
3b) Install Docker from link below
Docker - [Version 18.09.3 build 77a1f4eee](https://www.docker.com/products/docker-desktop) (later versions not tested)

This section is performed within the Docker command line.
4b) Before starting, make sure you have the most recent version and are operating on the deployed node for the current version.
```
$ git pull https://gitlab.msu.edu/naultran/zed.git
$ git checkout FAIRTox_V4.0_live
```

5b) Run bash script
```
./docker_container_rebuild.sh
```

6b) To ensure container is running
```
$ docker ps
```

7b) To access the ShinyApp, navigate to dockerhost:80 in your web browser of choice.

## Built With
* [R Shiny](http://shiny.rstudio.com)
* [Docker](https://www.docker.com/)

Live version IP: 35.10.112.105 (requires connection to MSU campus wifi or VPN)

## Versioning
All version control is handled through this Gitlab page. Previous versions can be provided upon request.

## Authors
* Rance Nault - [email] naultran@msu.edu
* Jack Dodson - [email] dodsonj3@msu.edu

## License
No licensing information available at this time

## Acknowldgements
* Special thanks to Björn Bos - Dockerizing a ShinyApp

--------
