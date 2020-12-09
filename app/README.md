# Zacharewski Lab Expression DataViewer (ZED) ShinyApp
<img src="https://img.shields.io/badge/-draft-blue.svg"> 
Version 3.1

Initiatives to promote data FAIRness - data that is findable, accessible, 
interoperable, and reusable - are becoming a priority to funding agencies. 
These principles allow new knowledge to be derived from existing data with 
minimal additional investment. We have developed FAIRTox, an open source 
data exploration and visualization application for large gene expression 
datasets to enhance their FAIRness. FAIRTox is built using a R Shiny framework, 
chosen due to its wide use as a tool for omics data analysis, and incorporates 
Tableau visualizations for heightened user interaction. Additionally, targeted 
SQL queries enable the management and navigation of data maintained in the 
dbZach databases. With this application users can compare the expression of 
specific genes in response to various experimental factors such as circadian 
time as well as dose and duration of exposure to environmental contaminants. 
Furthermore, these transcriptomic changes can be evaluated across diverse datasets, 
enabling the development of novel hypotheses. Processing and accession of such 
large sets can prove difficult for investigators and was a major source of 
troubleshooting during development. Therefore, the implementation of new tools 
to simplify the process of database navigation, such as optimized querying and 
loading strategies, was crucial. Furthermore, careful attention was paid to the 
user experience. For example, page layout, input methods, responsiveness, 
and documentation were designed with the end user in mind. The result is an 
intuitive interface that serves bioinformaticians and novices alike.

## Getting Started -- Local Version
Pull ZED to a fresh repsitory. Everything you need is provided, including 
the full database. If lacking permission, contact either Rance Nault or 
Jack Dodson; their information is below.

Repo link: https://gitlab.msu.edu/naultran/zed.git

* Future versions will have support for database access through Michigan State University servers.

### Prerequisites
The following software is required and is available for download following 
the corresponding links.

R - [Version 3.5.1](https://www.r-project.org/)
```
Packages:
-dplyr
-ggnewscale
-ggm
-ggplot2
-plotly
-plyr
-reshape2
-RSQLite
-shiny
-shinycssloaders
-shinyjs
-stringr
-tibble
-UpSetR
```

In 'app' folder in repository, open the server.R file in RStudio, and click Run App.

## Getting Started -- Docker (Live) Version
* NOT required if running locally

Docker - [Version 18.09.3 build 77a1f4eee](https://www.docker.com/products/docker-desktop) (later versions not tested)

This section is performed within the Docker command line.

Before starting, make sure you have the most recent version and are operating on the deployed node for the current version.
```
$ git pull https://gitlab.msu.edu/naultran/zed.git
$ git checkout FAIRTox.v.1
```

Viewing any active instances
```
$ docker ps
```

If any instances are currently running, make sure to end them using
```
$ docker stop [instance_name]
$ docker rm [instance_name]
```

Building docker image:
```
$ docker build --no-cache -t shinyserver .
```
Note the period in the previous command

Running the image:
```
$ docker run -d -p 80:80 --name zed shinyserver
```

To access the ShinyApp, navigate to dockerhost:80 in your web browser of choice.

## Built With
* [R Shiny](http://shiny.rstudio.com)
* [Docker](https://www.docker.com/)

Live version 1.0 available at 35.10.112.105 (must be connected to MSU network or MSU VPN)

## Contributing
For issues or suggestions for feature updates, contact us at [email] notify+naultran-zed-7367-issue-@gitlab.msu.edu

## Versioning
All version control is handled through the Gitlab page. Previous versions can be provided upon request.

## Authors
* Rance Nault - [email] naultran@msu.edu
* Jack Dodson - [email] dodsonj3@msu.edu

## License
No licensing information available at this time

## Acknowldgements
* Special thanks to Björn Bos - Dockerizing a ShinyApp

--------
