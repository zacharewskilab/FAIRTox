# Zacharewski Lab Expression DataViewer (ZED) ShinyApp Changelog
All notable changes will be documented in this file. 

Current up-to-date packages and versions - 2020-08-19
- askpass 1.1
- assertthat 0.2.1
- backports 1.1.8
- base64enc 0.1-3
- BH 1.72.0-3
- BiocGenerics 0.32.0
- BiocManager 1.30.10
- BiocVersion 3.10.1
- bit 4.0.4
- bit64 4.0.2
- blob 1.2.1
- callr 3.4.3
- cli 2.0.2
- colorspace 1.4-1
- commonmark 1.7
- scales 1.1.1
- shiny 1.5.0
- shinycssloaders 1.0.0
- shinyjs 1.1
- snow 0.4-3
- sourcetools 0.1.7
- stringi 1.4.6
- stringr 1.4.0
- sys 3.4
- testthat 2.3.2
- tibble 3.0.3
- tidyr 1.1.1
- tidyselect 1.1.0
- UpSetR 1.4.0
- utf8 1.1.4
- V8 3.2.0
- vctrs 0.3.2
- viridisLite 0.3.0
- withr 2.2.0
- xfun 0.16
- xtable 1.8-4
- yaml 2.2.1
- littler 0.3.8

## [Unreleased]
[3.2.0] - 2020-09-14
### Added
- Table in Dose-response tab now includes "QC" score for dataset completeness
- Table in Time-course tab now includes "QC" score for dataset completeness
- New table added to circadian-regulated tab to include "QC" tab for dataset completeness
- New major tab: Principal component analysis (PCA) for LaPres/Kaminski/Zacharewski counts files (currently on NAS)

[3.1.0] - 2020-08-18
### Added
- Changelog for all versions up to and including V3.1.0

### Changed
- Minor naming and labelling details app-wide
- Minor color updates for a more consistent UI theme

[3.0.0] - 2020-08-14
### Added
- shinyProxy implementation that allows for password protection
- Gene Information table containing external links to repos
- Confirmed use with up to 20 consecutive users

### Changed
- Updated Dockerfile to allow successful builing on office linux server

### Removed
- Dependency on VennDiagram package (no longer in use)

[2.0.0] - 2020-07-20
### Added
- groupBy functionality for GSEA. Works for species and dose
- GSEA tab now displays NES barplot and GSEA 4-layer plots for each subset
- V1.0 of completed GSEA tab up and running'
- Comments, formatting to GSEA code
- Added legends to DR and TC with plans of updating

### Changed
- Fixed circCounts and BE plots
- Updated circadian SQL command
- Updated time course SQL command

[1.9.0] - 2020-06-10
### Added
- GSEA now accepting uploaded files (tsv only)
- Heatmap now accepting uploaded files (tsv only)
- Asterisks now indicate significance in GSEA plot
- Comments to heatmap code
- All metadata filters for heatmap tab 

### Changed
- Refactored heatmap and GSEA code
- Heatmaps now only created for selected projects to speed up runtime
- Updated heatmap UI

[1.8.0] - 2020-05-25
### Added
- GSEA framework
- Fully functional GSEA for tzlab datasets

### Removed
- Testing file used for GSEA debugging

[1.7.0] - 2020-05-18
### Added
- Full support for UpSet plots with circadian data, filterable by ZT
- Gene List DataViewer UI framework created
- Heatmap in GLD now working with project 150 (code from Rance's tzheatmap package)

### Changed
- Greatly streamlines UpSet plotting code
- Issue with metadata filtering resolved
- Refactored new code since V1.5.0

### Removed
- Unnecessary print statements in server and heatmap testing files

[1.6.0] - 2020-05-04
### Added
- Global filter function, placed in separate file
- New source file containing all dataframe building functions

### Changed
- UpSet update: only show organs that apply to project(s) selected
- Updated SQL commands to pull more data and be more efficient

[1.5.0] - 2020-04-29
### Added
- Comments app-wide
- Organ metadata now included in UpSet tab

### Changed
- Plot button on UpSet tab now floats
- UI unique to UpSet tab
- Fold change filter now a text input box
- Each project in UpSet tab now has its own p1t filter
- Plot export code now much more streamlined

[1.4.0] - 2020-04-13
### Added
- Modular projects divs in UpSet plot with framework still working
- Hide/show filter panel when UpSet tab is selected

### Changed
- Now displaying only chemical metadata present in ExpressionChange table
- New framework for UpSet plots for better modularity and efficiency

[1.3.0] - 2020-04-03
### Changed
- Global code review across all files and functions
- Refactored function/variable names

[1.2.0] - 2020-03-23
### Added
- Settings tab: will be used to modify plot appearances without cluttering 
  other tabs
- Functionality to show/hide metadata filters in left side panel
- All plots are now resizable and exportable through settings tab

[1.1.0] - 2020-03-15
### Changed
- Converted venn diagrams to UpSet plots
- Resized plots app-wide to better fit browser window

### Removed
- Unnecessary csv output files used for debugging

[1.0.0] - 2020-03-12
### Added
- First working commit up and running on remote linux server inside Docker container
- IP address: 35.10.112.105:80

### Changed
- Dockerfile to include most up to date packages and dependencies
