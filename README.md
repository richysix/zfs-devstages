Notes to recreate the analysis in [White et al. 2017](https://doi.org/10.7554/eLife.30860)

### Contents
* GO-biolayout.md    Notes on running topGO on the gene clusters created using
[BioLayout Express](http://www.biolayout.org/)
* zfs-plots.Rmd      Rmarkdown file to reproduce the plots in the paper

### Counts Files

The zfs-plots Rmarkdown file expects two data files which do not exist in the repository for reasons of size.
To run the code you will need to first download the files below

* RNA-seq count data  
Download from [Figshare](https://doi.org/10.6084/m9.figshare.4704529)  
This needs to be named zfs-rnaseq-grcz10.tsv and put in the dataFiles/rnaseq/ directory

* DeTCT count data  
Download from [Figshare](https://doi.org/10.6084/m9.figshare.4622311)  
This needs to be named zfs-detct-grcz10.tsv and put in the dataFiles/detct/ directory

## Update

The RNA-seq data has been remapped to GRCz11 and counted against the Ensembl v111 annotation. Count files are available for download at [Figshare](https://doi.org/10.6084/m9.figshare.25858966)
