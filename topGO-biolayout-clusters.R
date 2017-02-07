#!/usr/bin/env Rscript

library(optparse)

# parse options
option_list <- list(
  make_option(c("-d", "--directory"), type="character", default='cwd',
              help="Working directory [default %default]" ),
  make_option(c("-O", "--OutputDirectory"), type="character", default='results',
              help="Output directory [default %default]" ),
  make_option("--clusterName", type="character", default='All',
              help="The name of a cluster to run topGO on. Defaults to all levels of cluster within input file"),
  make_option("--ensembl_version", type="integer", default=85,
              help="The Ensembl version used for the annotation of genes to GO terms. [default %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'topGO-biolayout-clusters.R',
    usage = "Usage: %prog [options] clusterFile" ),
  positional_arguments = 1
)

# load packages
packages <- c("topGO", "grid", "Rgraphviz" )
for( package in packages ){
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

#cmd_line_args <- list(
#   options = list(
#       directory="/lustre/scratch109/sanger/is1/rnaseq-zfs-stages-grcz10/topgo",
#       OutputDirectory='0.94',
#       clusterName = 'Cluster001',
#       verbose=TRUE ),
#   args = c('change_with_stage_r-0.94_pearson.cluster.info.txt')
#)

clusterDataFile <- cmd_line_args$args[1]
clusterData <- read.table(clusterDataFile, header=TRUE, sep="\t")

# get list of whole gene set
allGenes <- as.character(clusterData$gene_id)

# prep topGO
# read in all 3 ontologies and keep in a list
code2ontology <- list(
  "BP" = 'biological_process',
  "CC" = 'cellular_component',
  "MF" = 'molecular_function'
)
# Read in GO annotation
id2GOList <- list()
for( ontology in c('BP', 'CC', 'MF') ){
  mapFile <- file.path( 
    cmd_line_args$options[['directory']], 
    paste(code2ontology[[ontology]],
          paste0( 'e', cmd_line_args$options[['ensembl_version']] ),
           'map', sep="." ) 
  )
  id2GOList[[ontology]] <- readMappings(file = mapFile)
}

# function to return the significant genes using a 'pvalue'
sigLevel <- 0.05
topDiffGenes <- function(allScore) {
  return(allScore < sigLevel )
}

# runTopGO
# function to run topGO for a given ontology on a subset of genes
# clusterInfo = data.frame. subset of genes (i.e. a cluster)
runTopGO <- function( clusterInfo, ontologyClass ){
  # get sig genes from cluster
  clusterName <- clusterInfo$cluster[1]
  sigGenes <- rep(0.95, length(allGenes))
  names(sigGenes) <- allGenes
  sigGenes[ as.character(clusterInfo$gene_id) ] <- 0.01

  # Make topGOdata object
  # nodeSize=10 : prune GO hierarchy of terms associated with < 10 genes
  GOdata <- new("topGOdata", ontology=ontologyClass, allGenes=sigGenes,
                geneSel=topDiffGenes, annot=annFUN.gene2GO, 
                gene2GO=id2GOList[[ontologyClass]], nodeSize=10)
  
  # Run topGO
  resultKS.elim <- runTest(GOdata, algorithm="elim", statistic="ks")
  nodecount <- length(score(resultKS.elim))
  allRes <- GenTable(GOdata, elimKS=resultKS.elim, topNodes=nodecount)
  # Horrible way to get all the genes associated with each term
  allRes$Genes <- sapply(allRes$GO.ID,
      function(x) gsub('[c()" \n]', '', genesInTerm(GOdata, x)))
  sigRes <- allRes[suppressWarnings(as.numeric(allRes$elimKS)) < sigLevel,]
  
  outputPrefix <- file.path( cmd_line_args$options[['directory']], 
                        cmd_line_args$options[['OutputDirectory']],
                        paste0(clusterName, '.GO.', ontologyClass)  )
  write.table( allRes, file=paste0(outputPrefix, ".all.tsv"), quote=FALSE,
      row.names=FALSE, sep="\t" )
  
  # Write PDF
  pdf(paste0(outputPrefix, ".pdf"))
  try(suppressWarnings(showSigOfNodes(GOdata, score(resultKS.elim),
      firstSigNodes=5, useInfo="all")), silent=TRUE)
  try(suppressWarnings(showSigOfNodes(GOdata, score(resultKS.elim),
      firstSigNodes=10, useInfo="all")), silent=TRUE)
  try(suppressWarnings(showSigOfNodes(GOdata, score(resultKS.elim),
      firstSigNodes=nrow(sigRes), useInfo="all")), silent=TRUE)
  try(suppressWarnings(lapply(sigRes[,1],
      function(x) showGroupDensity(GOdata, x))), silent=TRUE)
  dev.off()
}

# run topGO for each domain (BP, MF, CC)
topGOResList <- vector( "list", length = 3 )
if( cmd_line_args$options[['clusterName']] == 'All' ){
  for( ontology in c('BP', 'MF', 'CC' ) ){
    topGOResList[[ontology]] <- by(clusterData, clusterData$cluster, runTopGO, ontology )
  }
} else{
  # subset data to specified cluster
  clusterInfo <- clusterData[ clusterData$cluster == cmd_line_args$options[['clusterName']] &
                                !is.na( clusterData$cluster ), ]
  for( ontology in c('BP', 'MF', 'CC' ) ){
    runTopGO( clusterInfo, ontology )
  }
}

