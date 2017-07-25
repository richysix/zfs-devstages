#!/usr/bin/env Rscript

library(optparse)

# parse options
option_list <- list(
  make_option(c("-d", "--directory"), type="character", default='cwd',
              help="Working directory [default %default]" ),
  make_option(c("-O", "--outputDirectory"), type="character", default='results',
              help="Output directory [default %default]" ),
  make_option("--ontology", type="character", default='BP',
              help="The ontology to use, BP, MF, CC or All [default %default]"),
  make_option("--significanceLevel", type="numeric", default=0.05,
              help="Significance level for adjusted pvalue for sig.tsv file [default %default]"),
  make_option("--mtc", type="character", default='None',
              help="Multiple testing correction method [default %default]"),
  make_option("--test", type="character", default='ks',
              help="Type of test to do. Should be one of 'ks' or 'fisher' [default %default]"),
  make_option("--minClusterSize", type="integer", default=10,
              help="Significance level for adjusted pvalue for sig.tsv file [default %default]"),
  make_option("--enrichmentOnly", action="store_true", default=FALSE,
              help="Only test for enrichments"),
  make_option("--depletionOnly", action="store_true", default=FALSE,
              help="Only test for depletions"),
  make_option("--clusterName", type="character", default='All',
              help="The name of a cluster to run topGO on. Defaults to all levels of cluster within input file"),
  make_option("--clusterNameColumn", type="character", default='cluster',
              help="The column name containing the cluster ids. [default %default]"),
  make_option("--ensembl_version", type="integer", default=85,
              help="The Ensembl version used for the annotation of genes to GO terms. [default %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'topGO-multiple-clusters.R',
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
#       directory="output",
#       outputDirectory='topgo/0.94-e85',
#       ontology = 'MF',
#       significanceLevel = 0.05,
#       minClusterSize = 10,
#       enrichmentOnly = TRUE,
#       depletionOnly = FALSE,
#       clusterName = 'Cluster001',
#       clusterNameColumn = 'cluster',
#       ensembl_version = 85,
#       test = 'fisher'
#       verbose=TRUE ),
#   args = c('/lustre/scratch117/maz/team31/projects/rnaseq-zfs-stages-grcz10/dataFiles/zfs-grcz10.tpm_r-0.94_pearson.cluster.info.txt')
#)

# go through options
# working directory
if( cmd_line_args$options[['directory']] == 'cwd' ){
  cmd_line_args$options[['directory']] <- getwd()
}else{
  setwd( cmd_line_args$options[['directory']] )
}

# ontology
if (cmd_line_args$options[['ontology']] == 'All') {
  cmd_line_args$options[['ontology']] <- c('BP', 'MF', 'CC')
}

# enrichment vs depletion
if (cmd_line_args$options[['enrichmentOnly']] &
      cmd_line_args$options[['depletionOnly']] ) {
  stop("Options enrichmentOnly and depletionOnly cannot be specified together.\nTo test for both enrichment and depletion don't add either option to the command")
}

if( cmd_line_args$options[['verbose']] ){
  cat( "Working directory:", cmd_line_args$options[['directory']], "\n", sep=" " )
  cat( "Ontology:", cmd_line_args$options[['ontology']], "\n", sep=" " )
  cat( "Significance Level:", cmd_line_args$options[['significanceLevel']], "\n", sep=" " )
  cat( "Multiple Testing Correction:", cmd_line_args$options[['mtc']], "\n", sep=" " )
  cat( "Test statistic:", cmd_line_args$options[['test']], "\n", sep=" " )
  cat( "Minimum Cluster Size:", cmd_line_args$options[['minClusterSize']], "\n", sep=" " )
  cat( "Enrichment Only:", cmd_line_args$options[['enrichmentOnly']], "\n", sep=" " )
  cat( "Depletion Only:", cmd_line_args$options[['depletionOnly']], "\n", sep=" " )
  cat( "Cluster Name:", cmd_line_args$options[['clusterName']], "\n", sep=" " )
  cat( "Cluster Name Column:", cmd_line_args$options[['clusterNameColumn']], "\n", sep=" " )
  cat( "Ensembl Version:", cmd_line_args$options[['ensembl_version']], "\n", sep=" " )
}

clusterDataFile <- cmd_line_args$args[1]
cluster_data <- read.table(clusterDataFile, header=TRUE, sep="\t", row.names = 2)

# get list of whole gene set
allGenes <- rownames(cluster_data)

# prep topGO
# read in all 3 ontologies and keep in a list
code2ontology <- list(
  "BP" = 'biological_process',
  "CC" = 'cellular_component',
  "MF" = 'molecular_function'
)
# Read in GO annotation
id2GOList <- list()
for( ontology in cmd_line_args$options[['ontology']] ){
  mapFile <- file.path( 
    paste(code2ontology[[ontology]],
          paste0( 'e', cmd_line_args$options[['ensembl_version']] ),
           'map', sep="." )
  )
  id2GOList[[ontology]] <- readMappings(file = mapFile)
}

################################################################################
# FUNCTIONS

#' sig_genes_for_cluster
#' creates a factor indicating which genes are in the cluster
#' 
#' Takes the cluster data and a cluster name and creates a integer factor
#'    with 0/1 to indicate which genes are in the cluster
#'    
#' @param cluster_data data.frame, rownames should be gene ids, 
#'        must contain a column named cluster_number
#' @param cluster_name integer
#'
#' @return sigGenes factor
sig_genes_for_cluster <- function(cluster_data, cluster_name, test){
  if( test == 'fisher' ){
    sigGenes <- rep(0, nrow(cluster_data))
    names(sigGenes) <- rownames(cluster_data)
    sigGenes[ cluster_data$cluster == cluster_name ] <- 1
    return( factor(sigGenes) )
  } else {
    sigGenes <- rep(0.95, length(allGenes))
    names(sigGenes) <- allGenes
    sigGenes[ cluster_data$cluster == cluster_name ] <- 0.01
    return( sigGenes )
  }
}

extract_subtrees <- function(graph, nodes){
  # create induced graph from the supplied nodes
  inducedG <- topGO::inducedGraph(graph, nodes)
  totalNodes <- numNodes(inducedG)

  # build levels of induced graph
  inducedGLevels <- topGO::buildLevels(inducedG)

  # create reversed graph of induced graph (same nodes, reversed direction of edges)
  reversedG <- topGO::reverseArch(inducedG)

  # create an environment to keep a track of which nodes we have visited
  nodesVisited <- new.env(hash = T, parent = emptyenv())

  # start at the root of the induced graph and continue until all nodes have been seen
  subtreeList <- vector('list', length = length(nodes))
  i <- 1
  for( level in as.character(seq_len(inducedGLevels$noOfLevels)) ){
    # break out of the loop if we have seen all the supplied nodes already
    if( length(ls(nodesVisited)) == totalNodes ){
      break
    } else{
      for (node in inducedGLevels$level2nodes[[level]]){
        # cat(node, "\n")
        if ( length(ls(nodesVisited)) == totalNodes ){
          break
        } else if (exists(node, envir = nodesVisited)) {
          next
        } else {
          # if the node is in the nodes list then get the induced graph from that node to the leaves using the reversed graph
          # if not carry on
          if( node %in% nodes ){
            subtree <- topGO::inducedGraph(reversedG, node)
            # add any nodes in this induced graph to the nodes visited environment
            nodesInSubtree <- topGO::nodesInInducedGraph(reversedG, node)
            invisible( lapply(nodesInSubtree, function(node){ nodesVisited[[node]] <- TRUE }) )
            
            # if induced subtree has more than one node reverse it and add to subtrees list
            if (numNodes(subtree) > 1){
              subtreeList[[i]] <- topGO::reverseArch(subtree)
            } else {
              subtreeList[[i]] <- subtree
            }
            i <- i + 1
          } else {
            nodesVisited[[node]] <- TRUE
          }
        }
      }
    }
  }
  # tidy up and return subtrees list
  numTrees <- sum( sapply(subtreeList, function(x){ !is.null(x) } ) )
  subTrees <- vector('list', length = numTrees)
  for( i in seq_len(numTrees) ){
    subTrees[[i]] <- subtreeList[[i]]
  }
  return(subTrees)
}

## function to return the significant genes using a 'pvalue'
sigLevel <- 0.05
topDiffGenes <- function(allScore) {
  return(allScore < sigLevel )
}

#' runTopGO
#' function to run GO enrichment using topGO
#' 
#' Takes the counts list split by stage and a stage index and runs GO enrichment
#'    using topGO on the gene set in the counts matrix.
#'    Outputs a file of results of the pvalue associated with each GO Term and
#'    pdfs showing the position of the most significant terms within the GO graph
#'    
#' @param cluster_number integer
#' @param cluster_data data.frame, rownames should be gene ids, 
#'        must contain a column named cluster_number
#' @param GOdata topGO object
#'
#' @return allRes data.frame
#'
runTopGO <- function( cluster_name, cluster_data, GOdata, test ){
  # update topGO object with gene list for current cluster number
  sigGenes <- sig_genes_for_cluster(cluster_data, cluster_name, test)
  
  # Run topGO
  if( test == 'ks' ){
    GOdata <- updateGenes(GOdata, sigGenes, geneSelFun = topDiffGenes)
    result <- runTest(GOdata, algorithm="elim", statistic="ks")
  } else {
    GOdata <- updateGenes(GOdata, sigGenes)
    result <- runTest(GOdata, algorithm="weight01", statistic="fisher")
  }
  nodecount <- length(score(result))
  allRes <- GenTable(GOdata, pvalue=result, topNodes=nodecount)
  # add cluster name and ontology domain
  allRes$cluster <- rep(cluster_name, nrow(allRes))
  allRes$domain <- rep(code2ontology[[ontology(GOdata)]], nrow(allRes))
  
  # adjust p values
  # first convert '< 1e-30' to a number
  allRes$pvalue[ allRes$pvalue == '< 1e-30' ] <- 1e-30
  if( cmd_line_args$options[['mtc']] == 'Bonferroni' ){
    allRes$padj <- as.numeric(allRes$pvalue) * length(clusterNames) # convert from character to numeric and adjust for number of clusters tested
    allRes$padj[ allRes$padj > 1 ] <- 1 # set anything greater than 1 to 1
  } else{
    allRes$padj <- as.numeric(allRes$pvalue) # convert from character to numeric
  }
  # Horrible way to get all the genes associated with each term
  allRes$Genes <- sapply(allRes$GO.ID,
                         function(x) gsub('[c()" \n]', '', genesInTerm(GOdata, x)))
  # reorder columns
  allRes <- allRes[, c('cluster', 'GO.ID', 'Term', 'domain', 'Annotated', 'Expected', 'Significant', 'pvalue', 'padj', 'Genes')]
  colnames(allRes) <- c('cluster', 'GO.ID', 'Term', 'domain', 'Annotated', 'Expected', 'Observed', 'pvalue', 'padj', 'Genes')
  
  # make sig list
  # remove terms where there is a depletion. only interested in enrichments
  if (cmd_line_args$options[['enrichmentOnly']]){
    sigRes <- allRes[ allRes$Observed/allRes$Expected > 1 & !is.nan(allRes$Observed/allRes$Expected),
                     !grepl('Genes', colnames(allRes)) ]
  } else if (cmd_line_args$options[['depletionOnly']]){
    sigRes <- allRes[ allRes$Observed/allRes$Expected < 1 & !is.nan(allRes$Observed/allRes$Expected),
                     !grepl('Genes', colnames(allRes)) ]
  }
  sigRes <- sigRes[ sigRes$padj < cmd_line_args$options[['significanceLevel']], ]
  
  # get subtrees
  sigNodes <- as.character(sigRes[['GO.ID']])
  if ( length(sigNodes) > 0 ) {
    sigResSubtreesList <- extract_subtrees(graph(GOdata), sigNodes )
    
    inducedG <- topGO::inducedGraph(graph(GOdata), sigNodes)
    inducedGLevels <- topGO::buildLevels(inducedG)
    # go through substrees and get levels from induced graph of all terms
    terms <- character( length = length(sigNodes) )
    subtree_numbers <- integer( length = length(sigNodes) )
    levels <- integer( length = length(sigNodes) )
    i <- 1
    for( tree_number in seq_len(length(sigResSubtreesList)) ){
      # get subgraph
      subgraph <- sigResSubtreesList[[tree_number]]
      
      for( node in nodes(subgraph) ){
        terms[i] <- node
        subtree_numbers[i] <- tree_number
        levels[i] <- inducedGLevels$nodes2level[[node]]
        i <- i + 1
      }
    }
    
    levels_df <- data.frame(
        'GO.ID' = terms,
        'SubtreeNumber' = subtree_numbers,
        'Level' = levels
    )
    sigResults <- merge(sigRes, levels_df)
    # order by subtree number and level
    sigResults <- sigResults[ order( sigResults$SubtreeNumber, sigResults$Level ), ]
    # reorder columns
    sigResults <- sigResults[, c('cluster', 'GO.ID', 'Term', 'domain', 'Annotated', 'Expected', 'Observed', 'pvalue', 'padj', 'SubtreeNumber', 'Level')]
    
    if (cmd_line_args$options[['enrichmentOnly']]) {
      sigFile <- file.path( cmd_line_args$options[['outputDirectory']],
                          paste0(cluster_name, '.GO.', ontology, '.enriched.sig.tsv') )
    } else if (cmd_line_args$options[['depletionOnly']]){
      sigFile <- file.path( cmd_line_args$options[['outputDirectory']],
                      paste0(cluster_name, '.GO.', ontology, '.depleted.sig.tsv') )
    } else{
      sigFile <- file.path( cmd_line_args$options[['outputDirectory']],
                      paste0(cluster_name, '.GO.', ontology, '.sig.tsv') )
    }
    write.table( sigResults, file=sigFile, quote=FALSE,
        row.names=FALSE, sep="\t" )
    
    pdf(file.path( cmd_line_args$options[['outputDirectory']],
                      paste0(cluster_name, '.GO.', ontology, '.pdf') ))
    try(suppressWarnings(showSigOfNodes(GOdata, score(result),
        firstSigNodes=5, useInfo="all")), silent=TRUE)
    try(suppressWarnings(showSigOfNodes(GOdata, score(result),
        firstSigNodes=10, useInfo="all")), silent=TRUE)
    try(suppressWarnings(showSigOfNodes(GOdata, score(result),
        firstSigNodes=nrow(sigRes), useInfo="all")), silent=TRUE)
    try(suppressWarnings(lapply(sigRes[,1],
        function(x) showGroupDensity(GOdata, x))), silent=TRUE)
    dev.off()
  }
  
  allFile <- file.path( cmd_line_args$options[['outputDirectory']],
                        paste0(cluster_name, '.GO.', ontology, '.all.tsv') )
  write.table(allRes, file = allFile, quote = FALSE,
      row.names = FALSE, sep = "\t" )

  return(allRes)
}
  
################################################################################

# run topGO for domain specified in options (BP, MF, CC)
if ( cmd_line_args$options[['clusterName']] == 'All' ){
  for ( ontology in cmd_line_args$options[['ontology']] ){
    cluster_name <- levels(cluster_data$cluster)[1]
    sigGenes <- sig_genes_for_cluster(cluster_data, cluster_name, cmd_line_args$options[['test']])
    
    # work out number of clusters above min_cluster_size for
    clusterSizes <- table(cluster_data$cluster) >= cmd_line_args$options[['minClusterSize']]
    clusterNames <- names(clusterSizes)[ clusterSizes ]
    
    # Make topGOdata object
    # nodeSize=10 : prune GO hierarchy of terms associated with < 10 genes
    if( cmd_line_args$options[['test']] == 'ks'){
      GOdata <- new("topGOdata", ontology=ontology, allGenes=sigGenes,
                    geneSel=topDiffGenes, annot=annFUN.gene2GO, 
                    gene2GO=id2GOList[[ontology]], nodeSize=10)
    } else {
      GOdata <- new("topGOdata", ontology=ontology, allGenes=sigGenes,
                    annot=annFUN.gene2GO, 
                    gene2GO=id2GOList[[ontology]], nodeSize=10)
    }
    
    # go through clusters
    topGOResList <- lapply(clusterNames, runTopGO, cluster_data, GOdata, cmd_line_args$options[['test']] )
  }
} else{
  for ( ontology in cmd_line_args$options[['ontology']] ){
    cluster_name <- cmd_line_args$options[['clusterName']]
    sigGenes <- sig_genes_for_cluster(cluster_data, cluster_name, cmd_line_args$options[['test']])
    
    # work out number of clusters above min_cluster_size for
    clusterSizes <- table(cluster_data$cluster) >= cmd_line_args$options[['minClusterSize']]
    clusterNames <- names(clusterSizes)[ clusterSizes ]
    if ( !(cluster_name %in% clusterNames) ) {
      stop( sprintf('Cluster is smaller than minimum cluster size: %d\n', cmd_line_args$options[['minClusterSize']]) )
    }
    
    # Make topGOdata object
    # nodeSize=10 : prune GO hierarchy of terms associated with < 10 genes
    if( cmd_line_args$options[['test']] == 'ks'){
      GOdata <- new("topGOdata", ontology=ontology, allGenes=sigGenes,
                    geneSel=topDiffGenes, annot=annFUN.gene2GO, 
                    gene2GO=id2GOList[[ontology]], nodeSize=10)
    } else {
      GOdata <- new("topGOdata", ontology=ontology, allGenes=sigGenes,
                    annot=annFUN.gene2GO, 
                    gene2GO=id2GOList[[ontology]], nodeSize=10)
    }
    topGORes <- runTopGO(cluster_name, cluster_data, GOdata, cmd_line_args$options[['test']])
  }
}
