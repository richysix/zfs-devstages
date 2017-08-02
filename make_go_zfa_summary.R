#!/usr/bin/env Rscript

# load packages
packages <- c("topGO", "grid", "Rgraphviz", "R2HTML" )
for( package in packages ){
  suppressPackageStartupMessages( suppressWarnings( library(package, character.only = TRUE) ) )
}

rootPath <- getwd()
for( dir in c(file.path('output', 'biolayout-clusters-files') ) ){
  dirPath <- file.path( rootPath, dir )
  if( !dir.exists(dirPath) ){
    dir.create( dirPath )
  }
}

ensVersion <- 85

gene_annotation_file <- 'dataFiles/annotation.txt'
gene_annotation <- read.table(gene_annotation_file, sep="\t", header = FALSE, quote = '"')

# make named vector of gene names for later
geneNames <- as.character( gene_annotation[,7] )
names(geneNames) <- as.character( gene_annotation[,1] )

code2ontology <- list(
  "BP" = 'biological_process',
  "CC" = 'cellular_component',
  "MF" = 'molecular_function'
)
# Read in GO annotation
# see GO-biolayout.md for how these map files are produced
id2GOList <- list()
for( ontology in c('BP', 'CC', 'MF') ){
  mapFile <- file.path( 
    rootPath,
    'output',
    paste0(code2ontology[[ontology]], '.e', ensVersion, '.map' )
  )
  id2GOList[[ontology]] <- readMappings(file = mapFile)
}

# read in cluster data
clusterFile <- file.path(rootPath, 'dataFiles', 'zfs-grcz10.tpm_r-0.94_pearson.cluster.info.txt')
cluster_data <- read.table(clusterFile, sep="\t", header=TRUE)
rownames(cluster_data) <- as.character( cluster_data$gene_id )
# get cluster info
# split cluster_data by cluster
clusters <- split(cluster_data, cluster_data$cluster)

# read in topGO results
topgo_results_file <- file.path(rootPath, 'output/topgo/0.94-e85-fisher/topGOResults.sig.tsv')
topgo_results <- read.table(topgo_results_file, sep = "\t", header = TRUE)

# collapse GO terms down to highest-level term that is unique at that level
#collapse_GO <- function(cluster_go){
#  subtree_sizes <- table(cluster_go$SubtreeNumber)
#  top_unique_levels <- vector('list', length = length(subtree_sizes))
#  i <- 1
#  for( subtree_number in names(subtree_sizes)){
#    subtree <- cluster_go[cluster_go$SubtreeNumber == subtree_number, ]
#    if (subtree_sizes[subtree_number] == 1) {
#      top_unique_levels[[i]] <- subtree
#    }
#    tree_levels <- table(subtree$Level)
#    # print(tree_levels)
#    prev_level <- NULL
#    top_unique_level <- NULL
#    for( level in names(tree_levels) ){
#      if(tree_levels[[level]] > 1){
#        if(is.null(prev_level)){
#          cat(sprintf('Cluster: %s, Domain: %s, Subtree: %s, Level: %s\n', 
#                      cluster_go$cluster[1], cluster_go$domain[1], 
#                      subtree_number, tree_levels))
#          stop('Gone too far')
#        } else {
#          top_unique_level <- prev_level
#          break
#        }
#      } else {
#        prev_level <- level
#      }
#    }
#    if(is.null(top_unique_level)){
#      top_unique_level <- prev_level
#    }
#    # cat(sprintf('Cluster: %s Domain: %s Subtree: %s TopLevel: %s\n', 
#    #                   cluster_go$cluster[1], cluster_go$domain[1], 
#    #                   subtree_number, top_unique_level))
#    top_unique_levels[[i]] <- subtree[subtree$Level == top_unique_level, ]
#    i <- i + 1
#  }
#  return(do.call(rbind, top_unique_levels))
#}
#
#collapsed_go_results <- 
#  do.call(rbind,
#          lapply(split(topgo_results, topgo_results$cluster),
#            function(cluster_subset){
#                do.call(rbind,
#                        lapply(split(cluster_subset, cluster_subset$domain),
#                        collapse_GO))
#            }
#          )
#  )
#

zfa_file <- file.path(rootPath, 'output', 'zfa', '0.94-e85', 'zfa.e85.sig.tsv')
biolayout_zfa_data <- read.table(zfa_file, sep = "\t", header=TRUE, quote = "\"")

# get list of whole gene set
allGenes <- rownames(cluster_data)

# construct a factor for the sig genes
sigGenes <- rep(0, length(allGenes))
names(sigGenes) <- allGenes
sigGenes[ rownames(clusters[[1]]) ] <- 1
sigGenes <- factor(sigGenes)

# Make topGOdata object
# nodeSize=10 : prune GO hierarchy of terms associated with < 10 genes
BP_GOdata <- new("topGOdata", ontology='BP', allGenes=sigGenes,
                  annot=annFUN.gene2GO, 
                  gene2GO=id2GOList[['BP']], nodeSize=10)
MF_GOdata <- new("topGOdata", ontology='MF', allGenes=sigGenes,
              geneSel=topDiffGenes, annot=annFUN.gene2GO, 
              gene2GO=id2GOList[['MF']], nodeSize=10)
CC_GOdata <- new("topGOdata", ontology='CC', allGenes=sigGenes,
              geneSel=topDiffGenes, annot=annFUN.gene2GO, 
              gene2GO=id2GOList[['CC']], nodeSize=10)

GO_data_list <- list(
  'biological_process' = BP_GOdata,
  'cellular_component' = CC_GOdata,
  'molecular_function' = MF_GOdata
)

# create mapping of ZFA term to ENSDARG ids
zfa_2_zdb_file <- file.path(rootPath, 'output', 'zfa', 'annotations.txt')
zfa_2_zdb <- read.table(zfa_2_zdb_file, sep = "\t", header=FALSE, quote = "\"")
zfa_2_zdb <- zfa_2_zdb[, c(2,5)]
colnames(zfa_2_zdb) <- c('zdb_id', 'zfa_id')

zdb_2_ens_file <- file.path(rootPath, 'output', 'zfa', 'zdb-gene-2-ens-gene.e85.txt')
zdb_2_ens <- read.table(zdb_2_ens_file, sep = "\t", header=FALSE, quote = "\"")
colnames(zdb_2_ens) <- c('zdb_id', 'ensembl_id')

zfa_2_ens <- merge(zfa_2_zdb, zdb_2_ens)

cluster_info_list <- vector('list', length = length(clusters))
for( cluster_num in seq_len(length(clusters))){
  cluster_name <- sprintf('Cluster%03d', cluster_num)
  cat(cluster_name, "\n")
  cluster_size <- nrow(clusters[[cluster_num]])
  go_results <- collapsed_go_results[ collapsed_go_results$cluster == cluster_name, ]
  cluster_info <- clusters[[cluster_num]]
  genes_ids_in_cluster <- new.env()
  for( gene in rownames(cluster_info)){
    genes_ids_in_cluster[[gene]] <- TRUE
  }
  
  # make a page for each cluster that is just the list of genes assigned to that cluster
  cluster_genes <- cluster_info[, c('gene_id', 'name', 'chr', 'start', 'end', 'biotype')]
  cluster_genes$gene_id <- 
    paste0('<a href="http://www.ensembl.org/Danio_rerio/Gene/Summary?db=core;g=', cluster_genes$gene_id,
           '" target="_blank" title="* This link opens in a new window">', cluster_genes$gene_id,
           '</a>')
  colnames(cluster_genes) <- c('Ensembl ID', 'Gene Name', 'Chr', 'Start', 'End', 'Biotype')
  cluster_genes_file <- file.path(rootPath, 'output', 'biolayout-clusters-files', paste0(cluster_name, '-genes.html') )
  cluster_genes_relative <- paste0('biolayout-clusters-files/', cluster_name, '-genes.html')
  cluster_size_link <- paste0('<a href="', cluster_genes_relative, '">',
                              cluster_size, '</a>')
  css_file <- '"biolayout-clusters.css"'
  cat('<html>\n', paste0('<head><title>', cluster_name, '</title>\n'),
      '<meta charset="UTF-8"></head>\n',
      '<body>\n', file = cluster_genes_file, append = FALSE )
  HTMLCSS(file = cluster_genes_file, CSSfile = css_file)
  HTML.title(paste0(cluster_name, ': Genes'), file = cluster_genes_file, HR = 1)
  HTML(cluster_genes, file = cluster_genes_file, row.names = FALSE, Border = 0)
  cat('</body></html>\n', file = cluster_genes_file, append = TRUE )
  
  if(nrow(go_results) == 0){
    go_terms <- 'NA'
    cluster_link <- cluster_name
    cluster_table_file <- NULL
  } else {
    go_terms <- paste(go_results$Term, paste0('(', go_results$GO.ID, ')'), 
                      sep = ' ', collapse = '<br>')
    
    gene_ids_list <- lapply(seq_len(nrow(go_results)), 
                            function(x, go_results){ 
                              GO_obj <- GO_data_list[[ as.character(go_results$domain[x]) ]]
                              return(genesInTerm(GO_obj, as.character(go_results$GO.ID[x]))[[1]])
                            },
                            go_results
    )
    gene_ids_list <- lapply(gene_ids_list,
                            function(gene_ids, genes_ids_in_cluster){
                              # checks which genes in term exist in cluster and return char vector
                              unlist(lapply(gene_ids,
                                     function(gene){
                                       if(exists(gene, genes_ids_in_cluster)){
                                         return(gene)
                                       }
                                      }) )
                            },
                            genes_ids_in_cluster)
    
    gene_names_list <- lapply(gene_ids_list,
                              function(gene_ids){
                                return(geneNames[gene_ids])
                              }
    )
    gene_names <- sapply(gene_names_list,
                            function(gene_names){
                              paste0(gene_names, sep = '', collapse = '<br>') }
    )
    gene_ids <- sapply(gene_ids_list,
                            function(gene_ids){
                              paste0('<a href="http://www.ensembl.org/Danio_rerio/Gene/Summary?db=core;g=',
                                      gene_ids, '" target="_blank" title="* This link opens in a new window">', gene_ids, '</a>', collapse = '<br>')
                              }
    )
    
    
    cluster_detail <- go_results[ , c('GO.ID', 'Term', 'domain', 
                                      'Annotated', 'Expected', 'Observed', 'padj')]
    # make GO.ID a link
    cluster_detail$GO.ID <- paste0('<a href="http://www.ebi.ac.uk/QuickGO-Beta/term/', cluster_detail$GO.ID,
                                   '">', cluster_detail$GO.ID, '</a>')
    
    cluster_detail$genes <- gene_names
    cluster_detail$gene_ids <- gene_ids
    colnames(cluster_detail) <- c('GO ID', 'Description', 'Domain', 'Annotated', 'Expected', 'Observed',
                                  'Adjusted p-value', 'Genes', 'Ensembl IDs')
    
    # make cluster specific page (linked to from cluster name on index page)
    # don't finish page in case there are ZFA enrichments
    cluster_table_file <- file.path(rootPath, 'output', 'biolayout-clusters-files', paste0(cluster_name, '.html') )
    css_file <- '"biolayout-clusters.css"'
    cat('<html>\n', paste0('<head><title>', cluster_name, '</title>\n'),
      '<meta charset="UTF-8"></head>\n',
      '<body>\n', file = cluster_table_file, append = FALSE )
    HTMLCSS(file = cluster_table_file, CSSfile = css_file)
    HTML.title(paste0(cluster_name, ': Detail'), file = cluster_table_file, HR = 1)
    cat('<a id="GO" href="#ZFA"><h3>Go to ZFA detail</h3></a>', file = cluster_table_file, append = TRUE)
    cat('<h2>GO</h2>', file = cluster_table_file, append = TRUE)
    HTML(cluster_detail, file = cluster_table_file, row.names = FALSE, Border = 0)
    
    cluster_link <- paste0('<a href="biolayout-clusters-files/', cluster_name,
                     '.html">', cluster_name, '</a>')
  }
  
    
  # ZFA
  zfa_for_cluster <- biolayout_zfa_data[ biolayout_zfa_data$Cluster == cluster_name, ]
  if(nrow(zfa_for_cluster) == 0){
    num_zfa_terms <- 0
    zfa_terms <- 'NA'
    if( !is.null(cluster_table_file) ){
      cat('</body></html>\n', file = cluster_table_file, append = TRUE )
    }
  } else {
    num_zfa_terms <- nrow(zfa_for_cluster)
    zfa_terms <- paste(zfa_for_cluster$name, paste0('(', zfa_for_cluster$ID, ')'), 
                      sep = ' ', collapse = '<br>')
    # get genes for ZFA terms
    gene_ids_list <- lapply(
      zfa_for_cluster$ID,
      function(zfa_term){
        as.character(zfa_2_ens$ensembl_id[ zfa_2_ens$zfa_id == as.character(zfa_term) ])
      }
    )
    gene_ids_list <- lapply(gene_ids_list,
                            function(gene_ids, genes_ids_in_cluster){
                              # checks which genes in term exist in cluster and return char vector
                              unlist(lapply(gene_ids,
                                            function(gene){
                                              if(exists(gene, genes_ids_in_cluster)){
                                                return(gene)
                                              }
                                            }) )
                            },
                            genes_ids_in_cluster)
    
    gene_names_list <- lapply(gene_ids_list,
                              function(gene_ids){
                                return(geneNames[gene_ids])
                              }
    )
    
    gene_names <- sapply(gene_names_list,
                         function(gene_names){
                           paste0(gene_names, sep = '', collapse = '<br>') }
    )
    gene_ids <- sapply(gene_ids_list,
                            function(gene_ids){
                              paste0('<a href="http://www.ensembl.org/Danio_rerio/Gene/Summary?db=core;g=',
                                      gene_ids, '" target="_blank" title="* This link opens in a new window">', gene_ids, '</a>', collapse = '<br>')
                              }
    )
    
    # make data frame to output as a table
    zfa_detail <- zfa_for_cluster[, c('ID', 'name', 'TotalInTerm', 'Expected', 'Observed', 'FoldEnrichment', 'p.adjusted')]
    zfa_detail$genes <- gene_names
    zfa_detail$gene_ids <- gene_ids
    colnames(zfa_detail) <- c('ZFA ID', 'Description', 'Annotated', 'Expected', 'Observed', 
                              'Fold Enrichment', 'Adjusted p-value', 'Genes', 'Ensembl IDs')
    
    if( is.null(cluster_table_file) ){
      # make page and finish
      cluster_table_file <- file.path(rootPath, 'output', 'biolayout-clusters-files', paste0(cluster_name, '.html') )
      css_file <- '"biolayout-clusters.css"'
      cat('<html>\n', paste0('<head><title>', cluster_name, '</title>\n'),
        '<meta charset="UTF-8"></head>\n',
        '<body>\n', file = cluster_table_file, append = FALSE )
      HTMLCSS(file = cluster_table_file, CSSfile = css_file)
      HTML.title(paste0(cluster_name, ': Detail'), file = cluster_table_file)
      cat('<a id="ZFA" href="#GO"><h3>Go to GO detail<h3></a>', file = cluster_table_file, append = TRUE)
      cat('<h2>ZFA</h2>', file = cluster_table_file, append = TRUE)
      HTML(zfa_detail, file = cluster_table_file, row.names = FALSE, Border = 0)
      cat('</body></html>\n', file = cluster_table_file, append = TRUE )
      
      # set cluster link
      cluster_link <- paste0('<a href="biolayout-clusters-files/', cluster_name,
                     '.html">', cluster_name, '</a>')
    } else {
      # finish page
      cat('<a id="ZFA" href="#GO"><h3>Go to GO detail<h3></a>', file = cluster_table_file, append = TRUE)
      cat('<h2>ZFA</h2>', file = cluster_table_file, append = TRUE)
      HTML(zfa_detail, file = cluster_table_file, row.names = FALSE, Border = 0)
      cat('</body></html>\n', file = cluster_table_file, append = TRUE )
    }
  }
  
  cluster_profile_tag <- 
    paste0('<img src="',
            file.path('biolayout-clusters-files',
              sprintf('Cluster%03d.png', cluster_num) ),
            '">'
    )

  cluster_info_list[[cluster_num]] <- data.frame(
    Cluster = cluster_link,
    Size = cluster_size_link,
    Num_GO = nrow(go_results),
    GO_Terms = go_terms,
    Num_ZFA = num_zfa_terms,
    ZFA_Terms = zfa_terms,
    Profile = cluster_profile_tag
  )
}

cluster_info <- do.call(rbind, cluster_info_list)
colnames(cluster_info) <- c('Cluster', 'Size', 'Number GO Terms enriched',
    'Enriched GO Terms', 'Number ZFA Terms enriched', 'Enriched ZFA Terms',
    'Expression Profile')

html_table_file <- file.path(rootPath, 'output', 'biolayout-clusters.html')
css_file <- paste0('"', file.path('biolayout-clusters-files', 'biolayout-clusters.css'), '"' )
cat('<html>\n', '<head><title>BioLayout Clusters</title>\n',
  '<meta charset="UTF-8"></head>\n',
  '<body>\n', file = html_table_file, append = FALSE )
HTMLCSS(file = html_table_file, CSSfile = css_file)
HTML.title("Table of BioLayout Cluster Information", file = html_table_file)
HTML('Click on the cluster links for more detail on GO/ZFA enrichments.<br>For the full list of genes assigned to the cluster, click on the link in the size column.', file = html_table_file)
HTML(cluster_info, file = html_table_file, row.names = FALSE, Border = 0)
cat('</body></html>\n', file = html_table_file, append = TRUE )
