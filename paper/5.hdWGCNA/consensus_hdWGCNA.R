# single-cell analysis package
library(Seurat)
# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
# network analysis & visualization package:
library(igraph)
# using the cowplot theme for ggplot
theme_set(theme_cowplot())
# set random seed for reproducibility
set.seed(123)

setwd('/slurm/home/yrd/liaolab/nieshuyang/AHBA_sc')
seurat_obj <- readRDS('/scRNA/adult_ctx.rds')
DefaultAssay(seurat_obj) <- 'RNA'
# Remove useless information to avoid potential error
seurat_obj@assays[["sketch"]]<- NULL
seurat_obj@reductions[["pca"]] <- NULL
seurat_obj@reductions[["integrated.rpca"]] <- NULL
seurat_obj@reductions[["umap"]] <- NULL
seurat_obj@commands[["ScaleData.sketch"]] <- NULL
seurat_obj@graphs <- list()
# Scale data
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
# Add rownames and colnames
rownames(seurat_obj@assays[["RNA"]]@layers[["data"]]) <- rownames(seurat_obj)
colnames(seurat_obj@assays[["RNA"]]@layers[["data"]]) <- colnames(seurat_obj)
rownames(seurat_obj@assays[["RNA"]]@layers[["counts"]]) <- rownames(seurat_obj)
colnames(seurat_obj@assays[["RNA"]]@layers[["counts"]]) <- colnames(seurat_obj)
rownames(seurat_obj@assays[["RNA"]]@layers[["scale.data"]]) <- rownames(seurat_obj)
colnames(seurat_obj@assays[["RNA"]]@layers[["scale.data"]]) <- colnames(seurat_obj)
# Setting up and running metacells
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = 'Consensus_hdWGCNA'
)
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("Cluster", "Dataset"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'integrated.rpca.full', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'Clustrer', # set the Idents of the metacell seurat object
  min_cells = 200
)
seurat_obj <- NormalizeMetacells(seurat_obj)

# Customed hdWGCNA function to further filter datasets with few meta-cells
SetMultiExpr_custom <- function(
    seurat_obj, group_name, use_metacells=TRUE, group.by=NULL, 
    multi.group.by = NULL, multi_groups = NULL, assay=NULL,
    slot='data', layer = 'data', mat=NULL, mat_group_delim=3, 
    wgcna_name=NULL, ...){
  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)
  # get the WGCNA genes:
  params <- GetWGCNAParams(seurat_obj, wgcna_name)
  gene_names <- GetWGCNAGenes(seurat_obj, wgcna_name)
  s_obj <- seurat_obj
  # get assay
  if(is.null(assay)){
    assay <- DefaultAssay(s_obj)
    warning(paste0('assay not specified, trying to use assay ', assay))
  }
  # get the different groups present if not specified by the user:
  if(is.null(multi_groups)){
    multi_groups <- as.character(unique(s_obj@meta.data[[multi.group.by]]))
  } else{
    seurat_groups <- as.character(unique(s_obj@meta.data[[multi.group.by]]))
    if(sum(multi_groups %in% seurat_groups) != length(multi_groups)){
      stop('Some or all groups specified in multi_groups not found in seurat_obj@meta.data[,multi.group.by]')
    }
  }
  # was a matrix supplied?
  if(is.null(mat)){
    # use metacells or whole seurat object?
    if(use_metacells){
      s_obj <- GetMetacellObject(seurat_obj, wgcna_name)
    } else{
      s_obj <- seurat_obj
    }
    # get the datExpr for each group
    datExpr_list <- lapply(multi_groups, function(cur_group){
      cur_datExpr <- SetDatExpr(
        seurat_obj,
        group_name = group_name,
        group.by = group.by,
        multi.group.by = multi.group.by,
        multi_group_name = cur_group,
        return_seurat = FALSE,
        use_metacells = use_metacells,
        wgcna_name = wgcna_name,
        assay = assay,
        slot = slot,
        layer = layer
      ) 
      as.matrix(cur_datExpr)
    })
  } else{
    sample_groups <- do.call(rbind, strsplit(rownames(mat), ':'))[,mat_group_delim]
    datExpr_list <- list()
    for(cur_group in multi_groups){
      cur_datExpr <- as.data.frame(mat[which(sample_groups == cur_group),])
      datExpr_list[[cur_group]]<- cur_datExpr
    }
  }
  # convert to multiExpr, get good genes:
  multiExpr <- WGCNA::list2multiData(datExpr_list)
  # Exclude Dataset with less than 100 metacells
  metacell_num <- unlist(lapply(1:length(multiExpr), function(i) {
    return(dim(multiExpr[[i]][['data']])[1])}))
  valid_dataset <- metacell_num > 100
  multiExpr <- multiExpr[valid_dataset]
  multi_groups <- multi_groups[valid_dataset]
  genes_use <- WGCNA::goodGenesMS(multiExpr)
  gene_names <- gene_names[genes_use]
  # subset the multiExpr by the good genes::
  datExpr_list <- lapply(1:length(multiExpr), function(i){
    multiExpr[[i]]$data[,genes_use]
  })
  multiExpr <- WGCNA::list2multiData(datExpr_list)
  names(multiExpr) <- multi_groups
  # update the WGCNA gene list:
  seurat_obj <- SetWGCNAGenes(seurat_obj, gene_names, wgcna_name)
  # set the multiExpr in the Seurat object
  seurat_obj@misc[[wgcna_name]]$multiExpr <- multiExpr
  seurat_obj
}

# Extract metadata for meta-cells
metacell_meta <- seurat_obj@misc[["Consensus_hdWGCNA"]][["wgcna_metacell_obj"]]@meta.data
group = "L5/6 IT Car3"
# Set up for hdWGCNA
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = group,
  metacell_location = 'Consensus_hdWGCNA'
)
# Choose a subset of datasets
metacell_sub <- metacell_meta[metacell_meta$Cluster==group,]
multi_groups <- unique(metacell_sub$Dataset)
# Setup expression matrix
seurat_obj <- SetMultiExpr_custom(
  seurat_obj,
  group_name = group,
  group.by = "Cluster",
  multi.group.by ="Dataset",
  multi_groups = multi_groups
)
# Filter datasets with metacell number less than 100
multi_expr <- seurat_obj@misc[[group]][["multiExpr"]]
metacell_num <- unlist(lapply(1:length(multi_expr), function(i) {
  return(dim(multi_expr[[i]][['data']])[1])}))
valid_dataset <- metacell_num > 100
multi_expr <- multi_expr[valid_dataset]
seurat_obj@misc[[group]][["multiExpr"]] <- multi_expr
# Test best soft power
seurat_obj <- TestSoftPowersConsensus(seurat_obj)
# Pick optimal soft power
power_tables <- GetPowerTable(seurat_obj) %>% dplyr::group_split(group)
soft_power <- sapply(power_tables, function(power_table){
  power_table %>% subset(SFT.R.sq >= 0.8 & Power > 3) %>% .$Power %>% min
})
# Avoid infinite soft power value
inf_idx <- which(is.infinite(soft_power))
if (length(inf_idx) > 0) {
  soft_power <- soft_power[-inf_idx]
  seurat_obj@misc[[group]][["multiExpr"]] <- seurat_obj@misc[[group]][["multiExpr"]][-inf_idx]
}
# Construct network
seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power=soft_power,
  consensus=TRUE,
  tom_name = "Car3_consensus",
)
# Calculate MEs and kMEs
seurat_obj <- ModuleEigengenes(seurat_obj, group.by.vars="Dataset", wgcna_name = group)
seurat_obj <- ModuleConnectivity(seurat_obj, group.by='Cluster', group_name = group, wgcna_name = group)
# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = paste0(group, '-M'),
  wgcna_name = group
)
saveRDS(seurat_obj@misc, file="scRNA/Car3_consensus.rds")

## Module preservation analysis
# Split whole dataset into reference and test dataset
seurat_obj@misc[["active_wgcna"]] <- group
ref_datasets <- names(seurat_obj@misc[[cur_group]][["multiExpr"]])
cells_ref <- colnames(seurat_obj)[seurat_obj@meta.data$Dataset %in% ref_datasets]
cells_query <- setdiff(colnames(seurat_obj), cells_ref)
seurat_ref <- subset(seurat_obj, cells = cells_ref)
seurat_query <- subset(seurat_obj, cells = cells_query)
seurat_query@misc <- list()
# Project modules from a reference to a query dataset
seurat_query <- ProjectModules(
  seurat_obj = seurat_query,
  seurat_ref = seurat_ref,
  group.by.vars = 'Dataset',
  wgcna_name = group,
  wgcna_name_proj = group,
  assay = "RNA" # assay for query dataset
)
# Set expression matrix for reference dataset
seurat_ref <- SetDatExpr(
  seurat_ref,
  group_name = group,
  group.by = "Cluster"
)
# Set expression matrix for query dataset:
seurat_query <- SetDatExpr(
  seurat_query,
  group_name = group,
  group.by = "Cluster"
)
# Run module preservation function
seurat_query <- ModulePreservation(
  seurat_query,
  seurat_ref = seurat_ref,
  parallel = FALSE,
  name = group,
  verbose = 3,
  n_permutations = 200,
  seed = 123,
)
saveRDS(seurat_query, "scRNA/ModulePreservation.rds")

# Pathway enrichment analysis
dbs <- c(
  'GO_Biological_Process_2025',
  'GO_Cellular_Component_2025',
  'GO_Molecular_Function_2025',
  'GWAS_Catalog_2025'
  )
# Perform enrichment tests
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs,
  max_genes = 50
)
# Retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)
write.csv(enrich_df, "scRNA/Module_enrichr.csv")