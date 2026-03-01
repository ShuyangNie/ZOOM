## Null simulation: VAM and UCell
library(Seurat)
library(VAM)
library(progress)
library(UCell)

setwd("/slurm/home/yrd/liaolab/nieshuyang/AHBA_sc")
seurat_obj <- readRDS("scRNA/adult_ctx.rds")
rownames(seurat_obj@assays[["RNA"]]@layers[["data"]]) <- rownames(seurat_obj)
colnames(seurat_obj@assays[["RNA"]]@layers[["data"]]) <- colnames(seurat_obj)
rownames(seurat_obj@assays[["RNA"]]@layers[["counts"]]) <- rownames(seurat_obj)
colnames(seurat_obj@assays[["RNA"]]@layers[["counts"]]) <- colnames(seurat_obj)
# Sample 100,000 cells
cell_sub <- read.csv("scRNA/benchmarking_100000cells.csv",row.names=1)
cell_sub <- rownames(cell_sub)
seurat_obj <- subset(seurat_obj, cells = cell_sub)
rownames(seurat_obj@assays[["RNA"]]@layers[["data"]]) <- rownames(seurat_obj)
colnames(seurat_obj@assays[["RNA"]]@layers[["data"]]) <- colnames(seurat_obj)
rownames(seurat_obj@assays[["RNA"]]@layers[["counts"]]) <- rownames(seurat_obj)
colnames(seurat_obj@assays[["RNA"]]@layers[["counts"]]) <- colnames(seurat_obj)
gc()

# This is a copy of VAM function, we excluded some potential error
VAM_Seurat <- function (seurat.data, gene.set.collection, 
                        center = FALSE, gamma = TRUE) {
  # Extract normalized data expression
  normalized.counts = seurat.data@assays$RNA@layers$data
  colnames(normalized.counts) = colnames(seurat.data)
  
  # Estimate technical variation and assign binary weights
  tech.var.prop = getTechVarPropV5(seurat.data)
  names(tech.var.prop) = rownames(seurat.data)
  gene.weights = rep(1, nrow(normalized.counts))
  names(gene.weights) = rownames(seurat.data)
  
  # Transpose
  gene.expr = Matrix::t(normalized.counts)
  rm(normalized.counts)
  p = ncol(gene.expr)
  n = nrow(gene.expr)
  
  # Extract information
  cell.ids = rownames(gene.expr)
  num.sets = length(gene.set.collection)
  set.names = names(gene.set.collection)
  gene.ids = colnames(gene.expr)
  set.sizes = unlist(lapply(gene.set.collection, length))
  min.set.size = min(set.sizes)
  median.set.size = median(set.sizes)
  results = list()
  results$distance.sq = matrix(0, nrow = n, ncol = num.sets, 
                               dimnames = list(cell.ids, set.names))
  results$cdf.value = matrix(0, nrow = n, ncol = num.sets, 
                             dimnames = list(cell.ids, set.names))
  message("Computing VAM distances for ", num.sets, " gene sets, ", 
          n, " cells and ", p, " genes.")
  message("Min set size: ", min.set.size, ", median size: ", 
          median.set.size)
  
  pb <- progress_bar$new(total = num.sets, format = "  [:bar] :percent :elapsed")
  for (i in 1:num.sets) {
    set.members = gene.set.collection[[i]]
    set.size = set.sizes[i]
    set.exprs = gene.expr[, set.members]
    set.weights = gene.weights[set.members]
    names(tech.var.prop) = colnames(gene.expr)
    vam.results = vam(gene.expr = set.exprs, tech.var.prop = tech.var.prop[set.members], 
                      gene.weights = set.weights, center = center, 
                      gamma = gamma)
    
    results$distance.sq[, i] = vam.results$distance.sq
    results$cdf.value[, i] = vam.results$cdf.value
    pb$tick()
  }
  return(results)
}
getTechVarPropV5 <- function (seurat.data) 
{
  if (!is(seurat.data@assays$RNA, "Assay5")) {
    stop("Trying to call Assay5 logic on non-V5 object!")
  }
  p = nrow(seurat.data@assays$RNA@layers$data)
  meta.data.cols = colnames(seurat.data@assays$RNA@meta.data)
  if (length(meta.data.cols) <= 3 | meta.data.cols[2] != "vf_vst_counts_variance" | 
      meta.data.cols[3] != "vf_vst_counts_variance.expected") {
    message("Did not find vst variance decomposition, setting technical variance proportion to 1")
    return(rep(1, p))
  }
  tech.var.prop = seurat.data@assays$RNA@meta.data[, 3]/seurat.data@assays$RNA@meta.data[, 
                                                                                         2]
  tech.var.prop[which(is.nan(tech.var.prop))] = 1
  return(tech.var.prop)
}

gene_nums <- c(30,50,75,100,150,200,250,300,400,500)
for (gene_num in gene_nums) {
  signatures <- vector("list", 100)
  names(signatures) <- paste0("Null_", 1:100)
  for (n in 1:100) {
    fname <- sprintf("Benchmark/Null_simulation/Null_SBP/null%d_PLS_report.csv", n)
    df <- read.csv(fname, row.names = 1)
    df <- df[df$Sign>0,]
    df <- df[1:gene_num,]
    signatures[[n]] <- rownames(df)
  }
  # VAM
  VAM_score <- VAM_Seurat(seurat.data=seurat_obj, gene.set.collection=signatures)
  VAM_score_df <- as.data.frame(VAM_score[["cdf.value"]])
  pvals_df <- 1 - VAM_score_df
  write.csv(pvals_df, file=paste0("Benchmark/Null_simulation/VAM_pvals_",gene_num,"genes.csv"))
  # UCell
  seurat_obj_UCell <- AddModuleScore_UCell(seurat_obj, features=signatures, name=NULL)
  # Compute p-values
  df <- seurat_obj_UCell@meta.data
  null_cols <- paste0("Null_", 1:100)
  norm_scores <- list()
  norm_pvals <- list()
  for (col in null_cols) {
    v_score <- df[[col]]
    v_norm_score <- (v_score - mean(v_score, na.rm = TRUE)) / sd(v_score, na.rm = TRUE)
    v_norm_p <- 1 - pnorm(v_norm_score)
    norm_scores[[col]] <- v_norm_score
    norm_pvals[[col]] <- v_norm_p
  }
  norm_pvals_df <- as.data.frame(norm_pvals)
  rownames(norm_pvals_df) <- cell_sub
  colnames(norm_pvals_df) <- paste0("pval_", 1:100)
  write.csv(norm_pvals_df, file=paste0("Benchmark/Null_simulation/UCell_pvals_",gene_num,"genes.csv"))
}