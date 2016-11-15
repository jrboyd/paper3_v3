library(xlsx)
library(magrittr)
source("H:/R_workspace/jrb_R_scripts/heatmap.3-split.R")
source("H:/R_workspace/jrb_R_scripts/heatmap.3-kmeans_wrapper.R")

replot_GREAT_matrix = function(working_directory, pdfname, n_low = 1, n_high = 1, cexRow = .5, is_enhancer = F){
  start_wd = getwd()
  setwd(working_directory)
  try(
    detach("package:xlsx", unload = T), silent = T)
  library(openxlsx)
  ctrl_comparisons_files = dir(pattern = ".xlsx")
  ctrl_comparisons = list()
  for(f in ctrl_comparisons_files){
    print(f)
    res = read.xlsx(f, sheet = "MSigDB Oncogenic Signatures")  
    key = basename(f) %>% sub(pattern = ".xlsx", replacement = "")
    ctrl_comparisons[[key]] = res
  }
  cnames = sapply(strsplit(names(ctrl_comparisons), "[- ]"), function(x)paste(x[-2], collapse = " "))
  if(is_enhancer) cnames = sapply(strsplit(names(ctrl_comparisons), "[- ]"), function(x)paste(x[c(-1,-3)], collapse = " "))
  pval_comparisons = lapply(ctrl_comparisons, function(x){
    return(x[,c("ID", "Hyper_Raw_PValue")])
  })
  pvals = pval_comparisons[[1]]
  for(i in 2:length(pval_comparisons)){
    pvals = suppressWarnings(merge(pvals, pval_comparisons[[i]], by = "ID", all = T))
  }
  colnames(pvals)[2:ncol(pvals)] = cnames
  pvals_mat = as.matrix(pvals[,2:ncol(pvals)])
  rownames(pvals_mat) = pvals[,1]
  pvals_mat = pvals_mat * nrow(pvals_mat)
  pvals_min = apply(pvals_mat, 1, min)
  kept = names(pvals_min)[pvals_min < .05]
  
  fe_comparisons = lapply(ctrl_comparisons, function(x){
    return(x[,c("ID", "Hyper_Fold_Enrichment")])
  })
  
  x = fe_comparisons[[1]]
  for(i in 2:length(fe_comparisons)){
    x = suppressWarnings(merge(x, fe_comparisons[[i]], by = "ID", all = T))
  }
  colnames(x)[2:ncol(x)] = cnames
  o = order(x[,2], decreasing = T)
  x = x[o,]
  mat = as.matrix(x[,2:ncol(x)])
  rownames(mat) = x[,1]
  source("H:/R_workspace/jrb_R_scripts/heatmap.3-split.R")
  source("H:/R_workspace/jrb_R_scripts/heatmap.3-kmeans_wrapper.R")
  # heatmap.3_kmeans_wrapper(mat, nclust = 6, margins = c(1,1))
  hmap_col = rgb(colorRamp(c(rep("white", n_low), "beige", "red", rep("magenta", n_high)))(0:100/100)/255)
  
  pdf(pdfname)
  heatmap.2(mat[kept,], margins = c(14,15), trace = "n", col = hmap_col, cexRow = cexRow)
  dev.off()
  setwd(start_wd)
}

#working dir list should be named
working_directory_list = list("MCF7" = "output/GREAT_comparisons_enhancers_from_ctrl/H3K27AC/", "MCF10A" = "output/GREAT_comparisons_MCF10A_enhancers_from_ctrl/H3K27AC/")
replot_GREAT_matrix_combine_dirs = function(working_directory_list, pdfname, n_low = 1, n_high = 1, cexRow = .5, is_enhancer = F, ...){
  start_wd = getwd()
  ctrl_comparisons = list()
  all_pvals = list()
  for(prefix in names(working_directory_list)){
    working_directory = working_directory_list[[prefix]]
    setwd(working_directory)
    try(
      detach("package:xlsx", unload = T), silent = T)
    library(openxlsx)
    ctrl_comparisons_files = dir(pattern = ".xlsx")
    
    for(f in ctrl_comparisons_files){
      print(f)
      res = read.xlsx(f, sheet = "MSigDB Oncogenic Signatures")  
      key = basename(f) %>% sub(pattern = ".xlsx", replacement = "")
      key = paste(prefix, key)
      ctrl_comparisons[[key]] = res
    }
    
    setwd(start_wd)
  }

  cnames = sapply(strsplit(names(ctrl_comparisons), "[- ]"), function(x)paste(x[-3], collapse = " "))
  if(is_enhancer) cnames = sapply(strsplit(names(ctrl_comparisons), "[- ]"), function(x)paste(x[c(-2,-4)], collapse = " "))
  # cnames = paste(prefix, cnames)
  pval_comparisons = lapply(ctrl_comparisons, function(x){
    return(x[,c("ID", "Hyper_Raw_PValue")])
  })
  pvals = pval_comparisons[[1]]
  for(i in 2:length(pval_comparisons)){
    pvals = suppressWarnings(merge(pvals, pval_comparisons[[i]], by = "ID", all = T))
  }
  colnames(pvals)[2:ncol(pvals)] = cnames
  # all_pvals[[prefix]] = pvals
  # pvals = merge(all_pvals[[1]], all_pvals[[2]], by = "ID")
  
  
  pvals_mat = as.matrix(pvals[,2:ncol(pvals)])
  rownames(pvals_mat) = pvals[,1]
  pvals_mat = pvals_mat * nrow(pvals_mat)
  pvals_min = apply(pvals_mat, 1, min)
  kept = names(pvals_min)[pvals_min < .05]
  
  fe_comparisons = lapply(ctrl_comparisons, function(x){
    return(x[,c("ID", "Hyper_Fold_Enrichment")])
  })
  
  x = fe_comparisons[[1]]
  for(i in 2:length(fe_comparisons)){
    x = suppressWarnings(merge(x, fe_comparisons[[i]], by = "ID", all = T))
  }
  colnames(x)[2:ncol(x)] = cnames
  o = order(x[,2], decreasing = T)
  x = x[o,]
  mat = as.matrix(x[,2:ncol(x)])
  rownames(mat) = x[,1]
  
  over_ctrl = grepl("no ctrl", colnames(mat))
  colnames(mat)[over_ctrl] = sub(" no ctrl", "+", colnames(mat)[over_ctrl])
  colnames(mat)[!over_ctrl] = sub(" ctrl no", "", colnames(mat)[!over_ctrl])
  colnames(mat)[!over_ctrl] = paste(colnames(mat)[!over_ctrl], "-")
  

  hmap_col = rgb(colorRamp(c(rep("white", n_low), "beige", "red", rep("magenta", n_high)))(0:100/100)/255)
  o = order(colnames(mat))
  mat = mat[,o]
  o = order(grepl("\\+", colnames(mat)))
  mat = mat[,o]
  setwd(start_wd)
  pdf(pdfname)
  heatmap.2(mat[kept,], margins = c(14,15), trace = "n", col = hmap_col, cexRow = cexRow, ...)
  dev.off()
  setwd(start_wd)
}