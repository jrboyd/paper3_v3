library(magrittr)
library(openxlsx)
library(xlsx)
source("setup.R")
source("module_STEIN_ChIPseq.R")
source("module_enhancers_Carroll.R")
source("functions_rGREAT.R")
get_res = function(res = bdgdiff_res, cell = "", mark = "", from = "", to = "", direction = ""){
  if(from != "") from = paste0("_", from, "_")
  if(to != "") to = paste0("_", to, "_")
  keep = grepl(cell, names(res)) &
    grepl(mark, names(res)) &
    grepl(from, names(res)) &
    grepl(to, names(res)) &
    grepl(direction, names(res))
  return(res[keep])
}


marks = c("H3K4ME3", "H3K27ME3", "H3K4AC", "H3K27AC")
# marks = marks[-1]
for(m in marks){
  setwd(wd)
  STEIN_SUBSET = get_res(res = STEIN_ChIPseq, cell = "MCF7", mark = m)
  
  #filter peak sets down to smallest set
  MIN_PEAKS = min(sapply(STEIN_SUBSET, length))
  if(MIN_PEAKS > 10^5) MIN_PEAKS = 10^5
  STEIN_SUBSET = lapply(STEIN_SUBSET, function(x){
    o = order(x$adj_p, decreasing = T)
    pos = 1:length(o)
    names(pos) = o
    pos = pos[as.character(1:length(pos))]
    keep = pos <= MIN_PEAKS
    x[keep]
  })
  
  compare_drugs = function(subject, query, with = T, esr1_only = F, enh_only = F){
    tmp = get_res(STEIN_SUBSET, from = subject)
    if(length(tmp) != 1) stop("bad subject match")
    subject_gr = tmp[[1]]
    tmp = get_res(STEIN_SUBSET, from = query)
    if(length(tmp) != 1) stop("bad query match")
    query_gr = tmp[[1]]
    is_cov = GRanges(coverage(GRangesList(subject_gr, query_gr)) > 0)
    subject_or_query = is_cov[is_cov$score]
    subject_name = paste(subject, query, sep = "-")
    if(enh_only){
      subject_name = paste("enh", subject_name, sep = "-")
      olaps = findOverlaps(query = query_gr, subject = GSE40129$MCF7_any_enhancer)
      query_gr = query_gr[queryHits(olaps)]
      olaps = findOverlaps(query = subject_or_query, subject = GSE40129$MCF7_any_enhancer)
      subject_or_query = subject_or_query[queryHits(olaps)]
    }
    if(esr1_only){
      subject_name = paste("ESR1", subject_name, sep = "-")
      olaps = findOverlaps(query = query_gr, subject = GSE40129$MCF7_e2_ESR1)
      query_gr = query_gr[queryHits(olaps)]
      olaps = findOverlaps(query = subject_or_query, subject = GSE40129$MCF7_e2_ESR1)
      subject_or_query = subject_or_query[queryHits(olaps)]
    }
    
    
    out = wrap_great(subject38 = subject_or_query, query38 = query_gr, 
                     subject_name = subject_name, query_name = query, with = with)
    return(out)
  }
  
  setwd(output_dir)
  sub_dir = paste0("GREAT_comparisons_from_ctrl_losses/", m)
  dir.create(sub_dir, recursive = T)
  setwd(sub_dir)
  
  try(detach("package:openxlsx", unload = T), silent = T)
  library(xlsx)
  
  # a = compare_drugs("bza", "ctrl", with = F)
  # a = compare_drugs("e2", "ctrl", with = F)
  # a = compare_drugs("e2bza", "ctrl", with = F)
  # a = compare_drugs("gc10", "ctrl", with = F)
  # a = compare_drugs("gc10bza", "ctrl", with = F)
  
  a = compare_drugs("ctrl", "bza", with = F)
  a = compare_drugs("ctrl", "e2", with = F)
  a = compare_drugs("ctrl", "e2bza", with = F)
  a = compare_drugs("ctrl", "gc10", with = F)
  a = compare_drugs("ctrl", "gc10bza", with = F)
  
  
  detach("package:xlsx", unload = T)
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
  hmap_col = rgb(colorRamp(c("white", "beige", "red", "magenta", "magenta", "magenta"))(0:100/100)/255)
  
  pdf(paste0(m, "_from_ctrl_msig_enrichments_loss.pdf"))
  heatmap.2(mat[kept,], margins = c(14,15), trace = "n", col = hmap_col, cexRow = .45)
  dev.off()
}