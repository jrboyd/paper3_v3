library(magrittr)
# library(openxlsx)
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

setwd(wd)
STEIN_SUBSET = get_res(res = STEIN_ChIPseq, cell = "MCF7", mark = "H3K27ME3")

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

# cov = coverage(GRangesList(STEIN_SUBSET))
# cov1 = coverage(STEIN_SUBSET[[1]]) + coverage(STEIN_SUBSET[[2]])
# olaps = findOverlaps(STEIN_SUBSET$MCF7_bza_H3K27AC, STEIN_SUBSET$MCF7_e2bza_H3K27AC)
# STEIN_SUBSET$MCF7_bza_H3K27AC[-queryHits(olaps)][500:520]
# 
# # vp_ctrl_v_e2 = venn_peaks(list(ctrl = STEIN_SUBSET$MCF7_ctrl_H3K27AC, e2 =STEIN_SUBSET$MCF7_e2_H3K27AC))
# is_cov = GRanges(coverage(GRangesList(STEIN_SUBSET[2:3])) > 0)
# ctrl_or_e2 = is_cov[is_cov$score]
# wrap_great(subject38 = ctrl_or_e2, query38 = STEIN_SUBSET$MCF7_e2_H3K27AC, 
#            subject_name = "ctrl-e2", query_name = "e2", with = T)
# 
# is_cov = GRanges(coverage(GRangesList(list(STEIN_SUBSET$MCF7_e2_H3K27AC, STEIN_SUBSET$MCF7_bza_H3K27AC))) > 0)
# e2_or_bza = is_cov[is_cov$score]
# wrap_great(subject38 = e2_or_bza, query38 = STEIN_SUBSET$MCF7_bza_H3K27AC, 
#            subject_name = "e2-bza", query_name = "bza", with = F)


library(xlsx)
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
sub_dir = "GREAT_comparisons_k27me3/to_ctrl_truncated"
dir.create(sub_dir, recursive = T)
setwd(sub_dir)
# e2bza_induced_over_e2 = compare_drugs("e2bza", "e2", with = F)
# e2_induced_over_e2bza = compare_drugs("e2", "e2bza", with = F)
# 
# only_e2_resp_res = compare_drugs("e2", "ctrl", with = F)
# 
# not_e2_resp_res = compare_drugs("ctrl", "e2", with = F)
# in_e2_resp_res = compare_drugs("ctrl", "e2", with = T)
# 
# bza_inhib_res = compare_drugs("e2", "bza", with = F)
# bza_notinhib_res = compare_drugs("e2", "bza", with = T)
# 
# only_e2_resp_res = compare_drugs("e2", "ctrl", with = F)
# ESR1_only_e2_resp_res = compare_drugs("e2", "ctrl", with = F, esr1_only = T)
# enh_only_e2_resp_res = compare_drugs("e2", "ctrl", with = F, enh_only = T)
# ESR1enh_only_e2_resp_res = compare_drugs("e2", "ctrl", with = F, esr1_only = T, enh_only = T)
# 
# compare_drugs("gc10", "e2", with = F)
# compare_drugs("e2", "gc10", with = F)

# a = compare_drugs("ctrl", "bza", with = F)
# a = compare_drugs("ctrl", "e2", with = F)
# a = compare_drugs("ctrl", "e2bza", with = F)
# a = compare_drugs("ctrl", "gc10", with = F)
# a = compare_drugs("ctrl", "gc10bza", with = F)
detach("package:openxlsx", unload = T)
library(xlsx)
a = compare_drugs(subject = "bza", query = "ctrl", with = F)
a = compare_drugs("e2", "ctrl", with = F)
a = compare_drugs("e2bza", "ctrl", with = F)
a = compare_drugs("gc10", "ctrl", with = F)
a = compare_drugs("gc10bza", "ctrl", with = F)

detach("package:xlsx", unload = T)
library(openxlsx)
ctrl_comparisons_files = dir(pattern = "ctrl no ctrl")
ctrl_comparisons = list()
for(f in ctrl_comparisons_files){
  print(f)
  res = read.xlsx(f, sheet = "MSigDB Oncogenic Signatures")  
  key = basename(f) %>% sub(pattern = ".xlsx", replacement = "")
  ctrl_comparisons[[key]] = res
}

pval_comparisons = lapply(ctrl_comparisons, function(x){
  return(x[,c("ID", "Hyper_Raw_PValue")])
})
pvals = pval_comparisons[[1]]
for(i in 2:length(pval_comparisons)){
  pvals = merge(pvals, pval_comparisons[[i]], by = "ID", all = T)
}
colnames(pvals)[2:6] = sapply(strsplit(names(fe_comparisons), "-"), function(x)x[1])
pvals_mat = as.matrix(pvals[,2:6])
rownames(pvals_mat) = pvals[,1]
pvals_mat = pvals_mat * nrow(pvals_mat)
pvals_min = apply(pvals_mat, 1, min)
kept = names(pvals_min)[pvals_min < .05]

fe_comparisons = lapply(ctrl_comparisons, function(x){
  return(x[,c("ID", "Hyper_Fold_Enrichment")])
})

x = fe_comparisons[[1]]
for(i in 2:length(fe_comparisons)){
  x = merge(x, fe_comparisons[[i]], by = "ID", all = T)
}
colnames(x)[2:6] = sapply(strsplit(names(fe_comparisons), "-"), function(x)x[1])
o = order(x$e2, decreasing = T)
x = x[o,]
mat = as.matrix(x[,2:6])
rownames(mat) = x[,1]
source("H:/R_workspace/jrb_R_scripts/heatmap.3-split.R")
source("H:/R_workspace/jrb_R_scripts/heatmap.3-kmeans_wrapper.R")
heatmap.3_kmeans_wrapper(mat, nclust = 6, margins = c(1,1))
hmap_col = rgb(colorRamp(c("beige", "red", "magenta"))(0:100/100)/255)

pdf("k27me3_msig_enrichments.pdf")
heatmap.2(mat[kept,], margins = c(8,15), trace = "n", col = hmap_col, cexRow = .65)
dev.off()
