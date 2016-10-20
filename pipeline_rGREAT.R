source("setup.R")
source("configure.R")
# source("module_DE_STEIN_RNAseq.R")
source("module_ESR1_and_enhancers.R")
source("functions_rGREAT.R")
# source("module_MCF7_enhancers.R")

# GSE40129$MCF7_e2_H3K4ME1
# GSE40129$MCF7_ctrl_H3K4ME1
# MCF7_H3K4ME1_venn = venn_peaks(list(ctrl = GSE40129$MCF7_ctrl_H3K4ME1, 
#                                     e2 = GSE40129$MCF7_e2_H3K4ME1))
# MCF7_H3K4ME1_and_ESR1_venn = venn_peaks(list(ESR1 = GSE40129$MCF7_e2_ESR1, 
#                                              H3K4ME1_e2 = GSE40129$MCF7_e2_H3K4ME1))

e2_olaps = venn_peaks(list(MCF7_e2_H3K4ME1 = GSE40129$MCF7_e2_H3K4ME1, MCF7_e2_H3K27AC = GSE40129$MCF7_e2_H3K27AC))
ctrl_olaps = venn_peaks(list(MCF7_ctrl_H3K4ME1 = GSE40129$MCF7_ctrl_H3K4ME1, MCF7_ctrl_H3K27AC = GSE40129$MCF7_ctrl_H3K27AC))

e2_enh = reduce(c(e2_olaps$`MCF7_e2_H3K4ME1 with MCF7_e2_H3K27AC` , e2_olaps$`MCF7_e2_H3K27AC with MCF7_e2_H3K4ME1`))
ctrl_enh = reduce(c(ctrl_olaps$`MCF7_ctrl_H3K4ME1 with MCF7_ctrl_H3K27AC` , ctrl_olaps$`MCF7_ctrl_H3K27AC with MCF7_ctrl_H3K4ME1`))

#combine enh sets, find middle and extend
any_enh = reduce(c(e2_enh, ctrl_enh))
widths = width(any_enh)
starts = start(any_enh)
mids = as.integer(starts + widths / 2)
ext = 5000
new_starts = mids - ext
new_ends = mids + ext
chrm = as.character(seqnames(any_enh))
enh_name = paste("enhancer", 1:length(chrm), sep = "_")
enh_bed = cbind(chrm, new_starts, new_ends, enh_name, rep(0, length(chrm)), rep(".", length(chrm)))
write.table(enh_bed, file = "Carroll_any_enhancers.bed", sep = "\t", row.names = F, col.names = F, quote = F)

enh_change = venn_peaks(list(e2_enh = e2_enh, ctrl_enh = ctrl_enh))
enh_gained = enh_change$`e2_enh no ctrl_enh`
enh_gained_bg = e2_enh
enh_lost = enh_change$`ctrl_enh no e2_enh`
enh_lost_bg = ctrl_enh

ESR1_peaks = GSE40129$MCF7_e2_ESR1
ESR1_peaks$enrichment = NULL #need to remove extra columns to do reduction
ESR1_peaks$adj_p = NULL
enh_set = enh_change[c(2,4)]
enh_set[["common_enh"]] = reduce(c(enh_change$`e2_enh with ctrl_enh`, enh_change$`ctrl_enh with e2_enh`))
i = 1
enh_ESR1_olaps = lapply(enh_set, function(x){
  input = list(x, ESR1_peaks)
  names(input) = c(names(enh_set)[i], "e2_ESR1")
  vp = venn_peaks(input)
  out = vp[c(2,4)]
  out[[paste(names(enh_set)[i], "with e2_ESR1")]] = reduce(c(vp[[1]], vp[[3]]))
  i <<- i + 1
  return(out)
})
enh_ESR1_olaps = unlist(enh_ESR1_olaps)
enh_ESR1_olaps = enh_ESR1_olaps[c(1,3,4,6,7,9)]
#e-e2 enh, c-ctrl enh, b-both ctrl and e2 enh
#E-ESR1 bound, n-not ESR1 bound
names(enh_ESR1_olaps) = c("en", "eE", "cn", "cE", "bn", "bE") 



# enh_gained_job = wrap_great(subject38 = e2_enh, 
#                             query38 = ctrl_enh, 
#                             subject_name = "e2_enh", 
#                             query_name = "ctrl_enh", 
#                             with = F)
# 
# enh_lost_job = wrap_great(subject38 = ctrl_enh, 
#                           query38 = e2_enh, 
#                           subject_name = "ctrl_enh", 
#                           query_name = "e2_enh", 
#                           with = F)
# 
# enh_ESR1_bound_job = wrap_great(subject38 = e2_enh, 
#                             query38 = GSE40129$MCF7_e2_ESR1, 
#                             subject_name = "e2_enh", 
#                             query_name = "e2_ESR1", 
#                             with = T)
# 
# enh_ESR1_absent_job = wrap_great(subject38 = e2_enh, 
#                                 query38 = GSE40129$MCF7_e2_ESR1, 
#                                 subject_name = "e2_enh", 
#                                 query_name = "e2_ESR1", 
#                                 with = F)
# 
# ESR1_at_enh_job = wrap_great(subject38 = GSE40129$MCF7_e2_ESR1, 
#                                  query38 = e2_enh, 
#                                  subject_name = "ESR1", 
#                                  query_name = "e2_enh", 
#                                  with = T)
# 
# ESR1_not_at_enh_job = wrap_great(subject38 = GSE40129$MCF7_e2_ESR1, 
#                              query38 = e2_enh, 
#                              subject_name = "ESR1", 
#                              query_name = "e2_enh", 
#                              with = F)
# 
# par(mfrow = c(1, 3))
# res_gained = plotRegionGeneAssociationGraphs(enh_gained_job, )
# res_lost = plotRegionGeneAssociationGraphs(enh_lost_job)
