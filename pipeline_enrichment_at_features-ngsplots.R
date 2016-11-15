library(magrittr)
source('setup.R')
# source("module_STEIN_ChIPseq.R")
source("module_enhancers_Carroll.R")
source("function_plot_ngsheatmap_extraData.R")
#consider enrichment by ngsplot and FE at promoters and enhancers
setwd(input_dir)
all_prof_files = dir("STEIN_ChIPseq_at_features/ngsplots/", full.names = T) %>% dir(full.names = T) %>% dir(full.names = T, patter = "RData")
enhancer_prof_files =  all_prof_files[grepl("Carroll_enhancers", all_prof_files)]
prof_names = basename(enhancer_prof_files) %>% 
  sub(pattern = "-E2", replacement = "Crl_e2") %>% 
  sub(pattern = "-Veh", replacement = "Crl_ctrl") %>% 
  sub(pattern = "_pooled.RData", replacement = "")
is_Carroll = grepl("Crl", prof_names)
Carroll_enhancer_prof_files =  enhancer_prof_files[is_Carroll]
names(Carroll_enhancer_prof_files) = prof_names[is_Carroll]
Carroll_enhancer_profs = lapply(Carroll_enhancer_prof_files, function(file){
  paste("load:", file) %>% print()
  load(file)
  prof = enrichList[[1]]
  rownames(prof) = sapply(strsplit(rownames(prof), ":"), function(x)x[1])
  return(prof)
})
# setwd(wd)

norm_prof = lapply(Carroll_enhancer_profs[c(1,4,2,5,3)], function(x){
  x = x[,41:61]
  # x = 2^x
  MIN = quantile(x, .2)
  MAX = quantile(x, .90)
  x = ifelse(x > MAX, MAX, x)
  x = ifelse(x < MIN, MIN, x)
  
  return(x)
})

lab = sub("MCF7Crl_", "", names(norm_prof)) %>% sub(pattern = "_", replacement = "\n")
# 
setwd(output_dir)
sel_ctrl = setdiff(GSE40129$MCF7_enhancer_sets$ctrl_only, GSE40129$MCF7_enhancer_sets$k4me3)
ctrl_enh_prof = lapply(norm_prof, function(x)x[sel_ctrl,])
pdf("Carroll_enhancers.final.ctrl.pdf")
res = heatmap.ngsplots_extraData2(ctrl_enh_prof, sideplot_lwd = 1, labels_below = lab, cex.col = .7, sidePlot_smoothing = 3,
                                  side_plots = list(1, 2:3, 4:5),
                                  profile_colors = c("black", "blue", "purple", "blue", "purple"))
dev.off()

sel_e2 = setdiff(GSE40129$MCF7_enhancer_sets$e2_only, GSE40129$MCF7_enhancer_sets$k4me3)
e2_enh_prof = lapply(norm_prof, function(x)x[sel_e2,])
pdf("Carroll_enhancers.final.e2.pdf")
res = heatmap.ngsplots_extraData2(e2_enh_prof, sideplot_lwd = 1, labels_below = lab, cex.col = .7, sidePlot_smoothing = 3,
                                  side_plots = list(1, 2:3, 4:5),
                                  profile_colors = c("black", "blue", "purple", "blue", "purple"))
dev.off()

sel_both = setdiff(GSE40129$MCF7_enhancer_sets$common, GSE40129$MCF7_enhancer_sets$k4me3) %>%
  setdiff(y = GSE40129$MCF7_enhancer_sets$e2_only) %>%
  setdiff(y = GSE40129$MCF7_enhancer_sets$ctrl_only)
both_enh_prof = lapply(norm_prof, function(x)x[sel_both,])
pdf("Carroll_enhancers.final.both.pdf")
res = heatmap.ngsplots_extraData2(both_enh_prof, sideplot_lwd = 1, labels_below = lab, cex.col = .7, sidePlot_smoothing = 3,
                                  side_plots = list(1, 2:3, 4:5),
                                  profile_colors = c("black", "blue", "purple", "blue", "purple"))
dev.off()

k4me3_free = setdiff(names(GSE40129$MCF7_any_enhancer), GSE40129$MCF7_enhancer_sets$k4me3)
