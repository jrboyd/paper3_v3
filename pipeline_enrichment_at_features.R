library(magrittr)
source('setup.R')
source("functions_ngsheatmaps.R")
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
  return(enrichList[[1]])
})
setwd(wd)


sel = order(runif(nrow(Carroll_enhancer_profs[[1]])))[1:5000]
sel_prof = lapply(Carroll_enhancer_profs[c(1,4,2,5,3)], function(x){
  MIN = quantile(x, .1)
  MAX = quantile(x, .98)
  x = ifelse(x > MAX, MAX, x)
  x = ifelse(x < MIN, MIN, x)
  return(x[sel,])
})
pdf("Carroll_enhancers.5k.pdf")
res = heatmap.ngsplots(sel_prof)
dev.off()











count_types = dir("STEIN_ChIPseq_at_features/counts/drugs_with_merged_inputs/", full.names = T)
names(count_types) = basename(count_types)
count_files = lapply(count_types, function(x){
  dir(x, pattern = "pool", full.names = T)
})

# tmp = lapply(count_files, function(x){
#   x = basename(x)
#   sapply(strsplit(x, "\\."), function(y)y[1])
# })
count_data = lapply(count_files, function(x){
  n = basename(x)
  n = sapply(strsplit(n, "\\."), function(y)y[1])
  names(x) = n
  lapply(x, function(y){
    print(paste("load:", y))
    read.table(y, row.names = 1)
  })
})
