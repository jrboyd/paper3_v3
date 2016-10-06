files = dir("H:/projects/Terri/paper3_v1/ngsplots_gb_no_prom_v3/", pattern = "MCF7", full.names = T)
keep = !grepl("MCF7_H", files); files = files[keep]
DE_ensg = rownames(DE_norm_counts)
all_DE_ngsprof = list()
for(f in files){
  print(f)
  name = basename(f) %>% sub("_pooled.RData", "", .)
  load(f)
  prof = enrichList[[1]]
  rownames(prof) = sapply(strsplit(rownames(prof), ":"), function(x)x[1])
  common = intersect(DE_ensg, rownames(prof))
  prof = prof[DE_ensg,]
  
  all_DE_ngsprof
}