source("setup.R")
setwd(input_dir)
files = dir(path = "ngsplot_BRCA_enhancers_no_gb//", full.names = T)
attributes = matrix("", nrow = length(files), ncol = 3)
all_profs = list()
i = 1
for(f in files){
  print(f)
  attribs = strsplit(basename(f), "_")[[1]][1:3]
  cell = attribs[1]
  drug = attribs[2]
  mark = attribs[3]
  load(f)
  prof = enrichList[[1]]
  if(!exists("grps")){
    grps = sapply(strsplit(rownames(prof), "_"), function(x)paste(x[-length(x)], collapse = "_"))
  }
  rownames(prof) = paste(grps, 1:nrow(prof))
  all_profs[[paste(cell, drug, mark, sep = "_")]] = prof
  attributes[i, ] = c(cell, drug, mark)
  i = i + 1
}
attributes = as.data.frame(attributes)
colnames(attributes) = c("cell", "drug", "mark")
labels = apply(attributes, 1, function(x)paste(x, collapse = "\n"))

norm_profs = lapply(all_profs, function(x){
  MIN = quantile(x, .2)
  MAX = quantile(x, .98)
  x = ifelse(x < MIN, MIN, x)
  x = ifelse(x > MAX, MAX, x)
  x = (x - mean(x)) / sd(x)
})

source("H:/R_workspace/jrb_R_scripts/heatmap.ngsplots_kmeans_with_sideplot.R")
# heatmap.ngsplots(all_profs[1:4], labels_below = labels[1:4])
pdf("BRCA_enhancer_heatmaps_no_gb_12clust_sp.pdf", width = 8, height = 12)
for(cl in unique(attributes$cell)){
  for(m in unique(attributes$mark)){
    main_title = paste(cl, m)
    print(main_title)
    keep = attributes$cell == cl & attributes$mark == m
    heatmap.ngsplots(norm_profs[keep], labels_below = as.character(attributes$drug[keep]), main_title = main_title, nclust = 10)
  }
}
dev.off()
