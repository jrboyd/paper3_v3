library(magrittr)
source('setup.R')
# source("module_STEIN_ChIPseq.R")
source("module_enhancers_Carroll.R")
source("function_plot_ngsheatmap_extraData.R")
#consider enrichment by ngsplot and FE at promoters and enhancers
setwd(input_dir)

setwd(input_dir)
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

bed_files = dir("STEIN_ChIPseq_at_features/", pattern = ".bed", full.names = T)
bed_dfs = list()
bed_grs = list()
for(f in bed_files){
  key = f %>% basename %>% 
    sub(pattern = ".bed", replacement = "")
  bed = read.table(f)
  colnames(bed) = c("chrm", "start", "end", "id", "score", "strand")
  bed$ucsc = paste0(bed$chrm, ":", bed$start, "-", bed$end)
  strand = bed$strand
  strand = sub(".", "*", strand)
  gr = GRanges(seqnames = bed[,1], ranges = IRanges(start = bed[,2], end = bed[,2]), strand = strand)
  gr$ucsc = paste0(bed[,1], ":", bed[,2], "-", bed[,3])
  names(gr) = bed[,4]
  rownames(bed) = as.character(bed$id)
  bed_dfs[[key]] = bed
  bed_grs[[key]] = gr
}


blacklist_files = character()#dir("STEIN_ChIPseq_at_features/", pattern = ".txt", full.names = T)


#load blacklists and remove
for(f in blacklist_files){
  key = f %>% basename %>% 
    sub(pattern = "blacklist_", replacement = "", f) %>% 
    sub(pattern = ".txt", replacement = "")
  blacklist = read.table(f, stringsAsFactors = F)[,1]
  whitelist = setdiff(rownames(count_data[[key]][[1]]), blacklist)
  count_data[[key]] = lapply(count_data[[key]], function(x){
    return(x[whitelist,, drop = F])
  })
}

count_raw = lapply(count_data, function(x){
  mat = matrix(0, nrow = nrow(x[[1]]), ncol = length(x))
  rownames(mat) = rownames(x[[1]])
  colnames(mat) = names(x)
  for(n in names(x)){
    print(n)
    dat = x[[n]]
    mat[rownames(dat), n] = dat[,1]
  }
  colnames(mat)  = colnames(mat) %>% sub(pattern = "_pooled", replacement = "")
  return(mat)
})



count_norm = lapply(count_data, function(x){
  mat = matrix(0, nrow = nrow(x[[1]]), ncol = length(x))
  rownames(mat) = rownames(x[[1]])
  colnames(mat) = names(x)
  for(n in names(x)){
    print(n)
    dat = x[[n]]
    mat[rownames(dat), n] = dat[,1]
  }
  colnames(mat)  = colnames(mat) %>% sub(pattern = "_pooled", replacement = "")
  norm = apply(mat, 2, function(col){
    col = col / sum(col) * 10^6
  })
  return(norm)
})

fe = lapply(count_norm, function(x){
  cells = strsplit(colnames(x), "_") %>% sapply(function(x){x[1]})
  drugs = strsplit(colnames(x), "_") %>% sapply(function(x){x[2]})
  marks = strsplit(colnames(x), "_") %>% sapply(function(x){x[3]})
  is_mark = marks != "input"
  mat = matrix(0, nrow = nrow(x), ncol = sum(is_mark))
  rownames(mat) = rownames(x)
  colnames(mat) = colnames(x)[is_mark]
  pseudo = .001
  for(n in colnames(mat)){
    treat = x[,n] + pseudo
    mark = strsplit(n, "_")[[1]][3]
    control_n = sub(mark, "input", n)
    control = x[,control_n] + pseudo
    mat[,n] = treat / control
  }
  return(mat)
})

source("H:/R_workspace/jrb_R_scripts/PCA.R")

data = fe$k4me3_promoters
data = data[,grepl("K27M", colnames(data)) & grepl("MCF10A", colnames(data))]
cells = sapply(strsplit(colnames(data), "_"), function(x)x[1])
drugs = sapply(strsplit(colnames(data), "_"), function(x)x[2])
marks = sapply(strsplit(colnames(data), "_"), function(x)x[3])
plotPCA(data, conditionLabels = drugs, replicateLabels = drugs, secondaryConditionLabels = drugs, n = 4, displayReps = T)
