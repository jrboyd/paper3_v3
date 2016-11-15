#calculates DE_res and saves the result
#uses RNA_MIN_FC and RNA_MIN_ADJP from config file
library(magrittr)
env_DE_res = new.env()
evalq(envir = env_DE_res, {
  files.10A = dir("Z:/RNA_seq_data/Breast/MCF10A_Pfizer/DESeq2/", pattern = "Significant_v21.txt", full.names = T)
  files.7 = dir("Z:/RNA_seq_data/Breast/MCF7_Pfizer/DESeq2/", pattern = "Significant_v21.txt", full.names = T)
  DE_res = list()
  cleanup_filename = function(f){
    sapply(basename(f) %>%
             sub(pattern = "BZAE2", replacement = "e2bza", x = .) %>%
             sub("BZAGC10", "gc10bza", x = .) %>%
             tolower() %>%
             strsplit("_"), function(x){
               paste("from", x[4], "to", x[2], sep = " ")
             })
  }
  filter_DE_res = function(f){
    dat = read.table(f)
    is_up = dat$log2FoldChange < log2(1/RNA_MIN_FC)
    is_down = dat$log2FoldChange > log2(RNA_MIN_FC)
    is_sig = dat$padj < RNA_MIN_ADJP
    res_down = rownames(dat)[is_up & is_sig]
    res_up = rownames(dat)[is_down & is_sig]
    return(list(up = res_up, down = res_down))
  }
  for(f in files.10A){
    name = cleanup_filename(f)
    name = paste("MCF10A_drugs_RNA", name)
    res = filter_DE_res(f)
    DE_res[[paste(name, "down")]] = res$down
    DE_res[[paste(name, "up")]] = res$up
  }
  for(f in files.7){
    name = cleanup_filename(f)
    name = paste("MCF7_drugs_RNA", name)
    res = filter_DE_res(f)
    DE_res[[paste(name, "down")]] = res$down
    DE_res[[paste(name, "up")]] = res$up
  }
  
  setwd(input_dir)
  DE_res = DE_res
  save(DE_res, file = paste0("DE_res.", RNA_MIN_FC, "FC.", RNA_MIN_ADJP, "adjP.save"))
  setwd(wd)
  
})

DE_res = env_DE_res$DE_res
