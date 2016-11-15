#filter by DE genes
#establish types of ESR1 binding (none, +ctrl, +e2-ctrl)
#promoters bound by ESR1 types
#establish enhancers
#enhancers bound by ESR1 types

#discover gene promoters that are bound by ESR1

# library(GenomicRanges)
# source("setup.R")
env_STEIN_ChIPseq = new.env()
evalq(envir = env_STEIN_ChIPseq, {
  
  file2gr = function(np_file){
    print(paste("load", np_file))
    #load narrowPeak
    np = read.table(np_file)
    #convert narrowPeak to GRanges
    gr = np2granges(np)
    return(gr)
  }
  
  setwd(input_dir)
  #from GSE14664 #these are not good data
  #switch to STEIN_ChIPseq
  STEIN_ChIPseq = list()
  files = dir("STEIN_ChIPseq_nomodel/", full.names = T)
  attribs = t(sapply(strsplit(basename(files), "[/_]"), function(x)x[1:3]))
  cells = unique(attribs[,1])
  drugs = unique(attribs[,2])
  marks = unique(attribs[,3])
  
  for(cl in cells){
    for(d in drugs){
      for(m in marks){
        keep = grepl(cl, files) & grepl(paste0("_", d, "_"), files) & grepl(m, files)
        if(sum(keep) != 1) stop("bad!")
        f = files[keep]
        STEIN_ChIPseq[[paste(cl, d, m, sep = "_")]] = file2gr(f)
      }
    }
  }
  setwd(wd)
})

STEIN_ChIPseq = env_STEIN_ChIPseq$STEIN_ChIPseq
# ESR1_binding = env_STEIN_ChIPseq$ESR1_binding
