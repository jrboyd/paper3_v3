
env_bdgdiff = new.env()
# source("H:/projects/Apps/diff_marks_v2/scripts/enrichment_testing.R")
evalq(envir = env_bdgdiff, {
  setwd(input_dir)
  files = dir(path = "bdgdiff_short_len/", pattern = "common", full.names = T)
  all_bdgdiff_tables = list()
  
  wrap.read.table = function(file){
    test = read.table(file, sep = "\n", nrows = 2)
    if(nrow(test) == 1){
      warning(paste("no results for:", file))
      blank_gr = GRanges(lr = numeric())
      return(blank_gr)
    } 
    
    tab = read.table(file, skip = 1)
    gr = GRanges(seqnames = tab[,1], ranges = IRanges(start = tab[,2], end = tab[,3]), lr = tab[,5])
    return(gr)
  }
  
  for(f in files){
    name = basename(f)
    name = sub("diff_", "", name)
    name = sub("_c3.0_common.bed", "", name)
    print(name)
    cl = strsplit(name, "_")[[1]][1]
    d1 = strsplit(name, "_")[[1]][2]
    d2 = strsplit(name, "_")[[1]][6]
    m = strsplit(name, "_")[[1]][3]
    
    top_name = paste(cl, m, sep = "_")
    bot_name = paste("from", d1, "to", d2, sep = "_")
    rbot_name = paste("from", d2, "to", d1, sep = "_")
    
    common_res = wrap.read.table(f)
    d1_res = wrap.read.table(sub("common.bed", "cond1.bed", f))
    d2_res = wrap.read.table(sub("common.bed", "cond2.bed", f))
    
    if(is.null(all_bdgdiff_tables[[top_name]])) all_bdgdiff_tables[[top_name]] = list()
    all_bdgdiff_tables[[top_name]][[bot_name]] = list("common" = common_res, "down" = d1_res, "up" = d2_res)
    all_bdgdiff_tables[[top_name]][[rbot_name]] = list("common" = common_res, "down" = d2_res, "up" = d1_res)
  }
  
  #sort names for neatness
  all_bdgdiff_tables = lapply(all_bdgdiff_tables, function(x){
    o = order(names(x))
    return(x[o])
  })
  
  #filter out regions that aren't wide enough
  all_bdgdiff_tables = lapply(all_bdgdiff_tables, function(group){
    lapply(group, function(pair){
      lapply(pair, function(direction){
        keep = width(direction) > 250
        return(direction[keep])
      })
    })
  })
  
  #create bg+up/down peaks by combining up/down with bg
  all_bdgdiff_tables = lapply(all_bdgdiff_tables, function(group){
    lapply(group, function(pair){
      pair$up_bg = c(pair$common, pair$up)
      pair$down_bg = c(pair$common, pair$down)
      return(pair)
    })
  })
  
  setwd(wd)
})
#compare up/down to common+up/down
all_bdgdiff_tables = env_bdgdiff$all_bdgdiff_tables