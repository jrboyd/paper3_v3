library(magrittr)

env_bdgdiff = new.env()
evalq(envir = env_bdgdiff, {
  
  setwd(input_dir)
  data_dir = "bdgdiff_nomodel"
  bdgdiff_files = dir(data_dir, full.names = T, pattern = "cond")
  bdgdiff_res = list()
  default_df = data.frame(chrm = factor(), start = integer(), end = integer(), id = character(), llr = numeric())
  for(f in bdgdiff_files){
    print(f)
    #check for no results
    if(nrow(read.table(f,  sep = "\n", nrows = 2)) < 2){
      df = default_df
    }else{
      df = read.table(f,  skip = 1)
      colnames(df) = colnames(default_df)
    }
    attrib = (f %>% basename %>% strsplit(split = "_"))[[1]]
    cell = attrib[2]
    mark = attrib[4]
    drug1 = attrib[3]
    drug2 = attrib[7]
    up = attrib[10]
    desc = paste("from", drug1, "to", drug2, sep = "_")
    if(up == "cond1.bed"){#drug 1 is up
      direction = "down"
    }else{
      direction = "up"
    }
    full_name = paste(cell, mark, desc, direction, sep = "_")
    bdgdiff_res[[full_name]] = df
    
    desc = paste("from", drug2, "to", drug1, sep = "_")#compliment
    if(up == "cond1.bed"){#drug 1 is up
      direction = "up"
    }else{
      direction = "down"
    }
    full_name = paste(cell, mark, desc, direction, sep = "_")
    bdgdiff_res[[full_name]] = df
  }
  
  get_bdgdiff_res = function(res = bdgdiff_res, cell = "", mark = "", from = "", to = "", direction = ""){
    if(from != "") from = paste0("from_", from, "_")
    if(to != "") to = paste0("to_", to, "_")
    keep = grepl(cell, names(res)) &
      grepl(mark, names(res)) &
      grepl(from, names(res)) &
      grepl(to, names(res)) &
      grepl(direction, names(res))
    return(res[keep])
  }
  
  bdgdiff_res.strict = lapply(get_bdgdiff_res(from = "ctrl"), function(x){
    x[x$llr >9,]
  })
  
  bdgdiff_res.strict.gr = lapply(bdgdiff_res.strict, function(x){
    gr = GRanges(seqnames = x$chrm, ranges = IRanges(start = x$start, end = x$end), id = x$id, llr = x$llr)
    return(gr)
  })
  
  MCF7_enhancers_intergenic = GSE40129$MCF7_any_enhancer[GSE40129$MCF7_enhancer_sets$k4me3_free]
  
  bdgdiff_res.strict.enhancers = lapply(bdgdiff_res.strict.gr, function(x){
    olaps = suppressWarnings(findOverlaps(query = MCF7_enhancers_intergenic, subject = x))
    enh_ids = names(MCF7_enhancers_intergenic)[queryHits(olaps)]
    return(enh_ids)
  })
  
  # get_bdgdiff_res(bdgdiff_res.strict.enhancers, cell = "MCF7", mark = "H3K27AC")
  
})

env_bdgdiff$bdg