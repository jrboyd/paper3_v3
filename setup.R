if(!exists("wd")){
  input_dir = "H:/projects/Terri/paper3_data/"
  wd = getwd()
  output_dir = paste0(wd, "/output")
  dir.create(output_dir, showWarnings = F)
}
setwd(wd)
library(GenomicRanges)
np2granges = function(np, top = NA){
  if(is.na(top)){
    top = nrow(np)
  }
  ref = 1:nrow(np)
  o = order(np[,9], decreasing = T)
  np = np[o,][1:top,]
  ref = ref[o][1:top]
  o = order(ref)
  np = np[o,]
  GRanges(seqnames = np[,1], ranges = IRanges(start = np[,2], end = np[,3]), enrichment = np[,7], adj_p = np[,9])
}

file2granges = function(file, min_adjp = 0){
  np = read.table(file)
  gr = np2granges(np)
  k = gr$adj_p > min_adjp
  return(gr[k])
}

venn_peaks = function(gr_pairList){
  if(length(gr_pairList) != 2) stop("length must be 2")
  a = gr_pairList[[1]]
  b = gr_pairList[[2]]
  a_name = names(gr_pairList)[1]
  b_name = names(gr_pairList)[2]
  olaps = findOverlaps(query = a, subject = b)
  qH = unique(queryHits(olaps))
  sH = unique(subjectHits(olaps))
  a_wb = a[qH]
  a_nb = a[-qH]
  b_wa = b[sH]
  b_na = b[-sH]
  out = list(
    a_wb,
    a_nb,
    b_wa,
    b_na
  )
  names(out) = c(
    paste(a_name, "with", b_name),
    paste(a_name, "no", b_name),
    paste(b_name, "with", a_name),
    paste(b_name, "no", a_name)
  )
  return(out)
}