#combine k27ac and k4me1 peaks

#list files
k27ac_files = dir("narrowPeaks_with_merged_inputs/", pattern = "MCF7.+H3K27AC_pooled_peaks.narrow", full.names = T)
k4me1_files = dir("narrowPeaks_e2_response/", pattern = "H3K4ME1_pooled", full.names = T)

#load as narrowPeak
names(k27ac_files) = sub("_pooled_peaks.narrowPeak", "", basename(k27ac_files))
read_np = function(x){
  np = read.table(x)
  return(np)
}
k27ac_nps = lapply(k27ac_files, read_np)
names(k4me1_files) = sub("_pooled_peaks.narrowPeak", "", basename(k4me1_files))
k4me1_nps = lapply(k4me1_files, read_np)

#convert narrowPeak to GRanges
np2gr = function(x){
  keep = x$V7 > 4
  np2granges(x[keep,])
}
k27ac_gr = lapply(k27ac_nps, np2gr)
k4me1_gr = lapply(k4me1_nps, np2gr)

venn_peaks = function(gr_pairList){
  if(length(gr_pairList) != 2) stop("length must be 2")
  a = gr_pairList[[1]]
  b = gr_pairList[[2]]
  a_name = names(gr_pairList)[1]
  b_name = names(gr_pairList)[2]
  olaps = findOverlaps(a, b)
  a_wb = a[queryHits(olaps)]
  a_nb = a[-queryHits(olaps)]
  b_wa = b[subjectHits(olaps)]
  b_na = b[-subjectHits(olaps)]
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

write_beds = function(gr_pairList){
  gr_venn = venn_peaks(gr_pairList)
  gr_venn_randbed = lapply(gr_venn, function(x){
    o = order(runif(length(x)))
    x = x[o][1:10000]
    as_bed = as.data.frame(x)[,1:3]
    as_bed = cbind(as_bed, paste("rand", 1:nrow(as_bed), sep = "_"), rep(0, nrow(as_bed)), rep(".", nrow(as_bed)))
    colnames(as_bed) = NULL
    o = order(as_bed[,2])
    as_bed = as_bed[o,]
    o = order(as.numeric(sub("chr", "", as_bed[,1])))
    as_bed = as_bed[o,]
    return(as_bed)
  })
  
  hidden = sapply(names(gr_venn_randbed), function(bed_nam){
    bed = gr_venn_randbed[[bed_nam]]
    bed = t(apply(bed, 1, function(x){
      avg = round(mean(as.numeric(x[2:3])))
      x[2] = avg - 2500
      x[3] = avg + 2500
      return(x)
    }))
    write.table(bed, file = paste0(bed_nam, ".bed"), quote = F, col.names = F, row.names = F, sep = '\t')
  })
}

write_beds(list(MCF7_ctrl_ESR1 = ESR1_ctrl_gr_all, MCF7_e2_ESR1 = ESR1_e2_gr_all))
write_beds(k4me1_gr)
write_beds(k27ac_gr[c(2,3)])

a = out$`MCF7_ctrl_H3K4ME1 with MCF7_e2_H3K4ME1`
b = out$`MCF7_e2_H3K4ME1 with MCF7_ctrl_H3K4ME1`
plot(log10(a$enrichment), log10(a$adj_p))
plot(log10(b$enrichment), log10(b$adj_p))
plot(log10(a$enrichment), log10(b$enrichment))
plot(log10(a$adj_p), log10(b$adj_p))


tmp = lapply(k27ac_gr, function(k27){
  lapply(k4me1_gr, function(k4){
    olaps = findOverlaps(query = k27, subject = k4)
    c(length(k27), length(unique(queryHits(olaps))), length(k4), length(unique(subjectHits(olaps))))
  })
})

x = k27ac_nps$MCF7_bza_H3K27AC
keep = x$V7 > 4
hist((x$V7[keep]))
