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
  files = dir("STEIN_ChIPseq/", full.names = T)
  attribs = t(sapply(strsplit(files, "[/_]"), function(x)x[3:5]))
  cells = unique(attribs[,1])
  drugs = unique(attribs[,2])
  marks = unique(attribs[,3])
  
  for(cl in cells){
    # STEIN_ChIPseq[[cl]] = list()
    for(d in drugs){
      # STEIN_ChIPseq[[cl]][[d]] = list()
      for(m in marks){
        keep = grepl(cl, files) & grepl(paste0("_", d, "_"), files) & grepl(m, files)
        if(sum(keep) != 1) stop("bad!")
        f = files[keep]
        # print(f)
        STEIN_ChIPseq[[paste(cl, d, m, sep = "_")]] = file2gr(f)
      }
    }
  }
  # 
  # 
  # library(ChIPpeakAnno)
  # load("ref/ensg_ref.save")
  # load("ref/enst_ref.save")
  # # enst_ref = enst_ref[enst_ref$tag == "basic"]
  # ref2gr = function(ref){
  #   gr = GRanges(seqnames = ref$chrm, ranges = IRanges(start = ref$start, end = ref$end), strand = ref$strand, 
  #           transcript_id = rownames(ref), gene_id = ref$gene_id)
  #   names(gr) = ref$gene_id
  #   return(gr)
  # }
  # ensg_gr = ref2gr(ensg_ref)
  # 
  # STEIN_ChIPseq_annot = lapply(STEIN_ChIPseq, function(x){
  #   annotatePeakInBatch(x, AnnotationData = ensg_gr)
  # }) 
  # 
  # STEIN_ChIPseq_pie = lapply(names(STEIN_ChIPseq_annot), function(name){
  #   x = STEIN_ChIPseq_annot[[name]]
  #   pie1(table(x$insideFeature), main = name)
  # })
  # 
  # STEIN_ChIPseq_upstream = lapply(STEIN_ChIPseq_annot, function(x){
  #   x = x[x$insideFeature == "upstream",]
  #   out = sapply(strsplit(names(x), "\\."), function(str){
  #     return(str[2])
  #   })
  #   out = unique(out)
  #   return(out)
  # })
  # 
  # STEIN_ChIPseq_downstream = lapply(STEIN_ChIPseq_annot, function(x){
  #   x = x[x$insideFeature == "downstream",]
  #   out = sapply(strsplit(names(x), "\\."), function(str){
  #     return(str[2])
  #   })
  #   out = unique(out)
  #   return(out)
  # })
  # 
  # sapply(STEIN_ChIPseq_upstream, function(x){
  #   sapply(STEIN_ChIPseq_upstream, function(y){
  #     return(length(intersect(x, y)) / length(x))
  #   })
  # })
  # 
  # sapply(STEIN_ChIPseq_downstream, function(x){
  #   sapply(STEIN_ChIPseq_downstream, function(y){
  #     return(length(intersect(x, y)) / length(x))
  #   })
  # })
  # 
  # library(rtracklayer)
  # 
  # plot_peaks_intersects = function(x, y, xlab = "x", ylab = "y"){
  #   y_all = np2granges(y)
  #   MIN = 50
  #   MAX = 10000
  #   xs = (MAX - MIN) * 0:20/20 + MIN
  #   xs = xs[xs < nrow(x)]
  #   
  #   hits = sapply(xs, function(top){
  #     gr = np2granges(x, top)
  #     olaps = findOverlaps(y_all, gr)
  #     subHits = unique(subjectHits(olaps))
  #     queHits = unique(queryHits(olaps))
  #     return(c(length(subHits), length(queHits)))
  #   })
  #   plot(xs, hits[1,], xlab = xlab, ylab = ylab)
  # }
  # 
  # #GSE32222 pooled from SRR1021787 SRR1021788 SRR1021789
  # esr1_binding = read.table("narrowPeaks_diff_ESR1//MCF7_ESR1_pooled_peaks.narrowPeak")
  # #plot 
  # plot_peaks_intersects(ESR1_ctrl_np, esr1_binding)
  # plot_peaks_intersects(ESR1_e2_np, esr1_binding)
  
  # keep = esr1_binding$V9 > 6
  # esr1_binding = esr1_binding[keep,]
  # o = order(esr1_binding$V9, decreasing = T)
  # esr1_binding = esr1_binding[o,]
  # o = order(as.numeric(rownames(esr1_binding)))
  # esr1_binding = esr1_binding[o,]
  # ESR1_diffpaper_gr_all = GRanges(seqnames = esr1_binding[,1], ranges = IRanges(start = esr1_binding[,2], end = esr1_binding[,3]))
  # ESR1_binding = list(diffpaper = ESR1_diffpaper_gr_all,
  #                     e2_induced = ESR1_e2_gr_all,
  #                     non_induced = ESR1_ctrl_gr_all)
  # save(ESR1_binding, file = "module_ESR1_binding.save")
})

STEIN_ChIPseq = env_STEIN_ChIPseq$STEIN_ChIPseq
# ESR1_binding = env_STEIN_ChIPseq$ESR1_binding
