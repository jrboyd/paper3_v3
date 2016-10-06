#filter by DE genes
#establish types of ESR1 binding (none, +ctrl, +e2-ctrl)
#promoters bound by ESR1 types
#establish enhancers
#enhancers bound by ESR1 types

#discover gene promoters that are bound by ESR1

library(GenomicRanges)
source("setup.R")

setwd(input_dir)
load("DE_res.save")
load("ref/ensg_ref.save")
ESR1_ctrl_file = dir("narrowPeaks_e2_response/", pattern = "control-ESR1_pool", full.names = T)
ESR1_e2_file = dir("narrowPeaks_e2_response/", pattern = "e2-ESR1_pool", full.names = T)
ESR1_ctrl_np = read.table(ESR1_ctrl_file)
ESR1_e2_np = read.table(ESR1_e2_file)

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
  GRanges(seqnames = np[,1], ranges = IRanges(start = np[,2], end = np[,3]))
}

ESR1_ctrl_gr_all = np2granges(ESR1_ctrl_np)
ESR1_e2_gr_all = np2granges(ESR1_e2_np)

x = ESR1_ctrl_np
y = ESR1_e2_np
plot_peaks_intersects = function(x, y){
  y_all = np2granges(y)
  MIN = 50
  MAX = 10000
  xs = (MAX - MIN) * 0:20/20 + MIN
  xs = xs[xs < nrow(x)]
  
  hits = sapply(xs, function(top){
    gr = np2granges(x, top)
    olaps = findOverlaps(y_all, gr)
    subHits = unique(subjectHits(olaps))
    queHits = unique(queryHits(olaps))
    return(c(length(subHits), length(queHits)))
  })
  plot(xs, hits[1,])
}



olaps = findOverlaps(ESR1_ctrl_gr, ESR1_e2_gr)
subHits = unique(subjectHits(olaps))
queHits = unique(queryHits(olaps))
length(subHits)
length(queHits)

esr1_binding = read.table("narrowPeaks_diff_ESR1//MCF7_ESR1_pooled_peaks.narrowPeak")
plot_peaks_intersects(ESR1_ctrl_np, esr1_binding)
plot_peaks_intersects(ESR1_e2_np, esr1_binding)

# keep = esr1_binding$V9 > 6
# esr1_binding = esr1_binding[keep,]
# o = order(esr1_binding$V9, decreasing = T)
# esr1_binding = esr1_binding[o,]
# o = order(as.numeric(rownames(esr1_binding)))
# esr1_binding = esr1_binding[o,]
esr1_binding = GRanges(seqnames = esr1_binding[,1], ranges = IRanges(start = esr1_binding[,2], end = esr1_binding[,3]))


load("DM_lists.ensg.padj0.01.fc2.save")
setwd(wd)
setwd(output_dir)

keep = ensg_ref$gene_type == "protein_coding"
ensg_ref = ensg_ref[keep,]



ref2prom_gr = function(ref){
  tsses = apply(ref, 1, function(y){
    if(y[6] == "+"){
      out = y[4]
    }else{
      out = y[5]
    }
    return(as.numeric(out))
  })
  ref$start = tsses - 2000
  ref$end = tsses + 2000
  gr = GRanges(seqnames = ref$chrm, 
               ranges = IRanges(start = ref$start, end = ref$end), 
               strand = ref$strand,
               gene_id = ref$gene_id, 
               gene_name = ref$gene_name)
  return(gr)
}
gr2dists = function(gr){
  nearest = distanceToNearest(x = gr, subject = esr1_binding)
  dists = as.data.frame(nearest)[,3]
  return(dists)
}


all_promoters_GR = ref2prom_gr(ensg_ref)
all_promoters_dists = gr2dists(all_promoters_GR)
keep = all_promoters_dists < 500
keep = ifelse(is.na(keep), F, keep)
all_ESR1_bound_promoters_GR = all_promoters_GR[keep]
ESR1_symbols = all_ESR1_bound_promoters_GR$gene_name

ESR1_filter = all_ESR1_bound_promoters_GR$gene_id
DE_esr1_res = lapply(DE_res, function(x){
  intersect(x, ESR1_filter)
})
DE_esr1_res_sym = lapply(DE_esr1_res, function(x){
  ensg_ref[x,]$gene_name
  
})
hist(log10(all_promoters_dists), main = "log10 distance of ESR1 binding to promoters")
cbind(sapply(DE_res, length), sapply(DE_esr1_res, length))

keep = !grepl("!", names(all_diffbind_tables), fixed = T)
all_diffbind_tables = all_diffbind_tables[keep]
keep = sapply(all_diffbind_tables, nrow) > 0
all_diffbind_tables = all_diffbind_tables[keep]

all_diffbind_gr_up = lapply(all_diffbind_tables, function(x){
  is_up = x[,8] > x[,7]
  up = x[is_up,]
  down = x[!is_up,]
  print(nrow(up))
  down_gr = GRanges(seqnames = up$seqnames, ranges = IRanges(start = up$start, end = up$end))
  return(down_gr)
})

all_diffbind_gr_down = lapply(all_diffbind_tables, function(x){
  is_up = x[,8] > x[,7]
  up = x[is_up,]
  down = x[!is_up,, drop = F]
  print(nrow(down))
  down_gr = GRanges(seqnames = down$seqnames, ranges = IRanges(start = down$start, end = down$end))
  return(down_gr)
})

all_esr1_diffbind_gr_up = t(sapply(all_diffbind_gr_up, function(x){
  
  olaps = findOverlaps(query = x, subject = esr1_binding)
  return(c(length(olaps) , (length(x))))
  # return(x[queryHits(olaps)])
}))

all_esr1_diffbind_gr_down = t(sapply(all_diffbind_gr_down, function(x){
  
  olaps = findOverlaps(query = x, subject = esr1_binding)
  return(c(length(olaps) , (length(x))))
  # return(x[queryHits(olaps)])
}))

xs = all_esr1_diffbind_gr_up[,2]
ys = all_esr1_diffbind_gr_up[,1] / (ifelse(all_esr1_diffbind_gr_up[,2] < 1, 1, all_esr1_diffbind_gr_up[,2]))
is_MCF7 = grepl("MCF7", names(xs))
colors = ifelse(is_MCF7, "red", "blue")
plot(xs, ys, col = colors, xlab = "number of DM peaks", ylab = "fraction bound by ESR1")
legend("topright", fill= c("red", "blue"), legend = c("MCF7", "MCF10A"))

xs = all_esr1_diffbind_gr_down[,2]
ys = all_esr1_diffbind_gr_down[,1] / (ifelse(all_esr1_diffbind_gr_down[,2] < 1, 1, all_esr1_diffbind_gr_down[,2]))
is_MCF7 = grepl("MCF7", names(xs))
colors = ifelse(is_MCF7, "red", "blue")
plot(xs, ys, col = colors, xlab = "number of DM peaks", ylab = "fraction bound by ESR1")
legend("topright", fill= c("red", "blue"), legend = c("MCF7", "MCF10A"))