#filter by DE genes
#establish types of ESR1 binding (none, +ctrl, +e2-ctrl)
#promoters bound by ESR1 types
#establish enhancers
#enhancers bound by ESR1 types

#discover gene promoters that are bound by ESR1

#enh588 from Korkmaz 2016 chr11:69,513,645-69,517,457, hits Carroll enhancer_5686

# library(GenomicRanges)
# source("setup.R")
env_ESR1_binding = new.env()
evalq(envir = env_ESR1_binding, {
  
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
  #switch to GSE40129
  GSE40129 = list()
  GSE40129$MCF7_e2_ESR1 = file2gr("GSE40129_MCF7_ESR1_enhancers/MCF7-E2_ESR1_pooled_peaks.narrowPeak")
  GSE40129$MCF7_e2_ESR1_idr = file2gr("GSE40129_MCF7_ESR1_enhancers/MCF7-E2_ESR1_pooled_peaks_passIDR.05.narrowPeak")
  GSE40129$MCF7_e2_H3K4ME1 = file2gr("GSE40129_MCF7_ESR1_enhancers/MCF7-E2_H3K4ME1_R1_peaks.narrowPeak")
  GSE40129$MCF7_ctrl_H3K4ME1 = file2gr("GSE40129_MCF7_ESR1_enhancers/MCF7-Veh_H3K4ME1_R1_peaks.narrowPeak")
  GSE40129$MCF7_e2_H3K27AC = file2gr("GSE40129_MCF7_ESR1_enhancers/MCF7-E2_H3K27AC_R1_peaks.narrowPeak")
  GSE40129$MCF7_ctrl_H3K27AC = file2gr("GSE40129_MCF7_ESR1_enhancers/MCF7-Veh_H3K27AC_R1_peaks.narrowPeak")
  
  #load k4me3 promoters
  promoters_bed = read.table("STEIN_ChIPseq_at_features/k4me3_promoters.bed")
  promoters_gr = GRanges(seqnames = promoters_bed[,1], 
                         ranges = IRanges(start = promoters_bed[,2], end = promoters_bed[,3]),
                         strand = promoters_bed[,6])
  
  #identify overlaps between k4me1 and k27ac in any combination
  #write bed file of enhancers present in ctrl or e2 conditions
  e2_olaps = venn_peaks(list(MCF7_e2_H3K4ME1 = GSE40129$MCF7_e2_H3K4ME1, MCF7_e2_H3K27AC = GSE40129$MCF7_e2_H3K27AC))
  ctrl_olaps = venn_peaks(list(MCF7_ctrl_H3K4ME1 = GSE40129$MCF7_ctrl_H3K4ME1, MCF7_ctrl_H3K27AC = GSE40129$MCF7_ctrl_H3K27AC))
  e2_enh = reduce(c(e2_olaps$`MCF7_e2_H3K4ME1 with MCF7_e2_H3K27AC` , e2_olaps$`MCF7_e2_H3K27AC with MCF7_e2_H3K4ME1`))
  ctrl_enh = reduce(c(ctrl_olaps$`MCF7_ctrl_H3K4ME1 with MCF7_ctrl_H3K27AC` , ctrl_olaps$`MCF7_ctrl_H3K27AC with MCF7_ctrl_H3K4ME1`))
  
  # #remove k4me3 promoters
  # e2_phits = findOverlaps(e2_enh, promoters_gr) %>% queryHits() %>% unique()
  # ctrl_phits = findOverlaps(ctrl_enh, promoters_gr) %>% queryHits() %>% unique()
  # e2_enh = e2_enh[-e2_phits]
  # ctrl_enh = ctrl_enh[-ctrl_phits]
  
  any_enh = reduce(c(e2_enh, ctrl_enh))
  starts = start(any_enh)
  widths = width(any_enh)
  mids = as.integer(starts + widths / 2)
  ext = 2000
  new_starts = mids - ext
  new_ends = mids + ext
  chrm = as.character(seqnames(any_enh))
  enh_name = paste("enhancer", 1:length(chrm), sep = "_")
  enh_bed = cbind(chrm, new_starts, new_ends, enh_name, rep(0, length(chrm)), rep(".", length(chrm)))
  # write.table(enh_bed, file = "Carroll_any_enhancers.bed", sep = "\t", row.names = F, col.names = F, quote = F)
  
  MCF7_any_enhancer = GRanges(seqnames = chrm, ranges = IRanges(new_starts, new_ends))
  
  #need to reappaly reduction after extension
  any_enh = reduce(MCF7_any_enhancer)
  starts = start(any_enh)
  widths = width(any_enh)
  mids = as.integer(starts + widths / 2)
  new_starts = mids - ext
  new_ends = mids + ext
  chrm = as.character(seqnames(any_enh))
  enh_name = paste("enhancer", 1:length(chrm), sep = "_")
  enh_bed = cbind(chrm, new_starts, new_ends, enh_name, rep(0, length(chrm)), rep(".", length(chrm)))
  
  
  MCF7_any_enhancer = GRanges(seqnames = chrm, ranges = IRanges(new_starts, new_ends))
  names(MCF7_any_enhancer) = enh_name
  
  write.table(enh_bed, file = "Carroll_any_enhancers.bed", sep = "\t", row.names = F, col.names = F, quote = F)
  GSE40129$MCF7_any_enhancer = MCF7_any_enhancer
  
  #divide into e2 unique and ctrl unique and common sets
  e2_hits = queryHits(findOverlaps(query = MCF7_any_enhancer, subject = e2_enh))
  ctrl_hits = queryHits(findOverlaps(query = MCF7_any_enhancer, subject = ctrl_enh))
  common = intersect(e2_hits, ctrl_hits)
  e2_only = setdiff(e2_hits, ctrl_hits)
  ctrl_only = setdiff(ctrl_hits, e2_hits)
  #set that intersects with k4me3 bound promoters
  k4me3 = names(MCF7_any_enhancer)[findOverlaps(query = MCF7_any_enhancer, subject = promoters_gr) %>% queryHits() %>% unique()]
  k4me3_free = setdiff(names(MCF7_any_enhancer), k4me3)
  #set bound by esr1
  esr1 = names(MCF7_any_enhancer)[findOverlaps(query = MCF7_any_enhancer, subject = GSE40129$MCF7_e2_ESR1) %>% queryHits() %>% unique()]
  esr1_free = setdiff(names(MCF7_any_enhancer), esr1)
  MCF7_enhancer_sets = list(common = enh_name[common], 
                            ctrl_only = enh_name[ctrl_only], 
                            e2_only = enh_name[e2_only], 
                            k4me3 = k4me3,
                            k4me3_free = k4me3_free,
                            esr1 = esr1,
                            esr1_free = esr1_free)
  GSE40129$MCF7_enhancer_sets = MCF7_enhancer_sets
  # tmp = findOverlaps(MCF7_any_enhancer[MCF7_enhancer_sets$ctrl_only], MCF7_any_enhancer[MCF7_enhancer_sets$e2_only])
  # MCF7_any_enhancer[MCF7_enhancer_sets$ctrl_only][queryHits(tmp)]
  # MCF7_any_enhancer[MCF7_enhancer_sets$e2_only][subjectHits(tmp)]
  setwd(wd)
})

GSE40129 = env_ESR1_binding$GSE40129
