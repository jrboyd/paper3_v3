# source("H:/R_workspace/jrb_R_scripts/parse_gtf.R")
# ensg_ref = parse_gtf("C:/Users/jrboyd/Downloads/gencode.v21.annotation_genes.gtf", rownames_attrib = "gene_id", feature_type = "gene", additional_attrib = "gene_type")
#filter by DE genes
#establish types of ESR1 binding (none, +ctrl, +e2-ctrl)
#promoters bound by ESR1 types
#establish enhancers
#enhancers bound by ESR1 types

#discover gene promoters that are bound by ESR1
# 
# library(GenomicRanges)
# source("setup.R")

setwd(input_dir)
load("DE_res.2FC.0.05adjP.save")
load("ref/ensg_ref.save")
load("module_ESR1_binding.save")
is_protein = ensg_ref$gene_type == "protein_coding"

ref2prom_gr = function(ref){
  tsses = apply(ref, 1, function(y){
    if(y[6] == "+"){
      out = y[4]
    }else{
      out = y[5]
    }
    return(as.numeric(out))
  })
  ref$start = tsses - 1
  ref$end = tsses + 1
  gr = GRanges(seqnames = ref$chrm, 
               ranges = IRanges(start = ref$start, end = ref$end), 
               strand = ref$strand,
               gene_id = ref$gene_id, 
               gene_name = ref$gene_name)
  return(gr)
}
gr2dists = function(from_gr, to_gr){
  nearest = distanceToNearest(x = from_gr, subject = to_gr)
  dists = as.data.frame(nearest)[,3]
  return(dists)
}

olaps = findOverlaps(ESR1_ctrl_gr_all, ESR1_e2_gr_all)
ESR1_ctrl_gr_wE2 = ESR1_ctrl_gr_all[queryHits(olaps)]
ESR1_ctrl_gr_nE2 = ESR1_ctrl_gr_all[-queryHits(olaps)]

olaps = findOverlaps(ESR1_e2_gr_all, ESR1_ctrl_gr_all)
ESR1_e2_gr_wCtrl = ESR1_e2_gr_all[queryHits(olaps)]
ESR1_e2_gr_nCtrl = ESR1_e2_gr_all[-queryHits(olaps)]

ESR1_sets = list(
  ESR1_ctrl_gr_all = ESR1_ctrl_gr_all, 
  ESR1_e2_gr_all = ESR1_e2_gr_all, 
  ESR1_ctrl_gr_wE2 = ESR1_ctrl_gr_wE2,
  ESR1_ctrl_gr_nE2 = ESR1_ctrl_gr_nE2,
  ESR1_e2_gr_wCtrl = ESR1_e2_gr_wCtrl,
  ESR1_e2_gr_nCtrl = ESR1_e2_gr_nCtrl,
  ESR1_diffpaper_gr_all = ESR1_diffpaper_gr_all)

#DE genes between e2 and ctrl
DEG = union(DE_res$`MCF7_drugs_RNA from e2 to ctrl down`, DE_res$`MCF7_drugs_RNA from e2 to ctrl up`)

#distance from DEG promoter to ESR1 peak
for(ER_name in names(ESR1_sets)){
  ER = ESR1_sets[[ER_name]]
  ref = ensg_ref[is_protein, ]
  ref = ref[intersect(rownames(ref), DEG), ]
  all_promoters_GR = ref2prom_gr(ref)
  all_promoters_dists = gr2dists(from_gr = all_promoters_GR, to_gr = ER)
  keep = all_promoters_dists < 5000
  all_promoters_dists = all_promoters_dists[keep]
  hist(all_promoters_dists, main = paste("distance of", ER_name, "to DE genes"))
}
for(ER_name in names(ESR1_sets)){
  ER = ESR1_sets[[ER_name]]
  all_promoters_GR = ref2prom_gr(ensg_ref[is_protein, ])
  all_promoters_dists = gr2dists(from_gr = all_promoters_GR, to_gr = ER)
  keep = all_promoters_dists < 5000
  all_promoters_dists = all_promoters_dists[keep]
  hist(all_promoters_dists, main = paste("distance of", ER_name, "to genes"))
}
