if(!exists("rGREAT_ready")) rGREAT_ready = F
if(!rGREAT_ready){
  library(magrittr)
  library(openxlsx)
  library(xlsx)
  source("setup.R")
  source("module_STEIN_ChIPseq.R")
  source("module_enhancers_Carroll.R")
  source("functions_rGREAT.R")
  source("functions_replot_GREAT_matrix.R")
  setwd(input_dir)
  k4me3_promoters = read.table("STEIN_k4me3_promoters_4ngsplot.bed") %>% bed2granges()
  k4ac_promoters = read.table("STEIN_k4ac_promoters_4ngsplot.bed") %>% bed2granges()
  setwd(wd)
  rGREAT_ready = T
}

#for all marks and cells, compare base_drug to treat_drug and inverse.  
#perform a GREAT analysis for each with unique compared to union of peaks
#peaks may be limited to promoters, enhancers, and/or esr1
#promoters and enhancers are mutually exclusive, so don't do that
#all output will appear in parent_dir with sub directories by cell then mark
# parent_dir = "output/GREAT_comparisons_test"; marks = c("H3K4AC", "H3K4ME3"); cells = c("MCF10A", "MCF7")
# base_drug = "ctrl"; treat_drug = "e2"; filter_to_promoters = T; filter_to_enhancers = F; filter_to_esr1 = T
run_rGREAT = function(parent_dir, marks, cells, base_drug, treat_drug, 
                      filter_to_promoters = F, filter_to_enhancers = F, filter_to_esr1 = F){
  start_dir = getwd()
  dir.create(parent_dir, recursive = T, showWarnings = F)
  parent_dir = normalizePath(parent_dir)
  setwd(parent_dir)
  for(m in marks){
    for(cl in cells){
      
      STEIN_SUBSET = get_res(res = STEIN_ChIPseq, cell = cl, mark = m)
      
      #filter peak sets down to smallest set
      MIN_PEAKS = min(sapply(STEIN_SUBSET, length))
      if(MIN_PEAKS > 10^5) MIN_PEAKS = 10^5
      STEIN_SUBSET = lapply(STEIN_SUBSET, function(x){
        o = order(x$adj_p, decreasing = T)
        pos = 1:length(o)
        names(pos) = o
        pos = pos[as.character(1:length(pos))]
        keep = pos <= MIN_PEAKS
        x[keep]
      })
      #create directory for output
      sub_dir = paste0(cl, "/", m)
      print(paste("working in:", sub_dir))
      dir.create(sub_dir, recursive = T, showWarnings = F)
      setwd(sub_dir) #return to parent dir at end of each loop
      
      try(detach("package:openxlsx", unload = T), silent = T)
      library(xlsx)
      #compare a to b, then b to a      
      base_comparison_results = compare_drugs(grs = STEIN_SUBSET, base_drug, treat_drug, with = F, 
                                              enh_only = filter_to_enhancers, prom_only = filter_to_promoters, esr1_only = filter_to_esr1)
      treat_comparison_results = compare_drugs(grs = STEIN_SUBSET, treat_drug, base_drug, with = F, 
                                               enh_only = filter_to_enhancers, prom_only = filter_to_promoters, esr1_only = filter_to_esr1)
      
      #save gene region associations to xlsx
      base_gene_regions = base_comparison_results$job %>% 
        plotRegionGeneAssociationGraphs(type = NULL) %>% 
        as.data.frame()
      treat_gene_regions = treat_comparison_results$job %>% 
        plotRegionGeneAssociationGraphs(type = NULL) %>% 
        as.data.frame()
      
      wb = createWorkbook()
      wb_name = paste("gene_regions", base_drug, "and", treat_drug, "combined.xlsx", sep = "_")
      if(filter_to_promoters) wb_name = paste0("prom-", wb_name)
      if(filter_to_enhancers) wb_name = paste0("end-", wb_name)
      if(filter_to_esr1) wb_name = paste0("ESR1-", wb_name)
      #these would ideally be lifted back over to hg38, ugh!
      base_sheet = createSheet(wb, paste(base_drug, "hg19"))
      addDataFrame(base_gene_regions, base_sheet, row.names = F)
      treat_sheet = createSheet(wb, paste(treat_drug, "hg19"))
      addDataFrame(treat_gene_regions, treat_sheet, row.names = F)
      
      saveWorkbook(wb, file = wb_name)
      setwd(parent_dir)
    }
  }
  setwd(start_dir)
}

get_res = function(res = bdgdiff_res, cell = "", mark = "", from = "", to = "", direction = ""){
  if(from != "") from = paste0("_", from, "_")
  if(to != "") to = paste0("_", to, "_")
  keep = grepl(cell, names(res)) &
    grepl(mark, names(res)) &
    grepl(from, names(res)) &
    grepl(to, names(res)) &
    grepl(direction, names(res))
  return(res[keep])
}

compare_drugs = function(grs, subject, query, with = T, esr1_only = F, enh_only = F, prom_only = F){
  tmp = get_res(grs, from = subject)
  if(length(tmp) != 1) stop("bad subject match")
  subject_gr = tmp[[1]]
  tmp = get_res(grs, from = query)
  if(length(tmp) != 1) stop("bad query match")
  query_gr = tmp[[1]]
  is_cov = GRanges(coverage(GRangesList(subject_gr, query_gr)) > 0)
  subject_or_query = is_cov[is_cov$score]
  subject_name = paste(subject, query, sep = "-")
  if(prom_only){
    prom = reduce(c(k4ac_promoters, k4me3_promoters))
    subject_name = paste("prom", subject_name, sep = "-")
    olaps = findOverlaps(query = query_gr, subject = prom)
    query_gr = query_gr[queryHits(olaps)]
    olaps = findOverlaps(query = subject_or_query, subject = prom)
    subject_or_query = subject_or_query[queryHits(olaps)]
  }
  if(enh_only){
    CARROLL_ENH = GSE40129$MCF7_any_enhancer[GSE40129$MCF7_enhancer_sets$k4me3_free]
    subject_name = paste("enh", subject_name, sep = "-")
    olaps = findOverlaps(query = query_gr, subject = CARROLL_ENH)
    query_gr = query_gr[queryHits(olaps)]
    olaps = findOverlaps(query = subject_or_query, subject = CARROLL_ENH)
    subject_or_query = subject_or_query[queryHits(olaps)]
  }
  if(esr1_only){
    CARROLL_ESR1 = GSE40129$MCF7_e2_ESR1
    subject_name = paste("ESR1", subject_name, sep = "-")
    olaps = findOverlaps(query = query_gr, subject = CARROLL_ESR1)
    query_gr = query_gr[queryHits(olaps)]
    olaps = findOverlaps(query = subject_or_query, subject = CARROLL_ESR1)
    subject_or_query = subject_or_query[queryHits(olaps)]
  }
  
  
  out = wrap_great(subject38 = subject_or_query, query38 = query_gr, 
                   subject_name = subject_name, query_name = query, with = with)
  return(out)
}