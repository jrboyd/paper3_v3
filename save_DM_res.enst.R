library("DiffBind")
library("xlsx")

source("H:/projects/Apps/diff_marks_v2/scripts/enrichment_testing.R")
diffbind_res = list()
files = dir(path = "input/diffbind_drugs_r2/", pattern = "diffbind.+.save", full.names = T)
for(f in files){
  print(f)
  load(f)
  diffbind_res[[basename(f)]] = ss
}

load("H:/projects/Terri/paper3_v1/ref/ref_granges.save")
load("H:/projects/Terri/paper3_v1/ref/supp_enst.save")
load("H:/projects/Apps/diff_marks_v2/ref/gene_ontolgies.save")
load("H:/projects/Apps/diff_marks_v2/ref/msigdb_byCategory.save")
load("H:/projects/Apps/diff_marks_v2/ref/ensg_dicts.save")
load("H:/projects/Apps/diff_marks_v2/ref/go_table.save")
load("H:/projects/Terri/paper3_v1/ref/enst_ref.save")
load("H:/projects/Terri/paper3_v1/ref/ensg_ref.save")

ensg_prom = enst_ref
k = ensg_prom$strand == "+"
ensg_prom[k,"end"] = ensg_prom[k,"start"] + 1000
ensg_prom[k,"start"] = ensg_prom[k,"start"] - 1000
k = !k
ensg_prom[k,"start"] = ensg_prom[k,"end"] - 1000
ensg_prom[k,"end"] = ensg_prom[k,"end"] + 1000

ensg_prom_granges = GRanges(ensg_prom$chrm, 
                            IRanges(ensg_prom$start, ensg_prom$end), 
                            transcript_id = rownames(ensg_prom),
                            gene_id = ensg_prom$gene_id, 
                            gene_name = ensg_prom$gene_name)
names(ensg_prom_granges) = rownames(ensg_prom)
my_ref = ensg_prom_granges

# dir.create("promoters_stringent", showWarnings = F)
save_filename = "input/DM_res.save"
min_fdr=.05
min_FC = 1

assign_ensg = function(peak_granges, ref_granges){
  olaps = findOverlaps(peak_granges, ref_granges, minoverlap = )
  out = my_ref[subjectHits(olaps)]
  return(out$transcript_id)
}

all_ensg_lists = list()
all_gene_lists = list()
all_diffbind_tables = list()
print("annotate diffbind peaks to genes...")
for(j in  1:length(diffbind_res)){
  exp_nam =  names(diffbind_res)[j]
  # if(grepl("compareReps", exp_nam)) next
  print(exp_nam)
  ss = diffbind_res[[exp_nam]]
  for(i in 1:length(ss$contrasts)){
    rprt = dba.report(ss, contrast = i)
    cont_nam = paste("from", paste(ss$contrasts[[i]][3:4], collapse = " to "))
    
    tab = as.data.frame(rprt)
    keep = tab$FDR < min_fdr
    nres = sum(keep)
    tab = tab[keep,, drop = F]
    if(nres > 0){
      keep = abs(tab$Fold) > min_FC
      nres = sum(keep)
      tab = tab[keep,, drop = F]
    }
    all_diffbind_tables[[paste(exp_nam, cont_nam)]] = tab
    print(paste0(cont_nam, " = ", nres))
    if(nres < 1){
      all_ensg_lists[[paste(exp_nam, cont_nam, "up")]] = character()
      all_ensg_lists[[paste(exp_nam, cont_nam, "down")]] = character()
      next
    }
    
    ranges = IRanges(tab$start, tab$end)
    tab_granges = GRanges(seqnames = tab$seqnames, ranges = ranges, strand = tab$strand,
                          FDR = tab$FDR, FC = tab$Fold, FE1 = tab[,7], FE2 = tab[,8])
    
    tab_granges_up = tab_granges[tab_granges$FE2 > tab_granges$FE1]
    tab_granges_down = tab_granges[tab_granges$FE2 < tab_granges$FE1]
    
    if(nrow(tab) < 1) next
    colnames(tab)[1] = "chrm"
    if(exists("assign_ensg")){
      ensg = assign_ensg(tab_granges_up, my_ref)
      all_ensg_lists[[paste(exp_nam, cont_nam, "up")]] = ensg
      
      ensg = assign_ensg(tab_granges_down, my_ref)
      all_ensg_lists[[paste(exp_nam, cont_nam, "down")]] = ensg
    }
  }
}

all_gene_lists = lapply(all_ensg_lists, function(x){
  out = unique(enst_ref[x,]$gene_name)
  return(out)
})
for(cl in c("MCF10A", "MCF7")){
  sapply(all_gene_lists, length)
  is_drug = grepl("compareDrugs", names(all_gene_lists))
  is_mcf7 = grepl(paste0(cl, "_drugs"), names(all_gene_lists))
  is_pw = !grepl("!", names(all_gene_lists))
  gene_lists = all_gene_lists[is_drug & is_mcf7 & is_pw]
  names(gene_lists) = sub("diffbind_", "", names(gene_lists))
  names(gene_lists) = sub("_drugs_", " ", names(gene_lists))
  names(gene_lists) = sub("_compareDrugs.save", "", names(gene_lists))
  hidden = sapply(names(gene_lists), function(x){
    xsplit = strsplit(x, " ")[[1]]
    direct = xsplit[7]
    from = xsplit[4]
    to = xsplit[6]
    newsplit = xsplit
    newsplit[4] = to
    newsplit[6] = from
    newsplit[7] = switch (direct,
                          up = "down",
                          down = "up"
    )
    gene_lists[[paste(newsplit, collapse = " ")]] <<- gene_lists[[x]]
    
  }); remove("hidden")
  gene_lists = gene_lists[sort(names(gene_lists))]
  
  
  marks = c("H3K4ME3", "H3K4AC", "H3K27AC", "H3K27ME3")
  for(m in marks){
    k = grepl(m, names(gene_lists))
    m_gene_lists = gene_lists[k]
    nondir_names = sort(unique(sapply(strsplit(names(m_gene_lists), " "), function(x)paste(x[-length(x)], collapse = " "))))
    summary_mat = matrix(0, nrow = length(nondir_names), ncol = 2)
    colnames(summary_mat) = c("up", "down")
    rownames(summary_mat) = nondir_names
    
    for(nam in nondir_names){
      up_list = m_gene_lists[[paste(nam, "up")]]
      down_list = m_gene_lists[[paste(nam, "down")]]
      summary_mat[nam, "up"] = length(up_list)
      summary_mat[nam, "down"] = length(down_list)
      
    }
    
    wb = createWorkbook()
    summary_sheet = createSheet(wb, "summary")
    summary_cb = CellBlock(summary_sheet, startRow = 1, startColumn = 1, noRows = nrow(summary_mat) + 1, noColumns = 4)
    # body_cb = CellBlock(summary_sheet, startRow = 2, startColumn = 2, noRows = nrow(summary_mat), noColumns = 2)
    # left_cb = CellBlock(summary_sheet, startRow = 2, startColumn = 1, noRows = nrow(summary_mat), noColumns = 1)
    froms = sapply(strsplit(rownames(summary_mat), " "), function(x)x[4])
    tos = sapply(strsplit(rownames(summary_mat), " "), function(x)x[6])
    CB.setRowData(cellBlock = summary_cb, x = c("from", "to", colnames(summary_mat)), rowIndex = 1)
    CB.setColData(cellBlock = summary_cb, x = froms, colIndex = 1, rowOffset = 1)
    CB.setColData(cellBlock = summary_cb, x = tos, colIndex = 2, rowOffset = 1)
    CB.setMatrixData(cellBlock = summary_cb, x = summary_mat, startRow = 2, startColumn = 3)
    for(i in 1:4){
      autoSizeColumn(sheet = summary_sheet, colIndex = i)
    }
    
    for(nam in nondir_names){
      up_list = m_gene_lists[[paste(nam, "up")]]
      down_list = m_gene_lists[[paste(nam, "down")]]
      list_mat = matrix(NA, nrow = max(length(up_list), length(down_list)), ncol = 2)
      colnames(list_mat) = c("up", "down")
      if(length(up_list) > 0) list_mat[1:length(up_list), "up"] = up_list
      if(length(down_list) > 0) list_mat[1:length(down_list), "down"] = down_list
      
      nam_split = strsplit(nam, " ")[[1]]
      sheet_name = paste(nam_split[4], "->", nam_split[6])
      print(sheet_name)
      list_sheet = createSheet(wb, sheet_name)
      list_cb = CellBlock(list_sheet, startRow = 1, startColumn = 1, noRows = nrow(list_mat) + 1, noColumns = 2)
      
      CB.setRowData(cellBlock = list_cb, x = colnames(list_mat), rowIndex = 1)
      if(nrow(list_mat) > 0 ) CB.setMatrixData(cellBlock = list_cb, x = list_mat, startRow = 2, startColumn = 1)
      autoSizeColumn(sheet = list_sheet, colIndex = 1)
      autoSizeColumn(sheet = list_sheet, colIndex = 2)
    }
    
    saveWorkbook(wb, file = paste(cl, m, "DM_lists.enst.xlsx", sep = "_"))
  }
  
}

save(all_gene_lists, all_diffbind_tables, all_ensg_lists, file = "input/DM_res.enst.save")
