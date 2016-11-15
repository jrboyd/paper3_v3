library("DiffBind")
library("xlsx")

source("H:/projects/Apps/diff_marks_v2/scripts/enrichment_testing.R")
diffbind_res = list()
files = dir(path = "input/diffbind_drugs_r2/", pattern = "diffbind.+.save", full.names = T)
keep = grepl("compareDrugs", files)
files = files[keep]
for(f in files){
  print(f)
  load(f)
  diffbind_res[[basename(f)]] = ss
}

min_fdr=10^-2
min_FC = 2
save_filename = paste0("input/", paste0("DM_tables.ensg.padj", min_fdr, ".fc", min_FC, ".save"))

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
    
  }
}

save(all_diffbind_tables, file = save_filename)
