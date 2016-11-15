library("DiffBind")
library("xlsx")

# source("H:/projects/Apps/diff_marks_v2/scripts/enrichment_testing.R")
diffbind_res = list()
setwd(input_dir)
files = dir(path = "diffbind_STEIN_merged_inputs/", pattern = "diffbind.+.save", full.names = T)
keep = grepl("compareDrugs", files)
files = files[keep]
for(f in files){
  print(f)
  load(f)
  diffbind_res[[basename(f)]] = ss
}
#set in configure.R
# CHIP_MIN_FDR = 10^-2
# CHIP_MIN_FC = 2
# save_filename = paste0("input/", paste0("DM_tables.ensg.padj", CHIP_MIN_FDR, ".fc", CHIP_MIN_FC, ".save"))

all_diffbind_tables = list()
print("annotate diffbind peaks to genes...")
for(j in  1:length(diffbind_res)){
  exp_nam =  names(diffbind_res)[j]
  # if(grepl("compareReps", exp_nam)) next
  print(exp_nam)
  ss = diffbind_res[[exp_nam]]
  for(i in 1:length(ss$contrasts)){
    cont_nam = paste("from", paste(ss$contrasts[[i]][3:4], collapse = " to "))
    
    rprt = dba.report(ss, contrast = i, fold = 1, th = 1)
    
    
    tab = as.data.frame(rprt)
    keep = tab$FDR < CHIP_MIN_FDR
    nres = sum(keep)
    tab = tab[keep,, drop = F]
    if(nres > 0){
      keep = abs(tab$Fold) > CHIP_MIN_FC
      nres = sum(keep)
      tab = tab[keep,, drop = F]
    }
    all_diffbind_tables[[paste(exp_nam, cont_nam)]] = tab
    print(paste0(cont_nam, " = ", nres))
    
  }
}

# save(all_diffbind_tables, file = save_filename)
