source("functions_rGREAT_runner.R")

all_marks = c("H3K4AC", "H3K4ME3", "H3K27AC", "H3K27ME3")
all_cells = c("MCF10A", "MCF7")
base = "ctrl"
all_treat = c("bza", "e2", "e2bza", "gc10", "gc10bza")

for(treat in all_treat){
  print(paste("***", treat, "***"))
  run_rGREAT(parent_dir = "output/GREAT_comparison_round2/anywhere", 
             marks = all_marks, 
             cells = all_cells, 
             base_drug = base, 
             treat_drug = treat)
}

for(treat in all_treat){
  print(paste("***", treat, "***"))
  run_rGREAT(parent_dir = "output/GREAT_comparison_round2/promoters", 
             marks = all_marks, 
             cells = all_cells, 
             base_drug = base, 
             treat_drug = treat, 
             filter_to_promoters = T)
}

for(treat in all_treat){
  print(paste("***", treat, "***"))
  run_rGREAT(parent_dir = "output/GREAT_comparison_round2/enhancers", 
             marks = all_marks, 
             cells = all_cells, 
             base_drug = base, 
             treat_drug = treat, 
             filter_to_enhancers = T)
}

for(treat in all_treat){
  print(paste("***", treat, "***"))
  run_rGREAT(parent_dir = "output/GREAT_comparison_round2/ESR1-enhancers", 
             marks = all_marks, 
             cells = all_cells, 
             base_drug = base, 
             treat_drug = treat, 
             filter_to_enhancers = T,
             filter_to_esr1 = T)
}