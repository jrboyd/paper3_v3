source('heatmap.ngsplots_kmeans_with_sideplot.R')

my_ngs_heatmap = function(my_prof, ensg, main_title = "",
                          cells = "MCF7", 
                          drugs = c("ctrl", "bza", "e2bza", "e2"), by_drug = F,
                          marks = c("H3K4ME3", "H3K4AC", "H3K27AC", "H3K27ME3"),
                          pdfname = NULL, nclust = 6, height = 8, vert_lines = c(.2*5/7, 1 - .2*5/7)){
  id = gsub(':','', paste(strsplit(date(), ' ')[[1]][4], collapse = '_'))
  # browser()
  if(is.null(pdfname)) pdfname = paste0("validate_", id, ".pdf")
  
  sel_prof = lapply(my_prof, function(x){
    return(x[ensg,])
  })
  
  is_cell = sapply(strsplit(names(sel_prof), "_"), function(x){
    any(x[1] == cells)
  })
  is_drug = sapply(strsplit(names(sel_prof), "_"), function(x){
    any(x[2] == drugs)
  })
  is_mark = sapply(strsplit(names(sel_prof), "_"), function(x){
    any(x[3] == marks)
  })
  sel_prof = sel_prof[is_cell & is_drug & is_mark]
  #sort by drugs
  for(d in rev(drugs)){
    o = order(sapply(strsplit(names(sel_prof), "_"), function(x){
      x[2] == d
    }), decreasing = T)
    sel_prof = sel_prof[o]
  }
  #sort by marks
  for(m in rev(marks)){
    o = order(sapply(strsplit(names(sel_prof), "_"), function(x){
      x[3] == m
    }), decreasing = T)
    sel_prof = sel_prof[o]
  }
  
  # pdf(pdfname, width = 10, height = height)
  if(by_drug){
    n_per_mark = length(cells) * length(marks)
  }else{
    n_per_mark = length(cells) * length(drugs)  
  }
  
  
  
  
  nc = 7 + n_per_mark
  lmat_custom = matrix(0, nrow = 4 + nclust, ncol = nc)
  lmat_custom[4 + nclust,7:(nc-1)] = 1
  lmat_custom[3 + nclust,7:(nc-1)] = 2:(n_per_mark+1)
  if(main_title == ""){
    main_title = paste(sep = "\n", 
                       paste(cells, collapse = ", "),
                       paste(drugs, collapse = ", "),
                       paste(marks, collapse = ", "))
  }
  lab_below = sapply(strsplit(names(sel_prof), "_"), function(x){paste(x[1:2], collapse = "_")})
  #assign colors
  mark_colors = c("H3K4AC" = "darkblue", "H3K27ME3" = "darkred", "H3K4ME3" = "darkgreen", "H3K27AC" = "darkorange")
  drug_colors = c("ctrl" = "darkgray", "bza" = "black", "e2bza" = "#fdae61", "e2" = "#d7191c", "gc10bza" = "#a6d96a", "gc10" = "#1a9641")
  if(by_drug){
    pc = drug_colors[sapply(strsplit(names(sel_prof), "_"), function(x){x[2]})]
    sp = lapply(marks, function(x){
      which(grepl(x, names(sel_prof)))
    })
  }else{
    pc = mark_colors[sapply(strsplit(names(sel_prof), "_"), function(x){x[3]})]
    sp = lapply(unique(lab_below), function(x){
      which(lab_below == x)
    })
  }
  #how should aggregate profiles be grouped in side plots
  
  
  res = heatmap.ngsplots(ngs_profiles = sel_prof, cex.main = 1,
                         nclust = nclust,
                         profile_colors = pc,
                         side_plots = sp,
                         labels_above = marks, 
                         labels_below = lab_below, 
                         cex.col = .6,
                         lmat_custom = lmat_custom,
                         sidePlot_smoothing = 5,
                         dashed_line_positions = vert_lines,
                         main_title = main_title,
                         labelWithCounts = T,
                         labels_right = "genes"
  )
  if(by_drug){
    plot0(); legend('center', legend = drugs, fill = drug_colors[drugs], horiz = T, bty = 'n')
    for(cl in unique(marks)){
      plot0(); text(.5,.5, cl)
    }
  }else{
    plot0(); legend('center', legend = marks, fill = mark_colors[marks], horiz = T, bty = 'n')
    for(cl in unique(lab_below)){
      plot0(); text(.5,.5, cl)
    }
  }
  
  
  # dev.off()                 
  #   source('function_write_myxlsx.R')
  #   write_myxlsx_clusters(res, filename = sub(".pdf", ".xlsx", pdfname))
  return(res)
}
