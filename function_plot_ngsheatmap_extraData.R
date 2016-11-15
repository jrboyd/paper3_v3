source("heatmap.ngsplots_kmeans_with_sideplot.R")
source("heatmap.ngsplots_kmeans_cluster_extraData.R")
plot_ngsheatmap_extraData = function(sel, ngs_prof, i, nclust, main_title, is_gb_prof = T, tsstes = c(.16, .84), ...){
  print(length(sel))
  if(length(sel) < 10){
    warning("too few to plot")
    layout(1)
    plot(0:1,0:1)
    text(.5,.5,'too few to plot')
    return()
  } 
  if(!exists("rna_dat")){
    load("ref/enst_ref.save")
    rna_dat = read.table('Z:/RNA_seq_data/Breast/MCF10a_vs_MCF7_vs_MDA_no_drug_treatment/DESeq2/using_MCF10A_replicates_1,2,4/DESeq2_MCF10a_MCF7_MDA_normalized_counts.txt', header = T, row.names = 1)
    rna_dat = t(apply(rna_dat, 1, function(x){
      return(c(mean(x[1:3]),mean(x[4:6]), mean(x[7:9])))
    }))
    rna_dat = rna_dat + 2
    colnames(rna_dat) = c('MCF10A', 'MCF7', 'MDA231')
    #add -1s for undetected
    ensg_needed = setdiff(unique(enst_ref$gene_id), rownames(rna_dat))
    rna_needed = matrix(2, nrow = length(ensg_needed), ncol = ncol(rna_dat))
    rownames(rna_needed) = ensg_needed
    rna_dat <<- rbind(rna_dat, rna_needed)
  }
  names(sel) = NULL
  sel_ensg = enst_ref[sel,]$gene_id
  sel_rna_dat = rna_dat[sel_ensg,]
  rownames(sel_rna_dat) = sel
#   
#   sel_ensg = enst_ref[sel,]$gene_id
#   sel_ensg = intersect(sel_ensg, rownames(rna_dat))
#   keep = sapply(sel, function(x){
#     any(enst_ref[x,]$gene_id == sel_ensg)
#   })
#   sel = sel[keep]
#   keep = !duplicated(enst_ref[sel,]$gene_id)
#   sel = sel[keep]
#   sel_ensg = enst_ref[sel,]$gene_id
  
  sel_prof = lapply(ngs_prof, function(x){
    return(x[sel,])
  })
  
  plot_rna = log2(sel_rna_dat)
  rownames(plot_rna) = sel
  exDat_linePlot = function(dat, ylim){
    as_FC = F
    xs = 1:ncol(dat)
    if(as_FC){
      dat = dat - dat[,1]
      ylim = c(-5,5)
      dat = ifelse(dat < -5, -5, dat)
      dat = ifelse(dat > 5, 5, dat)
    }
    plot(0, type = 'n', axes = F, xlim = c(.8,ncol(dat) + .2), ylim = ylim)
    hidden = apply(dat, 1, function(x){
      lines(xs, x, col = 'gray')
    })
    
    lines(xs, apply(dat, 2, median), col = "#f62626", lwd = 3)
    
  }
  cls = sapply(strsplit(names(sel_prof), '_'), function(x)x[1])
  hms = unique(sapply(strsplit(names(sel_prof), '_'), function(x)x[2]))
  half_n = length(sel_prof)/2
  side_plots = list(1:half_n, 1:half_n + half_n)
  profile_colors = RColorBrewer::brewer.pal(6, 'Dark2')[c(4:6,1:3)]
  names(profile_colors) = c("H1", "H7", "H9", "MCF10A", "MCF7", "MDA231")
  profile_colors["H7"] = 'gray'
  nr = 4 + nclust
  nc = 10
  lmat_custom = matrix(0, ncol = nc, nrow = nr)
  lmat_custom[2,nc + -3:-2] = 1
  lmat_custom[2,nc-1] = 2
  lmat_custom[nr,-3:-2 + nc] = 3
  lmat_custom[2,nc] = 4
  lmat_custom[nr-1,-3 + nc] = 5
  lmat_custom[nr-1,-2 + nc] = 6
  lmat_custom[nr-1,-1 + nc] = 7
  # for(i in 0:5){
  dashed_line_pos = .5
  if(is_gb_prof) dashed_line_pos = tsstes
  res = heatmap.ngsplots_extraData(sel_prof, profile_colors = profile_colors[cls], extraData = plot_rna, main_title = main_title,
                         doSidePlot = T, side_plots = side_plots, sidePlot_smoothing = 10, dashed_line_positions = dashed_line_pos,
                         extraData_plotFunction = exDat_linePlot, nclust = nclust, cex.col = 1.3, 
                         labelWithCounts = T, labels_right = 'genes', 
                         hmapColors = c('royalblue1','black', 'yellow'), 
                         labels_above = hms, labels_below = cls, lmat_custom = lmat_custom, ...)
  plot0();text(.5,.5, 'average profile')
  plot0();text(.5,.5, 'log2 gene\nexpression')
  plot0();legend('center', legend = unique(cls), fill = profile_colors[unique(cls)], horiz = T, bty = 'n')
  plot0();text(.5,.5, 'cluster size')
  plot0();text(.5,.5, 'H3K4me3')
  plot0();text(.5,.5, 'H3K27me3')
  plot0()
  text(.15,.9, 'MCF10A', srt = 270, adj = c(0,.5))
  text(.5,.9, 'MCF7', srt = 270, adj = c(0,.5))
  text(.85,.9, 'MDA231', srt = 270, adj = c(0,.5))
  return(res)
}