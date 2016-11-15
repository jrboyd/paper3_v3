#clustering is kmeans then sorted within clusters
#v2 adds 2 clustering modes and keep by profile as default
#extra uses supplied extra_data
#manual uses supplied cluster_members list
source("functions_movingAverage.R")
#source('heatmap.2.2.R')
library('png')
library('reshape')
library('ggplot2')
library('gplots')

plot0 = function(width = 1, height = 1){
  fudge = 0.037037037037
  plot(c(0+fudge*width, width-fudge * width), c(0+fudge*height, height-fudge * height), type = 'n', xlab = '', ylab = '', axes = F)
}


startRasterMode = function(width = 1, height = 1){
  png('tmp.png', units = 'px', width = width, height = height)
  par(mai = rep(0,4))
}

stopRasterMode = function(mai = NA){
  dev.off()
  if(length(mai) > 1 && !is.na(mai)){
    par(mai = mai)
  }
  plot0()
  rasterImage(readPNG("tmp.png", native = FALSE), xleft = 0, xright = 1, ytop = 1, ybottom = 0,
              interpolate = FALSE)
}

heatmap.ngsplots_extraData2 = function(ngs_profiles, 
                                       main_title = NULL, 
                                       lmat_custom = NULL,
                                       show_layout = F,
                                       profiles_to_plot = NA, 
                                       profile_colors = NA,
                                       side_plots = list(),#should be list of index array for side plots, if not list, will not be plotted.
                                       sideplot_lwd = 2.5,
                                       sidePlotSize = 1,
                                       nclust = 6, 
                                       labels_below = NA, 
                                       labels_above = NA, 
                                       fg_toPlot = character(), 
                                       labels_right = NA, 
                                       sortClustersByTotal = F,
                                       clusterWeights = NULL,
                                       hmapColors = c('royalblue1', 'black', 'yellow'), 
                                       labelWithCounts = F, 
                                       fg_label = NA, 
                                       label_clusters = T,
                                       key.lab = 'log2 FE',
                                       key.cex = 1.3,
                                       cex.main = 2.5,
                                       cex.row = 1.3,
                                       cex.col = 2,
                                       doSidePlot = T,
                                       sidePlot_smoothing = 1,
                                       cluster_type = c("profile", "extra", "manual"),
                                       manual_cluster_members = NULL,
                                       extraData = NULL,
                                       extraData_colors = NULL,
                                       extraData_plotFunction = NA,
                                       dashed_line_positions = .5,
                                       forPDF = T,
                                       globalScale = 1){
  if(length(profiles_to_plot) == 1 && is.na(profiles_to_plot)){
    profiles_to_plot = names(ngs_profiles)
  }
  if(length(labels_below) < 2 && is.na(labels_below)){
    labels_below = profiles_to_plot
  }
  if(length(profiles_to_plot) %% length(labels_below) != 0){
    stop('length of labels_below must divide profiles_to_plot evenly!')
  }
  if(length(profiles_to_plot) %% length(labels_above) != 0){
    stop('length of labels_above must divide profiles_to_plot evenly!')
  }
  doSidePlot = T
  if(!class(side_plots) == 'list'){#sidePlots omitted if list not supplied
    doSidePlot = F
  }
  if(length(side_plots) == 0){#if empty list (default), all profiles go in single side plot
    side_plots = list(1:length(profiles_to_plot))
  }
  
  prof = matrix(0, ncol = 0, nrow = nrow(ngs_profiles[[1]]))#assemble single matrix by joining selected profiles (matrices)
  hidden = lapply(profiles_to_plot, function(x){
    prof <<- cbind(prof, ngs_profiles[[x]])
  })
  
  cr = colorRamp(hmapColors)
  tmp = 0:100 / 100
  colors = rgb(cr(tmp)/255)
  
  len = length(labels_below)
  len_major = length(labels_above)
  cseps = 1:(len-1) * ncol(prof) / len #column separators are between joined profiles
  if(len == 1) cseps = -1
  cseps_major = 1:(len_major-1) * ncol(prof) / len_major
  if(len_major == 1) cseps_major = -1
  
  labels.above = rep('', ncol(prof))#column labels are centered above each profile
  if(!is.na(labels_above[1])){
    labels_loc = round((1:length(labels_above) -.499) * ncol(prof) / length(labels_above))
    labels.above[labels_loc] = labels_above
  }
  
  if(any(is.na(profile_colors))){
    profile_colors = RColorBrewer::brewer.pal(length(profiles_to_plot), 'Set1')
  }
  
  labels.below = rep('', ncol(prof))
  labels_loc = round((1:length(labels_below)-.499) * ncol(prof) / length(labels_below))
  labels.below[labels_loc] = labels_below
  
  set.seed(1)
  clust_prof = prof
  if(!is.null(clusterWeights)){
    if(length(clusterWeights) != ncol(clust_prof)) clusterWeights = rep(clusterWeights, length(ngs_profiles))
    profile_coefficient = rep(clusterWeights, length(ngs_profiles))
    clust_prof = t(apply(clust_prof, 1, function(x)x*clusterWeights))
  }
  if(cluster_type == "profile"){
    kmclust = kmeans(clust_prof, centers = nclust, iter.max = 10)
  }else if(cluster_type == "extra"){
    kmclust = kmeans(extraData - extraData[,1], centers = nclust, iter.max = 10)  
  }else if(cluster_type == "manual"){
    if(is.null(manual_cluster_members)){
      stop("must set manual_cluster_members if cluster_type is manual")
    }
    nclust = length(manual_cluster_members)
    kmclust = list()
    kmclust$size = sapply(manual_cluster_members, length)
    kmclust$cluster =  unlist(lapply(1:nclust, function(x)rep(x, kmclust$size[x])))
    names(kmclust$cluster) = unlist(manual_cluster_members)
    kmclust$centers = t(sapply(manual_cluster_members, function(x){
      colMeans(clust_prof[x,])
    }))
    rownames(kmclust$centers) = 1:nclust
  }else{
    stop(paste("invalid cluster type", cluster_type))
  }
  
  
  
  
  
  if(length(fg_toPlot) > 0){#extract fg_plot as special cluster
    #rseps = c(0, rseps)#add an
    kmclust$cluster[fg_toPlot] = 0#change cluster of fg
    kmclust$cluster = kmclust$cluster + 1
    
    
    nclust = nclust + 1
    kmclust$size = sapply(1:nclust, function(x){
      return(sum(kmclust$cluster == x))
    })
    kmclust$centers = rbind(colMeans(extraData[fg_toPlot,]), kmclust$centers)
    rownames(kmclust$centers) = 1:nclust
  }
  
  init_clust_data = function(clust, data){
    #initize cluster object
    #add data
    #sort by cluster id
    #add nclust
    o = order(clust$cluster)
    clust$data = data[o,, drop = F]
    names(clust$cluster) = NULL
    clust$cluster = clust$cluster[o]
    clust$nclust = length(clust$size)
    return(clust)
  }
  fetch_clust = function(clust, i){
    #returns the ith cluster
    if(length(clust$size) < i) stop('fetch_clust i greater than nclust')
    keep = clust$cluster == i
    return(clust$data[keep,, drop = F])
  }
  move_clust = function(clust, i_tomove, i_dest){
    #removes i_tomove and inserts it at i_dest
    #updates data, cluster, centers, and size
    if(clust$nclust < i_tomove) stop('move_clust i_tomove greater than nclust')
    if(clust$nclust < i_dest) stop('move_clust i_dest greater than nclust')
    curr_o = 1:clust$nclust
    new_o = curr_o[curr_o != i_tomove]
    new_o = c(new_o[1:length(new_o) < i_dest], i_tomove, new_o[1:length(new_o) > i_dest - 1]); new_o
    new_data = matrix(0, nrow = 0, ncol = ncol(clust$data))
    new_cluster = integer()
    hidden = sapply(new_o, function(x){
      keep = clust$cluster == x
      new_data <<- rbind(new_data, clust$data[keep,, drop = F])
      new_cluster <<- c(new_cluster, clust$cluster[keep])
    })
    clust$data = new_data
    clust$cluster = new_cluster
    clust$centers = clust$centers[new_o,]
    clust$size = clust$size[new_o]
    
    dict = 1:clust$nclust
    names(dict) = new_o
    
    clust$cluster = dict[as.character(clust$cluster)]
    rownames(clust$centers) = 1:clust$nclust
    return(clust)
  }
  replace_clust = function(clust, i, new_datai){
    #essentially removes cluster i and insert new_data in its place
    clust$size[i] = nrow(new_datai)
    clust$centers[i,] = colMeans(new_datai)
    before = clust$cluster < i
    after = clust$cluster > i
    data = rbind(clust$data[before,,drop = F], new_datai, clust$data[after,, drop = F])
    clust$data = data
    clusters = c(clust$cluster[before], rep(i, nrow(new_datai)) ,clust$cluster[after])
    clust$cluster = clusters
    return(clust)
  }
  
  kmclust = init_clust_data(kmclust, prof)
  for(i in 1:nclust){#sort within each cluster by nest kmclust
    dat = fetch_clust(kmclust, i)
    nsubclust = min(30, nrow(dat)/10)
    tmp = prof[rownames(dat),,drop = F]
    o = order(rowSums(tmp), decreasing = T)
    tmp = tmp[o,,drop = F]
    dat = dat[o,,drop = F]
    if(nsubclust > 2){
      iclust = kmeans(tmp, nsubclust)
      o = order(iclust$cluster)
      iclust$cluster = iclust$cluster[o]
      tmp = tmp[o,,drop = F]
      dat = dat[o,,drop = F]
      clust_o = hclust(dist(iclust$centers))$order
      new_o = integer()
      for(ci in clust_o){
        keep = iclust$cluster == ci
        new_o = c(new_o, (1:length(keep))[keep])
      }
      iclust$cluster = iclust$cluster[new_o]
      tmp = tmp[new_o,,drop = F]
      dat = dat[new_o,,drop = F]
    }
    kmclust = replace_clust(kmclust, i, dat)
  }
  
  kmclust_order = 1:nclust
  if(length(fg_toPlot) > 0){
    sortByDist = F
    if(sortByDist){
      km_dist = as.matrix(dist(kmclust$centers))
      o = order(km_dist[,1])
      kmclust_order = o
    }else{
      kmclust_order = order(apply(kmclust$centers, 1, sum), decreasing = T)
      not_1 = kmclust_order != 1#move 1 to beginning
      kmclust_order = c(1, kmclust_order[not_1])
    }
    
    #    
  }else{
    if(sortClustersByTotal){#sort clusters by total of their centers
      kmclust_order = order(apply(kmclust$centers, 1, sum), decreasing = T)
    }else{
      hiclust = hclust(dist(kmclust$centers))#use hierarchical clustering on centers
      kmclust_order = hiclust$order
    }
  }
  max_colors = 8
  rColorChoices = RColorBrewer::brewer.pal(min(max_colors, nclust), 'Pastel2')#set cluster id colors
  if(nclust > max_colors){
    rColorChoices = rColorChoices[(1:nclust-1) %% max_colors + 1]#repeat colors as needed
  }
  if(length(fg_toPlot) > 0){#correct for special selection color
    rColorChoices = c('white', rColorChoices[2:length(rColorChoices)-1])
  }
  
  if(cluster_type != "manual"){
    for(i in 1:length(kmclust_order)){#sort k means clusters 
      tomove = kmclust_order[i]
      dest = i
      kmclust = move_clust(kmclust, tomove, dest)
      kmclust_order = ifelse(kmclust_order <= i, kmclust_order + 1, kmclust_order)#position of clusters less than i must be increased by 1
      kmclust_order[1:i] = 1:i#up to i is sorted
      #     
      #     
      #     new_cluster = get_kmclust(kmclust_order[i])
      #     prof_ordered = rbind(prof_ordered, new_cluster)
      #     
      #     cluster_color = rColorChoices[i]
      #     rColors = c(rColors, rep(cluster_color, nrow(new_cluster)))
    }
  }
  rColors = rColorChoices[kmclust$cluster]
  names(rColors) = rownames(kmclust$data)
  
  
  rseps = 1:(nclust-1)#cacluated row seperations
  #kmclust$size = kmclust$size[kmclust_order]
  
  for(i in 1:nclust){
    size = kmclust$size
    start = 1
    if(i > 1){
      start = sum(kmclust$size[1:(i-1)]) + 1
    }
    end = sum(kmclust$size[1:(i)])
    if(i < nclust){
      rseps[i] = end
    }
  }
  
  if(labelWithCounts){
    if(is.null(labels_right)){
      labels_right = 1:nclust
    } else {
      labels_right = paste(kmclust$size, labels_right)
    }
  }else{
    if(length(labels_right) == 1){
      labels_right = rep(labels_right, nclust)
    }else if(length(labels_right) != nclust){
      stop('length of labels_right must be 1 or = to # of clusters (+1 if fg_toPlot is specified)')
    }
  }
  if(length(fg_toPlot) == 0 && !is.na(fg_label)){
    stop('cannot add fg_label, fg_toPlot is empty')
  }
  RowSideLabels = character()
  if(label_clusters){
    RowSideLabels = 1:nclust#paste('cluster', 1:nclust)
    if(length(fg_toPlot > 0)){
      if(!is.na(fg_label)){
        RowSideLabels = c(fg_label, 1:(nclust-1))# paste("cluster", 1:(nclust-1)))#special label for fg_toPlot
      }
    }
  }else if(!is.na(fg_label)){
    RowSideLabels = c(fg_label, rep('', nclust-1))#special label for fg_toPlot
  }
  extra_spacer = 0
  if(length(fg_toPlot) > 0){
    extra_spacer = 2
  }
  args = list(
    x = prof[rownames(kmclust$data),], key.lab = key.lab, key.cex = key.cex, 
    RowSideLabels = RowSideLabels,
    RowSideColors = rColors,
    labels.below = labels.below, cexRow = cex.row, cexCol = cex.col, col = colors, 
    rowsep.major = c(rep(rseps[1], extra_spacer), rseps), 
    colsep.minor = cseps, 
    colsep.major = cseps_major,
    sepwidth.minor = .01, 
    sepwidth.major = .05, 
    labels.above = labels.above, 
    na.color = 'red', 
    labels_rowsep = c(labels_right[1], rep('', extra_spacer), labels_right[2:length(labels_right)]), 
    main = main_title, cex.main = cex.main
  )
  new_centers = matrix(0, nrow = nrow(kmclust$centers), ncol = ncol(prof))
  for(i in unique(kmclust$cluster)){
    keep = kmclust$cluster == i
    rnames = rownames(kmclust$data)[keep]
    new_centers[i,] = colMeans(prof[rnames,])
  }
  kmclust$centers = new_centers
  kmclust <<- kmclust
  hidden = heatmap.2.2(args$x, key.lab = args$key.lab, args$key.cex, #key.par = list(mai = c(.5,0,0,0)),
                       clust = kmclust, globalScale = globalScale, forPDF = forPDF,
                       RowSideLabels = args$RowSideLabels,
                       RowSideColors = args$RowSideColors,
                       labels.below = args$labels.below, cexRow = args$cexRow, cexCol = args$cexCol, col = args$col, 
                       rowsep.major = args$rowsep.major, 
                       colsep.minor = args$colsep.minor, 
                       colsep.major = args$colsep.major,
                       sepwidth.minor = args$sepwidth.minor, 
                       sepwidth.major = args$sepwidth.major, 
                       labels.above = args$labels.above, 
                       dashed_line_positions = dashed_line_positions,
                       na.color = args$na.colors, 
                       labels_rowsep = args$labels_rowsep,
                       main = args$main,
                       cex.main = args$cex.main,
                       doSidePlot = doSidePlot, 
                       side_plots = side_plots, 
                       sideplot_lwd = sideplot_lwd,
                       sidePlotSize = sidePlotSize,
                       side_plot_colors = profile_colors,
                       sidePlot_smoothing = sidePlot_smoothing,
                       extraData = extraData,
                       extraData_plotFunction = extraData_plotFunction,
                       lmat_custom = lmat_custom, 
                       show_layout = show_layout,
                       extraData_colors = extraData_colors)
  cluster_members = sapply(1:kmclust$nclust, function(x){
    return(rownames(kmclust$data)[kmclust$cluster == x])
  }) 
  ngs_res = list(class_sizes = sapply(cluster_members, length), cluster_members = cluster_members, colors = unique(rColors), as_plotted = kmclust$data, kmclust = kmclust, args = args)
  return(ngs_res)
}

#dev.off()
#dummy out tracing, dummy out clustering
heatmap.2.2 = function (x,
                        clust,
                        col, 
                        lmat_custom = NULL,
                        show_layout = F,
                        globalScale = 1,
                        forPDF = T,
                        #dividing up the plot
                        colsep.minor = -5, 
                        colsep.major = -1, 
                        rowsep.minor = -1, 
                        rowsep.major = -1, 
                        sepwidth.minor = .01, 
                        sepwidth.major = .05, 
                        na.color = par("bg"), 
                        sidePlotSize = 1,
                        sideplot_lwd = 2.5,
                        #color and label for clusters
                        RowSideLabels = NA,
                        RowSideColors = NULL, 
                        
                        #title
                        main = NULL,
                        cex.main = 2,
                        
                        #labelling rows and columns
                        cexRow = 0.2 + 1/log10(nr), 
                        cexCol = 0.2 + 1/log10(nc),
                        labels.below = NULL, 
                        labels.above = NULL, 
                        labels_rowsep = NULL,
                        
                        #left margin size
                        left_mai = .2,
                        
                        #color key
                        key = T, 
                        key.height = 1.5,
                        key.lab = 'Color Key', 
                        key.cex = 1,
                        key.xtickfun = NULL, 
                        key.par = list(),
                        
                        #side plot of summary
                        doSidePlot = T,
                        side_plot_colors = NA,
                        side_plots = list(),
                        sidePlot_smoothing = 1,
                        dashed_line_positions = .5,
                        extraData = NULL,
                        extraData_colors = NULL,
                        extraData_plotFunction = NA) 
{
  
  x.original = x
  rowsep.minor = -1
  rowsep.minor.original = rowsep.minor
  rowsep.major.original = rowsep.major
  colsep.minor.original = colsep.minor
  colsep.major.original = colsep.major
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  
  dev.width = par('din')[1]
  dev.height = par('din')[2]
  nclust = length(clust$size)
  
  body_height = dev.height
  body_width = dev.width
  body_iy = 1:nclust
  body_ix = 1:2
  lmat = matrix(rep(1, nclust), ncol = length(body_ix), nrow = nclust)
  key_width = min(4,dev.width/3)
  lwid = c(key_width, body_width - key_width)
  lhei = rep(body_height /nclust, nclust)
  
  label.height = .5
  main.height = .8
  extraDataSize = 1
  divLinesSize = .3
  
  labRowSize = 1
  
  add_lmat_left = function(added_width){
    lmat <<- cbind(rep(0, nrow(lmat)), lmat)
    added_width = added_width * globalScale
    lwid[max(body_ix)] <<- lwid[max(body_ix)] - added_width
    lwid <<- c(added_width, lwid)
    body_ix <<- body_ix + 1
  }
  add_lmat_right = function(added_width, solid = T){
    if(solid){
      lmat <<- cbind(lmat, c(rep(0, min(body_iy)-1), rep(max(lmat)+1, nclust), rep(0, nrow(lmat)-max(body_iy))))
    }else{
      lmat <<- lmat <- cbind(lmat, c(rep(0, min(body_iy)-1), (max(lmat)+1):(max(lmat) + nclust), rep(0, nrow(lmat)-max(body_iy))))
    }
    added_width = added_width * globalScale
    lwid[max(body_ix)] <<- lwid[max(body_ix)] - added_width
    lwid <<- c(lwid, added_width)
  }
  
  add_lmat_top = function(added_height){
    lmat <<- rbind(max(lmat) + 1, lmat)
    added_height = added_height * globalScale
    lhei[body_iy] <<- lhei[body_iy] - added_height / nclust
    lhei <<- c(added_height, lhei)
    body_iy <<- body_iy + 1
  }
  add_lmat_bottom = function(added_height, body_xpos = -1){
    if(body_xpos > 0){
      new_row = rep(0, ncol(lmat))
      new_row[body_ix[body_xpos]] =  max(lmat) + 1
      lmat <<- rbind(lmat, new_row)
    }else{
      lmat <<- rbind(lmat, max(lmat) + 1)
    }
    added_height = added_height * globalScale
    lhei[body_iy] <<- lhei[body_iy] - added_height / nclust
    lhei <<- c(lhei, added_height)
  }
  if(!is.null(labels.above)){
    add_lmat_top(label.height)
  }
  if(!is.null(main)){
    add_lmat_top(main.height)
  }
  if(!is.null(labels.below)){
    add_lmat_bottom(label.height)
  }
  if(key){
    add_lmat_bottom(key.height, body_xpos = 1)
  }
  RowSideColors_size = .5
  if (!is.null(RowSideColors)) {
    if (!is.character(RowSideColors) || length(RowSideColors) != nrow(x)) 
      stop("'RowSideColors' must be a character vector of length nrow(x)")
    add_lmat_right(RowSideColors_size, solid = T)
  }
  if(doSidePlot | !is.null(extraData)){
    add_lmat_right(divLinesSize, solid = T)
  }
  for(i in 1:length(side_plots)){
    if(doSidePlot){
      add_lmat_right(sidePlotSize, solid = F)
    }
  }
  if(!is.null(extraData)){
    add_lmat_right(extraDataSize, solid = F)
  }
  
  if (!is.null(labels_rowsep)) {
    add_lmat_right(labRowSize, solid = T)
  }
  if(left_mai > 0){
    add_lmat_left(left_mai)
    lmat <- cbind(rep(0, nrow(lmat)), lmat)
    lwid[body_ix] = lwid[body_ix] - left_mai
    lwid <- c(left_mai, lwid)
    body_ix = body_ix + 1
  }
  if(!is.null(lmat_custom)){
    if(class(lmat_custom) != 'matrix'){
      stop('class of lmat_custom must be matrix')
    }
    if(any(dim(lmat_custom) != dim(lmat))){
      stop(paste("lmat_custom does not match lmat. dim was", 
                 paste(dim(lmat_custom), collapse = ', '), 
                 "dim should be", 
                 paste(dim(lmat), collapse = ', ')))
    }
    #print(dim(lmat_custom))
    #print(dim(lmat))
    lmat_custom = ifelse(lmat_custom > 0, lmat_custom + max(lmat), 0)
    lmat = lmat + lmat_custom
  }
  #print(lmat)
  #print(lhei)
  #print(lwid)
  if(show_layout){
    print(lmat)
    print(lhei)
    print(lwid)
  } 
  
  nf = layout(lmat, heights = lhei, widths = lwid)
  
  if(show_layout){
    layout.show(nf)
  } 
  
  breaks <- length(col) + 1
  symbreaks = any(x < 0, na.rm = TRUE)
  symkey = F#any(x < 0, na.rm = TRUE) || symbreaks
  
  if (length(breaks) == 1) {
    if (!symbreaks) 
      breaks <- seq(min(x, na.rm = T), max(x, na.rm = T), 
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function") 
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  
  x = x[nrow(x):1,]#reverse so heatmap top row is top row of input matrix
  colorClasses = unique(RowSideColors)
  colorClasses = colorClasses[(1:nclust-1) %% length(colorClasses) + 1]
  RowSideColors = rev(RowSideColors)
  nr <- nrow(x)
  nc <- ncol(x)
  
  insert_col = function(n, i){#insert n column following index i
    x <<- cbind(x[,1:i], matrix(NA, nrow = nrow(x), ncol = n),  x[,(i+1):ncol(x)])
    affected = colsep.minor > i
    colsep.minor[affected] <<- colsep.minor[affected] + n
    affected = colsep.major > i
    colsep.major[affected] <<- colsep.major[affected] + n
    if(!is.null(labels.below)) labels.below <<- c(labels.below[1:i], rep('', n), labels.below[(i + 1):length(labels.below)])
    if(!is.null(labels.above)) labels.above <<- c(labels.above[1:i], rep('', n), labels.above[(i + 1):length(labels.above)])
  }
  
  insert_row = function(n, i){#insert n rows following index i
    x <<- rbind(x[1:i,], matrix(NA, ncol = ncol(x), nrow = n),  x[(i+1):nrow(x),])
    affected = rowsep.minor >= i
    rowsep.minor[affected] <<- rowsep.minor[affected] + n
    affected = rowsep.major >= i
    rowsep.major[affected] <<- rowsep.major[affected] + n
    RowSideColors <<- c(RowSideColors[1:i], rep(NA, n), RowSideColors[(i+1):length(RowSideColors)])
  }
  
  do_col_major = colsep.major.original[1] > -1#test before value are changed
  if (colsep.minor.original[1] > -1){ #insert  columns in x and color
    for (i in 1:length(colsep.minor)){ 
      csep = colsep.minor[i]
      insert_col(round(sepwidth.minor*nc), csep)
    }
  }
  if (do_col_major){ #insert  columns in x and color
    for (i in 1:length(colsep.major)){ 
      csep = colsep.major[i]
      insert_col(round(sepwidth.major*nc), csep)
    }
  }
  x = x[nrow(x):1,]#reverse rows temporarily
  RowSideColors = rev(RowSideColors)
  do_row_major = rowsep.major[1] > -1#test before values are changed
  if (rowsep.minor.original[1] > -1){
    for (i in 1:length(rowsep.minor)){ 
      rsep = rowsep.minor[i]
      insert_row(round(sepwidth.minor*nr), rsep)
    }
  }
  if (do_row_major){
    for (i in 1:length(rowsep.major)){ 
      rsep = rowsep.major[i]
      insert_row(round(sepwidth.major*nr), rsep)
    }
  }
  
  #unreverse rows
  x = x[nrow(x):1,]
  RowSideColors = rev(RowSideColors)
  nc = ncol(x)
  nr = nrow(x)
  
  #plots the body of the heatap
  if(forPDF) startRasterMode(width = nc, height = nr)
  # save(x.original, x, nc, nr, breaks, col, file = "tmp.save")
  if(!forPDF) par(mai = rep(0,4))
  image(1:nc, 1:nr-1, t(x), xlim = c(0, nc), ylim = c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks)
  if(forPDF) stopRasterMode(mai = rep(0,4))
  
  #plots column labels above column
  srtCol = 0
  adjCol = c(.5,1)
  colCol = NULL
  
  apply_col_labels = function(col_labels, yshift = .5, yadj = .5){
    plot0(ncol(x))
    xpd.orig <- par("xpd")
    par(xpd = NA)
    xpos <- 1:nc-.5
    n = sum(col_labels != '')
    n_scale = 1 #/ n  #no longer necessary to resize text according to columns, plot gets wider
    
    text(x = xpos, y = yshift, labels = col_labels, adj = c(.5, yadj),
         cex = cexCol * globalScale * n_scale, srt = srtCol, col = colCol)
    par(xpd = xpd.orig)
  }
  if(!is.null(labels.above)){
    apply_col_labels(labels.above, .1, 0)
  }
  
  #draw title
  if(!is.null(main)){
    plot.new()
    par(xpd = NA)
    text(0,1, main, cex = cex.main * globalScale, adj = c(0,1))
  }
  
  #draws column labels below column
  if(!is.null(labels.below)){
    apply_col_labels(labels.below, .9, 1)
  }
  
  
  if (key) {
    mai <- c(.5, 0,.5,0) * globalScale
    par(mai = mai, cex = 0.75 * globalScale, mgp = c(2, 1, 0))
    if (length(key.par) > 0) 
      do.call(par, key.par)
    tmpbreaks <- breaks
    if (F) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)-.001
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)+.001
    }
    else {
      min.raw <- min(x, na.rm = T)
      max.raw <- max(x, na.rm = T)
    }
    z <- seq(min.raw, max.raw, by = min(diff(breaks)/4))
    
    if(forPDF) startRasterMode(width = length(z), height = 1)
    #draw the key color gradient
    if(!forPDF) par(mai = mai)
    image(z = matrix(rep(z,2), ncol = 2), xlim = c(0,1), ylim = c(0,1), col = col, breaks = tmpbreaks, 
          xaxt = "n", yaxt = "n")
    if(forPDF) stopRasterMode(mai = mai)
    if (is.null(key.xtickfun)) {
      lv <- -100:100 * 2
      keep = lv > min.raw & lv < max.raw
      lv = lv[keep]
      xmin = par('usr')[1]
      xmax = par('usr')[2]
      lvmin = min(lv)
      lvmax = max(lv)
      xv <- (lv - min.raw) / (max.raw - min.raw) * (xmax - xmin) + xmin
      xargs <- list(at = xv, labels = lv)
    }
    else {
      xargs <- key.xtickfun()
    }
    xargs$side <- 1
    par(xpd = NA)
    do.call(axis, xargs)
    #     MAX = max(abs(c(min(x), max(x))))
    #     par(usr = c(-MAX, MAX, 0, 1))
    #     axis(side = 1, )
    if (!is.na(key.lab)) {
      mtext(side = 3, key.lab, line = 2, padj = 1, 
            cex = key.cex * globalScale)
    }
  }
  if(!is.null(RowSideColors)){#plot colored block for each cluster
    tmp = as.factor(RowSideColors)
    vals = matrix(as.numeric(tmp), ncol = 1)
    navals = is.na(vals)[,1]
    naval_widths = integer(nclust - 1)#determine length of each block of NAs
    j = 1
    for(i in 1:length(naval_widths)){
      while(!navals[j]){
        j = j + 1
      }
      k = j
      while(navals[k]){
        k = k + 1
      }
      naval_widths[i] = k - j
      j = k
    }
    naval_widths = cumsum(c(0, naval_widths))
    lev = levels(tmp)
    mai = c(0,.05,0,.05) * globalScale
    if(forPDF) startRasterMode(height = nr, width = 2)
    #par(mai = c(0,.1,0,.1))
    if(!forPDF) par(mai = mai)
    image(0:1, 1:nr-1, t(cbind(vals, vals)), xlim = c(0, 1), ylim = c(0, nr), axes = FALSE, xlab = "", ylab = "", col = lev)
    if(forPDF) stopRasterMode(mai = mai)
    cluster_levels = 1:nclust
    # cluster_levels = cluster_levels[!is.na(cluster_levels)]
    cluster_ids = vals[,1]
    names(cluster_ids) = 1:length(cluster_ids)
    # cluster_ids = cluster_ids[!is.na(cluster_ids)]
    par(xpd = NA)
    # sep_added = (0:(nclust-1)) * sepwidth.minor
    cluster_starts = cumsum(c(0,rev(clust$size)[1:(nclust-1)]))+1 + naval_widths
    cluster_ends = cumsum(rev(clust$size)) + naval_widths
    for(i in 1:nclust){
      # save(clust, file = 'tmp.save')
      keep = (cluster_ids == cluster_levels[i])
      s = cluster_starts[i]
      e = cluster_ends[i]
      center = mean(c(s, e))
      if(forPDF) center = center / max(cluster_ends) #stupid bandaid correction for start/stopRaster model
      rowLab = rev(RowSideLabels)[i]
      text(.5,center, rowLab, adj = c(.5,.5), cex = cexRow * globalScale)
      
    }
    #par(xpd = NA)
    
  }
  if(doSidePlot | !is.null(extraData)){
    par(mai=rep(0,4))
    plot0()
    #     (x=c(0,1),frame.plot=FALSE, y=c(0,1), xaxt='n',yaxt='n', type="n", xlab="",
    #          ylab="",xlim=c(0,1),ylim=c(0,1))
    #no idea how these 'actual' margins are selected
    MIN=0#-.04
    MAX=1#1.04
    RNG=MAX-MIN
    
    ### draw the dotted lines connecting line plots to heatmap
    #line at top and bottom
    lines(c(0,1),c(0,0),lty=2)
    lines(c(0,1),c(1,1),lty=2)
    for(i in 1:(nclust-1)){
      #calculate number of rows covered so far as fraction of total
      hFrac_a = (cluster_ends[i]) / (max(cluster_ends))
      hFrac_b = (cluster_starts[i+1]-1) / (max(cluster_ends))
      hFrac_mean = mean(c(hFrac_a, hFrac_b))
      lplotFrac = i/nclust
      #x is always from min to max
      #y goes from variable fraction on heatmap to constant fraction of plotting column
      meet_in_middle = 0#fraction of x range where where 2 lines from hmap should join before meeting sideplots
      if(meet_in_middle > 0){
        MIM = 1 - meet_in_middle
        lines(c(MIN, MAX - MIM*RNG),c(RNG*hFrac_a,RNG*hFrac_mean),lty=2)
        lines(c(MIN, MAX - MIM*RNG),c(RNG*hFrac_b,RNG*hFrac_mean),lty=2)
      }
      
      lines(c(MIN + meet_in_middle * RNG, MAX), c(RNG*hFrac_mean,RNG*lplotFrac),lty=2)
    }
  }
  if(doSidePlot){
    for(i_side in 1:length(side_plots)){
      this_col = side_plots[[i_side]]
      avgA = matrix(0,nclust, ncol(clust$centers))
      for(i in 1:nrow(avgA)){
        avgA[i,] = clust$centers[i,]
      }
      
      #need to reverse order
      #o=o[nclust:1,]
      
      par(mar = c(0,0, 0,0))
      
      #draw line chart represented each class from clustering
      nsplits = length(colsep.minor)+1
      win = ncol(avgA) / nsplits
      prof_colors = side_plot_colors
      if(any(is.na(side_plot_colors))){
        side_plot_colors = RColorBrewer::brewer.pal(nsplits, 'Set1')
      }
      for (i in 1:nclust) { 
        xrange <- as.numeric(range(1:(ncol(avgA)/nsplits)))
        yrange <- range(avgA)
        xspace = (xrange[2] - xrange[1]) * .15
        yspace = (yrange[2] - yrange[1]) * .15
        # set up the plot 
        #op <- par(mar = rep(.01, 4))
        plot(xrange, yrange, xaxt='n',yaxt='n', type="n", xlab="",
             ylab="", ylim = c(min(yrange)-yspace, max(yrange)+yspace), xlim = c(min(xrange)-xspace, max(xrange)+xspace))  #c(minVal - .1*abs(minVal - maxVal),maxVal + .1*abs(minVal - maxVal))) 
        colors <- prof_colors[i] 
        linetype <- c(1:nclust) 
        plotchar <- seq(18,18+nclust,1)
        
        #   axis(side=1,tick=TRUE,at=days)
        vals <- avgA[i,]
        for(s in 1:nsplits){
          if(any(s == this_col)){
            start = (s - 1) * win + 1
            end = s * win
            xs = 1:(ncol(avgA)/nsplits)#first half of profile
            if(sidePlot_smoothing < 2){
              lines(xs, vals[start:end], type="l", lwd= sideplot_lwd,
                    lty=1, col=prof_colors[s], pch=16) 
            }else{
              lines(xs, movingAverage(vals[start:end], n = sidePlot_smoothing, centered = T), type="l", lwd=2.5,
                    lty=1, col=prof_colors[s], pch=16) 
            }}
        }
        for(dpos in dashed_line_positions){
          dx = (xrange[2] - xrange[1])*dpos + xrange[1]
          lines(rep(dx, 2), yrange, lty = 3)
        }
        
      }
    }
  }
  if(!is.null(extraData)){
    par(mai = rep(0, 4), xpd = T)
    YMIN = min(extraData)
    YMAX = max(extraData)
    colorClasses = extraData_colors
    for(i in 1:nclust){
      subset = rownames(clust$data)[clust$cluster == i]
      dat = extraData[subset,, drop = F]# + 64
      colnames(dat) = NULL
      if(is.function(extraData_plotFunction)){
        extraData_plotFunction(dat, c(YMIN,YMAX))
      }else{
        M = colMeans(dat)
        n = nrow(dat)
        s = apply(dat, 2, sd)
        SE = s / sqrt(n)
        E = 0
        if(n > 1) E = qt(.975, df=n - 1) * SE
        low = M - E
        mid = M
        high = M + E
        if(any(is.na(low))) low = mid
        if(any(is.na(high))) high = mid
        
        low = ifelse(low <= 0, 1, low)
        barplot2(mid, ci.l = low, ci.u = high, col = colorClasses, plot.ci = T, ylim = c(YMIN,YMAX), axes = F, space = 0, bty = 'o')
        for(pos in c(4,8,12,16)){
          lines(c(0,nclust+1), c(pos,pos), lty = 2) 
        }
        
        #axis(side = 2, at = c(250,500,1000,2000))
        box()
      }
    }
  }
  par(xpd = NA)
  if(!is.null(labels_rowsep)){
    par(mai = rep(0,4))
    plot0(height = nr)
    for(i in 1:length(labels_rowsep)){
      if(!doSidePlot && is.null(extraData)){
        rsep_prev = 1
        if(i > 1){
          rsep_prev = rowsep.major[i - 1] + sepwidth.major
        }
        rsep_curr = rowsep.major.original[i]
        if(i > 1){
          rsep_curr = rsep_prev - 1 + rowsep.major.original[i] - rowsep.major.original[i - 1]
        }
        if(i > length(rowsep.major.original)){
          rsep_curr = nrow(x)
        }
        rsep_mid = (mean(c(rsep_prev, rsep_curr)))
        ypos = nr - rsep_mid + 1
      }else{
        ypos = nr - (i - .5) / nclust * nr
      }
      
      text(.5, ypos, labels_rowsep[i], adj = c(.5,.5), cex = cexRow * globalScale )
    }
  }
}

do_example = F
if(do_example){
  
  nr = 100
  nc = 50
  nmods = 2
  nlines = 3
  
  prof_normed = lapply(1:(nmods * nlines),function(x){
    set.seed(x)
    out = matrix(runif(nr*nc)/10, nrow = nr, ncol = nc)
    affected = runif(nr) > runif(1)
    out[affected,] = out[affected,] * 10
    out = out * 100
    rownames(out) = paste('gene', 1:nr)
    return(out)
  })
  
  extraData = matrix(0, ncol = nlines, nrow = nr)
  rownames(extraData) = rownames(prof_normed[[1]])
  set.seed(1)
  for(i in 1:ncol(extraData)){
    
    out = runif(nr)
    affected = runif(nr) > runif(1)
    out[affected] = out[affected] + 5
    extraData[,i] = out
  }
  
  histone_mods = paste('mod', 1:nmods)
  cell_lines = paste('line', 1:nlines)
  prof_names = 1:(nmods*nlines)
  i = 1
  for(hm in histone_mods){
    for(cl in cell_lines){
      prof_names[i] = paste(cl, hm)
      i = i + 1
    }
  }
  names(prof_normed) = prof_names
  # pdf('compound_plots.pdf', width = 10, height = 8)
  nclust = 2
  lmat_nr = 3 + nclust
  lmat_nc = 11
  lmat_custom = matrix(0, ncol = lmat_nc, nrow = lmat_nr)
  lmat_custom[lmat_nr-1,lmat_nc-2] = 1
  lmat_custom[lmat_nr-1,lmat_nc-1] = 2
  lmat_custom[lmat_nr,-2:-1 + lmat_nc] = 3
  lmat_custom[1,lmat_nc] = 4
  # for(i in 0:5){
  main_title = 'simulated data'
  profile_colors = rep(RColorBrewer::brewer.pal(2, 'Set1')[1:2],3)
  profile_colors[6] = 'yellow'
  res = heatmap.ngsplots_extraData(ngs_profiles = prof_normed, nclust = nclust, extraData = extraData,
                                   profile_colors = profile_colors,side_plots = list(1:2,3:4, 5:6), labels_below = rep(cell_lines,2), 
                                   labels_above = histone_mods, lmat_custom = lmat_custom, labelWithCounts = T, labels_right = 'g')
  plot0();text(.5,.5, 'average profile')
  plot0();text(.5,.5, 'log gene expression')
  plot0();legend('center', legend = cell_lines, fill = RColorBrewer::brewer.pal(3, 'Set1'), horiz = T, bty = 'n')
  plot0();text(.5,.5, 'cluster size')
  
}
# dev.off()