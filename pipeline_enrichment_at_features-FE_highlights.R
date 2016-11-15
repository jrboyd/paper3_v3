library(magrittr)
source('setup.R')
# source("module_STEIN_ChIPseq.R")
source("module_enhancers_Carroll.R")
source("function_plot_ngsheatmap_extraData.R")
#consider enrichment by ngsplot and FE at promoters and enhancers
setwd(input_dir)

setwd(input_dir)
count_types = dir("STEIN_ChIPseq_at_features/counts/drugs_with_merged_inputs/", full.names = T)
names(count_types) = basename(count_types)
count_files = lapply(count_types, function(x){
  dir(x, pattern = "pool", full.names = T)
})



# tmp = lapply(count_files, function(x){
#   x = basename(x)
#   sapply(strsplit(x, "\\."), function(y)y[1])
# })
count_data = lapply(count_files, function(x){
  n = basename(x)
  n = sapply(strsplit(n, "\\."), function(y)y[1])
  names(x) = n
  lapply(x, function(y){
    print(paste("load:", y))
    read.table(y, row.names = 1)
  })
})

bed_files = dir("STEIN_ChIPseq_at_features/", pattern = ".bed", full.names = T)
bed_dfs = list()
bed_grs = list()
for(f in bed_files){
  key = f %>% basename %>% 
    sub(pattern = ".bed", replacement = "")
  bed = read.table(f)
  colnames(bed) = c("chrm", "start", "end", "id", "score", "strand")
  bed$ucsc = paste0(bed$chrm, ":", bed$start, "-", bed$end)
  strand = bed$strand
  strand = sub(".", "*", strand)
  gr = GRanges(seqnames = bed[,1], ranges = IRanges(start = bed[,2], end = bed[,2]), strand = strand)
  gr$ucsc = paste0(bed[,1], ":", bed[,2], "-", bed[,3])
  names(gr) = bed[,4]
  rownames(bed) = as.character(bed$id)
  bed_dfs[[key]] = bed
  bed_grs[[key]] = gr
}


blacklist_files = character()#dir("STEIN_ChIPseq_at_features/", pattern = ".txt", full.names = T)


#load blacklists and remove
for(f in blacklist_files){
  key = f %>% basename %>% 
    sub(pattern = "blacklist_", replacement = "", f) %>% 
    sub(pattern = ".txt", replacement = "")
  blacklist = read.table(f, stringsAsFactors = F)[,1]
  whitelist = setdiff(rownames(count_data[[key]][[1]]), blacklist)
  count_data[[key]] = lapply(count_data[[key]], function(x){
    return(x[whitelist,, drop = F])
  })
}

count_raw = lapply(count_data, function(x){
  mat = matrix(0, nrow = nrow(x[[1]]), ncol = length(x))
  rownames(mat) = rownames(x[[1]])
  colnames(mat) = names(x)
  for(n in names(x)){
    print(n)
    dat = x[[n]]
    mat[rownames(dat), n] = dat[,1]
  }
  colnames(mat)  = colnames(mat) %>% sub(pattern = "_pooled", replacement = "")
  return(mat)
})



count_norm = lapply(count_data, function(x){
  mat = matrix(0, nrow = nrow(x[[1]]), ncol = length(x))
  rownames(mat) = rownames(x[[1]])
  colnames(mat) = names(x)
  for(n in names(x)){
    print(n)
    dat = x[[n]]
    mat[rownames(dat), n] = dat[,1]
  }
  colnames(mat)  = colnames(mat) %>% sub(pattern = "_pooled", replacement = "")
  norm = apply(mat, 2, function(col){
    col = col / sum(col) * 10^6
  })
  return(norm)
})

fe = lapply(count_norm, function(x){
  cells = strsplit(colnames(x), "_") %>% sapply(function(x){x[1]})
  drugs = strsplit(colnames(x), "_") %>% sapply(function(x){x[2]})
  marks = strsplit(colnames(x), "_") %>% sapply(function(x){x[3]})
  is_mark = marks != "input"
  mat = matrix(0, nrow = nrow(x), ncol = sum(is_mark))
  rownames(mat) = rownames(x)
  colnames(mat) = colnames(x)[is_mark]
  pseudo = .001
  for(n in colnames(mat)){
    treat = x[,n] + pseudo
    mark = strsplit(n, "_")[[1]][3]
    control_n = sub(mark, "input", n)
    control = x[,control_n] + pseudo
    mat[,n] = treat / control
  }
  return(mat)
})



fe_tidy = lapply(fe, function(x){ #convert fe matrix to a tidy data frame
  id = factor(rep(rownames(x), ncol(x)))
  #cell_line = factor(unlist(lapply(strsplit(colnames(x), "_"), function(split_str){rep(split_str[1], nrow(x))})))
  cell_line = x %>% colnames %>% strsplit(split = "_") %>% #using pipes for fun
    lapply(FUN = function(split_str){
      rep(split_str[1], nrow(x))
    }) %>% unlist %>% factor
  
  drug = factor(unlist(lapply(strsplit(colnames(x), "_"), function(split_str){
    rep(split_str[2], nrow(x))
  })))
  mark = factor(unlist(lapply(strsplit(colnames(x), "_"), function(split_str){
    rep(split_str[3], nrow(x))
  })))
  enrichment = unlist(apply(x, 2, list))
  log2_enrichment = log2(enrichment)
  df = data.frame(id = id, 
                  cell_line = cell_line, 
                  drug = drug, 
                  mark = mark, 
                  enrichment = enrichment,
                  log2_enrichment = log2_enrichment)
  rownames(df) = NULL
  return(df)
})

library(ggplot2)
library(xlsx)
setwd(output_dir)
out = "fe_plots_highlights_common"
dir.create(out)
setwd(out)
# intersect_condition = GSE40129$MCF7_e2_ESR1_idr

##clean up enhancers by removing k4me3 bound
keep = GSE40129$MCF7_enhancer_sets$k4me3_free
k4me3_free_enhacers = merge(fe_tidy$Carroll_enhancers, data.frame(id = keep), by = "id")
fe_tidy$Carroll_enhancers = k4me3_free_enhacers

##set conditions that will be highlighted
#present in e2 only
# intersect_condition = GSE40129$MCF7_enhancer_sets$e2_only
##present in ctrl only
# intersect_condition = GSE40129$MCF7_enhancer_sets$ctrl_only
##present in both
intersect_condition = GSE40129$MCF7_enhancer_sets$common
intersection_gr = GSE40129$MCF7_any_enhancer[intersect_condition]


for(plot_cell_line in c("MCF7", "MCF10A")){
  plot_x_drug = "ctrl"
  # for(feature_type in names(fe_tidy)){
  for(feature_type in names(fe_tidy)[1]){
    print(paste(plot_cell_line, feature_type))
    for(plot_mark in c("H3K27AC", "H3K4ME3", "H3K4AC", "H3K27ME3")){
      print(paste0("  ", plot_mark))
      out_name = paste(plot_cell_line, feature_type, plot_mark, sep = "_")
      pdf(paste0(out_name, ".pdf"))
      diff_to_write = list()
      for(plot_y_drug in c("e2", "bza", "e2bza", "gc10", "gc10bza")){
        print(paste0("    ", plot_y_drug))
        
        dat = fe_tidy[[feature_type]]
        dat = subset(dat, cell_line == plot_cell_line & mark == plot_mark) # & (drug == "ctrl" | drug == "e2"), drop = T)
        
        plot_df = data.frame(x = subset(dat, drug == plot_x_drug)$log2_enrichment, 
                             y = subset(dat, drug == plot_y_drug)$log2_enrichment, 
                             id = subset(dat, drug == "ctrl")$id)
        MCF7_enhancers_intergenic
        rownames(plot_df) = plot_df$id
        plot_df$x = ifelse(plot_df$x < -2, -2, plot_df$x)
        plot_df$y = ifelse(plot_df$y < -2, -2, plot_df$y)
        keep = (plot_df$x + plot_df$y) > 1
        plot_df = plot_df[keep,]
        plot_gr = bed_grs[[feature_type]][plot_df$id]
        olaps = findOverlaps(query = plot_gr, subject = intersection_gr)
        bound = as.character(plot_df$id[queryHits(olaps)])
        # p = ggplot(data = plot_df, aes(x = x, y = y)) + geom_point()
        # # # print(p)
        # #
        # sp = p + geom_smooth(data = plot_df, method = "rlm", se = F)
        # print(sp)
        library(MASS)
        # plot_df$residy = resid(loess(y ~ x, data = plot_df))
        
        # plot_df$residx = resid(rlm(x ~ y, data = plot_df))
        
        # plot_df$mean = rowMeans(cbind(plot_df$x, plot_df$y))
        
        # p = ggplot(data = plot_df, aes(x = mean, y = residy)) + geom_point()
        # print(p)
        
        # plot_df$diff = abs(plot_df$residy) > 1 & (plot_df$x > 1 | plot_df$y > 1)
        plot_df$intersect = rep(F, nrow(plot_df))
        plot_df[bound,"intersect"] = T
        
        p = ggplot(data = plot_df, aes(x = x, y = y, col = intersect)) + 
          geom_point(data = subset(plot_df, !intersect), alpha = .05, color = "black") + 
          geom_point(data = subset(plot_df, intersect), alpha = 1, color = "red") +
          scale_color_manual(values =c("black", "red")) +
          labs(x = plot_x_drug, y = plot_y_drug, title = paste0(feature_type, ": ", plot_cell_line, " ", plot_mark))
        p = p + geom_smooth(method = "rlm", se = F)
        # p = p + facet_wrap(~intersect)
        
        startRasterMode(width = 800, height = 800)
        print(p)
        stopRasterMode()
        
        out_df = subset(plot_df, intersect)
        out_df$ucsc = bed_dfs[[feature_type]][as.character(out_df$id),"ucsc"]
        out_df$strand = bed_dfs[[feature_type]][as.character(out_df$id),]$strand
        diff_to_write[[plot_y_drug]] = out_df
        
      }
      dev.off()
      wb = createWorkbook()
      for(sname in names(diff_to_write)){
        s = createSheet(wb, sname)
        addDataFrame(diff_to_write[[sname]], s)
      }
      saveWorkbook(wb, paste0(out_name, ".xlsx"))
    }
  }
}

