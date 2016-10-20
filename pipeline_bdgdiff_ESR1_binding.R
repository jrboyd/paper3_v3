source('setup.R')
source("module_bdgdiff_STEIN_ChIPseq.R")
source("module_ESR1_and_enhancers.R")
source("functions_rGREAT.R")

clean = list()
for(i in 1:length(all_bdgdiff_tables)){
  a = names(all_bdgdiff_tables)[i]
  group = all_bdgdiff_tables[[i]]
  for(j in 1:length(group)){
    b = names(group)[j]
    pair = group[[j]]
    for(k in 1:3){
      c = names(pair)[k]
      gr = pair[[k]]
      clean[[paste(a,b,c)]] = gr
    }
  }
}

# lens = lapply(all_bdgdiff_tables, function(group){
#   t(sapply(group, function(pair){
#     sapply(pair[1:3], length)
#   }))
# })
lens = sapply(clean, length)


# olaps = lapply(all_bdgdiff_tables, function(group){
#   t(sapply(group, function(pair){
#     sapply(pair[1:3], function(x){
#       length(findOverlaps(subject = x, 
#                    query = GSE40129$MCF7_e2_ESR1))
#       
#     })
#   }))
# })
olaps = sapply(clean, function(x){
  olap = findOverlaps(subject = x,
                      query = GSE40129$MCF7_e2_ESR1)
  return(length(unique(subjectHits(olap))))
  
})
group_names = sapply(strsplit(names(lens), " "), function(x)x[1])
cells = sapply(strsplit(group_names, "_"), function(x)x[1])
marks = sapply(strsplit(group_names, "_"), function(x)x[2])
group_names = sub("_", "\n", group_names)
pair_names = sapply(strsplit(names(lens), " "), function(x)x[2])
direction_names = sapply(strsplit(names(lens), " "), function(x)x[3])
library(ggplot2)
df = data.frame(lens = lens, olaps = olaps, cell = cells, mark = marks, group = group_names, pair = pair_names, direction = direction_names)
gg = ggplot() +
  geom_jitter(data = df[grepl("from_ctrl_", pair_names), ],
              aes(x = group, 
                  y = olaps / lens, 
                  col = direction, 
                  shape = pair),
              width = .3, 
              height = 0)

print(gg)

olaps = sapply(clean, function(x){
  olap = findOverlaps(subject = x,
                      query = GSE40129$MCF7_e2_ESR1)
  return(length(unique(subjectHits(olap))))
  
})
group_names = sapply(strsplit(names(lens), " "), function(x)x[1])
cells = sapply(strsplit(group_names, "_"), function(x)x[1])
marks = sapply(strsplit(group_names, "_"), function(x)x[2])
group_names = sub("_", "\n", group_names)
pair_names = sapply(strsplit(names(lens), " "), function(x)x[2])
direction_names = sapply(strsplit(names(lens), " "), function(x)x[3])
library(ggplot2)
df = data.frame(lens = lens, olaps = olaps, cell = cells, mark = marks, group = group_names, pair = pair_names, direction = direction_names)
gg = ggplot() +
  geom_jitter(data = df[grepl("from_ctrl_", pair_names) & (grepl("_e2", pair_names) | grepl("_bza", pair_names)), ],
              aes(x = group, 
                  y = olaps / lens, 
                  col = direction, 
                  shape = pair),
              width = .3, 
              height = 0)

print(gg)

gg = ggplot() +
  geom_jitter(data = df[grepl("from_ctrl_", pair_names) & (grepl("_e2$", pair_names)), ],
              aes(x = group, 
                  y = olaps / lens, 
                  col = direction, 
                  shape = pair),
              width = .3, 
              height = 0)

print(gg)

gg = ggplot() +
  geom_jitter(data = df[grepl("from_ctrl_", pair_names) & (grepl("_e2bza$", pair_names)), ],
              aes(x = group, 
                  y = olaps / lens, 
                  col = direction, 
                  shape = pair),
              width = .3, 
              height = 0)

print(gg)
