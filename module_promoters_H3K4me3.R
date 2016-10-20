setwd(input_dir)
load("ref/enst_ref.save")
enst_ref$transcript_id = rownames(enst_ref)

files = dir("STEIN_ChIPseq_nomodel/", full.names = T)
STEIN_peaks = list()
for(f in files){
  name = sub("_peaks.narrowPeak", "", basename(f))
  print(name)
  gr = file2granges(f)
  STEIN_peaks[[name]] = gr
}

##get a random set of peaks
# x = STEIN_peaks$MCF7_e2_H3K4ME3
# n = 20
# o = order(runif(length(x)))[1:n]
# mat = as.data.frame(x[o])
# cbind(apply(mat, 1, function(x)paste(x[1:3], collapse = " ")))
ext = 1000
start = ifelse(enst_ref$strand == "+", enst_ref$start - ext, enst_ref$end - ext)
end = ifelse(enst_ref$strand == "+", enst_ref$start + ext, enst_ref$end + ext)
prom_ref = enst_ref
prom_ref$start = start
prom_ref$end = end

prom_gr = GRanges(seqnames = prom_ref$chrm, 
                  ranges = IRanges(start = prom_ref$start, end = prom_ref$end), 
                  strand = prom_ref$strand,
                  transcript_id = rownames(prom_ref),
                  gene_id = prom_ref$gene_id,
                  gene_name = prom_ref$gene_name, 
                  tag = prom_ref$tag,
                  transcript_support_level = prom_ref$transcript_support_level)
keep = prom_gr$tag == "basic" | prom_gr$tag == "alternative_5_UTR"
prom_gr = prom_gr[keep]
names(prom_gr) = prom_gr$transcript_id

is_k4me3 = grepl("H3K4ME3", names(STEIN_peaks))
STEIN_k4me3_proms = lapply(STEIN_peaks[is_k4me3], function(peaks){
  # peaks = STEIN_peaks$MCF7_e2_H3K4ME3
  
  olaps = findOverlaps(subject = peaks, query = prom_gr)
  hit_gr = prom_gr[unique(queryHits(olaps))]
  dupe_genes = unique(hit_gr$gene_name[duplicated(hit_gr$gene_name)])
  uniq_genes = setdiff(hit_gr$gene_name, dupe_genes)
  all_genes = c(dupe_genes, uniq_genes)
  all_tscripts = character(length = length(all_genes))
  names(all_tscripts) = all_genes
  for(i in 1:length(dupe_genes)){
    if(i %% 250 == 0) print(paste(i, "/", length(dupe_genes)))
    g = dupe_genes[i]
    df = data.frame(gene_name = g)
    m = merge(x = df, y = hit_gr, by = "gene_name")
    dupe_counts = sapply(unique(m$start), function(x)sum(m$start == x))
    if(max(dupe_counts) > 1){
      o = order(dupe_counts, decreasing = T)
      tscript = m$transcript_id[o][1]
    }else{
      o = order(m$transcript_support_level, decreasing = F, na.last = T)
      tscript = m$transcript_id[o][1]
    }
    # print(m)
    # print(tscript)
    all_tscripts[g] = tscript
  }
  for(i in 1:length(uniq_genes)){
    if(i %% 250 == 0) print(paste(i, "/", length(uniq_genes)))
    g = uniq_genes[i]
    df = data.frame(gene_name = g)
    m = merge(x = df, y = hit_gr, by = "gene_name")
    tscript = m$transcript_id
    all_tscripts[g] = tscript
  }
  
  names(hit_gr) = hit_gr$transcript_id
  hit_gr[all_tscripts]
})

STEIN_k4me3_genes = lapply(STEIN_k4me3_proms, function(x){
  x$gene_name
})

STEIN_tscripts = lapply(STEIN_k4me3_proms, function(x){
  x$transcript_id
})

all_genes = unique(unlist(STEIN_k4me3_genes))
all_mat = matrix(F, nrow = length(all_genes), ncol = length(STEIN_k4me3_genes))
colnames(all_mat) = names(STEIN_k4me3_genes)
rownames(all_mat) = all_genes
for(i in colnames(all_mat)){
  all_mat[STEIN_k4me3_genes[[i]], i] = T
}



all_tscripts = unique(unlist(STEIN_tscripts))
all_mat = matrix(F, nrow = length(all_tscripts), ncol = length(STEIN_tscripts))
colnames(all_mat) = names(STEIN_tscripts)
rownames(all_mat) = all_tscripts
for(i in colnames(all_mat)){
  all_mat[STEIN_tscripts[[i]], i] = T
}

# is_shared = rownames(all_mat)[rowSums(all_mat) == 12]
tscripts_gr = lapply(STEIN_tscripts, function(x){
  prom_gr[x]
})

best_ids = sapply(all_genes, function(gene){#pick one transcript to represent each gene
  print(gene)
  gene = data.frame(gene_name = gene)
  # gene = data.frame(gene_name = all_genes[i])
  matching = unlist(sapply(tscripts_gr, function(x){
    m = merge(x, gene, by = "gene_name")
    return(m$transcript_id)
  }))
  matching = enst_ref[matching,]
  o = order(matching$transcript_support_level, decreasing = F, na.last = T)
  matching = matching[o,]
  tscript_id = matching$transcript_id[1]
  return(tscript_id)
})

id2genes = names(best_ids)
names(id2genes) = best_ids

STEIN_k4me3_tscripts = lapply(STEIN_k4me3_genes, function(x){
  best_ids[x]
})

best_k4me3_tscripts = unique(unlist(STEIN_k4me3_tscripts))
k4me3_ref = prom_ref[best_k4me3_tscripts,]
k4me3_ref[,"name"] = paste(k4me3_ref$transcript_id, k4me3_ref$gene_id, k4me3_ref$gene_name, sep = ":")
k4me3_ref[,"score"] = rep(0, nrow(k4me3_ref))
k4me3_bed = as.matrix(k4me3_ref[,c("chrm", "start", "end", "name", "score", "strand" )])
k4me3_bed = gsub(" ", "", k4me3_bed)
write.table(k4me3_bed, file = "STEIN_k4me3_promoters_4ngsplot.bed", row.names = F, col.names = F, quote = F, sep = "\t")
