norm_counts = read.table("Z:/RNA_seq_data/Breast/MCF7_Pfizer/DESeq2/DESeq2_normalized_counts_v21.txt")
colnames(norm_counts) = sub(pattern = "BZAE2", replacement = "e2bza", x = colnames(norm_counts)) %>%
  sub("BZAGC10", "gc10bza", x = .) %>%
  tolower() %>% sub("mcf7", "MCF7", x = .) 


tmp = unlist(DE_res) %>% unique(); names(tmp) = NULL
DE_norm_counts = as.matrix(norm_counts[tmp,])

keep = rowSums(DE_norm_counts > MIN_COUNTS) >= 3
DE_norm_counts = DE_norm_counts[keep,]

save(DE_norm_counts, norm_counts, file = "input/norm_counts.save")