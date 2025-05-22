library(rtracklayer)
library(reshape2)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
gtf_file = args[1]
refbed_file = args[2]

df = import(gtf_file)
df = as.data.frame(df)

df = df[order(df$seqnames, df$start),]

start_end = df$start[df$strand == "-"] <= df$end[df$strand == "-"]

if("gene_id" %in% colnames(df)){
  df$gene_name = df$gene_id
}

if (unique(start_end) & length(unique(start_end)) == 1) {
  
  ft1_1 = df[df$strand != "-", ] %>%
    group_by(seqnames, transcript_id, gene_name, strand) %>%
    summarise(trans_start = min(start), trans_end = max(end), cds_start = min(start), cds_end = max(end))

  ft2_1 = df[df$strand != "-", ] %>%
    group_by(seqnames, transcript_id, gene_name, strand) %>%
    summarise(exon_start = paste(start, collapse = ","),
              exon_end = paste(end, collapse = ","))
  
  ft1_2 = df[df$strand == "-", ] %>%
    arrange(seqnames, desc(start)) %>%
    group_by(seqnames, transcript_id, gene_name, strand) %>%
    summarise(trans_start = max(end), trans_end = min(start), cds_start = max(end), cds_end = min(start))

  ft2_2 = df[df$strand == "-", ] %>%
    arrange(seqnames, desc(start)) %>%
    group_by(seqnames, transcript_id, gene_name, strand) %>%
    summarise(exon_start = paste(end, collapse = ","),
              exon_end = paste(start, collapse = ","))
  
  ft1 = rbind(ft1_1, ft1_2)
  ft2 = rbind(ft2_1, ft2_2)

} else {
  
  ft1 = df %>%
    group_by(seqnames, transcript_id, gene_name, strand) %>%
    summarise(trans_start = min(start), trans_end = max(end), cds_start = min(start), cds_end = max(end))

  ft2 = df %>%
    group_by(seqnames, transcript_id, gene_name, strand) %>%
    summarise(exon_start = paste(start, collapse = ","),
              exon_end = paste(end, collapse = ","))
}

ft = merge(ft1, ft2)

ft$type = "coding"

ft = ft[, c(1, 5:8, 4, 3, 2, 11, 9, 10)]
ft = ft[order(ft$seqnames, ft$trans_start),]

write.table(ft, refbed_file, sep = "\t", quote = F, col.names = F, row.names = F)
