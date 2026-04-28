#!/usr/bin/env Rscript
# ==============================================================================
# teprof2_refbed.R
# Generate TEProf2-specific refBed from candidates.getorf.fa and stringtie.refbed.
# This bridges TEProf2 ORF output with the coordinate mapping step in peptide screening.
# ==============================================================================

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  cat("Usage: Rscript teprof2_refbed.R <sample_name> <step10_csv> <stringtie_refbed> <getorf_fa> <output_refbed>\n")
  quit(status = 1)
}

sp1 = args[1]
step10_csv = args[2]
stringtie_refbed = args[3]
getorf_fa = args[4]
output_refbed = args[5]

# Read Step 10 statistics
dt = read.csv(step10_csv, stringsAsFactors = FALSE)
dt$id = paste(dt$Subfam, dt$Start.TE, dt$Gene, dt$Location.TE, dt$Gene, dt$Splice.Target, sep = "_")

# Handle intergenic entries
gl1 = grep("Intergenic", dt$Location.TE)
gl2 = grep("Intergenic", dt$Splice.Target)
gl_1 = intersect(gl1, gl2)
gl_2 = setdiff(gl1, gl2)
gl_3 = setdiff(gl2, gl1)

if (length(gl_1) != 0) {
  dt$id[gl_1] = paste(dt$Subfam[gl_1], dt$Start.TE[gl_1], "None", "None_None", "None", "None_None", sep = "_")
}
if (length(gl_2) != 0) {
  dt$id[gl_2] = paste(dt$Subfam[gl_2], dt$Start.TE[gl_2], "None", "None_None", dt$Gene[gl_2], dt$Splice.Target[gl_2], sep = "_")
}
if (length(gl_3) != 0) {
  dt$id[gl_3] = paste(dt$Subfam[gl_3], dt$Start.TE[gl_3], dt$Gene[gl_3], dt$Location.TE[gl_3], "None", "None_None", sep = "_")
}

# Read refBed
df = read.csv(stringtie_refbed, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Read getorf FASTA headers
fa = readLines(getorf_fa)
fa = sub(">", "", fa[grep(">", fa)])
fa_split = strsplit(fa, split = "_")
fa_len = sapply(fa_split, length)

# Handle variable-length IDs from getorf
spl1 = which(fa_len == 11)
for (j in 1:length(spl1)) {
  fa_split[[spl1[j]]] = c(paste(fa_split[[spl1[j]]][1:2], collapse = "_"), fa_split[[spl1[j]]][-c(1:2)])
}

spl2 = which(fa_len == 12)
for (j in 1:length(spl2)) {
  fa_split[[spl2[j]]] = c(paste(fa_split[[spl2[j]]][1:3], collapse = "_"), fa_split[[spl2[j]]][-c(1:3)])
}

fa_split = do.call(rbind, fa_split)
fa1 = apply(fa_split[, 1:8], 1, paste, collapse = "_")
fa2 = apply(fa_split[, 1:9], 1, paste, collapse = "_")

fa_dt = data.frame(id = fa1, ID = fa2)
fa_dt = unique(fa_dt)

# Merge with TEProf2 statistics
fa_dt_id = merge(fa_dt, dt[, c("Transcript.Name", "id")])
fa_dt_id_df = merge(df, fa_dt_id, by.x = "V8", by.y = "Transcript.Name")

# Select and reorder columns
TEpro2 = fa_dt_id_df[, c(2:8, 1, 9:11, 13)]

# Write output
write.table(TEpro2, output_refbed, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
cat("[teprof2_refbed] TEProf2 refBed written to:", output_refbed, "\n")
