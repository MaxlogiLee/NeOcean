#!/usr/bin/env Rscript
# ==============================================================================
# ms_results.R
# Process mass spectrometry results for NeOcean peptide screening.
# Supports optional I/L transformation (Isoleucine <-> Leucine ambiguity).
# ==============================================================================

library(openxlsx)
library(foreach)
library(dplyr)
library(reshape2)
library(stringr)

# ------------------------------------------------------------------------------
# Parse command-line arguments
# ------------------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  cat("Usage: Rscript ms_results.R <sample_name> <ms_xlsx> <output_dir> <do_il_transform>\n")
  cat("  sample_name      : Sample identifier (e.g., MCJ)\n")
  cat("  ms_xlsx          : Path to MS result Excel file\n")
  cat("  output_dir       : Output directory (will create res/peptide/ under it)\n")
  cat("  do_il_transform  : TRUE/FALSE - whether to perform I/L transformation\n")
  quit(status = 1)
}

sample_name = args[1]
ms_xlsx = args[2]
output_dir = args[3]
do_il_transform = toupper(args[4]) == "TRUE"

# ------------------------------------------------------------------------------
# Validate inputs
# ------------------------------------------------------------------------------
if (!file.exists(ms_xlsx)) {
  stop(paste("MS file not found:", ms_xlsx))
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

peptide_dir = file.path(output_dir, "res", "peptide")
if (!dir.exists(peptide_dir)) {
  dir.create(peptide_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Read MS data
# ------------------------------------------------------------------------------
fl = basename(ms_xlsx)
dt = read.xlsx(ms_xlsx)

# Keep required columns
dt = dt[, c("Peptide", "Length", "Sample")]
dt$Sample = sample_name

# ------------------------------------------------------------------------------
# I/L transformation functions
# ------------------------------------------------------------------------------
I2L = function(pep) {
  outres = list()
  for (i in 1:length(pep)) {
    loc = str_locate_all(pep[i], "I")[[1]][, 1]
    IL2_new = rep(pep[i], length(loc))
    for (k in 1:length(loc)) {
      str_sub(IL2_new[k], loc[k], loc[k]) = "L"
    }
    outres[[i]] = IL2_new
  }
  return(unique(do.call(c, outres)))
}

L2I = function(pep) {
  outres = list()
  for (i in 1:length(pep)) {
    loc = str_locate_all(pep[i], "L")[[1]][, 1]
    IL2_new = rep(pep[i], length(loc))
    for (k in 1:length(loc)) {
      str_sub(IL2_new[k], loc[k], loc[k]) = "I"
    }
    outres[[i]] = IL2_new
  }
  return(unique(do.call(c, outres)))
}

# ------------------------------------------------------------------------------
# Process peptides
# ------------------------------------------------------------------------------
if (do_il_transform) {
  # I2L: peptides containing I
  IL2 = dt[grep("I", dt$Peptide), ]
  if (nrow(IL2) > 0) {
    IL2_pep = list()
    for (j in 1:nrow(IL2)) {
      pep = IL2$Peptide[j]
      loc = str_locate_all(pep, "I")[[1]][, 1]
      res = list()
      for (m in 1:length(loc)) {
        res[[m]] = I2L(pep)
        pep = res[[m]]
      }
      IL2_pep[[j]] = data.frame(Peptide = IL2$Peptide[j], ItoL = unlist(res))
    }
    IL2_replace = do.call(rbind, IL2_pep)
    IL2 = merge(IL2_replace, IL2, by = "Peptide")
  } else {
    IL2 = data.frame(Peptide = character(), ItoL = character(), Length = integer(), Sample = character())
  }

  # L2I: peptides containing L
  LI2 = dt[grep("L", dt$Peptide), ]
  if (nrow(LI2) > 0) {
    LI2_pep = list()
    for (j in 1:nrow(LI2)) {
      pep = LI2$Peptide[j]
      loc = str_locate_all(pep, "L")[[1]][, 1]
      res = list()
      for (m in 1:length(loc)) {
        res[[m]] = L2I(pep)
        pep = res[[m]]
      }
      LI2_pep[[j]] = data.frame(Peptide = LI2$Peptide[j], ItoL = unlist(res))
    }
    LI2_replace = do.call(rbind, LI2_pep)
    LI2 = merge(LI2_replace, LI2, by = "Peptide")
  } else {
    LI2 = data.frame(Peptide = character(), ItoL = character(), Length = integer(), Sample = character())
  }

  # Peptides without I or L
  IL1 = dt[!grepl("I|L", dt$Peptide), ]
  IL1$ItoL = IL1$Peptide

  final_data = rbind(
    IL1[, c("Peptide", "ItoL", "Length", "Sample")],
    IL2[, c("Peptide", "ItoL", "Length", "Sample")],
    LI2[, c("Peptide", "ItoL", "Length", "Sample")]
  )
} else {
  # No I/L transformation
  final_data = dt
  final_data$ItoL = final_data$Peptide
}

# ------------------------------------------------------------------------------
# Save peptide list
# ------------------------------------------------------------------------------
peptide_list_file = file.path(peptide_dir, paste0(sample_name, "_peptide_list.txt"))
write.table(
  unique(final_data$ItoL),
  peptide_list_file,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------------------------
# Save MS summary table
# ------------------------------------------------------------------------------
ms_result_file = file.path(output_dir, "res", paste0(sample_name, "_MS_result.xlsx"))
mt = final_data
mt$Sample = factor(mt$Sample, levels = unique(mt$Sample))
df = dcast(
  mt,
  as.formula(paste(paste0(colnames(mt)[c(1, 2, 3)], collapse = "+"), "~", "Sample", sep = "")),
  value.var = "Sample",
  fun.aggregate = length
)
colnames(df)[4:ncol(df)] = paste("MS_", colnames(df)[4:ncol(df)], sep = "")
colnames(df)[1:2] = c("Peptide_raw", "Peptide")
write.xlsx(unique(df), ms_result_file, rowNames = FALSE)

# ------------------------------------------------------------------------------
# Save sample and database lists (used by Integrate_result.R)
# ------------------------------------------------------------------------------
sp_file = file.path(output_dir, "res", "sp_file.txt")
db_file = file.path(output_dir, "res", "database_file.txt")

sp = sample_name
db = c("BN", "NCBI", "Uniport", "Mutation", "Fusion", "IR2loc", "TElocal2loc", "TEprof22loc", "DN2loc")

write.table(sp, sp_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(db, db_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("[ms_results] Completed for sample:", sample_name, "\n")
cat("  Peptide list:", peptide_list_file, "\n")
cat("  MS result:", ms_result_file, "\n")
cat("  sp_file:", sp_file, "\n")
cat("  db_file:", db_file, "\n")
