# 加载所需包
library(openxlsx)
library(dplyr)
library(reshape2)

# Setting the file name and sample name
fl <- "sample_HLA_I.xlsx"
nm <- "sample" 
dt <- read.xlsx(file.path("data/MSdata", fl))  # Reading Mass Spectrometry Files

dt$Sample <- nm

dt <- dt[, c("Peptide", "Length", "Sample")]
dt$ItoL <- dt$Peptide  

write.table(
  unique(dt$ItoL),
  paste0("res/peptide/", nm, "_peptide_list.txt"),
  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)
mt <- dt
mt$Sample <- factor(mt$Sample, levels = unique(mt$Sample))
df <- dcast(mt, Peptide + Length + Sample ~ Sample, value.var = "Sample", fun.aggregate = length)

colnames(df)[4] <- paste0("MS_", nm)
colnames(df)[1:3] <- c("Peptide_raw", "Length", "Sample")
df$Peptide <- df$Peptide_raw  
df <- df[, c("Peptide_raw", "Peptide", "Length", paste0("MS_", nm))]

write.xlsx(unique(df), paste0("res/", nm, "_MS_result.xlsx"))

# 输出样本和数据库列表
sp <- nm
db <- c("Uniport", "Mutation", "IR", "TElocal", "TEprof2", "DN", "Fusion", "NCBI", "BN")

write.table(sp, "res/sp_file.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(db, "res/database_file.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
