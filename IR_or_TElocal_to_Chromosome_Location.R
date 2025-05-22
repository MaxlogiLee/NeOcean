library(foreach)

args = commandArgs(trailingOnly = TRUE)
input_file = args[1]
output_file = args[2]
merge_loc = args[3] == "TRUE"

dt = read.table(input_file, sep = "\t", header = T, stringsAsFactors = F)
item = colnames(dt)[2]

if (item == "IR") {
  dt$IR = gsub("::", "chr", dt$IR)
  dt$IR = gsub("\\(REVERSE SENSE\\) ", "", dt$IR)
}

gl = which(dt[, item] != "")
peptide_loc = c()

for (j in 1:length(gl)) {
  loc = as.data.frame(do.call(rbind, strsplit(strsplit(dt[gl[j], item], split = ";")[[1]], split = "\\[|\\]")))
  loc = cbind(loc, do.call(rbind, strsplit(loc$V2, split = " - ")), do.call(rbind, strsplit(loc$V4, split = "-")))
  colnames(loc)[5:8] = c("N_pos1", "N_pos2", "TC_pos1", "TC_pos2")
  loc$pos1 = foreach(k = 1:nrow(loc), .combine = c) %do% (as.numeric(loc$N_pos1)[k]:as.numeric(loc$N_pos2)[k])[as.numeric(loc$TC_pos1)[k] * 3 - 2]
  loc$pos2 = foreach(k = 1:nrow(loc), .combine = c) %do% (as.numeric(loc$N_pos1)[k]:as.numeric(loc$N_pos2)[k])[as.numeric(loc$TC_pos2)[k] * 3]
  loc$id = sub("_[0-9]*.$", "", loc$V1)

  pep_loc = c()
  for (m in 1:nrow(loc)) {
    if (item == "IR") {
      rg_chr = gsub(".*-(chr[0-9XY]+).*", "\\1", loc$V3[m])
      rg_range = gsub(".*(chr[0-9XY]+:[0-9]+-[0-9]+).*", "\\1", loc$V3[m])
      
      rg_start_end = strsplit(rg_range, ":|\\-")[[1]]
      rg_start = as.numeric(rg_start_end[2])
      rg_end = as.numeric(rg_start_end[3])
      
      if (is.na(rg_start) || is.na(rg_end)) {
        stop(paste("Invalid coordinates in:", loc$V3[m]))
      }
      
      rg_start = rg_start + 1
      rg_end = rg_end + 1
    }
    
    pep_start = loc$pos1[m]
    pep_end = loc$pos2[m]

    loc_vector = rg_start:rg_end
    pep_loc[m] = paste(item, "-", rg_chr, ":", loc_vector[pep_start], "-", loc_vector[pep_end], sep = "")
  }

  if (merge_loc) {
    peptide_loc[j] = paste(unique(pep_loc), collapse = ";")
  } else {
    peptide_loc[j] = paste(pep_loc, collapse = ";")
  }
}

dt$loc = ""
dt$loc[gl] = peptide_loc

write.table(dt, output_file, row.names = F, quote = F, sep = '\t')
