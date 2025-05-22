library(foreach)

args = commandArgs(trailingOnly = TRUE)
input_file = args[1]  
refbed_file = args[2]  
output_file = args[3]  
merge_loc = args[4] == "TRUE"  

dt = read.table(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
refBed = read.table(refbed_file, sep = "\t", header = FALSE)

colnames(refBed) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11")

item = colnames(dt)[2]
if (item == "DN") {
  id = "V8"
} else if (item == "TEprof2") {
  id = "V12"
}

gl = which(dt[, item] != "")
peptide_loc = c()

for (j in 1:length(gl)) {
  loc = as.data.frame(do.call(rbind, strsplit(strsplit(dt[gl[j], item], split = ";")[[1]], split = "\\[|\\]")))
  loc = cbind(loc, do.call(rbind, strsplit(loc$V2, split = " - ")), do.call(rbind, strsplit(loc$V4, split = "-")))
  colnames(loc)[5:8] = c("N_pos1", "N_pos2", "TC_pos1", "TC_pos2")
  
  loc$pos1 = foreach(k = 1:nrow(loc), .combine = c) %do% (as.numeric(loc$N_pos1)[k]:as.numeric(loc$N_pos2)[k])[as.numeric(loc$TC_pos1)[k] * 3 - 2]
  loc$pos2 = foreach(k = 1:nrow(loc), .combine = c) %do% (as.numeric(loc$N_pos1)[k]:as.numeric(loc$N_pos2)[k])[as.numeric(loc$TC_pos2)[k] * 3]
  
  loc$tcons = sub("_[0-9]*.$", "", loc$V1)
  
  pep_loc = c()
  for (m in 1:nrow(loc)) {
    is_reverse = grepl("REVERSE SENSE", loc$V3[m])
    
    bed = refBed[refBed[, id] == loc$tcons[m], ]
    exon_start = as.numeric(unlist(strsplit(bed$V10, split = ",")))
    exon_end = as.numeric(unlist(strsplit(bed$V11, split = ",")))
    
    exon_length = cumsum(abs(exon_start - exon_end) + 1)
    bed$V12 = paste0(exon_length, collapse = ",")
    
    pep_start = loc$pos1[m]
    pep_end = loc$pos2[m]
    
    pep_start_exon = which(pep_start - exon_length <= 0)[1]
    pep_end_exon = which(pep_end - exon_length <= 0)[1]
    
    loc_vector = foreach(k = 1:length(exon_start), .combine = c) %do% exon_start[k]:exon_end[k]
    
    if (pep_start_exon == pep_end_exon) {
      pep_pos = paste(item, "-", bed$V1, ":", loc_vector[pep_start], "-", loc_vector[pep_end], "(", bed$V6, ")", sep = "")
    } else {
      if (is_reverse) {
        pep_pos = paste(item, "-", paste(bed$V1, ":", loc_vector[pep_start], "-", exon_start[pep_start_exon], "(", bed$V6, ")", sep = ""), ",",
                        paste(bed$V1, ":", loc_vector[exon_length[pep_end_exon]], "-", loc_vector[pep_end], "(", bed$V6, ")", sep = ""), sep = "")
      } else {
        pep_pos = paste(item, "-", paste(bed$V1, ":", loc_vector[pep_start], "-", loc_vector[exon_length[pep_start_exon]], "(", bed$V6, ")", sep = ""), ",",
                        paste(bed$V1, ":", exon_start[pep_end_exon], "-", loc_vector[pep_end], "(", bed$V6, ")", sep = ""), sep = "")
      }
    }
    pep_loc[m] = pep_pos
  }
  
  if (merge_loc) {
    peptide_loc[j] = paste(unique(pep_loc), collapse = ";")
  } else {
    peptide_loc[j] = paste(pep_loc, collapse = ";")
  }
}

dt$loc = ""
dt$loc[gl] = peptide_loc

write.table(dt, output_file, row.names = FALSE, quote = FALSE, sep = '\t')
