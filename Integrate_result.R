library(foreach)
library(openxlsx)
library(Biostrings)

args = commandArgs(trailingOnly = TRUE)
ms_file = args[1]
sp_file = args[2]
db_file = args[3]

if (!dir.exists("res")) {
  dir.create("res")
}

db_grep = function(database, sp) {
  dt = list()
  k = 0
  result_df = NULL  
  
  for (i in 1:length(sp)) {
    grep_file_path = paste("res/", sp[i], "/", database, "_grep.txt", sep = "")
    
    if (file.exists(grep_file_path)) {
      k = k + 1
      dt[[i]] = read.delim(grep_file_path, sep = "\t", stringsAsFactors = F, header = T)
      
      if (!is.data.frame(dt[[i]])) {
        stop(paste("Error: ", grep_file_path, " is not a data frame."))
      }
      
      print(paste("Successfully read:", grep_file_path))  

      is_loc = grepl("2loc", database)
      if (is_loc) {
        dt[[i]] = dt[[i]][, c(1, 3)]  
        colnames(dt[[i]])[2] = sub("2loc", "", database)  
      }
      
      colnames(dt[[i]])[2] = paste(colnames(dt[[i]])[2], sp[i], sep = "_")
      
      if (k == 1) {
        result_df = dt[[i]]
      } else {
        result_df = merge(result_df, dt[[i]], all = TRUE)  
      }
    } else {
      print(paste("File not found:", grep_file_path))  
    }
  }

  if (is.null(result_df) || nrow(result_df) == 0) {
    stop("Error: No valid data found to merge. Please check the input files and paths.")
  }

  if (!is.data.frame(result_df)) {
    stop("Error: After merging, 'result_df' is not a data frame.")
  }

  result_df[] = lapply(result_df, function(x) {
    x = as.character(x)
    ifelse(x == "", NA, x)
  })

  result_df = result_df[rowSums(is.na(result_df[,-1, drop = F])) != k, ]
  
  result_df = unique(result_df)

  return(result_df)
}


  if (is.null(result_df) || nrow(result_df) == 0) {
    stop("Error: No valid data was found to merge. Please check the input files.")
  }

  if (!is.data.frame(result_df)) {
    stop("Error: After merging, 'result_df' is not a data frame.")
  }

  result_df[] = lapply(result_df, function(x) {
    x = as.character(x)
    ifelse(x == "", NA, x)
  })

  result_df = result_df[rowSums(is.na(result_df[,-1, drop = F])) != k, ]
  
  result_df = unique(result_df)

  return(result_df)
}

  if (nrow(result_df) == 0) {
    stop("Error: No valid data was found to merge. Please check the input files.")
  }

  if (!is.data.frame(result_df)) {
    stop("Error: After merging, 'result_df' is not a data frame.")
  }

  result_df[] = lapply(result_df, function(x) {
    x = as.character(x)
    ifelse(x == "", NA, x)
  })

  result_df = result_df[rowSums(is.na(result_df[,-1, drop = F])) != k, ]
  
  result_df = unique(result_df)

  return(result_df)
}

  if (is.null(result_df)) {
    stop("Error: No valid data was found to merge.")
  }

  if (!is.data.frame(result_df)) {
    stop("Error: After merging, 'result_df' is not a data frame.")
  }

  result_df[] = lapply(result_df, function(x) {
    x = as.character(x)
    ifelse(x == "", NA, x)
  })

  result_df = result_df[rowSums(is.na(result_df[,-1, drop = F])) != k, ]
  result_df = unique(result_df)

  return(result_df)
}
  
  if (!is.data.frame(result_df)) {
    stop("Error: The variable 'result_df' is not a data frame after merging.")
  }
  
  result_df[] = lapply(result_df, function(x) {
    x = as.character(x)
    ifelse(x == "", NA, x)
  })
  
  result_df = result_df[rowSums(is.na(result_df[,-1, drop = F])) != k, ]
  
  result_df = unique(result_df)
  
  return(result_df)
}
  
  if (!is.data.frame(result_df)) {
    stop("The variable 'result_df' is not a data frame.")
  }

  result_df[] = lapply(result_df, function(x) {
    x = as.character(x)
    ifelse(x == "", NA, x)
  })
  
  result_df = result_df[rowSums(is.na(result_df[,-1, drop = F])) != k, ]
  
  result_df = unique(result_df)
  
  return(result_df)
}

ms_dt = read.xlsx(ms_file)
sp = read.table(sp_file, sep = "\t", header = F, stringsAsFactors = F)$V1
db = read.table(db_file, sep = "\t", header = F, stringsAsFactors = F)$V1

ft = list()
ms = ms_dt

for (j in 1:length(db)) {
  grep_res = db_grep(db[j], sp)  
  ft[[j]] = merge(ms, grep_res)  
  
  columns_to_keep = intersect(c("Peptide", "Peptide_raw", "Length", "Found.By", 
                                paste("MS_", sp, sep = ""), 
                                paste(sub("2loc", "", db[j]), "_", sp, sep = "")), 
                              colnames(ft[[j]]))
  
  ft[[j]] = ft[[j]][, columns_to_keep]
  
  ft[[j]][is.na(ft[[j]])] = ""
  
  ms = ms[!is.element(ms$Peptide, ft[[j]]$Peptide), ]
  
  output_file = paste("res/", j, "_", sub("2loc", "", db[j]), "_MS_", nrow(ft[[j]]), ".xls", sep = "")
  write.table(ft[[j]], output_file, sep = "\t", row.names = F, quote = F)
}

write.xlsx(ms, paste("res/0_", "ungreped", "_MS_", nrow(ms), ".xlsx", sep = ""), rowNames = F)
write.table(unique(ms$Peptide), paste("res/0_", "ungreped", "_MS_", length(ms$Peptide), "_peplist.txt", sep = ""), sep = "\t", row.names = F, col.names = F, quote = F)

pep = ms$Peptide
Seq = list()

for (i in 1:length(pep)) {
  Seq[[i]] = foreach(j = 1:(nchar(pep[i]) - 1), .combine = rbind) %do% c(pep[i], nchar(pep[i]), 1, j, subseq(pep[i], start = 1, end = j), j + 1, nchar(pep[i]), subseq(pep[i], start = j + 1, end = nchar(pep[i])))
}

all_seq = do.call(rbind, Seq)
all_seq = as.data.frame(all_seq)
colnames(all_seq) = c("Peptide", "Length", "Sequence1_start", "Sequence1_end", "Sequence1", "Sequence2_start", "Sequence2_end", "Sequence2")

write.table(unique(all_seq), paste("res", "/", "Uniport_cut_peptide.txt", sep = ""), sep = "\t", row.names = F, quote = F)

write.table(union(all_seq$Sequence1, all_seq$Sequence2), paste("res", "/", "Uniport_cut_peptidelist.txt", sep = ""), sep = "\t", row.names = F, col.names = F, quote = F)
