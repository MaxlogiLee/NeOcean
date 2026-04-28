library(foreach)
library(openxlsx)
library(Biostrings)

# Read command-line arguments
args = commandArgs(trailingOnly = TRUE)
ms_file = args[1]
sp_file = args[2]
db_file = args[3]

# Check and create output folder
if (!dir.exists("res")) {
  dir.create("res")
}

# db_grep function to extract data from each sample's grep file
db_grep = function(database, sp) {
  dt = list()
  k = 0
  result_df = NULL  # Initialize result_df to prevent undefined errors
  
  # Ensure result_df is initialized at the start of the function
  for (i in 1:length(sp)) {
    grep_file_path = paste("res/", sp[i], "/", database, "_grep.txt", sep = "")
    
    # Check if file exists
    if (file.exists(grep_file_path)) {
      k = k + 1
      dt[[i]] = read.delim(grep_file_path, sep = "\t", stringsAsFactors = F, header = T)
      
      # Check if the data read is a data frame
      if (!is.data.frame(dt[[i]])) {
        stop(paste("Error: ", grep_file_path, " is not a data frame."))
      }
      
      print(paste("Successfully read:", grep_file_path))  # Output the successfully read file path

      # If it's a 2loc type dataset, modify columns
      is_loc = grepl("2loc", database)
      if (is_loc) {
        dt[[i]] = dt[[i]][, c(1, 3)]  # Select the first two columns
        colnames(dt[[i]])[2] = sub("2loc", "", database)  # Modify column name
      }
      
      # Append the sample name as a suffix to the second column
      colnames(dt[[i]])[2] = paste(colnames(dt[[i]])[2], sp[i], sep = "_")
      
      # Merge results for each sample
      if (k == 1) {
        result_df = dt[[i]]
      } else {
        result_df = merge(result_df, dt[[i]], all = TRUE)  # Merge data
      }
    } else {
      print(paste("File not found:", grep_file_path))  # Output file path to confirm whether the file exists
    }
  }

  # Check if result_df is NULL or empty after merging
  if (is.null(result_df) || nrow(result_df) == 0) {
    stop("Error: No valid data found to merge. Please check the input files and paths.")
  }

  # Ensure result_df is a data frame
  if (!is.data.frame(result_df)) {
    stop("Error: After merging, 'result_df' is not a data frame.")
  }

  # Replace empty strings with NA in each column
  result_df[] = lapply(result_df, function(x) {
    x = as.character(x)
    ifelse(x == "", NA, x)
  })

  # Remove rows where all values (except the first column) are NA
  result_df = result_df[rowSums(is.na(result_df[,-1, drop = F])) != k, ]
  
  # Remove duplicates
  result_df = unique(result_df)

  return(result_df)
}

# Read input files
ms_dt = read.xlsx(ms_file)
sp = read.table(sp_file, sep = "\t", header = F, stringsAsFactors = F)$V1
db = read.table(db_file, sep = "\t", header = F, stringsAsFactors = F)$V1

# Initialize a list to store results for each database
ft = list()
ms = ms_dt

# Process each database
for (j in 1:length(db)) {
  # Call db_grep function to extract data
  grep_res = db_grep(db[j], sp)  # Extract data
  ft[[j]] = merge(ms, grep_res)  # Merge results
  
  # Select columns to keep
  columns_to_keep = intersect(c("Peptide", "Peptide_raw", "Length", "Found.By", 
                                paste("MS_", sp, sep = ""), 
                                paste(sub("2loc", "", db[j]), "_", sp, sep = "")), 
                              colnames(ft[[j]]))
  
  # Select required columns
  ft[[j]] = ft[[j]][, columns_to_keep]
  
  # Replace NA values with empty strings
  ft[[j]][is.na(ft[[j]])] = ""
  
  # Update ms dataframe, remove peptides already merged
  ms = ms[!is.element(ms$Peptide, ft[[j]]$Peptide), ]
  
  # Write merged results to file
  output_file = paste("res/", j, "_", sub("2loc", "", db[j]), "_MS_result_", nrow(ft[[j]]), ".xls", sep = "")
  write.table(ft[[j]], output_file, sep = "\t", row.names = F, quote = F)
}

# Write unmatched results to file
write.xlsx(ms, paste("res/0_", "Unmatched_MS_result_", nrow(ms), ".xlsx", sep = ""), rowNames = F)
write.table(unique(ms$Peptide), paste("res/0_", "Unmatched_MS_result_", length(ms$Peptide), "_peplist.txt", sep = ""), sep = "\t", row.names = F, col.names = F, quote = F)

# Input UniProt peptide data
pep = ms$Peptide
Seq = list()

# Process each peptide for slicing
for (i in 1:length(pep)) {
  pep_str = as.character(pep[i])  # Ensure pep[i] is character type
  
  Seq[[i]] = foreach(j = 1:(nchar(pep_str) - 1), .combine = rbind) %do% c(pep_str, nchar(pep_str), 1, j, subseq(pep_str, start = 1, end = j), j + 1, nchar(pep_str), subseq(pep_str, start = j + 1, end = nchar(pep_str)))
}

# Combine all peptide sequences into a dataframe
all_seq = do.call(rbind, Seq)
all_seq = as.data.frame(all_seq)
colnames(all_seq) = c("Peptide", "Length", "Sequence1_start", "Sequence1_end", "Sequence1", "Sequence2_start", "Sequence2_end", "Sequence2")

# Write sliced peptide sequences to file
write.table(unique(all_seq), paste("res", "/", "Uniport_cut_peptide.txt", sep = ""), sep = "\t", row.names = F, quote = F)

# Get unique sequences and write to file
write.table(union(all_seq$Sequence1, all_seq$Sequence2), paste("res", "/", "Uniport_cut_peptidelist.txt", sep = ""), sep = "\t", row.names = F, col.names = F, quote = F)
