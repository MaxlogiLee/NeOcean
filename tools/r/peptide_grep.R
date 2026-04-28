# 加载必要的库
library(foreach)
library(stringr)
library(openxlsx)

# 函数定义，grep_peptide1 用于查找肽段在数据库中的位置
grep_peptide1 = function(database_file, sequence_merge = FALSE, peptide, database, split_str = NULL) {
  # 如果不需要合并FA序列（即序列已经合并好）
  if (!sequence_merge) {
    fa = readLines(database_file)
  } else {
    # 如果需要合并FA序列（此处可以保持原样，假设在你的文件中已合并）
    fa <- system(paste("awk '{if($0 ~ />/) name = $0; else seq[name] = seq[name] $0;} END {for(i in seq) print i\"\\n\"seq[i]}'", database_file), intern = TRUE)
  }
  
  # 提取序列（去掉头部的行）
  db = fa[-grep(">", fa)]
  names(db) = fa[grep(">", fa)]
  
  # grep并提取header
  gl_res = foreach(j = 1:length(peptide)) %do% system(paste("LC_ALL=C grep -B1 ", peptide[j], " ", database_file, "|grep '^>'", sep = ""), intern = T)
  
  # 定义匹配到的位置
  is.match = foreach(j = 1:length(peptide), .combine = c) %do% length(gl_res[[j]])
  match = which(is.match != 0)
  
  match_res = rep("", length(peptide))
  
  for (j in match) {
    gl_loc = foreach(k = 1:length(gl_res[[j]]), .combine = c) %do% paste(
      paste(sub("^>", "", gl_res[[j]][k]), "[", paste(str_locate_all(db[gl_res[[j]][k]], peptide[j])[[1]][, 1], 
      str_locate_all(db[gl_res[[j]][k]], peptide[j])[[1]][, 2], sep = "-"), "]", sep = ""), collapse = ";")
    match_res[j] = paste(gl_loc, collapse = ";")
  }

  grep_res = data.frame(Peptide = peptide, match = match_res)
  colnames(grep_res)[2] = c(database)
  
  return(grep_res)
}

# 读取命令行参数
args = commandArgs(trailingOnly = TRUE)
sp = args[1]
peptide_file = args[2]
database = args[3]
database_file = args[4]
sequence_merge = args[5] == "TRUE"  # TRUE 或 FALSE
outdir = args[6]

# 如果输出目录不存在，创建它
if (!dir.exists(outdir)) {
  dir.create(outdir)
}

# 读取肽段列表
peptide = read.table(peptide_file, sep = "\t", stringsAsFactors = F)$V1

# 调用函数进行匹配
res = grep_peptide1(database_file, sequence_merge, peptide, database)

# 将结果保存为文件
write.table(res, file.path(outdir, paste(database, "_grep.txt", sep = "")), sep = "\t", row.names = F, quote = F)

# 打印完成消息
print(paste(database, "grep Done!", sep = " "))
