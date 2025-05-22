#sp=tab$Name
library(foreach)
library(stringr)
library(openxlsx)
grep_peptide1=function(database_file,sequence_merge=FALSE,peptide,database,split_str=NULL){
  if(sequence_merge){
    fa <- system(paste("awk '{if($0 ~ />/) name = $0; else seq[name] = seq[name] $0;} END {for(i in seq) print i\"\\n\"seq[i]}'", database_file), intern = TRUE)
  }else{
    fa=readLines(database_file)
  }
  db=fa[-grep(">",fa)]
  names(db)=fa[grep(">",fa)]
  gl_res=foreach(j=1:length(peptide)) %do% system(paste("LC_ALL=C grep -B1 ",peptide[j]," ",database_file,"|grep '^>'",sep=""),intern=T)
  is.match=foreach(j=1:length(peptide),.combine=c) %do% length(gl_res[[j]])
  match=which(is.match!=0)
  match_res=rep("",length(peptide))
  for(j in match){
     gl_loc=foreach(k=1:length(gl_res[[j]]),.combine=c) %do% paste(paste(sub("^>","",gl_res[[j]][k]),"[",paste(str_locate_all(db[gl_res[[j]][k]],peptide[j])[[1]][,1],str_locate_all(db[gl_res[[j]][k]],peptide[j])[[1]][,2],sep = "-"),
                                                                             "]",sep=""),collapse=";")
         match_res[j]=paste(gl_loc,collapse=";")
    }
  
  grep_res=data.frame(Peptide=peptide,match=match_res)
  colnames(grep_res)[2]=c(database)
  return(grep_res)
}

args=commandArgs(trailingOnly=TRUE)
sp=args[1]
peptide_file=args[2]
database=args[3]
database_file=args[4]
sequence_merge=args[5]=="TRUE"
outdir=args[6]
if(!dir.exists(outdir)){
  dir.create(outdir)
}

peptide=read.table(peptide_file,sep="\t",stringsAsFactors = F)$V1
res=grep_peptide1(database_file,sequence_merge,peptide,database)
write.table(res,file.path(outdir,paste(database,"_grep.txt",sep="")),sep="\t",row.names = F,quote = F)
print(paste(database,"grep Done!",sep=" "))
