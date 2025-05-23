library(foreach)
library(openxlsx)
args=commandArgs(trailingOnly=TRUE)
pep_file=args[1]
out_dir=args[2]
cutoff=args[3]
target_db=args[4]
decoy_db=args[5]
#write.table(union(all_seq$Sequence1,all_seq$Sequence2),paste(outdir,"/","Uniport_cut_peptidelist.txt",sep=""),sep="\t",row.names = F,col.names = F,quote = F)
seq=readLines(pep_file)

prt_res=foreach(j=1:length(seq)) %do% system(paste("LC_ALL=C grep -B1 ",seq[j]," ",target_db,"|grep '^>'",sep=""),intern=T)
prt_N=foreach(j=1:length(seq),.combine = c) %do% length(prt_res[[j]])

decoy_res=foreach(j=1:length(seq)) %do% system(paste("LC_ALL=C grep -B1 ",seq[j]," ",decoy_db,"|grep '^>'",sep=""),intern=T)
decoy_N=foreach(j=1:length(seq),.combine = c) %do% length(decoy_res[[j]])

target_decoy=data.frame(Sequence=seq,Target=prt_N,Decoy=decoy_N)
target_decoy$FDR=target_decoy$Decoy/target_decoy$Target

write.table(target_decoy,paste(out_dir,"/target_decoy_result.xls",sep=""),sep="\t",row.names = F,quote=F) 

sig_seq=target_decoy$Sequence[target_decoy$Target!=0&target_decoy$FDR<cutoff]
write.table(sig_seq,paste(out_dir,"/target_decoy_signif_sequences.txt",sep=""),sep="\t",row.names = F,col.names = F,quote=F)
