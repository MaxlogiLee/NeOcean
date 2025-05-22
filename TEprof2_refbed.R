library(foreach)
sp1=c("QLX-P2","SGSL-D2","CXML-R1","FM-P1","PG-P1","JN-P2")
sp2=c("QLXP2","SGSLM","CXMLM","FMM","PGT","JNT")

for(i in 1:length(sp1)){
  dt=read.csv(paste("data/5.TEprof2/step10/",sp2[i],"_All_TE-derived_Alternative_Isoforms_Statistics.csv",sep=""))
  dt$id=paste(dt$Subfam,dt$Start.TE,dt$Gene,dt$Location.TE,dt$Gene,dt$Splice.Target,sep="_")
  gl1=grep("Intergenic",dt$Location.TE)
  gl2=grep("Intergenic",dt$Splice.Target)
  gl_1=intersect(gl1,gl2)
  gl_2=setdiff(gl1,gl2)
  gl_3=setdiff(gl2,gl1)
  if(length(gl_1)!=0){
    dt$id[gl_1]=paste(dt$Subfam[gl_1],dt$Start.TE[gl_1],"None","None_None","None","None_None",sep="_")
  }
  if(length(gl_2)!=0){
    dt$id[gl_2]=paste(dt$Subfam[gl_2],dt$Start.TE[gl_2],"None","None_None",dt$Gene[gl_2],dt$Splice.Target[gl_2],sep="_")
  }
  if(length(gl_3)!=0){
    dt$id[gl_3]=paste(dt$Subfam[gl_3],dt$Start.TE[gl_3],dt$Gene[gl_3],dt$Location.TE[gl_3],"None","None_None",sep="_")
  }
  df=read.csv(paste("data/5.TEprof2/refBed/",sp2[i],".stringtie.refbed",sep=""),header = F,sep="\t",stringsAsFactors = F)
  
  fa=readLines(paste("data/5.TEprof2/",sp2[i],"_candidates.getorf.fa",sep=""))
  fa=sub(">","",fa[grep(">",fa)])
  fa_split=strsplit(fa,split = "_")
  fa_len=foreach(j=1:length(fa_split),.combine = c) %do% length(fa_split[[j]])
  spl1=which(fa_len==11)
  for(j in 1:length(spl1)){
    fa_split[[spl1[j]]]=c(paste(fa_split[[spl1[j]]][1:2],collapse = "_"),fa_split[[spl1[j]]][-c(1:2)])
  }
  spl2=which(fa_len==12)
  for(j in 1:length(spl2)){
    fa_split[[spl2[j]]]=c(paste(fa_split[[spl2[j]]][1:3],collapse = "_"),fa_split[[spl2[j]]][-c(1:3)])
  }
  fa_split=do.call(rbind,fa_split)
  fa1=foreach(j=1:nrow(fa_split),.combine = c) %do% paste(fa_split[j,1:8],collapse = "_")
  fa2=foreach(j=1:nrow(fa_split),.combine = c) %do% paste(fa_split[j,1:9],collapse = "_")
  
  fa_dt=data.frame(id=fa1,ID=fa2)
  fa_dt=unique(fa_dt)
  fa_dt_id=merge(fa_dt,dt[,c("Transcript.Name","id")])
  fa_dt_id_df=merge(df,fa_dt_id,by.x = "V8",by.y="Transcript.Name")
  TEpro2=fa_dt_id_df[,c(2:8,1,9:11,13)]
  write.table(TEpro2,paste("data/5.TEprof2/",sp1[i],"_TEprof2.refbed",sep=""),sep="\t",col.names = F,row.names = F,quote = F)
}
