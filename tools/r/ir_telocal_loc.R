library(foreach)
args=commandArgs(trailingOnly=TRUE)
input_file=args[1]
output_file=args[2]
merge_loc=args[3]=="TRUE"

dt=read.table(input_file,sep="\t",header = T,stringsAsFactors = F)

item=colnames(dt)[2]

gl=which(dt[,item]!="")
peptide_loc=c()
for(j in 1:length(gl)){
  loc=as.data.frame(do.call(rbind,strsplit(strsplit(dt[gl[j],item],split=";")[[1]],split="\\[|\\]")))
  loc=cbind(loc,do.call(rbind,strsplit(loc$V2,split=" - ")),do.call(rbind,strsplit(loc$V4,split="-")))
  colnames(loc)[5:8]=c("N_pos1","N_pos2","TC_pos1","TC_pos2")
  loc$pos1=foreach(k=1:nrow(loc),.combine = c) %do% (as.numeric(loc$N_pos1)[k]:as.numeric(loc$N_pos2)[k])[as.numeric(loc$TC_pos1)[k]*3-2]
  loc$pos2=foreach(k=1:nrow(loc),.combine = c) %do% (as.numeric(loc$N_pos1)[k]:as.numeric(loc$N_pos2)[k])[as.numeric(loc$TC_pos2)[k]*3]
  loc$id=sub("_[0-9]*.$","",loc$V1)
  #一个一个地遍历,结果保存在pep_loc中
  pep_loc=c()
  for(m in 1:nrow(loc)){
    if(item=="IR"){
    rg_chr=unlist(strsplit(loc$id[m],split=":|\\-"))[2]
    rg_start=as.numeric(unlist(strsplit(unlist(strsplit(loc$id[m],split=":"))[2],split = "\\-"))[1])+1
    rg_end=as.numeric(unlist(strsplit(unlist(strsplit(loc$id[m],split=":"))[2],split = "\\-"))[2])+1
    }
    if(item=="TElocal"){
      rg_chr=unlist(strsplit(loc$id[m],split=":|\\_"))[2]
      rg_start=as.numeric(unlist(strsplit(unlist(strsplit(loc$id[m],split=":"))[2],split = "\\-"))[1])+1
      rg_end=as.numeric(unlist(strsplit(unlist(strsplit(loc$id[m],split=":"))[2],split = "\\-"))[2])+1
      rg_id=sub("_[0-9]*$","",unlist(strsplit(loc$id[m],split=":| "))[1])
    }
    #判断peptide的核苷酸起点和终点落在哪一个exon中,且是不是同一个exon
    pep_start=loc$pos1[m]
    pep_end=loc$pos2[m]
    
    #把坐标向量化，也就是生成RNA的坐标，整个长度就是exon的总长度,根据RNA上核苷酸的位置，直接就能确定它的location
    loc_vector=rg_start:rg_end
    if(item=="IR"){
    pep_loc[m]=paste(item,"-",rg_chr,":",loc_vector[pep_start],"-",
                    loc_vector[pep_end],sep="")
    }
    if(item=="TElocal"){
      pep_loc[m]=paste(item,"-",rg_chr,":",loc_vector[pep_start],"-",
                       loc_vector[pep_end],"(",rg_id,")",sep="")
    }
    
  }
  if(merge_loc){
    peptide_loc[j]=paste(unique(pep_loc),collapse = ";")
  }else{
    peptide_loc[j]=paste(pep_loc,collapse = ";")
  }
}

dt$loc=""
dt$loc[gl]=peptide_loc

write.table(dt,output_file,row.names = F,quote = F,sep='\t')