library(foreach)
library(openxlsx)
library(reshape2)
library(dplyr)
args=commandArgs(trailingOnly=TRUE)
grep_file=args[1]
peptide_file=args[2]
target_decoy_file=args[3]
ms_file=args[4]
sp_file=args[5]

dt=read.delim(grep_file,sep="\t",stringsAsFactors = F,header = T)
prt=foreach(i=1:nrow(dt),.combine = rbind) %do% data.frame(Peptide=dt$Peptide[i],Uniport_cut=strsplit(dt$Uniport_cut[i],split=";")[[1]])
info=foreach(i=1:nrow(prt)) %do% strsplit(prt$Uniport_cut[i],split=" ")[[1]]
prt$Uniport=foreach(i=1:nrow(prt),.combine = c) %do% sub("uniprot-sp\\||uniprot-tr\\|","",info[[i]][1])
prt$Start=foreach(i=1:nrow(prt),.combine = c) %do% strsplit(info[[i]][length(info[[i]])],split="\\[|\\-")[[1]][2]
prt$End=foreach(i=1:nrow(prt),.combine = c) %do% strsplit(info[[i]][length(info[[i]])],split="\\-|\\]")[[1]][2]

df=read.table(peptide_file,sep="\t",stringsAsFactors = F,header = T)
target_decoy=read.table(target_decoy_file,sep="\t",stringsAsFactors = F,header = T)

ft=merge(df,prt[,-2],by.x="Sequence1",by.y="Peptide")
ft=merge(ft,prt[,-2],by.x="Sequence2",by.y="Peptide")
ft=merge(ft,target_decoy,by.x="Sequence1",by.y="Sequence")
ft=merge(ft,target_decoy,by.x="Sequence2",by.y="Sequence")
colnames(ft)[9:20]=c("Uniport1","Start1","End1","Uniport2","Start2","End2","Target1","Decoy1","FDR1","Target2","Decoy2","FDR2")

#删除Linear的结果
ft_full=ft[((ft$Uniport1==ft$Uniport2)&(as.numeric(ft$End1)==as.numeric(ft$Start2)-1)),]

cut_anot=ft[!((ft$Uniport1==ft$Uniport2)&(as.numeric(ft$End1)==as.numeric(ft$Start2)-1)),]
cut_anot$pep=paste(cut_anot$Sequence1,cut_anot$Sequence2,sep="|")


cis_anot=cut_anot[cut_anot$Uniport1==cut_anot$Uniport2,]
trans_anot=cut_anot[(cut_anot$Uniport1!=cut_anot$Uniport2)&
                      !is.element(cut_anot$Peptide,cis_anot$Peptide),]
trans_anot=cut_anot[(cut_anot$Uniport1!=cut_anot$Uniport2),]

sel_col=c("Peptide","Length","Sequence1","Sequence2","pep","Target1","Decoy1","FDR1",
          "Target2","Decoy2","FDR2","Uniport1","Start1","End1","Uniport2","Start2","End2")
if(nrow(cis_anot)!=0){
  write.xlsx(cis_anot[,sel_col],"res/Cis_Peptide_annotation.xlsx",rowNames=F)
}
if(nrow(trans_anot)!=0){
  write.xlsx(trans_anot[,sel_col],"res/Trans_Peptide_annotation.xlsx",rowNames=F)
}


#整理table
sp=read.table(sp_file,sep="\t",header = F,stringsAsFactors = F)$V1

dat=read.xlsx(ms_file)
dat_melt=melt(dat[,c("Peptide",colnames(dat)[grep("MS",colnames(dat))])])
colnames(dat_melt)[1]="Peptide"
dat_melt=dat_melt[dat_melt$value!=0,]

cis_anot$Uniport=paste(cis_anot$Uniport1,"[",cis_anot$Start1,"-",cis_anot$End1,"],",
                       cis_anot$Uniport2,"[",cis_anot$Start2,"-",cis_anot$End2,"]",sep="")
trans_anot$Uniport=paste(trans_anot$Uniport1,"[",trans_anot$Start1,"-",trans_anot$End1,"],",
                         trans_anot$Uniport2,"[",trans_anot$Start2,"-",trans_anot$End2,"]",sep="")
cis_tab=cis_anot%>%group_by(Peptide)%>%summarise(Uniport=paste0(Uniport,collapse = ";"))
trans_tab=trans_anot%>%group_by(Peptide)%>%summarise(Uniport=paste0(Uniport,collapse = ";"))

cis_dat=merge(dat_melt,cis_tab)
cis_dat$variable=sub("MS","Cis",cis_dat$variable)
cis_res=dcast(unique(cis_dat),Peptide~variable,value.var = "Uniport")
cis_res=merge(dat,cis_res)
cis_res[is.na(cis_res)]=""
cis_res=unique(cis_res)
cis_res=cis_res[,intersect(c("Peptide","Peptide_raw","Length","Found.By",paste("MS_",sp,sep=""),paste("Cis","_",sp,sep="")),colnames(cis_res))]

trans_dat=merge(dat_melt,trans_tab)
trans_dat$variable=sub("MS","Trans",trans_dat$variable)
trans_res=dcast(unique(trans_dat),Peptide~variable,value.var = "Uniport")
trans_res=merge(dat,trans_res)
trans_res[is.na(trans_res)]=""
trans_res=unique(trans_res)
trans_res=trans_res[,intersect(c("Peptide","Peptide_raw","Length","Found.By",paste("MS_",sp,sep=""),paste("Trans","_",sp,sep="")),colnames(trans_res))]

ms=dat[!is.element(dat$Peptide,union(cis_res$Peptide,trans_tab$Peptide)),]
write.table(cis_res,paste("res/10_","cis","_MS_",nrow(cis_res),".xls",sep=""),sep="\t",row.names = F,quote = F)
write.table(trans_res,paste("res/10_","trans","_MS_",nrow(trans_res),".xls",sep=""),sep="\t",row.names = F,quote = F)
write.table(ms,paste("res/11_","ungreped","_MS_",nrow(ms),"_peplist.xls",sep=""),sep="\t",row.names = F,quote = F)


