##output a list of pdb file paths and the motif to align them on
##
setwd("./VX22/AF3")
prompt2PDBs<-read.csv(file="../bestPDBs/prompt1sPDBs.txt",sep=" ",skip=1,header=F)
head(prompt2PDBs)
library(stringr)
prompt1PDBs$sampleNumber<-str_extract(prompt1PDBs$V11,"sample(\\d+)")
sampleNumbers2<-str_extract(prompt1PDBs$sampleNumber,"\\d+")
head(sampleNumbers2)
prompt1PDBs$sampleNumber2<-sampleNumbers2
head(prompt1PDBs)
prompt1PDBs$motifs<-prompt1TRB[as.numeric(prompt1PDBs$sampleNumber2)+1,2]
head(prompt1PDBs)
motifs2<-str_split(prompt1PDBs$motifs,",")
prompt1PDBs$motif2<-prompt1PDBs$motifs
head(motifs2)
for(i in c(1:length(motifs2))){
  print(i)
  stringOfMotifs<-paste0(motifs2[[i]][1],"-",str_trim(motifs2[[i]][2]),"+",str_trim(motifs2[[i]][3]),
                         "-",str_trim(motifs2[[i]][5]),"+",str_trim(motifs2[[i]][6]),"-",str_trim(motifs2[[i]][11]),
                         "+",str_trim(motifs2[[i]][12]),"-",str_trim(motifs2[[i]][13]))
  prompt1PDBs[i,15]<-stringOfMotifs
}
head(prompt1PDBs)
prompt1PDBs$prompt1FileName<-paste0("./VX22/bestPDBs/prompt1s/",prompt1PDBs$V11)
head(prompt1PDBs)
toPrint<-prompt1PDBs[,c("prompt1FileName","motif2")]
write.table(toPrint,file="../prompt1PDBFilesAndMotifs.csv",quote=F,sep=",",col.names = T,row.names = F)


prompt2PDBs<-read.csv(file="./VX22/bestPDBs/prompt2/prompt2PDB.txt",sep=" ",skip=1,header=F)
library(stringr)
prompt2PDBs$sampleNumber<-str_extract(prompt2PDBs$V1,"sample(\\d+)")
sampleNumbers2<-str_extract(prompt2PDBs$sampleNumber,"\\d+")
prompt2PDBs$sampleNumber2<-sampleNumbers2
head(prompt2PDBs)
prompt2PDBs$motifs<-prompt2TRB[as.numeric(prompt2PDBs$sampleNumber2)+1,2]
motifs2<-str_split(prompt2PDBs$motifs,",")
prompt2PDBs$motif2<-prompt2PDBs$motifs
for(i in c(1:length(motifs2))){
  print(i)
  stringOfMotifs<-paste0(motifs2[[i]][1],"-",str_trim(motifs2[[i]][2]),"+",str_trim(motifs2[[i]][3]),
                         "-",str_trim(motifs2[[i]][5]),"+",str_trim(motifs2[[i]][6]),"-",str_trim(motifs2[[i]][11]),
                         "+",str_trim(motifs2[[i]][12]),"-",str_trim(motifs2[[i]][13]))
  prompt2PDBs[i,4]<-stringOfMotifs
}
colnames(prompt2PDBs)[1]<-"prompt2FileName"
prompt2PDBs$motifResi<-paste0("resi ",prompt2PDBs$motifs)
prompt2PDBs$fileName<-paste0("./VX22/bestPDBs/prompt2/prompt2PDB/",prompt2PDBs$prompt2FileName)
toPrint<-prompt2PDBs[,c("fileName","motifResi")]
write.table(toPrint,file="./VX22/bestPDBs//prompt2PDBFilesAndMotifs.txt",quote=F,sep="\t",col.names = T,row.names = F)


prompt1iiPDBs<-read.csv(file="./VX22/bestPDBs/prompt1iis/prompt1iiPDBs.txt",sep=" ",skip=1,header=F)
prompt1iiPDBs$sampleNumber<-str_extract(prompt1iiPDBs$V1,"\\d+")



##prompt1iiPDB part 20 onwards
prompt1ii_20<-prompt1iiPDBs[c(20:54),]
sampleNumbers2<-str_extract_all(prompt1ii_20$V1,"\\d+")
for(i in c(1:length(sampleNumbers2))){
  prompt1ii_20[i,3]<-sampleNumbers2[[i]][2]
}

prompt1ii_20$motifs<-prompt1iiTRB[as.numeric(prompt1ii_20$sampleNumber2)+1,2]
motifs2<-str_split(prompt1ii_20$motifs,",")
prompt1ii_20$motif2<-prompt1ii_20$motifs
for(i in c(1:length(motifs2))){
  print(i)
  stringOfMotifs<-paste0(motifs2[[i]][1],"-",str_trim(motifs2[[i]][2]),"+",str_trim(motifs2[[i]][3]),
                         "-",str_trim(motifs2[[i]][5]),"+",str_trim(motifs2[[i]][6]),"-",str_trim(motifs2[[i]][11]),
                         "+",str_trim(motifs2[[i]][12]),"-",str_trim(motifs2[[i]][13]))
  prompt1ii_20[i,5]<-stringOfMotifs
}


colnames(prompt1ii_20)[1]<-"prompt1iiFileName"
prompt1ii_20$motifResi<-paste0("resi ",prompt1ii_20$motif2)
prompt1ii_20$fileName<-paste0("./VX22/bestPDBs/prompt1iis/prompt1iiPDB/",prompt1ii_20$prompt1iiFileName)
toPrint<-prompt1ii_20[,c("fileName","motifResi")]
write.table(toPrint,file="./VX22/bestPDBs/prompt1ii_20FilesAndMotifs.txt",quote=F,sep="\t",col.names = T,row.names = F)


prompt2iiPDBs<-read.csv(file="./VX22/bestPDBs/prompt2ii/prompt2iiPDB.txt",sep=" ",skip=1,header=F)
prompt2iiPDBs$sampleNumber<-str_extract(prompt2iiPDBs$V1,"sample(\\d+)")
sampleNumbers2<-str_extract(prompt2iPDBs$sampleNumber,"\\d+")
prompt2iiPDBs$sampleNumber2<-sampleNumbers2
sampleNumbers3<-str_extract_all(prompt2iiPDBs$V1,"\\d+")
NAs<-is.na(prompt2iiPDBs$sampleNumber2)
indices<-which(NAs)
for(i in c(1:length(indices))){
prompt2iiPDBs[i,2]<-sampleNumbers3[[i]][2]
}

prompt2iiPDBs$motifs<-prompt2iiTRB[as.numeric(prompt2iiPDBs$sampleNumber)+1,2]
motifs2<-str_split(prompt2iiPDBs$motifs,",")
prompt2iiPDBs$motif2<-prompt2iiPDBs$motifs
for(i in c(1:length(motifs2))){
  print(i)
  stringOfMotifs<-paste0(motifs2[[i]][1],"-",str_trim(motifs2[[i]][2]),"+",str_trim(motifs2[[i]][3]),
                         "-",str_trim(motifs2[[i]][5]),"+",str_trim(motifs2[[i]][6]),"-",str_trim(motifs2[[i]][11]),
                         "+",str_trim(motifs2[[i]][12]),"-",str_trim(motifs2[[i]][13]))
  prompt2iiPDBs[i,4]<-stringOfMotifs
}
colnames(prompt2iiPDBs)[1]<-"prompt2iiFileName"
prompt2iiPDBs$motifResi<-paste0("resi ",prompt2iiPDBs$motif2)
prompt2iiPDBs$fileName<-paste0("./VX22/bestPDBs/prompt2ii/prompt2iiPDB/",prompt2iiPDBs$prompt2iiFileName)
toPrint<-prompt2iiPDBs[,c("fileName","motifResi")]
write.table(toPrint,file="./VX22/bestPDBs//prompt2iiPDBsFilesAndMotifs.txt",quote=F,sep="\t",col.names = T,row.names = F)
