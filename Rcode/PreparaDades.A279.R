# Set package version in place with 
#if (!require("checkpoint")) {
#  install.packages("checkpoint")
#}
#library(checkpoint)
#checkpoint("2015-11-02")

basepath <- "/home/xavi/Estudis/2015-10-NuriaBarbarroja-IMIBIC-A279/"

#se prueba de utilizar biomart para conseguir anotaciones aunque creo que luego no se utilizará.
source("http://bioconductor.org/biocLite.R")
if (!require("biomaRt")) {
  biocLite("biomaRt")
}
library(biomaRt) 
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
listAttributes(mart)[1:100,]

miRNA <- getBM(c("mirbase_id", "ensembl_gene_id", "start_position", "chromosome_name"),
               filters = c("with_mirbase"), values = list(TRUE), mart = mart)
dim(miRNA)#1952 4
head(miRNA)

#archivo de EC con anotaciones (rma+anotaciones)
setwd(paste0(basepath, "dades"))
anotacion<-read.csv("rma with anotations.csv",sep="\t",header=TRUE)
dim(anotacion)#36137 6
head(anotacion)

#nos quedamos solo con los humanos
anotacion.hg<-anotacion[anotacion$Species_Scientific_Name == "Homo sapiens",]
dim(anotacion.hg)#6631 6
rownames(anotacion.hg)<-anotacion.hg$Probe_Set_ID
head(anotacion.hg)

#se comprueban los uniques de Transcript_ID y se eliminan
#rownames(anotacion.hg)<-anotacion.hg$Transcript_ID
a<-as.vector(anotacion.hg$Transcript_ID)
length(a)#6631
length(unique(a))
anotacion.hg<-anotacion.hg[!duplicated(anotacion.hg[,5]),]

#nuevo data set de anotaciones: Probe.Set.ID + Transcript.ID.Array.Design.
anota<-data.frame(anotacion.hg$Transcript_ID,anotacion.hg$Probe_Set_ID)
dim(anota)#6050 2
head(anota)
rownames(anota)<-anota$anotacion.hg.Probe_Set_ID
colnames(anota)<-c("Transcript_ID","Probe_Set_ID")
head(anota)

#preparación archivo Rda
setwd(paste0(basepath, "dades"))
x <-read.csv2("rma.A279.summary.csv", row.names=1, sep=";", dec=",")
head(x)
dim(x)#36249 55
targets <- read.csv2 ("targets.BRB.A279.csv",sep="\t")
sample.names <- as.character(targets$ShortName)
save(x, targets, anota, file=file.path(paste0(basepath, "dades"),"dades.BRB.A279.Rda"))


#merge x con anota#########################
expresEG <- merge(x,anota,by="row.names")
dim(expresEG)#6050 51
head(expresEG)
#rownames(expresEG)<-expresEG$Row.names
rownames(expresEG)<-expresEG$Transcript_ID
#expresEG<-expresEG[,-c(1,57,58)]
expresEG<-expresEG[,-c(1,50,51)]
dim(expresEG)#6050 48
head(expresEG)


save(expresEG, targets,anota, sample.names,  file=file.path(paste0(basepath, "dades"),"MicrosBRB_expressEG.Rda"))

##############################################
### VSN
##############################################
setwd(paste0(basepath, "celfiles"))
if (!require("makecdfenv")) {
  biocLite("makecdfenv")
}
if (!require("affxparser")) {
  biocLite("affxparser")
}
library(makecdfenv)    
library(affxparser)
convertCdf("miRNA-4_0-st-v1.cdf", "mirna40cdf", version=4, verbose=TRUE) 
pkgpath <- paste0(basepath, "dades")
make.cdf.package("mirna40cdf", version = packageDescription("makecdfenv", field = "Version"), species="", unlink=TRUE, compress=FALSE, package.path = pkgpath)

system(paste0("R CMD INSTALL \"", pkgpath, "/mirna40cdf", "\""))

library(affy)
require(affy)
x<- ReadAffy()
class(x)
if (!require("vsn")) {
  biocLite("vsn")
}
d_vsn = vsnrma(x)
meanSdPlot(d_vsn)

