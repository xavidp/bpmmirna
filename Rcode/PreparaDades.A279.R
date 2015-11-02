# Set package version in place with 
install.packages("checkpoint")
library(checkpoint)
checkpoint("2015-11-02")

basepath <- "/home/xavi/Estudis/2015-10-NuriaBarbarroja-IMIBIC-A279/"

#se prueba de utilizar biomart para conseguir anotaciones aunque creo que luego no se utilizará.
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
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
setwd("/home/rgonzalo/Documents/Estudis/2015-07-TeresaGarcia-VHIR-A/dades")
x <-read.csv2("rma.csv", row.names=1, sep="\t", dec=".")
head(x)
dim(x)#36249 55
targets <- read.csv2 ("targets.csv",sep="\t")
sample.names <- as.character(targets$ShortName)
save(x, targets, anota, file=file.path("/home/rgonzalo/Documents/Estudis/2015-07-TeresaGarcia-VHIR-A/dades","dades.Rda"))


#merge x con anota#########################
expresEG <- merge(x,anota,by="row.names")
dim(expresEG)#6050 58
head(expresEG)
#rownames(expresEG)<-expresEG$Row.names
rownames(expresEG)<-expresEG$Transcript_ID
expresEG<-expresEG[,-c(1,57,58)]
dim(expresEG)#6050 55
head(expresEG)


save(expresEG, targets,anota, sample.names,  file=file.path("/home/rgonzalo/Documents/Estudis/2015-07-TeresaGarcia-VHIR-A/dades","MicrosTGB_expressEG.Rda"))

##############################################
### VSN
##############################################
setwd("/home/rgonzalo/Documents/Estudis/2015-07-TeresaGarcia-VHIR-A/CELFiles")
biocLite("makecdfenv")
library(makecdfenv)    
library(affxparser)
convertCdf("miRNA-4_0-st-v1.cdf", "mirna40cdf", version=4, verbose=TRUE) 
pkgpath <- "/home/rgonzalo/Documents/Estudis/2015-07-TeresaGarcia-VHIR-A/dades"
make.cdf.package("mirna40cdf", version = packageDescription("makecdfenv", field = "Version"), species="", unlink=TRUE, compress=FALSE, package.path = pkgpath)



require(affy)
x<- ReadAffy()
class(x)
d_vsn = vsnrma(x)
meanSdPlot(d_vsn)

