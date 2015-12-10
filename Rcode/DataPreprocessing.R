# Set package version in place with 
#if (!require("checkpoint")) {
#  install.packages("checkpoint")
#}
#library(checkpoint)
#checkpoint("2015-11-02")

#basepath <- "/home/xavi/Estudis/YYYY-MM-XXX-YYY-ANNN" # Pentinella
#basepath <- "/mnt/magatzem02/tmp/YYYY-MM-XXX-YYY-ANNN" # MainHead
basepath <- baseDir

#se prueba de utilizar biomart para conseguir anotaciones aunque creo que luego no se utilizará.
source("http://bioconductor.org/biocLite.R")
if (!require("biomaRt")) {
  biocLite("biomaRt")
}
library(biomaRt) 
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensembl=useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
dataset="hsapiens_gene_ensembl"
mart=useDataset(dataset, mart=ensembl)
listAttributes(mart)[1:100,]

miRNA <- getBM(c("mirbase_id", "ensembl_gene_id", "start_position", "chromosome_name"),
               filters = c("with_mirbase"), values = list(TRUE), mart = mart)
dim(miRNA)#1952 4
head(miRNA)

#archivo de EC con anotaciones (rma+anotaciones)
setwd(file.path(basepath, dataRelDir))
anotacion<-read.csv("rma_affy_annotated_ref.txt",sep="\t",header=TRUE)
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
length(unique(a)) #6050
anotacion.hg<-anotacion.hg[!duplicated(anotacion.hg[,5]),]
head(anotacion.hg)

#nuevo data set de anotaciones: Probe.Set.ID + Transcript.ID.Array.Design.
anota<-data.frame(anotacion.hg$Transcript_ID,anotacion.hg$Probe_Set_ID)
dim(anota)#6050 2
head(anota)
rownames(anota)<-anota$anotacion.hg.Probe_Set_ID
colnames(anota)<-c("Transcript_ID","Probe_Set_ID")
head(anota)

#preparación archivo Rda
setwd(file.path(basepath, dataRelDir))
xx <-read.csv2(paste0("rma.", aID, ".summary.txt"), header = TRUE, row.names=1,
               sep="\t", dec=",", comment.char = "#")
head(xx)
dim(xx)#36249 48
targets <- read.csv2 (targetsFileName, sep="\t")
# Remove the rows containing the samples indicated in this param "samples2remove" 
# since we do not want to consider them in the analysis
row2remove.idx <- which(targets$SampleName %in% samples2remove)  
# Remove rows only if there is some sample to be removed
if (length(row2remove.idx) > 0) {
  targets <- targets[-row2remove.idx,]
}
sample.names <- as.character(targets$ShortName)
# Get rid of samples that are eliminated from the targets
# In this run, 41.CEL, 42.CEL, 44.CEL files are eliminated
# samples2remove is defined in the basicA.R file, with some content like:
# samples2remove <- c("41.CEL","42.CEL","44.CEL")
# samplesX2remove <- c("X41.CEL","X42.CEL","X44.CEL")
# Remove rows only if there is some sample to be removed
if (length(samples2remove) > 1) {
  samplesX2remove <- paste0("X", samples2remove)
  col2remove.idx <- which(colnames(xx) %in% samplesX2remove)  
  xx <- xx[-col2remove.idx]
}
dim(xx) #36249    45
save(xx, targets, anota, file=file.path(basepath, dataRelDir,
                                        paste0("data.", aID,".Rda")))


#merge xx con anota#########################
expresEG <- merge(xx,anota,by="row.names")
dim(expresEG)#6050 48
head(expresEG)
#rownames(expresEG)<-expresEG$Row.names
rownames(expresEG)<-expresEG$Transcript_ID
#expresEG<-expresEG[,-c(1,57,58)]
ncol(expresEG) # 48
# Remove the first column and the last 2
expresEG<-expresEG[,-c(1,ncol(expresEG)-1,ncol(expresEG))]
dim(expresEG)#6050 45
head(expresEG)


save(expresEG, targets,anota, sample.names,  
     file=file.path(basepath, dataRelDir,
                    paste0("Micros_", aID,"_expressEG.Rda")))

##############################################
### VSN
##############################################
setwd(file.path(basepath, celfilesRelDir))
if (!require("makecdfenv")) {
  biocLite("makecdfenv")
}
if (!require("affxparser")) {
  biocLite("affxparser")
}
library(makecdfenv)    
library(affxparser)
if (!require("mirna41cdf")) {
  convertCdf("miRNA-4_0-st-v1.cdf", "mirna41cdf", version=4, verbose=TRUE) 
  pkgpath <- file.path(basepath, dataRelDir)
  make.cdf.package("mirna41cdf", version = packageDescription("makecdfenv", field = "Version"), species="", unlink=TRUE, compress=FALSE, package.path = pkgpath)
  
  system(paste0("R CMD INSTALL \"", pkgpath, "/mirna41cdf", "\""))
}

library(affy)
require(affy)
# Before reading the celfiles through ReadAffy, limit the filenames to get rid of the ones in samples2remove
fns <- list.celfiles(path=file.path(basepath, celfilesRelDir),full.names=TRUE)
# Remove rows only if there is some sample to be removed
if (length(samples2remove) > 1) {
  cel2remove.idx <- match(unique(grep(paste(samples2remove,collapse="|"),
                                      fns, value=TRUE)), fns)
  fns <- fns[-cel2remove.idx]
}
cat("Reading files:\n",paste(fns,collapse="\n"),"\n")
cat("Removed from the analysis:\n",paste(samples2remove,collapse="\n"),"\n")
##read a binary celfile
xx <- ReadAffy(filenames=fns)
#class(xx)
#str(xx)
# exprs(xx)
if (!require("vsn")) {
  biocLite("vsn")
}
d_vsn = vsnrma(xx) # d_vsn is the expression set
class(d_vsn)
meanSdPlot(d_vsn)
dim(exprs(d_vsn))
head(exprs(d_vsn))
save(d_vsn, file = "expres.filtered.Rda")

if (!require("naturalsort")) {
  install.packages("naturalsort")
}
require(naturalsort)
# order(row.names(pData(d_vsn))) # Wrong sorting with "order"
row.names(pData(d_vsn))
colnames(exprs(d_vsn))
table(row.names(pData(d_vsn)) == colnames(exprs(d_vsn)))
row.names(pData(d_vsn)) == colnames(exprs(d_vsn))

ns.d_vsn <-naturalsort(row.names(pData(d_vsn))) # ns -> natural sort
idx.d_vsn <- pData(d_vsn)[ns.d_vsn,] # Indexes from d_vsn sorted with natural sort
d_vsn2 <- d_vsn[,idx.d_vsn]
rownames(targets) <- as.character(targets$SampleName)
targets2 <- targets[naturalsort(as.character(targets$SampleName)),]
#class(targets$SampleName)

table(colnames(exprs(d_vsn2)) == rownames(targets2))
table(colnames(exprs(d_vsn2)) == rownames(targets))
colnames(exprs(d_vsn2)) == rownames(targets2)

pData(d_vsn2) <- targets2

# > order(pData(d_vsn))
# > pData(d_vsn) <- targets[order(pData(d_vsn)[,1])]
# Error in `[.data.frame`(targets, order(pData(d_vsn)[, 1])) : 
#   undefined columns selected
# > order(pData(d_vsn)[,1])
# > pData(d_vsn)[,1]

setwd(basepath)
