# Set package version in place with 
#if (!require("checkpoint")) {
#  install.packages("checkpoint")
#}
#library(checkpoint)
#checkpoint("2015-11-02")

#basepath <- "/home/xavi/Estudis/YYYY-MM-XXX-YYY-ANNN" # Pentinella
#basepath <- "/mnt/magatzem02/tmp/YYYY-MM-XXX-YYY-ANNN" # MainHead
basepath <- baseDir

#################################################
# 0. Creacio de diccionari d'anotacions
################################################

# #se prueba de utilizar biomart para conseguir anotaciones aunque creo que luego no se utilizará.
# source("http://bioconductor.org/biocLite.R")
# if (!require("biomaRt")) {
#   biocLite("biomaRt")
# }
# library(biomaRt) 
# #mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# ensembl=useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
# dataset="hsapiens_gene_ensembl"
# mart=useDataset(dataset, mart=ensembl)
# listAttributes(mart)[1:100,]
# 
# miRNA <- getBM(c("mirbase_id", "ensembl_gene_id", "start_position", "chromosome_name"),
#                filters = c("with_mirbase"), values = list(TRUE), mart = mart)
# dim(miRNA)#1952 4
# head(miRNA)

#archivo de EC con anotaciones (rma+anotaciones)
setwd(file.path(basepath, dataRelDir))
annotation.affy<-read.csv("rma.affy.annotated.ref.txt",sep="\t",header=TRUE)
dim(annotation.affy)#36137 6
head(annotation.affy)

#nos quedamos solo con los humanos
annotation.affy.hg<-annotation.affy[annotation.affy$Species_Scientific_Name == "Homo sapiens",]
dim(annotation.affy.hg)#6631 6
rownames(annotation.affy.hg)<-annotation.affy.hg$Probe_Set_ID
head(annotation.affy.hg)

#se comprueban los uniques de Transcript_ID y se eliminan
#rownames(annotation.affy.hg)<-annotation.affy.hg$Transcript_ID
a<-as.vector(annotation.affy.hg$Transcript_ID)
length(a)#6631
length(unique(a)) #6050
annotation.affy.hg<-annotation.affy.hg[!duplicated(annotation.affy.hg[,5]),]
head(annotation.affy.hg)

#nuevo data set de anotaciones: Probe.Set.ID + Transcript.ID.Array.Design.
anota<-data.frame(annotation.affy.hg$Transcript_ID,annotation.affy.hg$Probe_Set_ID)
dim(anota)#6050 2
head(anota)
rownames(anota)<-anota$annotation.affy.hg.Probe_Set_ID
colnames(anota)<-c("Transcript_ID","Probe_Set_ID")
head(anota)

#################################################
# 1. Lectura de dades normalitzades amb Expression Console i anotació via affymetrix reference
#################################################

# 1.1 Lectura 
# -------------------------
setwd(file.path(basepath, dataRelDir))
data.rma <-read.csv2(paste0("rma.", aID, ".summary.txt"), header = TRUE, row.names=1,
               sep="\t", dec=",", comment.char = "#")
head(data.rma)
dim(data.rma)#36249 48
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
  col2remove.idx <- which(colnames(data.rma) %in% samplesX2remove)  
  data.rma <- data.rma[-col2remove.idx]
}
dim(data.rma) #36249    45
save(data.rma, targets, anota, file=file.path(basepath, dataRelDir,
                                        paste0("data.rma.", aID,".Rda")))


# 1.2 Anotacio amb  "anota" 
# -------------------------
data.rma.anotated <- merge(data.rma,anota,by="row.names")
dim(data.rma.anotated) #6050 48
head(data.rma.anotated)
#rownames(data.rma.anotated)<-data.rma.anotated$Row.names
rownames(data.rma.anotated)<-data.rma.anotated$Transcript_ID
#data.rma.anotated<-data.rma.anotated[,-c(1,57,58)]
ncol(data.rma.anotated) # 48
# Remove the first column and the last 2
data.rma.anotated<-data.rma.anotated[,-c(1,ncol(data.rma.anotated)-1,ncol(data.rma.anotated))]
dim(data.rma.anotated)#6050 45
head(data.rma.anotated)


save(data.rma.anotated, targets, anota, sample.names,  
     file=file.path(basepath, dataRelDir,
                    paste0("data.rma.annotated.", aID,".Rda")))

##############################################
### VSN
##############################################

# 1. Crear arxiu d'identificacio sondes-posicio al xip "cel definition" file cdf
# -------------------------
setwd(file.path(basepath, celfilesRelDir))
if (!require("makecdfenv")) {
  biocLite("makecdfenv")
}
if (!require("affxparser")) {
  biocLite("affxparser")
}
library(makecdfenv)    
library(affxparser)
if (!require("mirna40cdf")) {
  convertCdf("miRNA-4_0-st-v1.cdf", "mirna40cdf", version=4, verbose=TRUE) 
  pkgpath <- file.path(basepath, dataRelDir)
  make.cdf.package("mirna40cdf", version = packageDescription("makecdfenv", field = "Version"), species="", unlink=TRUE, compress=FALSE, package.path = pkgpath)
  
  system(paste0("R CMD INSTALL \"", pkgpath, "/mirna40cdf", "\""))
}

# 2. Lectura dels celfiles (que requereix el .cdf)
# -------------------------
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

rawData2 <- ReadAffy(filenames=fns)
#class(rawData2)
#str(rawData2)
# exprs(rawData2)

# 3. Normalitzacio amb VSN
# -------------------------
if (!require("vsn")) {
  biocLite("vsn")
}
data.vsn.rma = vsnrma(rawData2) # data.vsn.rma is the expression set
class(data.vsn.rma)
meanSdPlot(data.vsn.rma)
dim(exprs(data.vsn.rma))
head(exprs(data.vsn.rma))
save(data.vsn.rma, file = file.path(basepath, dataRelDir, paste0("data.vsn.rma.", aID,".Rda")))

if (!require("naturalsort")) {
  install.packages("naturalsort")
}
require(naturalsort)
# order(row.names(pData(data.vsn.rma))) # Wrong sorting with "order"
row.names(pData(data.vsn.rma))
colnames(exprs(data.vsn.rma))
table(row.names(pData(data.vsn.rma)) == colnames(exprs(data.vsn.rma)))
row.names(pData(data.vsn.rma)) == colnames(exprs(data.vsn.rma))

ns.data.vsn.rma <-naturalsort(row.names(pData(data.vsn.rma))) # ns -> natural sort
idx.data.vsn.rma <- pData(data.vsn.rma)[ns.data.vsn.rma,] # Indexes from data.vsn.rma sorted with natural sort
data.vsn.rma2 <- data.vsn.rma[,idx.data.vsn.rma]
rownames(targets) <- as.character(targets$SampleName)
targets2 <- targets[naturalsort(as.character(targets$SampleName)),]
#class(targets$SampleName)

table(colnames(exprs(data.vsn.rma2)) == rownames(targets2))
table(colnames(exprs(data.vsn.rma2)) == rownames(targets))
colnames(exprs(data.vsn.rma2)) == rownames(targets2)

pData(data.vsn.rma2) <- targets2

#--------------------------------
#############################
# FEATURE (GENE, miRNA, ...)  SELECTION
#############################
eset_norm <- data.vsn.rma2 # data.vsn.rma2 ja té les columnes i files (covariables) coincidint amb la info del targets
dim(exprs(eset_norm))
#[1] 36249    48
head(exprs(eset_norm))

## PART QUE ELIMINA LES FILES AMB IDENTICS VALORS D'EXPRESSIO (FERRAN, JUNY 2012)
## guardem com a valides les files de exprs(eset_norm) que NO estiguin duplicades
## (en cas de duplicats, es guarda la primera instancia de cada una com a bona)
repes <- duplicated(exprs(eset_norm), MARGIN=1)
#table(repes)
# repes
# FALSE  TRUE 
# 25111 11138
exprs(eset_norm) <- exprs(eset_norm)[!repes,]
#dim(eset_norm)
# Features  Samples 
# 25111       48


# XXX. BRB279. Repeteixo el process de filtrar d'abans, que sé que es quedan només amb els microRNA d'humans de l'expression set
# i ho torno a assignar al mateix objecte eset_norm, per a continuar amb el Basic Pipe estandard
eset_norm.hg <- eset_norm[rownames(exprs(eset_norm)) %in% as.character(annotation.affy.hg$Probe.Set.Name),]
#class(eset_norm.hg)
#dim(exprs(eset_norm.hg))
##[1] 5596   48
#[1] 5607   48
#head(exprs(eset_norm.hg))

# Eliminem els features que contenen la paraula "control", per que semblen ser controls d'Affymetrix o de l'experiment concret, i a més, un d'ells dona NA a totes les mostres
control.idx <- grep("control", rownames(exprs(eset_norm.hg)))
eset_norm.hg.nc <- eset_norm.hg[-control.idx,] # .nc stands for "No Controls"
#dim(exprs(eset_norm.hg.nc))
# 5595   48

# I per tant, creo a ma l'objecte data.eset.vsn.rma2.filtered, a partir dels valors d'expressió dels microRNA d'humans obtinguts en el pas anterior
data.eset.vsn.rma2.filtered <- eset_norm.hg.nc
# --------------------------------------



# Save Rda with normalized data
save(data.vsn.rma2, data.eset.vsn.rma2.filtered, file=file.path(basepath, dataRelDir,
                                        paste0("data.vsn.rma2.", aID,".Rda")))
save(annotation.affy.hg, anota, file=file.path(basepath, dataRelDir,
                                  paste0("data.annotation.affy.hg.Rda")))
save(data.rma, data.rma.anotated, file=file.path(basepath, dataRelDir,
                              paste0("data.rma.", aID,".Rda")))
save(rawData2, file=file.path(basepath, dataRelDir,
                              paste0("data.raw2.", aID,".Rda")))

# > order(pData(data.vsn.rma))
# > pData(data.vsn.rma) <- targets[order(pData(data.vsn.rma)[,1])]
# Error in `[.data.frame`(targets, order(pData(data.vsn.rma)[, 1])) : 
#   undefined columns selected
# > order(pData(data.vsn.rma)[,1])
# > pData(data.vsn.rma)[,1]

setwd(basepath)
