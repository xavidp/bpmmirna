## ----librerias, eval=TRUE------------------------------------------------
installifnot <- function (pckgName){
  if(!(require(pckgName, character.only=TRUE))){
    source("http://Bioconductor.org/biocLite.R")
    biocLite(pckgName)
  }
}
installifnot("Biobase")
installifnot("hgu133a.db")
installifnot("affy")
installifnot("affyPLM")
installifnot("arrayQualityMetrics")
installifnot("genefilter")
installifnot("limma")
installifnot("hgu133a.db")
installifnot("annotate")
installifnot("annaffy")
installifnot("hwriter")
installifnot("gplots")
installifnot("GOstats")
if(!require(SortableHTMLTables)) install.packages("SortableHTMLTables")

## ----preparaDirectorios, eval=TRUE---------------------------------------
baseDir <- "/home/xavi/Estudis/2015-10-NuriaBarbarroja-IMIBIC-A279"
workingDir <-setwd(baseDir)
dataDir <-file.path(workingDir, "dades")
resultsDir <- file.path(workingDir,"results")
celfilesDir <- file.path(workingDir,"celfiles")
setwd(workingDir)

## ----loadData------------------------------------------------------------
#load(file=file.path(resultsDir, "datos.normalizados.Rda"))
source("Rcode/PreparaDades.A279.R")

#############################
# FILTRATION
#############################
eset_norm <- d_vsn2 # d_vsn2 ja té les columnes i files (covariables) coincidint amb la info del targets
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

# # XXX BRB279: Aquí hauré de fer a ma probablement el filtrat dels microRNA que són controls d'Affymetrix per treure'ls.
# head(exprs(eset_norm))
# dim(exprs(eset_norm))
# #?subset
# ## Which miRNA affy probeset names are present in the df of miRNA for Humans created in the PreparaDades.A279.R file? (anotacion.hg)
# ## Summary of how many match the condition
# table(rownames(exprs(eset_norm)) %in% as.character(anotacion.hg$Probe.Set.Name)) 
# #
# #FALSE  TRUE 
# #19506  5596
# FALSE  TRUE 
# 19504  5607 
# ## Which ones match the condition? 
# # Get their indexes
# # set_norm.hg.idx <- which(rownames(exprs(eset_norm)) %in% as.character(anotacion.hg$Probe.Set.Name))
# # check that their are the same amount as in the table above
# length(set_norm.hg.idx)
# #[1] 5596
# dim(exprs(eset_norm))
# eset_norm.hg <- exprs(eset_norm)[set_norm.hg.idx,]
# dim(eset_norm.hg)
# class(eset_norm.hg)
# class(eset_norm)
# XXX. BRB279: Això ha funcionat per filtrar per humà però els filtrat el tinc com a matrix, i no com a expression set.

# XXX. BRB279. Repeteixo el process de filtrar d'abans, que sé que es quedan només amb els microRNA d'humans de l'expression set
# i ho torno a assignar al mateix objecte eset_norm, per a continuar amb el Basic Pipe estandard
eset_norm.hg <- eset_norm[rownames(exprs(eset_norm)) %in% as.character(anotacion.hg$Probe.Set.Name),]
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

# I per tant, creo a ma l'objecte exprs.filtered, a partir dels valors d'expressió dels microRNA d'humans obtinguts en el pas anterior
exprs.filtered <- eset_norm.hg.nc


## ----matDesign, eval=TRUE------------------------------------------------
#design<-matrix(
#            c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,
#              0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,
#              0,0,0,0,0,0,0,0,0,0,1,1,1,1,1),
#            nrow=15,byrow=F)
#colnames(design)<-c("A", "B", "L")
#rownames(design) <-  sampleNames 
#print(design)
setwd(baseDir)
targetsFileName <-"targets.BRB279.txt"
targets <- read.table(file = file.path(dataDir, targetsFileName),
                      header = TRUE, sep = "\t")

column2design <- 4   # Columna del ''targets'' en que es basa la matriu de disseny
# En dissenys d'un factor el nombre de grups = nombre de nivells
# En dissenys de més d'un factor nombre de grups = nivells(1)*nivells(2)*...

# XXX Si haguessin efectes batch, o altres efectes a tenir en compte com a factors, caldria afegir noves línies aquí.
lev <- targets[,column2design]  # nom de la columna del targets.txt que conte
# l'agrupacio per nivells (que ha de coincidir amb la de column2design 
# (i.e. no té per que ser $Grupo; veure comentaris més a munt on es defineix el column2design)
# CAL RECORDAR que aquesta forma de definir la matriu e disseny te sentit si es basa en UNA SOLA COLUMNA del targets
# Si s'hagues de basar en un factor fet a partir de dues o més columnes caldria fer-ho a ma (de moment)
# Podeu veure exemples al document "Com Crear Un Projecte basic Pipe"

### design estandar sense efecte batch:
design <- model.matrix( ~ 0 + lev)        
### si hi haguessin efectes batch, o altres efectes a tenir en compte com a factors, caldria afegir noves línies aqui, per exemeple:
#   batch1 <- as.factor(targets$AmplBatch)
#   batch2 <- as.factor(targets$HybrBatch)
# i el design inclouria aquests batchos despres del lev, aixi:
#   design <- model.matrix( ~0 +lev +batch1 +batch2)

#colnames(design) <- levels(lev)
### la definicio de colnames(design) pot ser molt mes complexa si hi ha efectes batch al design, com aquest exemple:
#   colnames(design)<- c(levels(lev),"batch12","batch13","batch14","batch15","batch16","batch17","batch22","batch23")


colnames(design) <- levels(lev)
# rownames(design) <- rownames(targets)
rownames(design) <- targets$ShortName  # Serveix per poder associar les files del targets i les de la matriu de disseny
# Si el targets no te rownames no podrem

numParameters <- ncol(design)

print(design) #comentar aquesta linia si no es vol visualitzar la matriu de disseny

## ----setContrasts, eval=TRUE---------------------------------------------
require(limma)
#cont.matrix <- makeContrasts (
#      AvsB = B-A,
#      AvsL = L-A,
#      BvsL = L-B,
#      levels=design)
#print(cont.matrix)
cont.matrix <- makeContrasts(                   ### RECORDAR QUE AIXÒ ES UN EXEMPLE!!!
  ## EstudiA279b:
  #                                 CANvsCTL = CAN - CTL, # Aquesta es fara mes endavant (columna diferent del targets, etc)
  ## EstudiA279a:
  CTL.RskvsCAN.Epid  = CTL.Rsk  - CAN.Epid,
  CTL.RskvsCAN.Ade   = CTL.Rsk  - CAN.Ade,
  CTL.RskvsCAN.Mic   = CTL.Rsk  - CAN.Mic,
  CTL.noRvsCAN.Epid  = CTL.noR  - CAN.Epid,
  CTL.noRvsCAN.Ade   = CTL.noR  - CAN.Ade,
  CTL.noRvsCAN.Mic   = CTL.noR  - CAN.Mic,
  CAN.EpidvsCAN.Ade  = CAN.Epid - CAN.Ade,
  CAN.EpidvsCAN.Mic  = CAN.Epid - CAN.Mic,
  CAN.AdevsCAN.Mic   = CAN.Ade  - CAN.Mic,
  CTL.RskvsCTL.noR   = CTL.Rsk  - CTL.noR,
  levels = design)

print(cont.matrix) #comentar aquesta linia si no es vol visualitzar la matriu de contrasts



## ----linearmodelfit,echo=F, eval=TRUE------------------------------------
require(limma)
fit<-lmFit(exprs.filtered, design)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)

compName <- c("G1.CTLRisk.vs.CANtype", "G2.CTLnoR.vs.CANtype", 
              "G3.CANTypes", "G4.CTLRisk.vs.CTLnoR") # si son comparacions multiples, fer tant noms com grups de comparacions (N) hi hagi
### 2a part de grups comparacions
### Correspon a "EstudiA279b"
#compName <- c("G5.CAN.vs.CTL")

wCont <- list(1:3, 4:6, 7:9, 10) # Relacionat amb la contrastsMatrix. 
# Llista amb N vectors, que defineixen els N conjunts (grups) de contrastos (comparacions)
# si N>1, cal indicar els rangs per separat
# e.g. list(1:8, 9:13, 14:17)
# Si hi ha només un grup de comparacions (p.e. Estres Termic), i amb dos nivells (estrés versus no estres), aquí es posaria com a:
# list(1:1)

pValCutOff <- c(0.01, 0.01, 0.01, 0.01) # si N>1, indicar el cut-off per cada conjunt de comparacions

# e.g. c(0.01, 0.05, 0.01) o bé c(rep(0.01,3))
# Com a màxim a la UEB es posa 0.25 i a adjMethod posar no ajustat ("none"). 
adjMethod <- c("none", "none", "none", "none")  # si N>1, indicar mètode per cada conjunt de comparacions
# e.g. c("none", "BH", "none") o bé c(rep("BH",3))

## Posar aqui els valors de minLogFC de cadascuna de les comparacions a fer
minLogFoldChange <- c(0, 0, 0, 0) # canviar aixo si s'ha decidit considerar sols els casos amb |logFC| >= que un valor minim
# e.g. c(1, 2, 1) o be c(rep(0, 3)) indicar minLogFC per cada grup de comparacions

## Controlem que el nombre d'elements dels parametres anteriors sigui igual
if(class(wCont)!="list") warning("L'objecte wCont no és una llista! Això pot provocar errors en els passos següents.")
if(length(wCont)!=length(compName)) warning("L'objecte wCont ha de tenir el mateix nombre d'elements que compName!")
if(length(pValCutOff)!=length(compName)) warning("L'objecte pValCutOff ha de tenir el mateix nombre d'elements que compName!")
if(length(adjMethod)!=length(compName)) warning("L'objecte adjMethod ha de tenir el mateix nombre d'elements que compName!")
if(length(minLogFoldChange)!=length(compName)) warning("L'objecte minLogFoldChange ha de tenir el mateix nombre d'elements que compName!")


#############################################################
# TOP TABLES
#############################################################

topTab <- list()
#class(topTab)
#str(topTab)

for (ii in 1:length(wCont)) { # ii is the index of the list with the multiple comparison group names
  #wCont[ii]
  for (jj in 1:length(wCont[[ii]])) { # jj is the index of the list with the single comparisons from within each group of comparisons
    topTab[[ wCont[[ii]][jj] ]] <-  topTable (fit.main, number=nrow(fit.main), 
                                              coef=colnames(cont.matrix)[ wCont[[ii]][jj] ], 
                                              adjust=adjMethod[ii], 
                                              lfc=abs(minLogFoldChange[ii]))
    # head(topTab[[ wCont[[ii]][jj] ]] )

    # Write the resulting topTable to disk
    write.csv2(topTab[[ wCont[[ii]][jj] ]], 
               file=file.path( resultsDir, paste("Selected.Genes.in.comparison.",
               colnames(cont.matrix)[ wCont[[ii]][jj] ], ".csv", sep="")) )
#     library(xtable)
#     print( xtable(topTab[[ wCont[[ii]][jj] ]], type="html"), 
#           file=file.path( resultsDir, paste("Selected.Genes.in.comparison.",
#                      colnames(cont.matrix)[ wCont[[ii]][jj] ],".html", sep="" ))
#           )
    outFile <- paste("Selected.Genes.in.comparison.",
                     colnames(cont.matrix)[ wCont[[ii]][jj] ], sep="" )
    outTitle <- paste("Selected.Genes.in.comparison: ", colnames(cont.matrix)[ wCont[[ii]][jj] ], sep="")
    # For some reason, the html produced doesn't contain the rownames, so we pre-pend them as the first column
    rm(topTab2)
    topTab2 <- cbind(rownames(topTab[[ wCont[[ii]][jj] ]]), 
                topTab[[ wCont[[ii]][jj] ]] )
    colnames(topTab2)[1] <- "featureName"
    #head(topTab2)
#     write.htmltable(x = topTab2, 
#                     file=file.path( resultsDir, outFile ),
#                     title =  outTitle,
#                     open = "wt")
    sortable.html.table(df = topTab2, 
                        output.file = paste0(outFile, "-sortable.html"),
                        output.directory = resultsDir,
                        page.title = outTitle )
    
    #head(topTab)
    #dim(topTab)
  }
}


# ## ----print=FALSE, echo=TRUE, eval=TRUE-----------------------------------
# topTab_AvsB <- topTable (fit.main, number=nrow(fit.main), coef="AvsB", adjust="fdr")
# topTab_AvsL <- topTable (fit.main, number=nrow(fit.main), coef="AvsL", adjust="fdr")
# topTab_BvsL  <- topTable (fit.main, number=nrow(fit.main) , coef="BvsL", adjust="fdr")

setwd(baseDir)
# Carrega numGeneChangedFC.R
source("numGeneChangedFC.R")

setwd(resultsDir)

if(!require(readr)) install.packages("readr")
require(readr)

numGeneChangedFC(filenames=grep("Selected.Genes.in.comparison.*.csv",dir(),value=TRUE),
                 comparisons= colnames(cont.matrix),
                 FC=0)


## ----volcanos, results=tex,echo=FALSE, eval=TRUE-------------------------
for(i in 1:ncol(cont.matrix)){
  compName <-colnames(cont.matrix)[i]
  file=paste("volcanoPlot", compName, ".pdf", sep="")
  pdf(file=file.path(workingDir, "images", file), paper="special", width=6, height=6)
  volcanoplot(fit.main, coef=i, highlight=10, names=fit.main$ID, 
            main=paste("Differentially expressed genes",compName, sep="\n"))
  abline(v=c(-1,1))
  dev.off()
  cat("\\includegraphics{", file, "}\n\n", sep="")
}

## ----CuantosGenes, echo=F, eval=TRUE-------------------------------------
cat("Numero de genes con un p--valor inferior a 0.05 en cada comparación:\n")
  cat ("En la comparación 'A vs B': ", sum(topTab_AvsB$adj.P.Val<=0.05),"\n")
  cat ("En la comparación 'A vs L': ", sum(topTab_AvsL$adj.P.Val<=0.05),"\n")
  cat ("En la comparación 'B vs L': ", sum(topTab_BvsL$adj.P.Val<=0.05),"\n")  


## ----decideTests.2, echo=F, eval=TRUE------------------------------------
res<-decideTests(fit.main, method="separate", adjust.method="fdr", p.value=0.05, lfc=1)
sum.res.rows<-apply(abs(res),1,sum)
res.selected<-res[sum.res.rows!=0,] 
print(summary(res))

## ----venn1,fig=T, eval=TRUE----------------------------------------------
vennDiagram (res.selected[,1:3], main="Genes in common #1", cex=0.9)


## ----htmlPages-----------------------------------------------------------
listOfTables <- list(AvsB = topTab_AvsB, AvsL = topTab_AvsL, BvsL = topTab_BvsL) 
for (i in 1:length(listOfTables)){
  # Seleccionamos la "topTable"
  topTab <- listOfTables[[i]]
  # Escogemos los grupos de sondas a incluir en la tabla
  whichGenes<-topTab["P.Value"]<0.05
  selectedIDs <- topTab$ID[whichGenes]
  # Los convertimos a identificadores Entrez ("EG") y a Gene Symbols
  genes<- getEG(selectedIDs, "hgu133a.db")
  simbols <-getSYMBOL(selectedIDs, "hgu133a.db")
  # Haremos la columna de Entrez sea hiperenlazable
  paraEnlace <- list (misgenes=genes)
  # Preparamos el data.frame con el que se creará el archivo de resultados
  otherNames = data.frame(selectedIDs, simbols, topTab[whichGenes,-1])
  names(otherNames) = c("Affy ID", "Gene Symbol", colnames(topTab)[-1])
  # Invocamos la función "htmlpage"
  comparison <- names(listOfTables)[i]
  htmlpage(paraEnlace, 
           filename =file.path(resultsDir, 
           paste("Selected Genes in comparison ",comparison,".html", sep="")) , 
           title = paste("Diff. expressed genes in comparison ", comparison, sep=""), 
           othernames = otherNames, 
           table.head = c("Entrez IDs", names(otherNames)),
           table.center = TRUE, 
           repository=list("en"))
}

## ----prepareData, eval=TRUE----------------------------------------------
probeNames<-rownames(res)
probeNames.selected<-probeNames[sum.res.rows!=0]
exprs2cluster <-exprs(eset_rma)[probeNames.selected,]
colnames(exprs2cluster)<-sampleNames
color.map <- function(grupo) { 
  if (grupo=="A"){
    c<- "yellow" 
  }else{ 
    if (grupo=="B"){
      c<- "red"
    }else{
      c<- "blue"
   }
  }
return(c)}

## ----plotHeatMap1, fig=T, eval=TRUE--------------------------------------
grupColors <- unlist(lapply(pData(eset_rma)$Group, color.map))
heatmap(exprs2cluster, col=rainbow(100), ColSideColors=grupColors, cexCol=0.9)

## ----plotHeatMap2, fig=T, eval=TRUE--------------------------------------
grupColors <- unlist(lapply(pData(eset_rma)$Group, color.map))
require("gplots")
heatmap.2(exprs2cluster, 
          col=bluered(75), scale="row",
          ColSideColors=grupColors, key=TRUE, symkey=FALSE, 
          density.info="none", trace="none", cexCol=1)


## ----listaArchivos,echo=FALSE,print=FALSE,results=tex, eval=TRUE---------
require(gdata)
listaArchivos <-read.table(file.path(resultsDir,"listaArchivos.txt"), head=TRUE, sep="\t") 
stopifnot(require(xtable))
x.big<-xtable(listaArchivos,
    label='listaArchivos',
    caption='Lista de archivos generados en este análisis')
print(x.big,tabular.environment='longtable',floating=FALSE)

## ----listaArchivos2html,echo=FALSE,print=FALSE, eval=TRUE----------------
require(hwriter)
hwrite(listaArchivos,file.path(resultsDir, "listaArchivos.html"))

