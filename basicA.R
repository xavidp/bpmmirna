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
if(!require(doMC)) install.packages("doMC")
if(!require(devtools)) install.packages("devtools")
if(!require(rCharts)) install_github('rCharts', 'ramnathv')
if(!require(plotly)) install.packages("plotly")
#if(!require(plotly)) devtools::install_github("ropensci/plotly")
if(!require(Nozzle.R1)) install.packages( "Nozzle.R1", type="source" );

#Load required libraries
# "Nozzle: a report generation toolkit for data analysis pipelines"
# http://bioinformatics.oxfordjournals.org/content/early/2013/02/17/bioinformatics.btt085
library(doMC)
require(devtools)
library(rCharts)
require( Nozzle.R1 )

## ----preparaDirectorios, eval=TRUE---------------------------------------
baseDir <- "/home/xavi/Estudis/2015-10-NuriaBarbarroja-IMIBIC-A279"
workingDir <- baseDir
dataDir <-file.path(workingDir, "dades")
resultsDir <- file.path(workingDir,"results")
celfilesDir <- file.path(workingDir,"celfiles")
setwd(workingDir)
## number of cores to use from your computer (with doMC package)
nCores <- 4 # 1 # In case of doubt, use just 1 core.

# Report with Nozzle.R1
report.filename <- "my_report_brba279"
# Remove any previous leftover
if (file.exists(paste0(report.filename, ".html"))) file.remove(paste0(report.filename, ".html"))
if (file.exists(paste0(report.filename, ".RData"))) file.remove(paste0(report.filename, ".RData"))

# Phase 1: create report elements
report.r <- newCustomReport( "My Report BRB A279" );
report.s <- newSection( "Results Generated (Files, Tables and Figures)" );

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

# Report
report.ss1 <- newSection( "Targets" );

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

report.p1 <- newParagraph( "Sample names, groups and color codes used in the analysis" );
report.t1 <- newTable( targets, "Targets file" ); # w/ caption

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
#  CANvsCTL           = CAN - CTL, # Aquesta es fara mes endavant (columna diferent del targets, etc)
  CANvsCTL           = ( (CAN.Epid*9/24) + (CAN.Ade*12/24) + (CAN.Mic*3/24) ) 
                        - ( (CTL.noR*16/24) + (CTL.Rsk*8/24) ), # Aquesta es fara mes endavant (columna diferent del targets, etc)
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

#print(cont.matrix) #comentar aquesta linia si no es vol visualitzar la matriu de contrasts
report.ss2 <- newSection( "Contrasts Matrix" );
report.p2 <- newParagraph( "Contrasts matrix: which sample types (groups) are used in each comparison" );
report.t2 <- newTable( cont.matrix, "Contrasts Matrix" ); # w/ caption

condSplitter <- "vs" # Condition splitter in the Comparison String. This param is used later on to obtain the condition names of the comparison

## ----linearmodelfit,echo=F, eval=TRUE------------------------------------
require(limma)
fit<-lmFit(exprs.filtered, design)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)

compGroupName <- c("G1.CAN.vs.CTL", "G2.CTLRisk.vs.CANtype", "G3.CTLnoR.vs.CANtype", 
              "G4.CANTypes", "G5.CTLRisk.vs.CTLnoR") # si son comparacions multiples, fer tant noms com grups de comparacions (N) hi hagi
### 2a part de grups comparacions
### Correspon a "EstudiA279b"
#compGroupName <- c("G5.CAN.vs.CTL")

wCont <- list(1:1, 2:4, 5:7, 8:10, 11) # Relacionat amb la contrastsMatrix. 
# Llista amb N vectors, que defineixen els N conjunts (grups) de contrastos (comparacions)
# si N>1, cal indicar els rangs per separat
# e.g. list(1:8, 9:13, 14:17)
# Si hi ha només un grup de comparacions (p.e. Estres Termic), i amb dos nivells (estrés versus no estres), aquí es posaria com a:
# list(1:1)

# make a simple list with multiple comparison names and the single comparisons that they contain
if (exists("compNamesAll")) rm(compNamesAll)
compNamesAll <- list()
for (ii in 1:length(wCont)) {
  compNamesAll[[ii]] <- colnames(cont.matrix)[ wCont[[ii]] ]
  names(compNamesAll)[ii] <- compGroupName[[ii]]  
}
report.ss3 <- newSection( "Comparisons performed" );
report.p3 <- newParagraph( "Multiple Comparison groups, with the comparisons they contain each" );
# convert the list of All comparison names into a data frame so that it can be printed easily in the html report with Nozzle
df.compNamesAll <- data.frame(t(sapply(compNamesAll,c)))
report.t3 <- newTable( df.compNamesAll, "Multiple comparison groups" ); # w/ caption

pValCutOff <- rep(0.01, 5) # c(0.01, 0.01, 0.01, 0.01, 0.01) # si N>1, indicar el cut-off per cada conjunt de comparacions

# e.g. c(0.01, 0.05, 0.01) o bé c(rep(0.01,3))
# Com a màxim a la UEB es posa 0.25 i a adjMethod posar no ajustat ("none"). 
adjMethod <- rep("none",5)  # si N>1, indicar mètode per cada conjunt de comparacions
# e.g. c("none", "BH", "none") o bé c(rep("BH",3))

## Posar aqui els valors de minLogFC de cadascuna de les comparacions a fer
minLogFoldChange <- rep(0, 5) #c(0, 0, 0, 0, 0) # canviar aixo si s'ha decidit considerar sols els casos amb |logFC| >= que un valor. Si no, només s'empra per als HeatMaps minim
# e.g. c(1, 2, 1) o be c(rep(0, 3)) indicar minLogFC per cada grup de comparacions

# Place the previous key parameters together in a single df so that they can be easily printed lataer on in the report with Nozzle
key.params <- rbind(pValCutOff, adjMethod, minLogFoldChange)
# Row names need to be added as the first column, in order to be shown in the Report with Nozzle
key.params <- cbind(rownames(key.params), key.params)
colnames(key.params) <- c("Parameter", colnames(df.compNamesAll))

report.ss4 <- newSection( "Filtering" );
report.p4 <- newParagraph( "Filtering was performed with the following parameters");
report.t4 <- newTable( key.params, "Key parameters used");

## Controlem que el nombre d'elements dels parametres anteriors sigui igual
if(class(wCont)!="list") warning("L'objecte wCont no és una llista! Això pot provocar errors en els passos següents.")
if(length(wCont)!=length(compGroupName)) warning("L'objecte wCont ha de tenir el mateix nombre d'elements que compGroupName!")
if(length(pValCutOff)!=length(compGroupName)) warning("L'objecte pValCutOff ha de tenir el mateix nombre d'elements que compGroupName!")
if(length(adjMethod)!=length(compGroupName)) warning("L'objecte adjMethod ha de tenir el mateix nombre d'elements que compGroupName!")
if(length(minLogFoldChange)!=length(compGroupName)) warning("L'objecte minLogFoldChange ha de tenir el mateix nombre d'elements que compGroupName!")

# Set the number of cores to use
registerDoMC(nCores)

#############################################################
# TOP TABLES
#############################################################
if (exists("topTab")) rm(topTab);topTab <- list()
if (exists("topTabExtended")) rm(topTabExtended);topTabExtended <- list()
#class(topTab)
#str(topTab)

# Temporary BD containing the expression values for the selected features
BD <- exprs(exprs.filtered)
if (exists("my.cels.idx")) rm(my.cels.idx); my.cels.idx <- list()
if (exists("my.cels")) rm(my.cels); my.cels <- list()

# Create topTab object in memory
# ----------------------
# Simple loop first out of parallel execution so that object topTab is kept in memory later on further down in the analysis
for (ii in 1:length(wCont)) { # ii is the index of the list with the multiple comparison group names
  #wCont[ii]
  for (jj in 1:length(wCont[[ii]])) { # jj is the index of the list with the single comparisons from within each group of comparisons

    topTab[[ wCont[[ii]][jj] ]] <-  topTable (fit.main, number=nrow(fit.main), 
                                          coef=colnames(cont.matrix)[ wCont[[ii]][jj] ], 
                                          adjust="fdr", 
                                          lfc=0)
  }
}

foreach (ii = 1:length(wCont)) %dopar% { # ii is the index of the list with the multiple comparison group names
  #wCont[ii]
  foreach (jj = 1:length(wCont[[ii]])) %dopar% { # jj is the index of the list with the single comparisons from within each group of comparisons

    my.compName <- colnames(cont.matrix)[ wCont[[ii]][jj] ]
    my.conds <- unlist(strsplit(my.compName, condSplitter)) 
    # my.cels.idx and my.cels below are the list of cel files to be appended to the TopTable object
    # before the csv generation, but avoided when printing the html version
    # We keep the values in a list since we will use this info later when creating heatmaps
    my.cels.idx[[ wCont[[ii]][jj] ]] <- grep( paste0(my.conds[1],"|",my.conds[2]), targets$Grupo)
    my.cels[[ wCont[[ii]][jj] ]] <- as.character(targets$SampleName[ my.cels.idx[[ wCont[[ii]][jj] ]] ] ) 
    write.csv2(cbind(my.cels.idx[[ wCont[[ii]][jj] ]], my.cels[[ wCont[[ii]][jj] ]]), 
               file=file.path( resultsDir, paste("my.cels.in.comparison.",
                                                 my.compName, ".csv", sep="")) )
    
    #head(topTab[[ wCont[[ii]][jj] ]] )
    topTabExtended[[ wCont[[ii]][jj] ]] <- cbind(topTab[[ wCont[[ii]][jj] ]],
                  BD[rownames( topTab[[ wCont[[ii]][jj] ]] ),
                     my.cels[[ wCont[[ii]][jj] ]] 
                     ])

    # Write the resulting topTable to disk
    write.csv2(topTabExtended[[ wCont[[ii]][jj] ]], 
               file=file.path( resultsDir, paste("Selected.Genes.in.comparison.",
               colnames(cont.matrix)[ wCont[[ii]][jj] ], ".csv", sep="")) )

    outFile <- paste("Selected.Genes.in.comparison.",
                     colnames(cont.matrix)[ wCont[[ii]][jj] ], sep="" )
    outTitle <- paste("Selected.Genes.in.comparison: ", colnames(cont.matrix)[ wCont[[ii]][jj] ], sep="")
    # For some reason, the html produced doesn't contain the rownames, so we pre-pend them as the first column
    if (exists("topTab.tmp")) rm(topTab.tmp)
    topTab.tmp <- cbind(rownames(topTab[[ wCont[[ii]][jj] ]]), 
                topTab[[ wCont[[ii]][jj] ]] )
    colnames(topTab.tmp)[1] <- "ID"
    #head(topTab2)

# Disabled in favor of the dTable created below
#     write.htmltable(x = topTab2, 
#                     file=file.path( resultsDir, outFile ),
#                     title =  outTitle,
#                     open = "wt")

# Disabled in favor of the dTable created below
#      sortable.html.table(df = topTab2, 
#                          output.file = paste0(outFile, "-sortable.html"),
#                          output.directory = resultsDir,
#                          page.title = outTitle )

    # When requested, create the dTable, filterable and sortable, etc.
    create.dTable <- T
    if (create.dTable == TRUE) {
      # Create a dTable, a filterable html table: sortable columns plus search box that filster records in real time
      # uses dTable from rCharts.
      # Only the first 7 are the needed ones for the html file generated. Example:
      # (rowname)                         ID      logFC  AveExpr         t      P.Value adj.P.Val  B
      # ENSG00000252190_st ENSG00000252190_st  0.1148340 4.995241  3.976190 0.0001824323 0.3810522  0.6738682
      filterable.dTable=dTable(topTab.tmp[,1:7], sPaginationType = "full_numbers")
      filterable.dTable$templates$script =  "http://timelyportfolio.github.io/rCharts_dataTable/chart_customsort.html" 
      for (cc in 1:7) { # XXXX Only the first 7 are the needed ones for the html file generated
        filterable.dTable$params$table$aoColumns[cc] =
          list( list(sType = "string_ignore_null", sTitle = colnames(topTab.tmp[cc])) )
      }
      filterable.dTable$save(file.path(resultsDir, 
                                       paste0(outFile, "-dTable.html")))    
    } # end of if create.dTable

  }
}


###################################################
## NumGeneChanged
###################################################
setwd(baseDir)
# Carrega numGeneChangedFC.R
source("numGeneChangedFC.R")

setwd(resultsDir)

if(!require(readr)) install.packages("readr")
require(readr)

numGeneChangedFC(filenames=grep("Selected.Genes.in.comparison.*.csv",dir(),value=TRUE),
                 comparisons= colnames(cont.matrix),
                 FC=0) # FC needs to be hardcoded to Zero at this step


###################################################
## Volcano Plots
###################################################
## ----volcanos, results=tex,echo=FALSE, eval=TRUE-------------------------
for (ii in 1:length(wCont)) { # ii is the index of the list with the multiple comparison group names
  #wCont[ii]
  for (jj in 1:length(wCont[[ii]])) { # jj is the index of the list with the single comparisons from within each group of comparisons
  # ii <- 1; jj<- 1
    my.compName <-colnames(cont.matrix)[ wCont[[ii]][jj] ]
    file=paste("volcanoPlot", my.compName, ".pdf", sep="")
    pdf(file=file.path(resultsDir, file), paper="special", width=6, height=6)
    # Set volcanoPointNames. 
    ## Recent versions of limma seem to not write the feature name as fitmai$ID anymore, 
    ## but just provide the feature name as the row.name  
    if (is.null(fit.main$ID)) {
      volcanoPointNames <- rownames(fit.main)
    } else {
      volcanoPointNames <- fit.main$ID
    }
    volcanoplot(fit.main, coef= wCont[[ii]][jj] , highlight=10, names=volcanoPointNames, 
                main=paste("Differentially expressed genes", my.compName, sep="\n"))
    abline(v=c(-1,1))
    dev.off()
    #cat("\\includegraphics{", file, "}\n\n", sep="")
    
  }
}

###################################################
## Venn Diagram
###################################################
if(!require(VennDiagram)) install.packages("VennDiagram")
require(VennDiagram)

# Re-set the needed lists to zero just in case
if (exists("fileVenn")) rm(fileVenn); fileVenn    <- list()
if (exists("listVenn")) rm(listVenn); listVenn    <- list()
if (exists("pValString")) rm(pValString); pValString  <- list()
if (exists("tmpVenn")) rm(tmpVenn)

for (ii in 1:length(wCont)) { # ii is the index of the list with the multiple comparison group names
  #wCont[ii]
  for (jj in 1:length(wCont[[ii]])) { # jj is the index of the list with the single comparisons from within each group of comparisons
    # ii <- 1; jj<- 1
    ## ------------------------------------------------
    ## Seleccio toptables i llistat de genes
    tmpVenn <- topTab[[ wCont[[ii]][jj] ]] 
    fileVenn[[ wCont[[ii]][jj] ]] <- tmpVenn
    #head(tmpVenn)
    #head(exprs(exprs.filtered))
    #targets
    if ( adjMethod[ii] == "none" ) {
      listVenn[[ wCont[[ii]][jj] ]] <- as.character(rownames(tmpVenn[tmpVenn$P.Value < pValCutOff[ii],]))
      pValString[ii]  <- "P.Value"
    } else {
      listVenn[[ wCont[[ii]][jj] ]] <- as.character(rownames(tmpVenn[tmpVenn$adj.P.Val < pValCutOff[ii],]))
      pValString[ii]  <- "Adj.P.Value"
    }
    
  } # end the loop of jj, to have all fileVenn and listVenn created for a multiple comparison

  # In case there are 2 or more comparisons, create a vennDiagram for them.
  if (length(wCont[[ii]]) > 1) {
    
    mainTitle <- paste0("Venn diagram for ", compGroupName[ii]," (", pValString[ii]," < ", pValCutOff[ii], ")") ## Titol
    
    ## Creació Venn Diagram
    venn.plot <- venn.diagram(listVenn[ wCont[[ii]] ], # The list of DE features in each comparison of each multiple comparison group
                              category.names = colnames(cont.matrix)[ wCont[[ii]]  ], ## Comparacions
                              fill = rainbow( length(wCont[[ii]]) ),
                              #fill = c("tomato", "orchid4", "turquoise3"),
                              alpha = 0.50,
                              resolution = 600,
                              cat.cex = 0.9,
                              main = mainTitle,
                              filename = NULL)
    pdf(file.path(resultsDir , paste( "vennDiagram", compGroupName[ii], 
                                      pValString[ii], pValCutOff[ii], "pdf", sep=".")))
    grid.draw(venn.plot)
    dev.off()
  } # end of venn.diagram generation (when appropriate) 

} # end of ii loop, the index of the list with the multiple comparison group names


# ## ----htmlPages-----------------------------------------------------------
# listOfTables <- list(AvsB = topTab_AvsB, AvsL = topTab_AvsL, BvsL = topTab_BvsL) 
# for (i in 1:length(listOfTables)){
#   # Seleccionamos la "topTable"
#   topTab <- listOfTables[[i]]
#   # Escogemos los grupos de sondas a incluir en la tabla
#   whichGenes<-topTab["P.Value"]<0.05
#   selectedIDs <- topTab$ID[whichGenes]
#   # Los convertimos a identificadores Entrez ("EG") y a Gene Symbols
#   genes<- getEG(selectedIDs, "hgu133a.db")
#   simbols <-getSYMBOL(selectedIDs, "hgu133a.db")
#   # Haremos la columna de Entrez sea hiperenlazable
#   paraEnlace <- list (misgenes=genes)
#   # Preparamos el data.frame con el que se creará el archivo de resultados
#   otherNames = data.frame(selectedIDs, simbols, topTab[whichGenes,-1])
#   names(otherNames) = c("Affy ID", "Gene Symbol", colnames(topTab)[-1])
#   # Invocamos la función "htmlpage"
#   comparison <- names(listOfTables)[i]
#   htmlpage(paraEnlace, 
#            filename =file.path(resultsDir, 
#            paste("Selected Genes in comparison ",comparison,".html", sep="")) , 
#            title = paste("Diff. expressed genes in comparison ", comparison, sep=""), 
#            othernames = otherNames, 
#            table.head = c("Entrez IDs", names(otherNames)),
#            table.center = TRUE, 
#            repository=list("en"))
# }

## ----prepareData, eval=TRUE----------------------------------------------
if (exists("exprs2cluster")) rm(exprs2cluster); exprs2cluster <- list()
if (exists("groupColors")) rm(groupColors); groupColors <- list()
if (exists("my.cels.df")) rm(my.cels.df); my.cels.df <- list()
if (exists("my.cels.idx")) rm(my.cels.idx); my.cels.idx <- list()
if (exists("my.cels")) rm(my.cels); my.cels <- list()
#str(listVenn)

for (ii in 1:length(wCont)) { # ii is the index of the list with the multiple comparison group names
  #wCont[ii]
  for (jj in 1:length(wCont[[ii]])) { # jj is the index of the list with the single comparisons from within each group of comparisons
    # ii <- 1; jj<- 1
    my.compName <-colnames(cont.matrix)[ wCont[[ii]][jj] ]
    # read mycels from the stored file on disk, created in the section related to TopTables
    my.cels.df[[ wCont[[ii]][jj] ]] <-read.table(file.path(resultsDir, 
                                                        paste0("my.cels.in.comparison.", my.compName, ".csv")
                                                        ),
                                              head=TRUE, sep=";") 
    my.cels.idx[[ wCont[[ii]][jj] ]] <- as.numeric(as.character(my.cels.df[[ wCont[[ii]][jj] ]] [,2]))
    my.cels[[ wCont[[ii]][jj] ]] <- as.character(my.cels.df[[ wCont[[ii]][jj] ]] [,3])
    
    exprs2cluster[[ wCont[[ii]][jj] ]] <- exprs(exprs.filtered)[listVenn[[ wCont[[ii]][jj] ]] ,
                                           my.cels[[ wCont[[ii]][jj] ]]   ]
    #str(exprs2cluster[[ wCont[[ii]][jj] ]])
  
#     ## ----plotHeatMap1, fig=T, eval=TRUE--------------------------------------
#     groupColors[[ wCont[[ii]][jj] ]] <-  as.character(pData(exprs.filtered)$ColoresTipoT[ my.cels.idx[[ wCont[[ii]][jj] ]] ])
#     heatmap(exprs2cluster[[ wCont[[ii]][jj] ]], col=rainbow(100),
#             ColSideColors=groupColors[[ wCont[[ii]][jj] ]], cexCol=0.9)
#     

    ## ----plotHeatMap2, fig=T, eval=TRUE--------------------------------------
    groupColors[[ wCont[[ii]][jj] ]] <-  as.character(pData(exprs.filtered)$ColoresTipoT[ my.cels.idx[[ wCont[[ii]][jj] ]] ])
    require("gplots")
    mainTitle <- paste0(my.compName) ## Titol
    # Save Heatmap to file on disk
    pdf(file.path(resultsDir , paste( "heatmap", my.compName, 
                                      pValString[ii], pValCutOff[ii], "pdf", sep=".")))
    heatmap.2(exprs2cluster[[ wCont[[ii]][jj] ]], 
              col=bluered(75), scale="row",
              ColSideColors=groupColors[[ wCont[[ii]][jj] ]], key=TRUE, symkey=FALSE, 
              density.info="none", trace="none", cexCol=1, main = mainTitle)    
    
    dev.off()

    ## ----plotHeatMap2, fig=T, eval=TRUE--------------------------------------
    require(plotly)
    # remove.packages("plotly")
    #install.packages("plotly")
    require(scales)
    py <- plot_ly(username='ueb', key='2gfg7ckkuz') # open plotly connection
    # See: https://plot.ly/r/getting-started/ for setting plotly variables for the R session
    Sys.setenv("plotly_username"="ueb") # it shouldn't be needed but it is in my machine "pentinella", for some reason
    Sys.setenv("plotly_api_key"="2gfg7ckkuz") # it shouldn't be needed but it is in my machine "pentinella", for some reason
    mainTitle <- paste0("Heatmap_dynamic_", my.compName) ## Titol
    # Save Heatmap to file online (it can't be saved on disk, it seems, other than in temporary session)
    # Testing Plotly
    # It works in the laptop with VHIR_Externa or MainHead, but not with the desktop (using ethernet cable behind proxy)
    py <- plot_ly(z = exprs2cluster[[ wCont[[ii]][jj] ]],
            x = colnames(exprs2cluster[[ wCont[[ii]][jj] ]]),
            y = rownames(exprs2cluster[[ wCont[[ii]][jj] ]]),
            type = "heatmap")
    plotly_POST(py, filename=mainTitle, world_readable=TRUE)
    # url of the files generated:
    ## Private image:
    # https://plot.ly/~ueb/55.embed (private image, non edited)
    ## Public image:
    #plotly_POST(py, filename=mainTitle, world_readable=TRUE)
    # https://plot.ly/~ueb/52.embed (public image, non edited)

# Success! Created a new plotly here -> https://plot.ly/~ueb/61
# Success! Created a new plotly here -> https://plot.ly/~ueb/63
# Success! Created a new plotly here -> https://plot.ly/~ueb/65
# Success! Created a new plotly here -> https://plot.ly/~ueb/67
# Success! Created a new plotly here -> https://plot.ly/~ueb/69
# Success! Created a new plotly here -> https://plot.ly/~ueb/71
# Success! Created a new plotly here -> https://plot.ly/~ueb/73
# Success! Created a new plotly here -> https://plot.ly/~ueb/75
# Success! Created a new plotly here -> https://plot.ly/~ueb/77
# Success! Created a new plotly here -> https://plot.ly/~ueb/79
# Success! Created a new plotly here -> https://plot.ly/~ueb/81

  } # end the loop of jj
} # end of ii loop, the index of the list with the multiple comparison group names


# ## ----listaArchivos,echo=FALSE,print=FALSE,results=tex, eval=TRUE---------
# require(gdata)
# listaArchivos <-read.table(file.path(resultsDir,"listaArchivos.txt"), head=TRUE, sep="\t") 
# stopifnot(require(xtable))
# x.big<-xtable(listaArchivos,
#     label='listaArchivos',
#     caption='Lista de archivos generados en este análisis')
# print(x.big,tabular.environment='longtable',floating=FALSE)
# 
# ## ----listaArchivos2html,echo=FALSE,print=FALSE, eval=TRUE----------------
# require(hwriter)
# hwrite(listaArchivos,file.path(resultsDir, "listaArchivos.html"))

# Report with Nozzle.R1
# Phase 2: assemble report structure bottom-up
report.ss1 <- addTo( report.ss1, report.p1, report.t1 )
report.ss2 <- addTo( report.ss2, report.p2, report.t2 )
report.ss3 <- addTo( report.ss3, report.p3, report.t3 )
report.ss4 <- addTo( report.ss4, report.p4, report.t4) # parent, child_1, ..., child_n 

report.s <- addTo( report.s, report.ss1, report.ss2, report.ss3, report.ss4 );
report.r <- addTo( report.r, report.s );

# Phase 3: render report to file
writeReport( report.r, filename=report.filename ); # w/o extension
#Two files called my_report.html and my_report.RData will be written to the current working directory.
