###################################################
## BasicA.R
##
## Basic Analysis of Microarrays when the standard 
## "BasicPipe" pipeline from the UEB-VHIR cannot be used
## See git hostory to see contributors
## to this script.
##
## UEB-VHIR (2015) http://ueb.vhir.org 
##
###################################################

#############################
# Package dependencies
#############################
## ----librerias, eval=TRUE------------------------------------------------
installifnot <- function (pckgName){
  if(!(require(pckgName, character.only=TRUE))){
    source("http://Bioconductor.org/biocLite.R")
    biocLite(pckgName)
  }
}
installifnot("Biobase")
installifnot("affy")
installifnot("arrayQualityMetrics")
installifnot("genefilter")
installifnot("limma")
installifnot("annotate")
installifnot("gplots")

if(!require(SortableHTMLTables)) install.packages("SortableHTMLTables")
if(!require(doMC)) install.packages("doMC")
if(!require(devtools)) install.packages("devtools")
if(!require(rCharts)) install_github('rCharts', 'ramnathv')
if(!require(plotly)) install.packages("plotly")
#if(!require(plotly)) devtools::install_github("ropensci/plotly")
if(!require(Nozzle.R1)) install.packages( "Nozzle.R1", type="source" );
if(!require(VennDiagram)) install.packages("VennDiagram")
if(!require(stringr)) install.packages("stringr")
if(!require(xml2)) install.packages("xml2")
if(!require(data.table)) install.packages("data.table") 
#For help type ?data.table or https://github.com/Rdatatable/data.table/wiki
#The fastest way to learn (by data.table authors): https://www.datacamp.com/courses/data-analysis-the-data-table-way

#if(!require(venneuler)) install.packages("venneuler")

#Load required libraries
# "Nozzle: a report generation toolkit for data analysis pipelines"
# http://bioinformatics.oxfordjournals.org/content/early/2013/02/17/bioinformatics.btt085
library(doMC)
require(devtools)
library(rCharts)
require( Nozzle.R1 )
require(stringr)

###################################
# Basic parameters for the script
###################################
analysisName <- "YYYY-MM-XXX-YYY-ANNN"            # Nom del directori del projecte (e.g. "aaaa-mm-NomCognomInvestigador-CENTRE-IdEstudi")
aID <- "LUQ320"                                               # Canviar XXXnnn per les 3 lletres de l'investigador/a + l'Id numeric de l'estudi
## ----preparaDirectorios, eval=TRUE---------------------------------------
baseDir <- paste0("/home/xavi/Estudis/", analysisName) # Pentinella
#baseDir <- "/mnt/magatzem02/tmp/YYYY-MM-XXX-YYY-ANNN" # MainHead
#baseDir <- "."
workingDir <- baseDir
setwd(workingDir)
# Relative Paths to subdirs
dataRelDir <- "dades"
docsRelDir <- "docs"
resultsRelDir <- "results"
celfilesRelDir <- "celfiles"
reportsRelDir <- "reports"
# Absolute Paths to subdirs
dataDir <-file.path(workingDir, dataRelDir)
docsDir <-file.path(workingDir, docsRelDir)
resultsDir <- file.path(workingDir, resultsRelDir)
celfilesDir <- file.path(workingDir, celfilesRelDir)
reportsDir <- file.path(workingDir, reportsRelDir)

# Sourcing the whole set of functions from the standard Basic Pipe: AnalysisFunctions2Pack.R 
source(file.path(baseDir,"Rcode", "AnalysisFunctions2Pack.R"))

## number of cores to use from your computer (with doMC package)
nCores <- 4 # 1 # In case of doubt, use just 1 core.

# Name of the targets file
targetsFileName <- paste0("targets.", aID, ".txt")
column2design <- 6   # Columna del ''targets'' en que es basa la matriu de disseny
# En dissenys d'un factor el nombre de grups = nombre de nivells
# En dissenys de més d'un factor nombre de grups = nivells(1)*nivells(2)*...

datacreacioTargets <- "YYYY-MM-DD1"             # format "aaa-mm-dd"
dataInici <- "YYYY-MM-DD2"

analysisTitle <- "Differentially expressed miRNA between ... and ..."   # Titol de l'estudi (en angles)
analysts <- "Xavier de Pedro, Ferran Brians&oacute; and Alex S&aacute;nchez"                               # Nom analista/es i Alex Sanchez

clientNames <- "Foo Bar"                                         # Nom investigador(a) responsable
lab_or_depart <- "Instituto Foo de Bar (FooBar)"                # Nom del lab o Institucio si es forani (e.g. "Neurovascular diseases")
contact_email <- "foobar@example.com"                                     # e-mail de l'investigador(a)

UEB <- TRUE                                                 # TRUE => capcalera d'UEB
# FALSE => capcalera EstBioinfo

# If some samples were too rare in the PCA from the ArrayQualityMetrics report (even at the raw report), consider
# running the script again removing those from the analysis. You can define which samples to be taken out of the 
# analysis through this parameter.
samples2remove <- "" # c("41.CEL","42.CEL","44.CEL")

# Report with Nozzle.R1
report.filename <- "ResultsFiles"
# Remove any previous leftover
if (file.exists(paste0(report.filename, ".html"))) file.remove(paste0(report.filename, ".html"))
if (file.exists(paste0(report.filename, ".RData"))) file.remove(paste0(report.filename, ".RData"))

# Phase 1: create report elements
report.r <- newCustomReport( "Results Files for Analysis ", aID );
#report.r <- setReportSubTitle( report.r, "Analysis of differentially expressed miRNA between ... cancer and ... control samples from Affymetrix miRNA 4.0 plate arrays");

report.s0a <- newSection( "Overview" );
report.s0a.p1 <- newParagraph( "Results Files, Version 1" );
report.s0a.p2 <- newParagraph( "This html file contains links to a series of files with results generated \
                    in the analyses performed in the Statistics and Bioinformatics Unit.");
report.s0a.p3 <- newParagraph( "Briefly stated the main goal of the study (indeed the first part of the study) \
                    is to select differentially expressed microRNAs associated with different subtypes of lung cancer. \
                    The samples are divided in two main groups: 'Control' with two subgroups \
                    ('With Risk Factor' and 'Without Risk Factor') \
                    and 'Cancer' with three subtypes ('Epidermoid', 'Adenocarcinoma', and 'Microcytic').");
report.s0a.p4 <- newParagraph( "The comparisons performed have been indicated by the researchers.");
report.s0a.p5 <- newParagraph( "Most results are stored as ", asStrong(".txt"), " or ", asStrong(".csv"),
                    ", as well as ", asStrong(".html"), " or ", asStrong(".pdf"), " files\
                    in order to avoid dependency of any specific software different than \
                    general purpose office suits, internet browsers or pdf readers." );
 
report.s0b <- newSection( "Main Documents" );
study.proposalFileName <- file.path(resultsRelDir, "YYYY-MM-XXX-YYY-ANNN-StudyBudget.pdf")
study.proposalRelFileName <- file.path(resultsRelDir, study.proposalFileName)
report.s0b.p1 <- newParagraph( "Study proposal: ", asLink(study.proposalRelFileName, study.proposalRelFileName ));
study.reportFileName <- file.path(resultsRelDir, "YYYY-MM-XXX-YYY-ANNN-MainReport.pdf")
study.reportRelFileName <- file.path(resultsRelDir, study.reportFileName)
report.s0b.p2 <- newParagraph( "Report with description of methodology and main results: ", asEmph("Not available yet")) ;
                               #, asLink(study.reportRelFileName, study.reportRelFileName ));

report.s1 <- newSection( "Base information" );

#############################
# DATA LOADING
#############################
# If the object with normalized data exists, load it. Otherwise, generate it again
if (!exists("data.eset.vsn.rma2.filtered")){
  # Object data.eset.vsn.rma2.filtered doesn't exist in memory
  # load it from disk (if present on disk)
  if (file.exists(file.path(dataRelDir,
                          paste0("data.vsn.rma2.", aID,".Rda")))) {
  ## ----loadData------------------------------------------------------------
  load(file=file.path(dataRelDir, paste0("data.vsn.rma2.", aID,".Rda")))
  } else {
  ## ----generateData------------------------------------------------------------
  source("Rcode/DataPreprocessing.R")
  }
}
#exprs.filtered <- exprs(data.eset.vsn.rma2.filtered)
eset.vsn.rma2.filtered <- data.eset.vsn.rma2.filtered
class(eset.vsn.rma2.filtered)


# Report
report.s1s1 <- newSection( "Samples and Groups" );

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
targets <- read.table(file = file.path(dataDir, targetsFileName),
                      header = TRUE, sep = "\t")

# Keep a copy of targets with all samples for the report (regardless of the fact that we may end up removing some samples from the targets file for further analysis)
targets.all <- targets

# Remove the rows containing the samples indicated in this param "samples2remove" 
# since we do not want to consider them in the analysis
row2remove.idx <- which(targets$SampleName %in% samples2remove)  
# Remove rows only if there is some sample to be removed
if (length(row2remove.idx) > 0) {
  targets <- targets[-row2remove.idx,]
}

report.s1p1 <- newParagraph( "The table below -usually known as the 'Targets' table- shows the classification \
                    of each sample, that is if the sample corresponds to an affected ('Cancer') or \
                    unaffected ('CTL') patient and to which subtype it has been assigned. \
                    The table also contains information about how the sample is labelled (to be shown in plots) \
                    or which colors are used in plots for each Group-Subgroup combination.");
report.s1t1 <- newTable( targets.all, "Targets file" ); # w/ caption
report.s1s1 <- addTo( report.s1s1, report.s1p1, report.s1t1 )


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

# Assign the first column names of the object "design" with the right names as indicated as levels of lev.
# If there is batch effect, there will be extra columns in object "design", that do not neet to be touched probably
# unless there are 2 or more columns with batch effect information
colnames(design)[1:length(levels(lev))] <- levels(lev)
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
cont.matrix <- makeContrasts(
  ### ---------------------------------------------------------------------------
  ### REMEMBER THAT THIS IS JUST AN EXAMPLE!!!
  ### AND THAT IF YOU REMOVE SAMPLES FOR PROCESSING AFTER QUALITY CONTROL
  ### YOU NEED TO RE-ADAPT THE RELATIVE WEIGHT OF EACH CASE IN THE CONTRASTS MATRIX
  ### (I.E. YOU NEED TO TWEAK THE NUMBERS CORRESPONDING TO THE NEW NUMBER OF CEL FILES
  ### FOR EACH CASE IN WHICH YOU REMOVED SAMPLES)
  ### ---------------------------------------------------------------------------
  CANvsCTL           =  ( (CAN.ht25.NoM*5/10) + (CAN.ht25.Met*5/10) + (CAN.lt25.NoM*5/10) + (CAN.lt25.Met*5/10) ) -
                        ( (CTL.ht25.NoM*5/10) + (CTL.ht25.Met*5/10) + (CTL.lt25.NoM*4/10) + (CTL.lt25.Met*6/10) ), 
  CAN.ht25vsCAN.lt25  = ( (CAN.ht25.NoM*5/10) + (CAN.ht25.Met*5/10) ) - ( (CAN.lt25.NoM*5/10) + (CAN.lt25.Met*5/10) ),
  CAN.MetvsCAN.NoM    = ( (CAN.ht25.Met*5/10) + (CAN.lt25.Met*5/10) ) - ( (CAN.ht25.NoM*5/10) + (CAN.lt25.NoM*5/10) ),
  CTL.ht25vsCTL.lt25  = ( (CTL.ht25.NoM*5/10) + (CTL.ht25.Met*5/10) ) - ( (CTL.lt25.NoM*4/10) + (CTL.lt25.Met*6/10) ),
  CTL.MetvsCTL.NoM    = ( (CTL.ht25.Met*5/10) + (CTL.lt25.Met*6/10) ) - ( (CTL.ht25.NoM*5/10) + (CTL.lt25.NoM*4/10) ),
  levels = design)

cont.matrix.file <- "contrasts.matrix.csv"
outFileNameRelPath <- file.path( resultsRelDir, cont.matrix.file )
write.csv2(cont.matrix, file=outFileNameRelPath )

splitterCond        <- "vs" # Condition splitter in the Comparison String
                            # (in CAN.V1.C1vsCTL.V2.C2 it would be "vs" to split the conditions into CAN.V1.C1 and CTL.V2.C2) 
                            # This param is used later on to obtain the condition names of the comparison
splitterIntracond   <- "."  # Intracondition splitter in the Comparison String
                            # (in CAN.V1.C1vsCTL.V2.C2 it would be "." to split the condition CAN.V1.C1 into the vars CAN V1 and C1 
                            #  and CTL.V2.C2 into the vars CTL V2 and C2) 
                            # This param is used later on to obtain the condition names of the comparison

compGroupName <- c("G1.CAN.vs.CTL", "G2.CAN.OBEvsCAN.nOB", "G3.CAN.Met.vs.CAN.NoM", 
                   "G4.CTL.OBEvsCTL.nOB", "G5.CTL.Met.vs.CTL.NoM") # si son comparacions multiples, fer tant noms com grups de comparacions (N) hi hagi


wCont <- list(1:1, 2:2, 3:3, 4:4, 5:5) # Relacionat amb la contrastsMatrix. 
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

#print(cont.matrix) #comentar aquesta linia si no es vol visualitzar la matriu de contrasts
report.s1s2 <- newSection( "Comparisons performed" );
report.s1p2a <- newParagraph( "The researchers have indicated which comparisons they were interested in. \
                            Some comparisons can be naturally grouped, for example because \
                            they are comparing all different cancer subtypes to a given control type, \
                            or because they are doing all possible comparisons between cancer subtypes. \
                            These cases are described as 'Groups of Comparisons'. \
                            This is not only an aestehical distinction because while a series of informations \
                            are provided for each individual comparison, some addititional results \
                            are provided for each group of comparisons (essentially for each group \
                            having more than one comparison)." );

report.s1p2b <- newParagraph( "For each individual comparison the following information is provided:");
report.s1l1 <- newList( isNumbered=FALSE,
                newParagraph( "A 'Top Table' containing the list of all features analyzed, sorted from smallest to greatest P-value." ),
                newParagraph( "A 'Volcano plot' showing how big -or how small- is the biological or the statistical effect." )
                ); # end of list

#   list <- newList( isNumbered=FALSE,
#                 newParagraph( "Nozzle also provides lists." ),
#                 newParagraph( "They can even be nested." ),
#                 newList( isNumbered=TRUE, 
#                          newParagraph( "They may contain other lists." ),
#                          newParagraph( "Or selected other elements." ),
#                          newParagraph( "Numbered." )  
#                          ),
#                 newParagraph( "Or unnumbered." )  
#                 ); # end of list		

report.s1p2c <- newParagraph( "See the text in sections 'Top Table' and 'Volcano Plot' \
                              for a more detailed description of their content and meaning. \
                              Besides the previous tables and plots, for each group of comparisons \
                              it is provided:");
report.s1l2 <- newList( isNumbered=FALSE,
                        newParagraph( "A 'Venn Diagram' showing the number of features that have been selected \
                                      in each comparison (according to a standard cutoff) and \
                                      allowing to determine how many of these features are in common \
                                      between the different comparisons in the same group. \
                                      Obviously these Venn Diagrams are only depicted \
                                      if there is more than one comparison in the group." ),
                        newParagraph( "A 'Heat Map' showing -according to a color code- \
                                      if each feature selected in this group of comparisons \
                                      is up or down regulated in each sample that has been included \
                                      in the comparison." )
                      ); # end of list	

report.s1s2 <- addTo( report.s1s2, report.s1p2a, report.s1p2b, report.s1l1, 
                      report.s1p2c, report.s1l2 )

## ----linearmodelfit,echo=F, eval=TRUE------------------------------------
require(limma)
fit<-lmFit(eset.vsn.rma2.filtered, design)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)

# convert the list of All comparison names into a data frame so that it can be printed easily in the html report with Nozzle
df.compNamesAll <- data.frame(t(sapply(compNamesAll,c)))
report.s1t2 <- newTable( df.compNamesAll, "Groups of comparisons" ); # w/ caption
report.s1s2 <- addTo( report.s1s2, report.s1t2 )

pValCutOff <- rep(0.01, length(compNamesAll) ) # c(0.01, 0.01, 0.01, 0.01, 0.01) # si N>1, indicar el cut-off per cada conjunt de comparacions

# e.g. c(0.01, 0.05, 0.01) o bé c(rep(0.01,3))
# Com a màxim a la UEB es posa 0.25 i a adjMethod posar no ajustat ("none"). 
adjMethod <- rep("none", length(compNamesAll) )  # si N>1, indicar mètode per cada conjunt de comparacions
# e.g. c("none", "BH", "none") o bé c(rep("BH",3))

## Posar aqui els valors de minLogFC de cadascuna de les comparacions a fer
minLogFoldChange <- rep(0, length(compNamesAll) ) #c(0, 0, 0, 0, 0) # canviar aixo si s'ha decidit considerar sols els casos amb |logFC| >= que un valor. Si no, només s'empra per als HeatMaps minim
# e.g. c(1, 2, 1) o be c(rep(0, 3)) indicar minLogFC per cada grup de comparacions

# Place the previous key parameters together in a single df so that they can be easily printed lataer on in the report with Nozzle
key.params <- rbind(pValCutOff, adjMethod, minLogFoldChange)
# Row names need to be added as the first column, in order to be shown in the Report with Nozzle
key.params <- cbind(rownames(key.params), key.params)
colnames(key.params) <- c("Parameter", colnames(df.compNamesAll))

#report.s1s3 <- newSection( "Section removed/refactored (it doesn't exist anymore)" );

#############################
# Quality Control (raw)
#############################
# Section borrowed from the the standard Basic Pipe code.
#rawData <- xx
# or ????
rawData <- rawData2

# phenoData(rawData)$Grupo<-targets$Grupo
# phenoData(rawData)$ShortName<-targets$ShortName
# phenoData(rawData)$Colores<-targets$Colores
# 
# arrayQualityMetrics(expressionset =rawData,
#                     outdir = file.path(resultsDir, "QCDir.raw"),
#                     force = TRUE,
#                     intgroup = "Grupo",
#                     do.logtransform = FALSE)


# Add the section to the report 
report.s1s4 <- newSection( "Quality control (raw data)" );
report.s1p4a1 <- newParagraph( "Different types of quality checks have been performed on the raw data \
                            before deciding that they were valid for the analysis. \
                            Indeed they have been repeated twice: once with raw data and \
                            another with normalized data (next section). A comprehensive report \
                            is provided for the raw data (QCDir.raw/index.html) \
                            and for the normalized data (QCDir.norm/index.html) to help the user to understand \
                            whether a particular array can be considered as an outlier.");
report.s1p4a2 <- newParagraph( "The check consists of a relatively high number of summaries and plots \
                               that can be seen below. Each plot is acompanied by a brief description \
                               of what this means and how it may be interpreted. The quality report \
                               also contains a summary table providing advice on which samples \
                               might be candidates for being considered as outliers and should perhaps \
                               be removed.");
samples2remove.collapsed <- paste(samples2remove, collapse=", ") 
report.s1p4a3 <- newParagraph( "After reviewing the quality checks on raw data it has been decided to remove three samples: ",
                               asStrong( samples2remove.collapsed ));
report.s1h4 <- newHtml( "<iframe src=\"", resultsRelDir, "/QCDir.raw/index.html\" frameborder=1 height=600 scrolling=auto width=\"900\"></iframe>", style="background-color: snow;" )
report.s1s4 <- addTo( report.s1s4, report.s1p4a1, report.s1p4a2, report.s1p4a3, report.s1h4 ) 


#############################
# Quality Control (norm)
#############################
# Section borrowed from the the standard Basic Pipe code.

phenoData(eset.vsn.rma2.filtered)$Grupo<-targets$Grupo
phenoData(eset.vsn.rma2.filtered)$ShortName<-targets$ShortName
phenoData(eset.vsn.rma2.filtered)$Colores<-targets$Colores

# arrayQualityMetrics(expressionset =eset.vsn.rma2.filtered,
#                     outdir = file.path(resultsDir, "QCDir.norm"),
#                     force = TRUE,
#                     intgroup = "Grupo",
#                     do.logtransform = FALSE)

# Add the section to the report 
report.s1s5 <- newSection( "Quality control (normalized data)" );
report.s1p5 <- newParagraph( "A similar quality check has been performed on the samples after they have been normalized. \
                              The results are consistent with those obtained with raw data confirming the decision \
                              to remove samples: ", asStrong( samples2remove.collapsed ));
report.s1h5 <- newHtml( "<iframe src=\"", resultsRelDir, "/QCDir.norm/index.html\" frameborder=1 height=600 scrolling=auto width=\"900\"></iframe>", style="background-color: snow;" )
report.s1s5 <- addTo( report.s1s5, report.s1p5, report.s1h5 ) 

#############################
# Display Feature selection Parameters
#############################
#report.s1s6 <- newSection( "Feature selection" );
#report.s1p6 <- newParagraph( "Feature selection was performed with the following parameters");
#report.s1t6 <- newTable( key.params, "Key parameters used");
#report.s1s6 <- addTo( report.s1s6, report.s1p6, report.s1t6 ) # parent, child_1, ..., child_n 

## Controlem que el nombre d'elements dels parametres anteriors sigui igual
if(class(wCont)!="list") warning("L'objecte wCont no és una llista! Això pot provocar errors en els passos següents.")
if(length(wCont)!=length(compGroupName)) warning("L'objecte wCont ha de tenir el mateix nombre d'elements que compGroupName!")
if(length(pValCutOff)!=length(compGroupName)) warning("L'objecte pValCutOff ha de tenir el mateix nombre d'elements que compGroupName!")
if(length(adjMethod)!=length(compGroupName)) warning("L'objecte adjMethod ha de tenir el mateix nombre d'elements que compGroupName!")
if(length(minLogFoldChange)!=length(compGroupName)) warning("L'objecte minLogFoldChange ha de tenir el mateix nombre d'elements que compGroupName!")

#############################
# Feature annotation reference
#############################
anotacion.hg.fileName.noext <- "featureAnotation" 
outFile <- anotacion.hg.fileName.noext # Name set in DataPreprocessing.R
outFileName <- paste0(outFile, ".csv")
outFileNameRelPath <- file.path( resultsRelDir, outFileName )
  
report.s1s7 <- newSection( "Feature annotation reference" );
report.s1p7 <- newParagraph( "List of Transcript ID for the correspoding Probeset names and database name where to look for more information.");
# Write the resulting files to the report
report.s1file7a.csv <- newHtml( "File (CSV): <a href=", outFileNameRelPath,">",
                                outFileNameRelPath, "</a>",
                                style="background-color: snow;" )

# When requested, create the dTable, filterable and sortable, etc.
create.dTable.featAnot <- T
if (create.dTable.featAnot == TRUE) {
  # Create a dTable, a filterable html table: sortable columns plus search box that filters records in real time
  # uses dTable from rCharts.
  # Only the first 6 are the needed ones for the html file generated. Example:
  #  Probe_Set_ID Species_Scientific_Name Sequence.Type Sequence.Source   Transcript_ID  Probe.Set.Name
  #  20500112            Homo sapiens         miRNA         miRBase   hsa-let-7a-5p MIMAT0000062_st
  featAnot.dTable=dTable(anotacion.hg[,1:ncol(anotacion.hg)], sPaginationType = "full_numbers", 
                         bScrollInfinite = T,
                         bScrollCollapse = T,
                         sScrollY = "300px",
                         aaSorting=list() ) # This list() should indicate to avoid any attempt of sorting client side at display time 
                         # aaSorting=list(c(4, "asc")) ) # This c(4, "asc") sorts ascending on the 5th column. 
  #str(featAnot.dTable)
  # featAnot.dTable$LIB$url
  # featAnot.dTable$html_assets # ...$js & ...$css
  # Disabled tweaking these params here since they change the url for some of the paths, not all (js and css from datatables, for instance, are not converted into relative paths)
  #featAnot.dTable$html_assets$js  <- "foobarjs"
  #featAnot.dTable$html_assets$css <- "foobarcss"
    #     <!doctype HTML>
    #       <meta charset = 'utf-8'>
    #       <html>
    #       <head>
    #       <link rel='stylesheet' href='/home/xavi/R/x86_64-pc-linux-gnu-library/3.2/rCharts/libraries/datatables/css/jquery.dataTables.css'>
    #       <link rel='stylesheet' href='/home/xavi/R/x86_64-pc-linux-gnu-library/3.2/rCharts/libraries/datatables/css/jquery.dataTables_themeroller.css'>
    #       <link rel='stylesheet' href='/home/xavi/R/x86_64-pc-linux-gnu-library/3.2/rCharts/libraries/datatables/css/demo_table.css'>
    #       <link rel='stylesheet' href='foobarjs'>
    #       
    #       <script src='/home/xavi/R/x86_64-pc-linux-gnu-library/3.2/rCharts/libraries/datatables/js/jquery.js' type='text/javascript'></script>
    #       <script src='/home/xavi/R/x86_64-pc-linux-gnu-library/3.2/rCharts/libraries/datatables/js/jquery.dataTables.js' type='text/javascript'></script>
    #       <script src='foobarcss' type='text/javascript'></script>
    #       (...)
    #     
  
  
  # If you want no initial sorting client side because it's been sorted server side already, 
  # you can disable this client initialization with
  # /* Disable initial sort */ From http://stackoverflow.com/a/4964423
  #      "aaSorting": [], which in R code might imply writing:
  #      "aaSorting": list(), which in R code might imply writing:
  #featAnot.dTable$templates$script =  "http://timelyportfolio.github.io/rCharts_dataTable/chart_customsort.html" 
  for (cc in 1:ncol(anotacion.hg)) { # XXXX Only the first 7 are the needed ones for the html file generated
    featAnot.dTable$params$table$aoColumns[cc] =
      list( list(sType = "string_ignore_null", sTitle = colnames(anotacion.hg[cc])) )
  }
  outFileName <- paste0(outFile, "-dTable.html")
  outFileNameRelPath <- file.path( resultsRelDir, outFileName )
  featAnot.dTable$save(outFileNameRelPath) 
} # end of if create.dTable

report.dTable.featAnot <- TRUE
if (report.dTable.featAnot == TRUE) {
  # Assign file names again just in case create.dTable.featAnot was FALSE but report.dTable.featAnot is TRUE
  outFileName <- paste0(outFile, "-dTable.html")
  outFileNameRelPath <- file.path( resultsRelDir, outFileName )
  
  # Write the resulting files to the report
  report.s1s7h1 <- newHtml( "<iframe src=\"", outFileNameRelPath, "\" frameborder=1 height=400 scrolling=auto width=\"1000\"></iframe>", style="background-color: snow;" )
  report.s1file7b <- newHtml( "File (HTML): <a href=\"", outFileNameRelPath,"\">", outFileNameRelPath, "</a>",
                              style="background-color: snow;" )
  report.s1s7 <- addTo( report.s1s7, report.s1p7, report.s1s7h1, report.s1file7a.csv, report.s1file7b)
} else { # if no report of dTable, add at least, the link to the csv file
  report.s1s7 <- addTo( report.s1s7, report.s1p7, report.s1file7a.csv)
} # end of if report dTable


# Set the number of cores to use
registerDoMC(nCores)

#############################################################
# TOP TABLES
#############################################################
# Add it to the report
report.s2 <- newSection( "Results (Files)" );
report.s2s1 <- newSection( "Top Tables" );
report.s2p1 <- newParagraph( "The analysis to select differentially expressed features (here a \"feature\" is a microRNA) \
                             has been based on adjusting a linear model with empirical bayes moderation \
                             of the variance. This is a technique similar to ANOVA which was specifically \
                             developed for microarray data analysis.");
report.s2p2 <- newParagraph( "Each comparison yields a list of features sorted from lowest to highest p-value \
                             which is equivalent to saying from most to least differentially expressed. \
                             This is called generically a '", asEmph("Top Table"), "' and the meaning \
                             of the columns it contains is the following:");
report.s2l1 <- newList( isNumbered=FALSE, 
                        newParagraph( asStrong("Probe.Set.Name"), " is the identifier of each feature in the array."),
                        newParagraph( asStrong("Transcript_ID"), " is the identifier of each feature in the reference database."),
                        newParagraph( asStrong("logFC"), " is the estimated fold change between conditions, which corresponds to the mean difference in log-scale."),
                        newParagraph( asStrong("t"), " is a 'moderated-t' similar to the usual student's t statistic."),
                        newParagraph( asStrong("P.Value"), " is the p-value corresponding to t."),
                        newParagraph( asStrong("adj.P.Val"), " is the adjusted p-value corresponding to t (see below)."),
                        newParagraph( asStrong("B"), " is the B-statistic. It roughly indicates the logarithm of the odds that a feature\
                                      is differentially expressed vs that it is not.\ 
                                      It may be interpreted as follows: When it is greater than 0 it is more likely that the feature \
                                      is differentially expressed than that it is not differentially expressed.")
                        );
report.s2p3 <- newParagraph( "If one wishes to have a statistically grounded criteria, the selection of \
                             the differentially expressed features should be based on adjusted p-values \
                             (less than 0.01) or B statistic (greater than 0).");
report.s2p4 <- newParagraph( "If these criteria yield too few features, table 5 indicate how many features \
                             would be selected by applying a less restrictive criteria such as calling \
                             differentially expressed features with, for instance, adjusted p-values less \
                             than 0.25 or unadjusted p-values smaller than 0.05.");
report.s2p5 <- newParagraph( "In this study the selection of differentially expressed features \
                             has been based on ", asStrong("unadjusted p-values under 0.01"));
report.s2p6 <- newParagraph( "In the sections below the top tables and volcano plots for each comparison and \
                             group of comparisons are depicted." );
report.s2p7 <- newParagraph( "Each section is labelled with two names: the name of the group of comparisons, \
                             and the name of the specific comparison." );
report.s2p8 <- newParagraph( "The tables and plots can be explored interactively on the window or \
                             can be downloaded in different formats for further exploration." );

report.s2s1 <- addTo( report.s2s1, report.s2p1, report.s2p2, report.s2l1, report.s2p3, 
                      report.s2p4, report.s2p5, report.s2p6, report.s2p7, report.s2p8)

if (exists("topTab")) rm(topTab);topTab <- list()
if (exists("topTabExtended")) rm(topTabExtended);topTabExtended <- list()
#class(topTab)
#str(topTab)

# Temporary BD containing the expression values for the selected features
BD <- exprs(eset.vsn.rma2.filtered)
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

topTabLoop <- foreach (ii = 1:length(wCont)) %do% { # ii is the index of the list with the multiple comparison group names
  #wCont[ii] ; ii <- 1; jj <- 1;
  foreach (jj = 1:length(wCont[[ii]])) %do% { # jj is the index of the list with the single comparisons from within each group of comparisons
    
    my.compName <- colnames(cont.matrix)[ wCont[[ii]][jj] ]
    my.conds <- unlist(strsplit(my.compName, splitterCond, fixed = TRUE)) 
    my.conds.units <- unlist(strsplit(my.conds, splitterIntracond, fixed = TRUE)) 
    
    # my.cels.idx and my.cels below are the list of cel files to be appended to the TopTable object
    # before the csv generation, but avoided when printing the html version
    # We keep the values in a list since we will use this info later when creating heatmaps
#    my.cels.idx[[ wCont[[ii]][jj] ]] <- grep( paste0(my.conds, collapse="|"), targets$Grupo)
    my.cels.idx[[ wCont[[ii]][jj] ]] <- grep( paste0(my.conds.units, collapse="|"), targets$Grupo)
    my.cels[[ wCont[[ii]][jj] ]] <- as.character(targets$SampleName[ my.cels.idx[[ wCont[[ii]][jj] ]] ] ) 
    write.csv2(cbind(my.cels.idx[[ wCont[[ii]][jj] ]], my.cels[[ wCont[[ii]][jj] ]]), 
               file=file.path( resultsDir, paste("celfiles.in.comparison.",
                                                 my.compName, ".csv", sep="")) )
    
    #head(topTab[[ wCont[[ii]][jj] ]] )
    topTabExtended[[ wCont[[ii]][jj] ]] <- cbind(topTab[[ wCont[[ii]][jj] ]],
                                                 BD[rownames( topTab[[ wCont[[ii]][jj] ]] ),
                                                    my.cels[[ wCont[[ii]][jj] ]] 
                                                    ])
    
    # Write the resulting topTable to disk
    outFile <- paste("Selected.Features.in.comparison.",
                     colnames(cont.matrix)[ wCont[[ii]][jj] ], sep="")
    outFileName <- paste0(outFile, ".csv")
    outFileNameRelPath <- file.path( resultsRelDir, outFileName )
    
    
    write.csv2(topTabExtended[[ wCont[[ii]][jj] ]], 
               file=outFileNameRelPath )
    # Write the resulting files to the report
    report.s2file1a.csv <- newHtml( "File (CSV): <a href=", outFileNameRelPath,">",
                                outFileNameRelPath, "</a>",
                                style="background-color: snow;" )

    
    outTitle <- paste("Selected.Features.in.comparison: ", colnames(cont.matrix)[ wCont[[ii]][jj] ], sep="")
    # For some reason, the html produced doesn't contain the rownames, so we pre-pend them as the first column
    if (exists("topTab.tmp")) rm(topTab.tmp)
    topTab.tmp <- cbind(rownames(topTab[[ wCont[[ii]][jj] ]]), 
                        topTab[[ wCont[[ii]][jj] ]] )
    colnames(topTab.tmp)[1] <- "Probe.Set.Name"
    #head(topTab.tmp)
    # Add column TranscriptID at this point, merging through the common column Probe.Set.Name
    #head(anotacion.hg)
    tid <- data.frame(anotacion.hg$Transcript_ID,anotacion.hg$Probe.Set.Name)
    colnames(tid) <- str_replace(colnames(tid), "anotacion.hg.", "")
    #head(tid)
    topTab.tmp <- merge(topTab.tmp, tid, by = "Probe.Set.Name")
    # head(topTab.tmp.tid)
    # dim(topTab.tmp)
    # dim(topTab.tmp.tid)

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
    create.dTable.topTab <- T
    if (create.dTable.topTab == TRUE) {
      ## Sort topTab.tmp before creating the dTable
      # No need to re-sort since it's already sorted by p-value
      #head(topTab.tmp[,1:ncol(topTab.tmp)])
      
      # Create a dTable, a filterable html table: sortable columns plus search box that filster records in real time
      # uses dTable from rCharts.
      # Only the first 8 are the needed ones for the html file generated. Example:
      # (rowname)                         Probe.Set.Name      logFC  AveExpr         t      P.Value adj.P.Val  B Transcript_ID
      # ENSG00000252190_st ENSG00000252190_st  0.1148340 4.995241  3.976190 0.0001824323 0.3810522  0.6738682   foobar
      filterable.dTable=dTable(topTab.tmp[,1:ncol(topTab.tmp)], sPaginationType = "full_numbers", 
                              # aaSorting=list() ) # This list() should indicate to avoid any attempt of sorting client side at display time 
                               aaSorting=list(c(4, "asc")) ) # This c(4, "asc") sorts ascending on the 5th column. 
                              # If you want no initial sorting client side because it's been sorted server side already, 
                              # you can disable this client initialization with
                              # /* Disable initial sort */ From http://stackoverflow.com/a/4964423
                              #      "aaSorting": [], which in R code might imply writing:
                              #      "aaSorting": list(), which in R code might imply writing:
      filterable.dTable$templates$script =  "http://timelyportfolio.github.io/rCharts_dataTable/chart_customsort.html" 
      for (cc in 1:ncol(topTab.tmp)) { # XXXX Only the first 8 are the needed ones for the html file generated
        filterable.dTable$params$table$aoColumns[cc] =
          list( list(sType = "string_ignore_null", sTitle = colnames(topTab.tmp[cc])) )
      }
      outFileName <- paste0(outFile, "-dTable.html")
      outFileNameRelPath <- file.path( resultsRelDir, outFileName )
      filterable.dTable$save(outFileNameRelPath) 
    } # end of if create.dTable
    
    report.dTable.topTab <- TRUE
    if (report.dTable.topTab == TRUE) {
      outFileName <- paste0(outFile, "-dTable.html")
      outFileNameRelPath <- file.path( resultsRelDir, outFileName )
      
      # Create a new subsection (hidden by default) to display the html as an iframe
      report.s2s1s1 <- newSection( names(compNamesAll)[ii], " | ", compNamesAll[[ii]][jj] );
      
      # Write the resulting topTable files to the report
      report.s2s1h1 <- newHtml( "<iframe src=\"", outFileNameRelPath, "\" frameborder=1 height=450 scrolling=auto width=\"1000\"></iframe>", style="background-color: snow;" )
      report.s2file1b <- newHtml( "File (HTML): <a href=\"", outFileNameRelPath,"\">", outFileNameRelPath, "</a>",
                                  style="background-color: snow;" )
      report.s2s1s1 <- addTo( report.s2s1s1, report.s2s1h1, report.s2file1a.csv, report.s2file1b)
      report.s2s1 <- addTo( report.s2s1, report.s2s1s1)
    } else { # if no report of dTable, add at least, the link to the csv file
      report.s2s1 <- addTo( report.s2s1, report.s2file1a.csv)
    } # end of if report dTable
    
  } # end of the jj loop
  return( c(report.s2s1) )  # only needed when the loop run in parallel, 
  # to preserve some objects created within the parallel loop for later reuse
  #return( c(topTabExtended) )
} # end of the ii loop

#str(topTabLoop)
#topTabLoop[[1]]$elements[[2]]$text
#rm(report.s2s1)
#report.s2s1

###################################################
## NumFeatureChanged (formarly named NumGeneChanged )
###################################################
# Add it to the report
report.s2s2 <- newSection( "Number of features changed in each case" );
report.s2p2 <- newParagraph( "Summary table with the number of features changed in each case for each pValue Type (adjusted or not) and Cutoff");
report.s2s2 <- addTo( report.s2s2, report.s2p2)

setwd(baseDir)
# Carrega numGeneChangedFC.R
source(file.path(baseDir,"Rcode", "numFeatureChangedFC.R"))

setwd(resultsDir)

if(!require(readr)) install.packages("readr")
require(readr)

numFeatureChangedFilenames <- paste0("Selected.Features.in.comparison.", 
                                  colnames(cont.matrix), 
                                  ".csv")
# Code derived from Miriam's code
numFeatureChangedFC(filenames=numFeatureChangedFilenames,
                 comparisons= colnames(cont.matrix),
                 FC=0) # FC needs to be hardcoded to Zero at this step

# ---
# Forma alternativa d'obtenir l'objecte numGenesChanged, a partir dels objectes de topTable en memoria,
# per comprovar que coincideix amb la formad e la Miriam de lapply a partir d'arxius en disk de genes de cada comparacio
# Conclusió: un cop arreglat el DataPreprocessing.R, i siguen coherents en la forma de procedir, 
# coincideixen els dos mètodes. Xavier. 18 Des'15.
# Codi de l'Alex per fer numGenesChanged - Pendent de concloure si cal canviar-ho o no
# https://github.com/alexsanchezpla/scripts/blob/master/numGenesChanged.R
#source("https://raw.githubusercontent.com/alexsanchezpla/scripts/master/numGenesChanged.R")
# # Code derived from Alex's code
# # It can be aggregated with a loop over all comparisons in colnames(cont.matrix)
# numFeatureChanged.df <- data.frame(matrix(nrow = 12))[,-1]
# for (ii in 1:ncol(cont.matrix)) {
#   colnames(cont.matrix)[ii]
#   numFeatureChanged.tmp <- numGenesChanged(topTab[[ii]],
#                                                  colnames(cont.matrix)[ii] )
#   numFeatureChanged.df <- cbind(numFeatureChanged.df, numFeatureChanged.tmp)
# }
# # Write it to disk
# write.csv2( numFeatureChanged.df, file=paste0("numFeaturesChangedAlex.", aID, ".csv") )
# ---


#class(numFeatureChanged.df); class(numFeatureChanged.tmp)
# head(numFeatureChanged.df)
outFileName <- paste("numFeaturesChangedFC0.csv",sep="")
outFileNameRelPath <- file.path( resultsRelDir, outFileName )
numFeatureChangedFC.df <- read.table(file = file.path(resultsDir, outFileName),
                                  header = TRUE, sep = ";")
# Write the resulting files to the report
report.s2t2a <- newTable( numFeatureChangedFC.df[,1:6], 
                         file=outFileNameRelPath,
                         "Number of features changed between comparisons for given p.value cutoffs and methods (adjusted p.value or not). First comparisons." ); # w/ caption
# report.s2t2b <- newTable( numFeatureChangedFC.df[,c(1,7:ncol(numFeatureChangedFC.df))], 
#                           file=outFileNameRelPath,
#                           "Number of features changed between comparisons for given p.value cutoffs and methods (adjusted p.value or not). Last comparisons." ); # w/ caption
# report.s2s2 <- addTo( report.s2s2, report.s2t2a, report.s2t2b)
report.s2s2 <- addTo( report.s2s2, report.s2t2a)
#report.s2file <- newParagraph( "File: <a href=\"", outFileNameRelPath,"\">",
#                                 outFileNameRelPath, "</a>")
#report.s2s2 <- addTo( report.s2s2, report.s2file)

###################################################
## Volcano Plots
###################################################
# Add it to the report
report.s2s3 <- newSection( "Volcano plots" );
report.s2p3 <- newParagraph( "Files with the volcano plot for each comparison");
report.s2s3 <- addTo( report.s2s3, report.s2p3)

## ----volcanos, results=tex,echo=FALSE, eval=TRUE-------------------------
for (ii in 1:length(wCont)) { # ii is the index of the list with the multiple comparison group names
  #wCont[ii]
  for (jj in 1:length(wCont[[ii]])) { # jj is the index of the list with the single comparisons from within each group of comparisons
    # ii <- 1; jj<- 1
    my.compName <-colnames(cont.matrix)[ wCont[[ii]][jj] ]
    # Compose the filenames
    outFileName.noext <- paste("volcanoPlot", my.compName, sep="")
    outFileNameRelPath.noext <- file.path( resultsRelDir, outFileName.noext )
    outFileName.pdf <- paste0(outFileName.noext, ".pdf")
    outFileName.png <- paste0(outFileName.noext, ".png")
    outFileNameRelPath.pdf <- file.path( resultsRelDir, outFileName.pdf )
    outFileNameRelPath.png <- file.path( resultsRelDir, outFileName.png )
    
    # Set volcanoPointNames. 
    ## Recent versions of limma seem to not write the feature name as fitmai$ID (Probe.Set.Name in this case) anymore, 
    ## but just provide the feature name as the row.name  
    if (is.null(fit.main$Probe.Set.Name)) {
      volcanoPointNames <- rownames(fit.main)
    } else {
      volcanoPointNames <- fit.main$Probe.Set.Name
    }

    # Rename volcanoPointNames to use TranscriptID instead of Probe.Set.Name. 
    # E.g.: "14q0" instead of "14q0_st". Or "hsa-let-7a-5p" instead of "MIMAT0000062_st"  
    #head(volcanoPointNames)
    #head(anotacion.hg)
    df1 <- data.frame(volcanoPointNames)
    colnames(df1) <-"Probe.Set.Name"
    df2 <- anotacion.hg
    volcanoPointNames_df <- merge(df1,df2, all.x=TRUE)
    volcanoPointNames_df <- volcanoPointNames_df$Transcript_ID
    #length(volcanoPointNames_df)
    
#     # ---------------------------------
#     # Attempt to reproduce the same behavior with package data.tables instead, which is k own to be more efficient
#     # However, I end up commenting out this chunk of code since ther results where not identical in some cases, 
#       due to the sorting performed internally in the data.table way, which seems to be needed for it's improvement in performance
#       and I don't know yet how to force to keep the same order in the first results, so I end up using the first methods which keeps 
#       the original order we need.
#     # ---------------------------------
#     require(data.table)
#     dt1 <- data.table(volcanoPointNames)
#     colnames(dt1) <-"Probe.Set.Name"
#     dt2 <- data.table(anotacion.hg)
#     #?merge.data.table
#     setkey(dt1,"Probe.Set.Name")
#     setkey(dt2,"Probe.Set.Name")
#     volcanoPointNames_dt <- dt2[dt1]
#     dim(volcanoPointNames_dt)
#     volcanoPointNames_dtid <- volcanoPointNames_dt$Transcript_ID
#     volcanoPointNames_dtidpn <- data.table(volcanoPointNames_dt$Probe.Set.Name, volcanoPointNames_dt$Transcript_ID)
#     colnames(volcanoPointNames_dtidpn) <- c("Probe.Set.Name", "Transcript_ID")
#     
#     table(volcanoPointNames_df == volcanoPointNames_dtid)
#     diff.idx <- which(volcanoPointNames_df != volcanoPointNames_dtid)
#     volcanoPointNames_df[diff.idx[1]-1]
#     volcanoPointNames_df[diff.idx[1:2]]
#     volcanoPointNames_df[c(diff.idx[1]-1, diff.idx[1:2])];
# 
#     volcanoPointNames_dtid[diff.idx[1]-1];
#     volcanoPointNames_dtid[diff.idx[1:2]];
#     volcanoPointNames_dtidpn[c(diff.idx[1]-1, diff.idx[1:2])];
#     
#     anotacion.hg[anotacion.hg$Transcript_ID=="ACA59",]
#     volcanoPointNames[grep("ACA59",volcanoPointNames)]
#     str(volcanoPointNames_df);str(volcanoPointNames_dtid);
#     head(volcanoPointNames_df);head(volcanoPointNames_dtid);
#     tail(volcanoPointNames_df);tail(volcanoPointNames_dtid);
#     table(sort(volcanoPointNames_df) == sort(volcanoPointNames_dtid))

                                
    # Generate the pdf
    pdf(file=outFileName.pdf, paper="special", width=6, height=6)
    volcanoplot(fit.main, coef= wCont[[ii]][jj] , highlight=10, names=volcanoPointNames_df, 
                main=paste("Differentially expressed features in ", my.compName, sep="\n"))
    abline(v=c(-1,1))
    dev.off()
    #cat("\\includegraphics{", file, "}\n\n", sep="")
    # Generate the png
    png(file=outFileName.png)
    volcanoplot(fit.main, coef= wCont[[ii]][jj] , highlight=10, names=volcanoPointNames_df, 
                main=paste("Differentially expressed features in ", my.compName, sep="\n"))
    abline(v=c(-1,1))
    dev.off()
    
    # Add it to the report as screenshot
    # figure file paths
    report.s2f3a <- newFigure( outFileNameRelPath.png, 
                               fileHighRes=outFileNameRelPath.pdf,
                               "Differentially expressed features for the comparison ", my.compName );
    report.s2s3 <- addTo( report.s2s3, report.s2f3a)

    # Write the resulting files to the report
    report.s2file <- newHtml( "File (PDF): <a href=\"", outFileNameRelPath.pdf,"\">",
                              outFileNameRelPath.pdf, "</a>",
                              style="background-color: snow;" )
    report.s2s3 <- addTo( report.s2s3, report.s2file)
  }
}

###################################################
## Venn Diagram
###################################################
# Add it to the report
report.s2s4 <- newSection( "Venn Diagrams" );
report.s2p4 <- newParagraph( "Files with the Venn Diagrams for the groups of comparisons (cases of more than one comparison in each group)");
report.s2s4 <- addTo( report.s2s4, report.s2p4)

if(!require(VennDiagram)) install.packages("VennDiagram")
require(VennDiagram)

# Re-set the needed lists to zero just in case
if (exists("listVenn.tid")) rm(listVenn.tid); listVenn.tid    <- list() # tid = Transcript_ID
if (exists("listVenn.psn")) rm(listVenn.psn); listVenn.psn    <- list() # psn = Probe.Set.Name
if (exists("pValString")) rm(pValString); pValString  <- list()
if (exists("tmpVenn")) rm(tmpVenn)

for (ii in 1:length(wCont)) { # ii is the index of the list with the multiple comparison group names
  #wCont[ii]
  for (jj in 1:length(wCont[[ii]])) { # jj is the index of the list with the single comparisons from within each group of comparisons
    # ii <- 1; jj<- 1
    ## ------------------------------------------------
    ## Seleccio toptables i llistat de Features
    tmpVenn <- topTab[[ wCont[[ii]][jj] ]] 
    tmpVenn <- cbind(rownames(tmpVenn), tmpVenn); colnames(tmpVenn)[1] <- "Probe.Set.Name"
    tmpVenn.tid <- merge(tmpVenn, tid, by = "Probe.Set.Name")
    # Reorder results by p.value again
    tmpVenn.tid <-  tmpVenn.tid[order(tmpVenn.tid$P.Value),]
    # Re-set rownames to Transcript_ID this time
    rownames(tmpVenn.tid) <- tmpVenn.tid$Transcript_ID
    
  
    # ListVenn objects:
    #listVenn.tid stands for Transcript_ID in the list for the VennDiagrams
    #listVenn.psn stands for Probe.Set.Name in the list for the VennDiagrams
    if ( adjMethod[ii] == "none" ) {
      listVenn.tid[[ wCont[[ii]][jj] ]] <- as.character(rownames(tmpVenn.tid[tmpVenn.tid$P.Value < pValCutOff[ii],]))
      listVenn.psn[[ wCont[[ii]][jj] ]] <- as.character(tmpVenn.tid[tmpVenn.tid$P.Value < pValCutOff[ii], "Probe.Set.Name"])
      pValString[ii]  <- "P.Value"
    } else {
      listVenn.tid[[ wCont[[ii]][jj] ]] <- as.character(rownames(tmpVenn.tid[tmpVenn.tid$adj.P.Val < pValCutOff[ii],]))
      listVenn.psn[[ wCont[[ii]][jj] ]] <- as.character(tmpVenn.tid[tmpVenn.tid$adj.P.Val < pValCutOff[ii], "Probe.Set.Name"])
      pValString[ii]  <- "Adj.P.Value"
    }
    
  } # end the loop of jj, to have all listVenn.tid created for a multiple comparison
  
  # In case there are 2 or more comparisons, create a vennDiagram for them.
  if (length(wCont[[ii]]) > 1) {
    
    mainTitle <- paste0("Venn diagram for ", compGroupName[ii]," (", pValString[ii]," < ", pValCutOff[ii], ")") ## Titol
    
    ## Creació Venn Diagram
    venn.plot <- venn.diagram(listVenn.tid[ wCont[[ii]] ], # The list of DE features in each comparison of each multiple comparison group
                              category.names = colnames(cont.matrix)[ wCont[[ii]]  ], ## Comparacions
                              fill = rainbow( length(wCont[[ii]]) ),
                              #fill = c("tomato", "orchid4", "turquoise3"),
                              alpha = 0.50,
                              resolution = 600,
                              cat.cex = 0.9,
                              main = mainTitle,
                              filename = NULL)
    # Compose the filenames
    outFileName.noext <- paste( "vennDiagram", compGroupName[ii], 
                          pValString[ii], pValCutOff[ii], sep=".")
    outFileNameRelPath.noext <- file.path( resultsRelDir, outFileName.noext )
    outFileName.pdf <- paste0(outFileName.noext, ".pdf")
    outFileName.png <- paste0(outFileName.noext, ".png")
    outFileNameRelPath.pdf <- file.path( resultsRelDir, outFileName.pdf )
    outFileNameRelPath.png <- file.path( resultsRelDir, outFileName.png )
    # Generate the pdf
    pdf(outFileName.pdf)
    grid.draw(venn.plot)
    dev.off()
    # Generate the png
    png(outFileName.png)
    grid.draw(venn.plot)
    dev.off()
    
    # Add it to the report as screenshot
    # figure file paths
    report.s2f4a <- newFigure( outFileNameRelPath.png, 
                               fileHighRes=outFileNameRelPath.pdf,
                               "Venn Diagram for the comparison group ", compGroupName[ii] );
    report.s2s4 <- addTo( report.s2s4, report.s2f4a)
    
    # Write the resulting files to the report
    report.s2file <- newHtml( "File (PDF): <a href=\"", outFileNameRelPath.pdf,"\">",
                              outFileNameRelPath.pdf, "</a>",
                              style="background-color: snow;" )
    report.s2s4 <- addTo( report.s2s4, report.s2file)
    
    ############################
    ## Save on disk the list of Features for each group of multiple comparisons
    ## Derived from 
    ## https://github.com/miriamMota/scripts/blob/master/Bioinf/VennDiagram.R
    ############################
    xx.1 <- listVenn.tid[ wCont[[ii]] ]
    names(xx.1) <- compNamesAll[[ii]]
    combs <-  unlist(lapply(1:length(xx.1), 
                            function(j) combn(names(xx.1), j, simplify = FALSE)),
                     recursive = FALSE)
    names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
    #str(combs)
    
    elements <- lapply(combs, function(i) Setdiff(xx.1[i], xx.1[setdiff(names(xx.1), i)]))
    n.elements <- sapply(elements, length)
    list_res <- list(elements= elements, n.elements=n.elements) 
    
    # Create a df with the list of features present in each area of the venn diagram sections
    seq.max <- seq_len(max(n.elements))
    mat <- sapply(elements, "[", i = seq.max)
    mat[is.na(mat)]<-""
    mat <- data.frame(mat)
    
    # Add the first line with the number of elements in each group of features
    # check that col names are in the same order before merging the rows from the 2 objects
    table(colnames(t(data.frame(n.elements)[1])) == colnames(mat)) 
    # merge the values of the two objects and store them inthe same mat object
    mat <- rbind( t(data.frame(n.elements)[1]) , mat)
    
    # Compose the filename
    outFileName <- paste( "vennDiagram", compGroupName[ii],
                          pValString[ii], pValCutOff[ii], "csv", sep=".")
    outFileNameRelPath <- file.path( resultsRelDir, outFileName )
    
    # Write results to disk
    write.csv2(mat, outFileName, row.names=FALSE)
    
    # Write the resulting files to the report
    report.s2file <- newHtml( "File (CSV): <a href=\"", outFileNameRelPath,"\">",
                              outFileNameRelPath, "</a>",
                              style="background-color: snow;" )
    report.s2s4 <- addTo( report.s2s4, report.s2file)
    
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

###################################################
## Heatmaps
###################################################
# Add it to the report
report.s2s5 <- newSection( "Heatmap plots" );
report.s2p5 <- newParagraph( "Files with the Heatmap plots for each comparison");
report.s2s5 <- addTo( report.s2s5, report.s2p5)

## ----prepareData, eval=TRUE----------------------------------------------
if (exists("exprs2cluster")) rm(exprs2cluster); exprs2cluster <- list()
if (exists("groupColors")) rm(groupColors); groupColors <- list()
if (exists("my.cels.df")) rm(my.cels.df); my.cels.df <- list()
if (exists("my.cels.idx")) rm(my.cels.idx); my.cels.idx <- list()
if (exists("my.cels")) rm(my.cels); my.cels <- list()
if (exists("heatmap.plotly.link")) rm(heatmap.plotly.link); heatmap.plotly.link <- list()

#str(listVenn.tid)

for (ii in 1:length(wCont)) { # ii is the index of the list with the multiple comparison group names
  #wCont[ii]
  for (jj in 1:length(wCont[[ii]])) { # jj is the index of the list with the single comparisons from within each group of comparisons
    # ii <- 1; jj<- 1
    my.compName <-colnames(cont.matrix)[ wCont[[ii]][jj] ]
    # read mycels from the stored file on disk, created in the section related to TopTables
    my.cels.df[[ wCont[[ii]][jj] ]] <-read.table(file.path(resultsDir, 
                                                           paste0("celfiles.in.comparison.", my.compName, ".csv")
    ),
    head=TRUE, sep=";") 
    my.cels.idx[[ wCont[[ii]][jj] ]] <- as.numeric(as.character(my.cels.df[[ wCont[[ii]][jj] ]] [,2]))
    my.cels[[ wCont[[ii]][jj] ]] <- as.character(my.cels.df[[ wCont[[ii]][jj] ]] [,3])
    
    # head(exprs(eset.vsn.rma2.filtered))
    exprs2cluster[[ wCont[[ii]][jj] ]] <- exprs(eset.vsn.rma2.filtered)[listVenn.psn[[ wCont[[ii]][jj] ]] ,
                                                                my.cels[[ wCont[[ii]][jj] ]]   ]
    #str(exprs2cluster[[ wCont[[ii]][jj] ]])
    #class(exprs2cluster[[ wCont[[ii]][jj] ]])
    #head(exprs2cluster[[ wCont[[ii]][jj] ]])
    exprs2cluster.tid <- data.frame(exprs2cluster[[ wCont[[ii]][jj] ]], rownames(exprs2cluster[[ wCont[[ii]][jj] ]]) )
    colnames(exprs2cluster.tid) <- c(colnames(exprs2cluster[[ wCont[[ii]][jj] ]]), "Probe.Set.Name") 
    #head(exprs2cluster.tid)
    #head(tid)
    rownames(tid) <- tid$Probe.Set.Name
    #mx <- exprs2cluster[[ wCont[[ii]][jj] ]]
    mx <- exprs2cluster.tid
    my <- tid
    # Merge both objects
    mxy <- merge(mx, my, by = intersect(names(mx), names(my)))
    head(mxy)
    # reassign rownames, but this time use Transcript_ID
    rownames(mxy) <- mxy$Transcript_ID
    # Remove columns of Transcript_ID and Probe.Set.Name after we confirmed by hand that the merge process did as expected 
    mxy.tid <- mxy[!names(mxy) %in% c("Probe.Set.Name", "Transcript_ID")]
    #head(mxy.tid)
    
    groupColors[[ wCont[[ii]][jj] ]] <-  as.character(pData(eset.vsn.rma2.filtered)$Colores[ my.cels.idx[[ wCont[[ii]][jj] ]] ])
    mainTitle <- paste0(my.compName) ## Titol
    # Compose the filename
    outFileName.noext <- paste( "heatmap", my.compName, 
                                pValString[ii], pValCutOff[ii], sep=".")
    outFileNameRelPath.noext <- file.path( resultsRelDir, outFileName.noext )
    
    #     ## ----plotHeatMap1, fig=T, eval=TRUE--------------------------------------
    # Save Heatmap to file on disk as PDF
    pdf(paste0(outFileName.noext, ".pdf"))
    heatmap(as.matrix(mxy.tid), col=bluered(75),
            ColSideColors=groupColors[[ wCont[[ii]][jj] ]], cexCol=0.9,
            main = mainTitle, xlab="Samples", ylab="")
            # main = mainTitle, xlab="Samples", ylab="Features")
    dev.off()

    #     ## ----plotHeatMap1, fig=T, eval=TRUE--------------------------------------
    # Save Heatmap to file on disk as PNG
    png(paste0(outFileName.noext, ".png"))
    heatmap(as.matrix(mxy.tid), col=bluered(75),
            ColSideColors=groupColors[[ wCont[[ii]][jj] ]], cexCol=0.9,
            main = mainTitle, xlab="Samples", ylab="")
            #main = mainTitle, xlab="Samples", ylab="Features")
    dev.off()
    
#     ## ----plotHeatMap2, fig=T, eval=TRUE--------------------------------------
#      require("gplots")
#      heatmap.2(as.matrix(mxy.tid), 
#                col=bluered(75), scale="row",
#                ColSideColors=groupColors[[ wCont[[ii]][jj] ]], key=TRUE, symkey=FALSE, 
#                density.info="none", trace="none", cexCol=1, main = mainTitle)    
# 
#     dev.off()
    
    # Add it to the report as screenshot
    # figure file paths
    figureFile1 <- paste0(outFileNameRelPath.noext, ".png");
    figureFileHighRes1a <- paste0(outFileNameRelPath.noext, ".pdf");
    report.s2f5a <- newFigure( figureFile1, fileHighRes=figureFileHighRes1a,
                               "Heatmap plot for the comparison ", my.compName );
    report.s2s5 <- addTo( report.s2s5, report.s2f5a)


    # Write the resulting files to the report
    report.s2file <- newHtml( "File (PDF): <a href=\"", outFileNameRelPath.noext,".pdf\">",
                              outFileNameRelPath.noext, ".pdf</a>",
                              style="background-color: snow;" )
    report.s2s5 <- addTo( report.s2s5, report.s2file)

    # When requested, create the plotly heatmap locally as png
    plotly.heatmaps.create <- FALSE # Disabled to avoid re-genearing them each time the script is rerun
    # When requested, create the dynamic plotly heatmap posted in the plot.ly server
    plotly.heatmaps.post <- FALSE # Disabled to avoid re-genearing them each time the script is rerun
    # When requested, report the dynamic plotly heatmap posted in the plot.ly server
    plotly.heatmaps.report <- FALSE # Disabled to avoid re-genearing them each time the script is rerun
    
    if (plotly.heatmaps.create == TRUE) {
      ## ----plotHeatMap2, fig=T, eval=TRUE--------------------------------------
      require(plotly)
      # remove.packages("plotly")
      #install.packages("plotly")
      require(scales)
      py <- plot_ly(username='ueb', key='2gfg7ckkuz') # open plotly connection
      # See: https://plot.ly/r/getting-started/ for setting plotly variables for the R session
      Sys.setenv("plotly_username"="ueb") # it shouldn't be needed but it is in my machine "pentinella", for some reason
      Sys.setenv("plotly_api_key"="2gfg7ckkuz") # it shouldn't be needed but it is in my machine "pentinella", for some reason
      mainTitle <- paste0("Heatmap_plotly", my.compName) ## Titol
      # Save Heatmap to file online (it can't be saved on disk, it seems, other than in temporary session)
      # Testing Plotly
      # It works in the laptop with VHIR_Externa or MainHead, but not with the desktop (using ethernet cable behind proxy)
      #Axis and legend default names seem to be the name of the variable containing them. Therefore we rename the vars here:
      Expressions <- mxy.tid
      Samples <- colnames(mxy.tid)
      Features <- rownames(mxy.tid)
      py <- plot_ly(z = Expressions,
                    x = Samples,
                    y = Features,
                    color=bluered(8),
                    type = "heatmap")
      # Create a new local plot.ly heatmap as png
      Png <- plotly_IMAGE(py, out_file = file.path(resultsDir, paste0(mainTitle, ".png")))
      
      # POST (create) a new plot.ly heatmap in plot.ly server only when requested
      if (plotly.heatmaps.post == TRUE) {
        plotly_POST(py, filename=mainTitle, world_readable=TRUE)
      }
      
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
      heatmap.plotly.link[[ wCont[[ii]][jj] ]] <- paste0("https://plot.ly/~ueb/",
                                                         seq(61, 81, by = 2)[ wCont[[ii]][jj] ],
                                                         ".embed") # 
    } # end of chunk to create plotly heatmaps when requested
    
    if (plotly.heatmaps.report == TRUE) {
      # Add it to the report
      # figure file paths
      figureFile1 <- file.path(resultsRelDir, paste0(mainTitle, ".png"));
      #figureFileHighRes1b <- "https://plot.ly/~ueb/61.embed";
      figureFileHighRes1b <- heatmap.plotly.link[[ wCont[[ii]][jj] ]]
      report.s2f5b <- newFigure( figureFile1, fileHighRes=figureFileHighRes1b, 
                                 "Dynamic version of the Heatmap for the comparison ", my.compName, 
                                 ". Hovering and zooming are possible." );
      
      #report.s2p6 <- newParagraph( "Heatmaps produced. See ", asReference( report.s2f5a ), "for instance");
      report.s2s5 <- addTo( report.s2s5, report.s2f5b)
    } # end of plotly.heatmaps.report
    
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

###################################################
## Acknowledgement
###################################################
# Add the common request to cite our service in their papers at the Acknowledgements section
report.s10 <- newSection( "Other considerations" );

report.s10.p1 <- newParagraph( "To evaluate the usefulness of our work we will be pleased if you inform us about the scientific products as proceedings, papers, thesis, etc) derived from this work. \
We will be very glad if you include in the Acknowledgements section of the paper the following sentence:" );

report.s10.p2 <- newParagraph( asCode( "\"Statistical analysis ( Bioinformatics ) has been carried out\
                                       in the Statistics and Bioinformatics Unit (UEB) at \
                                       Vall d'Hebron Research Institute (VHIR) - http://ueb.vhir.org \"") );

report.s10 <- addTo( report.s10, report.s10.p1, report.s10.p2);

###################################################
## Make the report
###################################################
# Report with Nozzle.R1
# Phase 2: assemble report structure bottom-up
report.s0a <- addTo( report.s0a, report.s0a.p1, report.s0a.p2, report.s0a.p3, report.s0a.p4, report.s0a.p5);
report.s0b <- addTo( report.s0b, report.s0b.p1, report.s0b.p2);
#report.s1 <- addTo( report.s1, report.s1s1, report.s1s2, report.s1s3, report.s1s4, report.s1s5,
#                    report.s1s6, report.s1s7);
report.s1 <- addTo( report.s1, report.s1s1, report.s1s2, report.s1s4, report.s1s5,
                    report.s1s7); # sections s1s3 and s1s6 have been removed from the report by Alex.
report.s2 <- addTo( report.s2, report.s2s1, report.s2s2, report.s2s3, report.s2s4, report.s2s5 );
report.r <- addTo( report.r, report.s0a, report.s0b, report.s1, report.s2 );
report.r <- addTo( report.r, report.s10 ); # "Other considerations" section, at the end.

# Ensure that the report is created at the baseDir, and not at resultsDir
setwd(baseDir)

# Settings
# set report maintainer information
report.r <- setMaintainerName( report.r, "UEB - VHIR" );
report.r <- setMaintainerEmail( report.r, "ueb@vhir.org" );
report.r <- setMaintainerAffiliation( report.r, "Statistics and Bioinformatics Unit - Vall d'Hebron Research Institute" );

# set the copyright notice for this report
report.r <- setCopyright( report.r, owner="UEB - VHIR", year=2015, statement="Some rights reserved.", url="http://ueb.vhir.org" ); 

# set contact information for error reports
report.r <- setContactInformation( report.r, email="ueb@vhir.org", subject="Problem with this Report", message="Hello!\n\nPlease describe the issue here.", label="Report an Issue" );

#report.r <- setCustomScreenCss( report.r, "paper.css" );

# Phase 3: render report to file
writeReport( report.r, filename=report.filename ); # w/o extension
#Two files called my_report.html and my_report.RData will be written to the current working directory.

# Phase 4: edit the report and other html files produced so that css and js libraries use relative paths and not absolute, so that they work when moved to another computer

# For html edition as simple text from R, we can use the base functions "readLines" and "cat" (to file)
base_url <- "/home/xavi/R/x86_64-pc-linux-gnu-library/3.2/rCharts/libraries/datatables/"
# Copy the css and js files from the reports folder to the results folder for later reuse within the ResultsFiles.html generated 
#file.link(from, to)
# Create folders for js and css in results Dir if they don't exist yet
folders2create <- c("js", "css")
files2create <- list(c("jquery.js", "jquery.dataTables.js"),
                     c("demo_table.css", "jquery.dataTables_themeroller.css", "jquery.dataTables.css"))
names(files2create) <- folders2create
# Start loop over folders
for (folder in folders2create) {
  # Check if folder exists. If not, create it.
  if (!dir.exists(file.path(resultsDir, folder))) dir.create(file.path(resultsDir, folder))
  
  # Start loop over files within each folder
  for (f2c in files2create[[folder]]) { # f2c = file to create
    # prepare to copy files if they do not exist yet in destination
    to    <- file.path(resultsDir, folder, f2c)
    if (!file.exists(to)) {
      from  <- file.path(reportsDir, folder, f2c)
      resultsFileLink <- file.link(from, to)
    } # end if clause for file copy 
  } # loop over files
} # loop over folders

# to be compressed and send elsewhere for the researcher with all the required css and js files in it. Copyright licenses should be respected (citation, etc).
filenames <- list.files(path = resultsDir, pattern = "html")
setwd(resultsDir)
## Replace base_url with nothing since these js and css folders have been copied to the resultsDir folder.
#foreach( ff in filenames ){
replaceResult <- foreach (ff = 1:length(filenames)) %dopar% {
  #ff <- 1
  xx <- readLines(filenames[ff])
  yy <- gsub( base_url, "", xx )
  cat(yy, file=filenames[ff], sep="\n")
  
}


# copy over the results folder the other main documents produced elsewhere
files2create <- c(study.proposalFileName, study.reportFileName)
for (f2c in files2create) { # f2c = file to create
  # prepare to copy files if they do not exist yet in destination
  to    <- file.path(resultsDir, f2c)
  if (!file.exists(to) && file.exists(from)) {
    from  <- file.path(docsDir, f2c)
    resultsFileLink <- file.link(from, to)
  } # end if clause for file copy 
} # loop over files

#   # For html edition from R, See: http://blog.rstudio.org/2015/04/21/xml2/
#   require(xml2)
#   base_url <- "/home/xavi/R/x86_64-pc-linux-gnu-library/3.2/rCharts/libraries/datatables/"
#   #file2edit <- file.path(resultsRelDir, paste0(report.filename, ".html"))
#   file2edit <- paste0(report.filename, ".html")
#   #f2e <- read_html(file2edit)
#   # xml_url(f2e)
#   # xml_contents(f2e)
#   # url_relative(f2r, basde_url)
