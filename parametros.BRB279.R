options(width = 170)
###################################################
### chunk number 1: variablesEstudi
###################################################
###line 31 "ComCrearUnProjecteParametres.rnw"

datacreacioTargets <- "2015-11-02"             # format "aaa-mm-dd"
         dataInici <- "2015-11-02"

   unitDisk <- "G"                                                    # Windows unit
windowsPath <- paste(unitDisk,"dades/Estudis/Microarrays", sep=":")   # Windows path

 linuxPath1 <- "/mydocs/UEB/estudis/microarrays/gene"                 # Linux path (e.g. "/media/TREKSTOR/Estudis/Microarrays")
 linuxPath2 <- "/home/xavi/Estudis"
  linuxPath <- linuxPath2

     estudi <- "2015-10-NuriaBarbarroja-IMIBIC-A279"            # Nom del directori del projecte (e.g. "aaaa-mm-NomCognomInvestigador-CENTRE-IdEstudi")
         eP <- "BRB279"                                               # Canviar XXXnnn per les 3 lletres de l'investigador/a + l'Id numeric de l'estudi

  titolEstudi <- "Differentially expressed miRNA between lung cancer and control samples"   # Titol de l'estudi (en angles)
    analistes <- "Xavier de Pedro and Alex S&aacute;nchez"                               # Nom analista/es i Alex Sanchez

   nomClients <- "Nuria Barbarroja Puerto"                                         # Nom investigador(a) responsable
 lab_o_depart <- "Instituto Maimónides de Inv. Biom. de Córdoba (IMIBIC)"                # Nom del lab o Institucio si es forani (e.g. "Neurovascular diseases")
contact_email <- "nuria.barbarroja.exts@juntadeandalucia.es"                                     # e-mail de l'investigador(a)

          UEB <- TRUE                                                 # TRUE => capcalera d'UEB
                                                                      # FALSE => capcalera EstBioinfo

### AIXO S'HAURIA DE POSAR AL CHUNK 10 (el de la design matrix)
column2design <- 4   # Columna del ''targets'' en que es basa la matriu de disseny
                     # En dissenys d'un factor el nombre de grups = nombre de nivells
                     # En dissenys de més d'un factor nombre de grups = nivells(1)*nivells(2)*...


  mylib.BioC1 <- "/mydocs/R/Bioconductor/mylibs/BioC-2.7"             
  mylib.BioC2 <- NULL
        mylib <- mylib.BioC2      # Ruta dels paquets instal·lats (p.e. "/mydocs/R/Bioconductor/mylibs/BioC")
      mylib.R <- .libPaths()[1]   # Ruta dels paquets de R instal·lats (p.e."/mydocs/R/packages")


###################################################
### chunk number 2: directories
###################################################
###line 57 "ComCrearUnProjecteParametres.rnw"

## Actualització 2011-1: 
## S'ha simplificat aqueta secció. Ara tot penja del workingDir

SO <- version[["os"]]

if (SO=="linux-gnu")
{
  workingDir <- file.path(linuxPath, estudi)
}else{
#  memory.limit(size = 4095)
  workingDir <- file.path(windowsPath, estudi) 
}

    dataDir <- file.path(workingDir, "dades")
celFilesDir <- file.path(workingDir, "celfiles")
 resultsDir <- file.path(workingDir, "results")
 reportsDir <- file.path(workingDir, "informes")
    tempDir <- file.path(workingDir, "temp")
    codeDir <- file.path(workingDir, "RCode")  

setwd(workingDir)



###################################################
### chunk number 3: annotations
###################################################
###line 100 "ComCrearUnProjecteParametres.rnw"

rawDataType <- "affyGeneArray"     # Els tipus previstos son: "text", "affy3primeIVT", "affyGeneArray", "affyExonArray", "illuminaBeadArray"

exonStudy <- TRUE                  # TRUE => estudi d'exons o genearrays

chipPackAvailable <- FALSE                  # si es TRUE s'ha de definir chipPackage (opcio per si els arrays son 3prime)
platformDesignPackAvailable <- TRUE         # si es TRUE cal definir platformDesignPackage (opcio per si son exon o genearrays)
textAnnotationsPackageAvailable <- FALSE    # si es TRUE definir textAnnotationsPackage (quan no sigui ni chip ni platform)

cdfPackage <- NULL   		   # Nomes si chipPackAvailable està TRUE
                                   # En els altres casos deixar-lo a NULL

# URL on cercar noms de paquets de bioconductor:
# http://www.bioconductor.org/packages/release/data/annotation/
chipPackage <- NULL                            # si chipPackAvailable es FALSE deixar-ho a NULL. Té format com bovine.db, que aquí es posaria com a "bovine"
textAnnotationsPackage <- NULL                 # si textAnnotationsPackageAvailable es FALSE deixar-ho a NULL
				    

  entrezTableFileName <- "Entrezs.Rda"
 symbolsTableFileName <- "Symbols.Rda"
controlsTableFileName <- "controls.Rda"
  annotationsFileName <- "Annotations"

## Els següents parametres si es deixen en NULL s'omplen sols quan s'analitzen Genearrays o arrays d'exons. Per la resta d'arrays encara s'ha de fer manual. Si et vols "saltar" parts del pipeline, s'han d'omplir manualment.

platformDesignPackage <- NULL		        # "pd.ragene.1.1.st.v1"
orgPackage <- NULL          	                # nom del paquet corresponent a l'organisme, sense .db al final. "org.Rn.eg"
specie <- NULL		                        # XXX CANVIAR EL NOM DE L'ESPECIE SI CAL
organisme <- NULL                               # CANVIAR L'ABREVIATURA DE L'ORGANISME SI CAL! 
					        #(e.g: mmu = Mus musculus; rno = Rattus norvegicus;



###################################################
### chunk number 4: experimentInformation
###################################################
###line 151 "ComCrearUnProjecteParametres.rnw"

htmlInfo <- list(To = nomClients,
                 Description = titolEstudi,
                 Analysts = analistes,
                 Contact = "Alex Sanchez (alex.sanchez@vhir.org)")



###################################################
### chunk number 5: FileNames
###################################################
###line 189 "ComCrearUnProjecteParametres.rnw"

targetsFileName <- paste("targets", eP, "txt", sep=".")                      # Nom arxiu targets
#targetsFileName <- paste("targets", eP, "chipsRemoved", "txt", sep=".")      # Nom arxiu targets sense els xips que s'han eliminat

  MainReportName <- "aaaa-mm-NomCognomInvestigador-CENTRE-999-MainReport.pdf"                                   # Nom de l'arxiu MainReport
  StudyProposal <- "aaaa-mm-NomCognomInvestigador-CENTRE-999-StudyProposal.pdf"        # Nom de l'arxiu de la proposta

    rawDataFileName <- "rawData.Rda"                          # Nom de l'arxiu de dades crues (llegides del .CEL)
expres.all.FileName <- "expres.Rda"                           # Nom de l'arxiu on es desen les expressions
      linksFileName <- paste("Links", eP, "txt", sep=".")     # Nom de l'arxiu on es desen els enllaços


###################################################
### chunk number 6: QCgrups
###################################################
###line 217 "ComCrearUnProjecteParametres.rnw"

perParts <- FALSE
 QCgrups <- NULL


###################################################
### chunk number 7: QCFileNames
###################################################
###line 260 "ComCrearUnProjecteParametres.rnw"

## Actualització 2011-4
## Si els arrays de l'estudi són d'exons no farem servir els controls de qualitat 
## implementats al Basic Pipe sino el paquet arrayQualityMetrics
## Aquesta elecció es pren al fer l'anàlisi tot i que es defineix aquí.

if(exonStudy)
{
  QCDir <- file.path(resultsDir, "QCDir")      # Per passar-ho al paquet arrayQualityMetrics
}else{
  QCDir <- file.path(resultsDir, "QCDir")      # Per passar-ho al paquet arrayQualityMetrics
  SSQCFileName <- "SSQCPlots.pdf"
  affyQCFileName <- "QCReport.pdf"
  PLMFileName <- "PLMQCPlots.pdf"
}

separatePlots <- TRUE


###################################################
### chunk number 8: normalizationParameters
###################################################
###line 302 "ComCrearUnProjecteParametres.rnw"

### Actualització 2011-5: S'ha afegit la possibilitat de fixar el mètode
### de normalització segons si treballem amb arrays d'exons o 3'IVT
### i, en el segon cas, segons si sumaritzem a nivell d'exó o de transcrit.
### Caldria que els exons tinguessin un nou paràmetre summarizationLevel?

if(exonStudy)
{  ### En el cas d'exons assumim RMA, i normMethod indica quin tipus de normalització (a nivell d'exon, de gen...)
  normMethod <- "core" ### Possibles Valors:
                       ###  -> a nivell d'exon: "probeset"
                       ###  -> a nivell de gen:     "core" : Valor per defecte. Només s'agafen les anotacions més fiables
                       ###                      "extended"
                       ###                          "full"
}else{  
  normMethod <- "RMA"
}

 normalized.eset.FileName <- "normalizedData.Rda"   
  normalized.all.FileName <- "normalized.all"
normalized.Plots.FileName <- "NormalizedPlots.pdf"

fileType <- "csv2"

 PCAFile <- "PCAVariability"
PCAPlots <- TRUE


###################################################
### chunk number 9: filterParameters
###################################################
###line 401 "ComCrearUnProjecteParametres.rnw"

### Actualització 2011-6: PENDENT- 
### Cal trobar el prefix dels transcrits de control per treure'ls en el filtratge

### Actualització 2011-7 (pendent)
### La definició dels grups per filtratge de senyal es fa per defecte 
### fent servir una columna del targets que s'ha de dir "grupo"
### Convindria que en algun punt això es defineixi a l'estil del \Rcode{column2design}
### es a dir es digui quina columna es fa servir per definir els grups i punt.
## Ja de pas recordem que quan hi ha més de d'un factor basar el filtratge en un d'ells no té massa sentit.

SignalFilter <- TRUE
signalThreshold <- 50
signalThreshold.as.percentage <- TRUE
targets <- read.table(file = file.path(dataDir, targetsFileName), header = TRUE, sep = "\t")
grups.For.SignalFilter <- targets["Grupo"]
signalFilter.Function <- "filt.by.Signal"

VarFilter <- TRUE
variabilityThreshold <- 50
variabilityThreshold.as.percentage <- TRUE
variability.Function <- "sdf"

pairing.Function <- NULL

expres.filtered.FileName <- "expres.filtered.Rda"
normalized.filtered.FileName <- "normalized.filtered"

### Actualització 2011-8: S'ha afegit la localització explicita de l'arxiu 
### amb els valors de filtratge (L'arxiu de filtratge ja existia pero no es mostrava)

FilteringReportFileName <- paste("FilteringReport", eP, "txt", sep = ".")


###################################################
### chunk number 10: genericMatrices
###################################################


###################################################
### chunk number 11: designMatrix
###################################################
###line 523 "ComCrearUnProjecteParametres.rnw"

if (!(require("limma", lib = mylib, character.only = TRUE))){
  biocLite("limma")
  require("limma", lib = mylib, character.only = TRUE)
}
      
targets <- read.table(file = file.path(dataDir, targetsFileName),
                      header = TRUE, sep = "\t")

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


###################################################
### chunk number 12: setContrasts
###################################################
###line 539 "ComCrearUnProjecteParametres.rnw"

if (!(require("limma", lib = mylib, character.only = TRUE))){
  biocLite("limma")
  require("limma", lib = mylib, character.only = TRUE)
}

### Matriu de contrasts a construir a partir de la proposta d'estudi i el targets

### A la funcio makeContrasts escriurem les comparacions com
###
###                  AvsD = A - D
###
### on  AvsD : es el nom de la comparacio
###    A - D : es la diferencia a fer entre condicions

contrastsMatrix <- makeContrasts(AvsD = A - D,  ### RECORDAR QUE AIXÒ ES UN EXEMPLE!!!
                                 BvsD = B - D,
                                 CvsD = C - D,
                                 DvsE = D - E,
                                 levels = design)

print(contrastsMatrix) #comentar aquesta linia si no es vol visualitzar la matriu de contrasts


###################################################
### chunk number 13: paremetersLists
###################################################
###line 564 "ComCrearUnProjecteParametres.rnw"

   lmParsList <- list()
   mcParsList <- list()
clustParsList <- list()
   GOParsList <- list()
 KEGGParsList <- list()

source(file.path(codeDir, "AnalysisFunctions2Pack.R"))


###################################################
### chunk number 14: lmParameters
###################################################
###line 598 "ComCrearUnProjecteParametres.rnw"

comparName <- "Estudi" # XXX Aquest paràmetre es pot editar per ajudar a reconèixer l'estudi en alguns dels pdf's i gràfics. Cal que sigui breu i sense espais.

Estudi <- list(dades = NULL,
               expresFileName = expres.filtered.FileName,
               targets = targets,
               designMat = design,
               contMat = contrastsMatrix,
               whichContrasts = 1:ncol(contrastsMatrix),
               anotPack = NULL,
               outputDir = resultsDir,
               ExpressionsAndTop = TRUE,
               showLmParams = TRUE, # posar a TRUE si es vol incloure els camps derivats als ExpressAndTop.csv
               use.dupCorr = FALSE,
               block = NULL,
               nDups = 1,
               comparisonName = comparName,  
               ENTREZIDs = "entrezTable",
               SYMBOLIDs = "symbolsTable",
               fileOfLinks = linksFileName,
               fitFileName = paste("fit", comparName, "Rda", sep = "."),
               csvType=fileType,
               rows2HTML = NULL,
               anotFilename = annotationsFileName
               )

lmParsList <- add2parsList(lmParsList, Estudi)


###################################################
### Chunk 15: Gene Annotations
###################################################
###no apareix a "ComCrearUnProjecteParametres.rnw"?

p <- Estudi

### Parametres:
###
###            specie : Nom de cientific de l'especie. Aquest parametre es necessari per a crear els links
###                     en els repositoris de la KEGG i de ENSEMBL. En l'actualitat nomes es contempla les
###                     seguents especies per a la KEGG: "Homo_sapiens", "Mus_musculus", "Rattus_norvegicus",
###                     "Bos_taurus", "Danio_rerio", "Sus_scrofa". El parametre per defecte esta a NULL
###         info2show : Tipus d'anotació que es vol extreure. En l'actualitat es pot extreure la seguent informacio
###
###                           "Affymetrix", "EntrezGene", "GenBank","KEGG", "GO", "RefSeq", "UniGene",
###                           "Ensembl", "UniProt", "PubMed", GeneSymbol", "GeneName", "Cytoband", "Enzyme"
###                     Per defecte es mostren "Affymetrix", "EntrezGene"
###
###   numGenesPerPage : Valor màxim de registres mostrats a cada fitxer d'anotacions.
###          maxGenes : Numero maxim de gens que es mostraran per arxiu HTML d'anotacions. En cas de superar aquest
###                     nombre es crea una nova pagina amb els seguents (i.e. si maxGenes = 100 i tenim 1000 gens a
###                     anotar creara 10 pagines). Per defecte pren maxGenes = NULL i mostra tots els gens

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### NOTA (feta el 31.05.2011)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Per ara les funcions associades a les anotacions no estan pensades per
### a fer un ParsList ja que nomes genera un full d'anotacions comu
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

anotList <- list(fitMain = NULL,
                 fitFileName = p$fitFileName,
                 my.IDs = p$ENTREZIDs,
                 anotPackage = orgPackage,
                 toHTML = TRUE,
                 outputDir = p$outputDir,
                 anotFilename = p$anotFilename,
                 titleAnotations = "Annotations for all genes analyzed",
		 specie = specie,
                 info2show = c( "Affymetrix", "EntrezGene", "GeneSymbol", "GeneName", "KEGG", "GO"), #"PubMed"), # XXXX Columnes a generar en el fitxer d'anotacions. Com menys columnes, més ràpid anirà l'estudi
                 linksFile = p$fileOfLinks,
                 numGenesPerPage = NULL
                 )


###################################################
### chunk number 16: mcParametersLoop
###################################################
###line 695 "ComCrearUnProjecteParametres.rnw"

p <- Estudi

## Definim els noms dels grups de comparacions i els contrastos associats a cada grup:
compName <- c("Group1", "Group2", "Group3") # si son comparacions multiples, fer tant noms com grups de comparacions (N) hi hagi
                        # e.g. c("Dose", "Cells", "Time")

wCont <- list(1:3, 4:6, 7:8) # Relacionat amb la contrastsMatrix. 
                   # Llista amb N vectors, que defineixen els N conjunts (grups) de contrastos (comparacions)
                   # si N>1, cal indicar els rangs per separat
                   # e.g. list(1:8, 9:13, 14:17)
                   # Si hi ha només un grup de comparacions (p.e. Estres Termic), i amb dos nivells (estrés versus no estres), aquí es posaria com a:
                   # list(1:1)

pValCutOff <- c(0.01, 0.01, 0.01) # si N>1, indicar el cut-off per cada conjunt de comparacions
                      # e.g. c(0.01, 0.05, 0.01) o bé c(rep(0.01,3))
                      # Com a màxim a la UEB es posa 0.25 i a adjMethod posar no ajustat ("none"). 
adjMethod <- c("none", "none", "none")  # si N>1, indicar mètode per cada conjunt de comparacions
                      # e.g. c("none", "BH", "none") o bé c(rep("BH",3))

## Posar aqui els valors de minLogFC de cadascuna de les comparacions a fer
minLogFoldChange <- c(1, 1, 1) # canviar aixo si s'ha decidit considerar sols els casos amb |logFC| >= que un valor minim
                               # e.g. c(1, 2, 1) o be c(rep(0, 3)) indicar minLogFC per cada grup de comparacions

## Controlem que el nombre d'elements dels parametres anteriors sigui igual
if(class(wCont)!="list") warning("L'objecte wCont no és una llista! Això pot provocar errors en els passos següents.")
if(length(wCont)!=length(compName)) warning("L'objecte wCont ha de tenir el mateix nombre d'elements que compName!")
if(length(pValCutOff)!=length(compName)) warning("L'objecte pValCutOff ha de tenir el mateix nombre d'elements que compName!")
if(length(adjMethod)!=length(compName)) warning("L'objecte adjMethod ha de tenir el mateix nombre d'elements que compName!")
if(length(minLogFoldChange)!=length(compName)) warning("L'objecte minLogFoldChange ha de tenir el mateix nombre d'elements que compName!")

## Creem tants objectes com grups de comparacions 
for(i in 1:length(compName))
{  
  mci <- list(fitMain = NULL,
              fitFileName = p$fitFileName,
              whichContrasts = wCont[[i]],
              comparisonName = compName[i],
              titleText = paste("for",
                                 ifelse(adjMethod[i]=="none","p-values","adj. p-values"),
                                 "<",
                                 pValCutOff[i], 
                                 "and |logFC| >",
                                 minLogFoldChange[i], sep = " "),
              anotPackage = p$anotPack,
              my.symbols = p$SYMBOLIDs,
              outputDir = p$outputDir,
              fileOfLinks = p$fileOfLinks,
              multCompMethod = "separate",
              adjustMethod = adjMethod[i], 
              selectionType = "any",
              P.Value.cutoff = pValCutOff[i],
              plotVenn = TRUE,
              colsVenn = NULL,
	      vennColors= c("red","yellow","green","blue","pink"),  ## Es poden modificar els colors del diagrama de venn
              cexVenn = 1,
              geneListFName = paste("geneList",
                                    compName[i],
                                    ifelse(adjMethod[i]=="none","pvalues","adj-pvalues"),
                                    "LT",
                                    pValCutOff[i],
                                    "Rda",
                                    sep = "."),
              minLogFC = minLogFoldChange[i],
              csvType = p$csvType)
  mcParsList <- add2parsList(mcParsList, mci)
}


###################################################
### chunk number 17: clustPars
###################################################
###line 793 "ComCrearUnProjecteParametres.rnw"

## Actualització 2011-11: No cal fer servir el marray per la paleta
## El paquet gplots és més natural

if(!require(gplots, lib = mylib.R)){
  install.packages("gplots", lib = mylib.R)
}

pal <- colorpanel(n = 32, low = "green", mid = "white", high = "magenta") # XXXX Si cal editar els colors del Heatmap, cal fer-ho aquí.  e.g. colorpanel(n = 32, low = "green", mid = "black", high = "magenta")

p <- Estudi

for(i in 1:length(wCont))
{
  s2clust <- which(as.logical(apply(design[,as.logical(apply(abs(as.matrix(contrastsMatrix[,wCont[[i]]])),1,sum))],1,sum)))
  
  clustPar <- list(expres = NULL,
                   expresFileName = p$expresFileName,
                   geneListFName = paste("geneList",
                                    compName[i],
                                    ifelse(adjMethod[i]=="none","pvalues","adj-pvalues"),
                                    "LT",
                                    pValCutOff[i],
                                    "Rda",
                                    sep = "."),
                   genes2cluster = NULL,
                   samples2cluster = s2clust,
                   sampleNames = as.character(p$targets$ShortName)[s2clust],
                   comparisonName = compName[i],
                   anotPackage = p$anotPack,
                   my.symbols = p$SYMBOLIDs,
                   outputDir = p$outputDir,
                   fileOfLinks = p$fileOfLinks,
                   numClusters = 2,
                   rowDistance = NULL,
                   colDistance = NULL,
                   RowVals = TRUE,
                   ColVals = FALSE,
                   escala = "row",
                   colorsSet = pal,
                   densityInfo = "density",
                   colsForGroups = as.character(p$targets$Colores)[s2clust],
                   cexForColumns = 0.8,
                   cexForRows = 0.8,
                   Title = paste(compName[i],
                                 "with",
                                 ifelse(adjMethod[i]=="none","pvalues","adj-pvalues"),
                                 "<",
                                 pValCutOff[i], ifelse(minLogFoldChange[i]==0, "", paste("\n and |logFC|>=", minLogFoldChange[i], sep=""))),
                   paste("Comparison:", compName[i], sep=" "),
                   csvType = p$csvType)
  
  clustParsList <- add2parsList(clustParsList, clustPar)
}


###################################################
### chunk number 18: GOPars
###################################################
###line 854 "ComCrearUnProjecteParametres.rnw"
 
p <- Estudi

GOPar <- list(fitFileName = p$fitFileName,
              whichContrasts = p$whichContrasts,        
              comparisonName = p$comparisonName,
              anotPackage =  orgPackage,
              my.IDs = p$ENTREZIDs,        
              addGeneNames = TRUE, 
              outputDir = p$outputDir, 
              fileOfLinks = p$fileOfLinks,
              fitMain = NULL,
              cutoffMethod = "adjusted",  # XXX pot ser "adjusted" o "none" 
              P.Value.cutoff = rep(0.05, length(p$whichContrasts)), # XXXX es podria voler adjustar si al final surten pocs elements a la llista, 
                                                                    # i ho podríem posar a 0.05 ó 0.25, per que acabin sortint algunes desenes o centenars de gens anotats a la GO.
              pvalGOterms = 0.01,
              min.count = 3,
              ontologias = c("MF", "BP", "CC"),
              testDirections = c("over", "under"),
              minNumGens = 0)
 
GOParsList <- add2parsList(GOParsList, GOPar)


###################################################
### chunk number 19: KEGGPars
###################################################

p <- Estudi

KEGGPar <- list(fitFileName = p$fitFileName,
                whichContrasts = p$whichContrasts,
                comparisonName = p$comparisonName,
                outputDir = p$outputDir,
                anotPackage = orgPackage, 
                organisme = organisme,    
                my.IDs = p$ENTREZIDs,
                addGeneNames = TRUE,
                outputDir = p$outputDir,
                fileOfLinks = p$fileOfLinks,
                fitMain = NULL,
                cutoffMethod = "adjusted", # XXX pot ser "adjusted" o "unadjusted"
                P.Value.cutoff = rep(0.05, length(p$whichContrasts)), # XXXX es podria voler adjustar si al final surten pocs elements a la llista, 
                                                                    # i ho podríem posar a 0.05 ó 0.25, per que acabin sortint algunes desenes o centenars de gens anotats a la KEGG.
                pvalKEGGterms = 0.05, # Potser d'ha de canviar per que s'agafin més coses.
                minNumGens = 0) 

KEGGParsList <- add2parsList (KEGGParsList, KEGGPar)


###################################################
### chunk number 19: PowerAnalysis
###################################################
## de moment, esta preparat per una sola comparacio entre 2 grups (estil A - B)

sampleA.size <- 3    # nombre de mostres (CEL files) del 1r grup a incloure a l'analisi de potencia (funcio pilotData)
sampleB.size <- 3    # nombre de mostres (CEL files) del 2n grup  "  "  "  "
powerAnalysisTitle <- "X" # XXX <- replace X with the variable name for this power analysis plot, with no spaces (it will be used for the pdf plot tile and pdf file name)

sample.sizes <- c(2, 3, 4, 5, 6, 7, 8, 9, 10) # llista de possibles mides mostrals a avaluar
 fdr.values <- c(0.01, 0.05) # llista de les FDR que es vulguin considerar en el calcul (0.01 i 0.05 per defecte)
                             # Podran ser mes grans o mes petits en funcio de la 'proportion of non-differentially expressed genes'
                             # (veure output text de ss i ssMod anteriors)
                             # No te sentit que aquests fdr siguin majors que la proporcio indicada.




###################################################
### chunk number 20: Session Info Parameters
###################################################

dataFi <- "2014-02-11"        # format "aaaa-mm-dd"

sessionInfoFileName <- "sessionInfo.txt"
outDir <- tempDir


###################################################
# chunk number 21: Memory Profiling
###################################################

profile_memory <- FALSE      # Param to profile memory usage



