###
# Basic-Pipe, as of AUgust 2011 is contained in the following-to-source file
###

source(file.path(codeDir, "AnalysisFunctions2Pack.R"))

###################################################
### chunk number 0: Memory Profiling (optional)
###################################################

profileMemory (profile_memory, eP)

###################################################
### chunk number 1: loadLibraries
###################################################

cat(format(Sys.time(), "%d/%m/%y %H:%M:%S - ") );
cat("Loading Functions and Packages...\n")


###
# Carrega dels paquets necessaris al llarg del codi d'anàlisi
# Important: nomes es carrega el paquet d'anotacions
#                         per al que la variable booleana es TRUE
# Semopre es carrega el orgPackage
###


loadPackages( TRYTOINSTALL,
              chipPackAvailable, platformDesignPackAvailable,
              chipPackage,
              runMultiCore)

##################################################
### chunk number 2: loadData
###################################################
#line 117 "ComCrearUnProjecteAnalisis.rnw"

### Creació de  l'arxiu de vincles a partir del qual es farà el 'Results_File'

cat(format(Sys.time(), "%d/%m/%y %H:%M:%S - ") ); cat("Loading Data...\n")

if(!file.exists(linksFileName)) createLinksFile(linksFileName)

my.targets <- read.AnnotatedDataFrame(file.path(dataDir, targetsFileName), header = TRUE, row.names = 1)
my.fileNames <- file.path(celFilesDir, rownames(pData(my.targets)))

addToLinksFile(linksFileName,
               StudyProposal,
               categ = 'INFO', 
               desc = "Study Proposal")

addToLinksFile(linksFileName,
               MainReportName,
               categ = 'INFO', 
               desc = "Report with description of the main results")

addToLinksFile(linksFileName,
               targetsFileName,
               categ = 'DATA', 
               desc = "Samples and covariates information file")

rawData <- readOrLoad.RawData(readCELS = readCELS,
                              phenoDat = my.targets, 
                              fileNames = my.fileNames,
                              dataFName = rawDataFileName,
                              outputDir = dataDir,
                              exonSt = exonStudy)
#dim(exprs(rawData))
#head(exprs(rawData))

if(rawDataType=="affy3primeIVT"){
  dates <-  pData(protocolData(rawData))
  targetsWithDatesCervell <- cbind(pData(my.targets), scandate=substr(dates$ScanDate,1,10))
  write.table(targetsWithDatesCervell, 
            file.path(dataDir,paste("withDates", targetsFileName, sep=".")), 
            sep="\t", row.names=F)
}

if(!exonStudy)
{
  annotation(rawData) <- old2db(chipPackage)
}else{
  annotation(rawData) <- platformDesignPackage
}

###################################################
# creaOCarregaAnotacions
##################################################

anotacions <- createOrLoadAnnotations (loadAnnotations=!createAnnotations,
                                       chipPackAvailable = chipPackAvailable,
                                       platformDesignPackAvailable = platformDesignPackAvailable,
                                       chipPackage = chipPackage,
                                       platformDesignPackage = platformDesignPackage,
                                       outputDir = dataDir,
                                       annotationsFileName = annotationsFileName,
                                       entrezTableFileName = entrezTableFileName,
                                       symbolsTableFileName = symbolsTableFileName,
                                       controlsTableFileName = controlsTableFileName
                                       )
  entrezTable <- anotacions$Entrez
 symbolsTable <-anotacions$Symbols
controlsTable <- anotacions$controls

 experimentInfo <- new("MIAME",
                       name = htmlInfo$To,
                       lab = lab_o_depart ,
                       contact = contact_email,
                       title = htmlInfo$Description,
                       abstract = "",
                       url = "",
                       other = htmlInfo)

experimentData(rawData) <- experimentInfo

###################################################
### chunk number 3: qualityControl
###################################################

if(preProcess)
{
  cat(format(Sys.time(), "%d/%m/%y %H:%M:%S - ") ); cat("Processing arrayQualityMetrics...\n")

  if(exonStudy)
  {  
    phenoData(rawData)$Grupo<-targets$Grupo
    phenoData(rawData)$ShortName<-targets$ShortName
    phenoData(rawData)$Colores<-targets$Colores
  }




  arrayQualityMetrics(expressionset =rawData,
                      outdir = paste(QCDir, "raw", sep = "."),
                      force = TRUE,
                      intgroup = "Grupo",
#                     grouprep = TRUE,
                      do.logtransform = FALSE)
#                     sN = as.character(targets$ShortName))
#                     split.plots = FALSE)

  addToLinksFile (linksFileName,
                 "QCDir.raw/index.html",  
                  categ = 'QC',
                  desc = "Quality Control Plots of raw data with ArrayQualityMetrics package")
}


###################################################
### chunk number 4: normalitzacio
###################################################

cat(format(Sys.time(), "%d/%m/%y %H:%M:%S - ") ); cat("Normalizing Data...\n")

eset_norm <- normalization(my.data = rawData,
                           method = normMethod,    
                           targetsinfo = my.targets,
                           inputDir = celFilesDir,
                           loadFile = !doNorm,
                           normalizedFName = normalized.eset.FileName,
                           outputDir = dataDir,
                           exonSt = exonStudy)

## XXXX Faig un bypass del eset_norm i li assigno el d_vsn que hem fet amb el PreparaDades.A279.R
## Tot i que aquest d_vsn té només un feature, i a dins té:
# assayData: 1 features, 48 samples 
#> dim(exprs(rawData))
#[1] 292681     48
# I després d'assignar "eset_norm <- d_vsn" més avall:
#> dim(exprs(eset_norm))
#[1] 36249    48

## i en canvi el rawData dels CEL files llegits pel basic pipe té: 
# assayData: 292681 features, 48 samples 

eset_norm <- d_vsn2 # d_vsn2 ja té les columnes i files (covariables) coincidint amb la info del targets
dim(exprs(eset_norm))
head(exprs(eset_norm))

## PART QUE ELIMINA LES FILES AMB IDENTICS VALORS D'EXPRESSIO (FERRAN, JUNY 2012)
## guardem com a valides les files de exprs(eset_norm) que NO estiguin duplicades
## (en cas de duplicats, es guarda la primera instancia de cada una com a bona)
repes <- duplicated(exprs(eset_norm), MARGIN=1)
#table(repes)
exprs(eset_norm) <- exprs(eset_norm)[!repes,]
#dim(eset_norm)


### Afegit per Alex l'Agost del 2012. AIXO NO HI HAVIA DE SER JA???

if(!(is.null(pData(eset_norm)$ShortName)))
  colnames(exprs(eset_norm))<- pData(eset_norm)$ShortName

### En el cas d'exons s'ha d'especificar aquestes variables:

if(exonStudy)
{  
  phenoData(eset_norm)$Grupo<-targets$Grupo
  phenoData(eset_norm)$ShortName<-targets$ShortName
  phenoData(eset_norm)$Colores<-targets$Colores
}


###################################################
### chunk number 5: QCNorm
###################################################
#line 247 "ComCrearUnProjecteAnalisis.rnw"

### QCDir de dades normalitzades amb arrayQualityMetrics
# En el cas d'exons s'ha d'especificar aquestes variables:
if(exonStudy)
{  
  phenoData(eset_norm)$Grupo<-targets$Grupo
  phenoData(eset_norm)$ShortName<-targets$ShortName
  phenoData(eset_norm)$Colores<-targets$Colores
}

 if(preProcess)
{
  arrayQualityMetrics(expressionset = (eset_norm),
                      outdir = paste(QCDir, "norm", sep = "."),
                      force = TRUE,
                      intgroup = "Grupo",
#                     grouprep = TRUE,
                      do.logtransform = FALSE)
#                     sN = as.character(targets$ShortName))
#                     split.plots = FALSE)

  addToLinksFile (linksFileName,
                 "QCDir.norm/index.html",  
                  categ = 'QC',
                  desc = "Quality Control Plots of normalized data with ArrayQualityMetrics package")

  
}

### QC amb els nostres grafics
if(normPlots)             
{
  cat(format(Sys.time(), "%d/%m/%y %H:%M:%S - ") ); cat("Generating NormPlots...\n")

  myCex <- 0.8

  ### normalized plots (boxplots, dendrogram, PCA)
 
  normplots2File(my.data = eset_norm,
                 sampleNames = pData(eset_norm)$ShortName,
                 my.colors = as.character(pData(eset_norm)$Colores), 
                 my.groups = pData(eset_norm)$Grupo,
                 my.method = "average",
                 my.cex = myCex,
                 posText = 2,
                 dim3 = FALSE,
                 fileName = normalized.Plots.FileName,
                 outputDir = resultsDir,
                 PCAPlots = TRUE,
                 csv = fileType)
              
  addToLinksFile (linksFileName,
                  normalized.Plots.FileName,
                  categ = 'QC',
                  desc = "Extra Quality Control Plots based on normalized data")
  
}


###################################################
### chunk number 6: removeChips eval=FALSE
###################################################
### #line 296 "ComCrearUnProjecteAnalisis.rnw"
###
### if(removeChips)
### {
###    chips2remove <- which(pData(rawData)$ShortName %in% c("HEM-2-99", "LAC-2-48", "LAC-2-119"))
###    chipsRemoved <- TRUE
###                     
###    rawData.ok <- rawData [, -chips2remove]
###    my.targets <- phenoData(rawData.ok)
### 
###    eset_norm <- normalization(rawData.ok,
###                               normMethod,                                  
###                               my.targets,
###                               inputDir = celFilesDir,
###                               loadFile = !doNorm,
###                               normalizedFName = normalized.eset.FileName,
###                               outputDir = dataDir)
###    source(file.path(dataDir, "parametros.chipsRemoved.XXX.R"))
### }


###################################################
### chunk number 7: filterData
###################################################
#line 346 "ComCrearUnProjecteAnalisis.rnw"

######################################################################
### Actualitzacio 2011 (PENDENT)
### Els fragment de codi d'aqui sota es redundant.
### S'ha d'arreglar treient el filterstring de la funcio filterData
######################################################################

if (FILTRAR){
  if(processaFILTRE)
    {
      cat(format(Sys.time(), "%d/%m/%y %H:%M:%S - ") ); cat("Applying filters...\n")
      exprs.filtered <- filterData(expres = exprs(eset_norm),                     
                               controls = names(controlsTable),
                               removeNAs = TRUE, 
                               entrezs = entrezTable,
                               bySignal = SignalFilter,
                               signalThr = signalThreshold,
                               grups = pData(eset_norm)$Grupo,
                               sigFun.Name = signalFilter.Function,
                               sigThr.as.perc = signalThreshold.as.percentage,
                              
                               byVar = VarFilter,
                               variabilityThr = variabilityThreshold,
                               varFun.Name = variability.Function,
                               varThr.as.perc = variabilityThreshold.as.percentage,
                              
                               pairingFun.Name = pairing.Function,
                               targets = my.targets,

                               doReport = doReport,
                               outputDir = tempDir,
                               filteringReportFName = FilteringReportFileName)
    }else{
      load(file.path(resultsDir,expres.filtered.FileName))
      exprs.filtered <- expres}
}else{
  exprs.filtered <- exprs (eset_norm)
}


###################################################
### chunk number 8: saveData
###################################################
#line 383 "ComCrearUnProjecteAnalisis.rnw"

#################################################################
# Actualitzacio 2011_2-11 (PENDENT)
# Originalment hi havia un camp anotPackage que es consultava 
# per treure els symbols.
# Això es treura pero de moment cal 
#       passar-li els symbolsTable i 
#       deixar l'anotpackage a NULL
#################################################################

cat(format(Sys.time(), "%d/%m/%y %H:%M:%S - ") ); cat("Saving Data...\n")

saveData(expres = exprs(eset_norm),
         expres.csv.FileName = normalized.all.FileName,
         csvType=fileType,
         description = "Normalized values for all genes",  
         anotPackage = NULL,
         symbolsVector = symbolsTable,
         SYMBOL = "SYMBOL", 
         expres.bin.FileName = expres.all.FileName, 
         linksFile = linksFileName,
         outputDir = resultsDir)

saveData(expres = exprs.filtered,
         expres.csv.FileName = normalized.filtered.FileName,
         csvType = fileType,
         description = "Normalized values for filtered genes",  
         anotPackage = NULL,
         symbolsVector = symbolsTable,
         SYMBOL = "SYMBOL", 
         expres.bin.FileName = expres.filtered.FileName, 
         linksFile = linksFileName,
         outputDir = resultsDir)


###################################################
### chunk number 9: lmAnalysis
###################################################
#line 456 "ComCrearUnProjecteAnalisis.rnw"

if (processaLM)
{
  cat(format(Sys.time(), "%d/%m/%y %H:%M:%S - ") ); cat("Processing lmAnalysis... \n")

  for(ix in 1:length(lmParsList))
  {
    fit.Main   <- doLmAnalysis(lmParsList[ix])
  }
}


###################################################
### chunk number 10: annotateResults
###################################################

#line 606 "ComCrearUnProjecteAnalisis.rnw"

#######################################################################
# Actualitzacio 2011_2-12 (PENDENT)
# En aquests moments no existeix un paràmetre propi per les anotacions 
# i s'està fent servir el metaparàmetre de les comparacions múltiples 
# ja que aquest conté tota la informació necessària. 
# Això s'ha de revisar ja que ara no volem fer servir els paquets de xip
# però tenim els d'organisme, que necessitarien l'ENTREZ no els SYMBOLS
# Comprovar si funciona be en tots els casos.
#######################################################################

#line 476 "ComCrearUnProjecteAnalisis.rnw"

if (annotateGenes)
{
   cat(format(Sys.time(), "%d/%m/%y %H:%M:%S - ") ); cat("Processing Annotation... \n")

   doGeneAnnotation(anotList)
}


###################################################
### chunk number 11: multipleCompAnalysis
###################################################
#line 491 "ComCrearUnProjecteAnalisis.rnw"

if (processaMultComp)
{
  cat(format(Sys.time(), "%d/%m/%y %H:%M:%S - ") ); cat("Processing MultComp... \n")

  for(ix in 1:length(mcParsList))
  {
    geneList.MCi <- doMultCompAnalysis(mcParsList[ix])   #geneList.MC cont?? la llista de gens seleccionats en alguna de les comparacions...
  }
}


###################################################
### chunk number 12: clusterAnalysis
###################################################
#line 508 "ComCrearUnProjecteAnalisis.rnw"

if (processaCluster)
{
  cat(format(Sys.time(), "%d/%m/%y %H:%M:%S - ") ); cat("Processing ClusterAnalysis... \n")

  for(ix in 1:length(clustParsList))
  {
    hm.Estudi <- doClusterAnalysis(clustParsList[ix])
  }
} 


###################################################
### chunk number 13: GOAnalysis
###################################################
#line 524 "ComCrearUnProjecteAnalisis.rnw"

if (processaGO)
{
  cat(format(Sys.time(), "%d/%m/%y %H:%M:%S - ") ); cat("Processing GOAnalysis... \n")

  if(runMulticore == 2 || runMulticore == 3){  # 0=none | 1=mclapply | 2=foreach | 3=mcapply & foreach
    foreach(ix = 1:length(GOParsList)) %dopar% {
      GOList <- doGOAnalysis (GOParsList[ix])
    }
  }else{
    for(ix in 1:length(GOParsList)) {
      GOList <- doGOAnalysis (GOParsList[ix])
    }
  }
}


###################################################
### chunk number 14: KEGGAnalysis
###################################################

if (processaKEGG)
{
  cat(format(Sys.time(), "%d/%m/%y %H:%M:%S - ") ); cat("Processing KEGG Analysis... \n")

  for(i in 1:length(KEGGParsList))
  {
    KEGGList <- doKEGGAnalysis(KEGGParsList[i])
  }
}
###################################################
### chunk number 15: powerAnalysis
###################################################
### based on SSPA package vignette
### by now, it's ready for a single comparison between two groups (such as A - B)

if (processaPower)
{
  cat(format(Sys.time(), "%d/%m/%y %H:%M:%S - ") ); cat("Processing power analysis... \n")

  ## load data from Rda file with fit.main (returning the fit.main object)
  load(file.path(resultsDir, Estudi$fitFileName))  # Estudi$fitFileName s'haura generat a partir de qualsevol dels chunks anteriors

  ## assign fit.main to fit for SSPA original code reusability
  fit <- fit.main # agafem el fit.main procedent de l'analisi pero com que 
                  #   aquest ja es moderat hem de recalcular els ordinary.t

  ## compute ordinary and assign moderated t values
  ordinary.t <- fit$coef/fit$stdev.unscaled/fit$sigma  # formula que dona el paquet SSPA per l'ordinary.t
  moderated.t <- fit$t # com que el fit.main procedent de l'analisi ja es moderat per ebayes
                       #   agafem els valors t directament d'alla per tenir els moderats

  ## assign degrees of freedom
  nu <- fit$df.residual[1]
  nu0 <- fit$df.prior

  ## create ordinary and moderated SSPA's data pilots
  pd <- pilotData(name = "ord.pd", statistics = ordinary.t[,1], 
                  samplesize=sqrt(1/((1/sampleA.size)+(1/sampleB.size))),
                  sampleSizeA = sampleA.size,  #sampleA.size s'estableix al parametros
                  sampleSizeB = sampleB.size,  #sampleB.size s'estableix al parametros
                  dof = nu, nullDist = "student")
  pd
  pdMod <- pilotData(name = "mod.pd", statistics = moderated.t[,1], 
                     samplesize=1,
                     sampleSizeA = sampleA.size + nu0/2,
                     sampleSizeB = sampleB.size + nu0/2, 
                     dof = nu + nu0, nullDist = "student")
  pdMod

  ## create ordinary and moderated SSPA's sample size objects
  ss <- sampleSize(pd)
  ss
  ssMod <- sampleSize(pdMod)
  ssMod

  ## create ordinary and moderated SSPA's power objects
  pwr <- predictpower(ss, plot = FALSE, samplesizes = sample.sizes, alpha = fdr.values)
  pwrMod <- predictpower(ssMod, plot = FALSE, samplesizes = sample.sizes + nu0/2, alpha = fdr.values)

  ## power plot with all data series, ordinary and moderated for each distinct fdr
  PowerPlotFileName <- paste0("PowerPlot", powerAnalysisTitle,".pdf")
  PowerPlotTitle <- paste0("Power Curves for ", powerAnalysisTitle)
  
  pdf(file.path(resultsDir, PowerPlotFileName))
    matplot(sample.sizes, pwr, ylim = c(0.0, 1.1), type = "b", col = 1, pch = 1, 
            ylab = "Power", xlab = "Sample size per group", main=PowerPlotTitle)
    matlines(sample.sizes, pwrMod, col = 2, type = "b", pch = 1)
    #legend("bottomleft", colnames(pwr), col = 1, lty = 1:length(fdr.values))
    legend("bottomright", c("ordinary", "moderated"), col = 1:2, lty = 1)
  dev.off()
  
  # Write the appropriate data into the links file, for the inclusion of those result files in the 
  # Resultfiles.html generated at the end.
  addToLinksFile(p$fileOfLinks, PowerPlotFileName, categ = 'POWER', desc = PowerPlotTitle)
  
  ### CODI ALTERNATIU A CONSIDERAR SI S'UTILITZA EL PACKET ssize ENLLOC DEL SSPA
  #########################
  ### ssize
  #########################
  #load("results/expres.filtered.Rda")
  #desvs <- apply(expres, 1, function(x) sd(x[1:3]))
  #n<- 3
  #fold.change<- 2^1.5
  #sig.level <-0.01
  #power<-0.8
  #meansdev<-mean(desvs)
  #
  # 2) What is necessary sample size for 80% power using 3 measurements/patient
  #    assuming Delta=1.0, Alpha=0.05 and Observed SDs?
  #
  #all.size  <- ssize(sd=desvs, delta=log2(fold.change),
  #                   sig.level=sig.level, power=power)
  #ssize.plot(all.size, lwd=2, col="magenta", xlim=c(1,20))
  #xmax <- par("usr")[2]-1; ymin <- par("usr")[3] + 0.05
  #legend(x=xmax, y=ymin,
  #       legend= strsplit( paste("fold change=",fold.change,",",
  #                              "alpha=", sig.level, ",",
  #                              "power=",power,",",
  #                              "# genes=", length(desvs), sep=''), "," )[[1]],
  #       xjust=1, yjust=0, cex=1.0)
  #title("Sample Size to Detect 2-Fold Change")
}

###################################################
### chunk number 16: writeReport
###################################################
#line 541 "ComCrearUnProjecteAnalisis.rnw"

if (processaReport)
{
  LinksFile2Html(linksFileName,
                 resultsDir,
                 htmlInfo,
                 IndexDir = "ResultFiles/",   # IndexDir : "ResultFiles/" (Paràmetre opcional)
                 UEB = UEB)                   # si UEB = TRUE => Capçalera UEB (indicat al parametros.XXX.R)
}

###################################################
### chunk number 17: writeSessionInfo
###################################################

if (processaSessionInfo)
{
  save.sessionInfo(fileName = sessionInfoFileName,
                   targetsdate = datacreacioTargets,
                   beginningDate = dataInici,
                   endingDate = dataFi,
                   outputDir = outDir
                   )
}


###################################################
### chunk number 18: memoryProfiling
###################################################
# Other examples of usage in: http://www.r-bloggers.com/examples-of-profiling-r-code/
if (profile_memory) # If monitor memmory is on, show a summary at the end
{
  cat(format(Sys.time(), "%d/%m/%y %H:%M:%S - ") ); cat("Profiling memory... \n")
  Rprof(NULL)
  print(summaryRprof(filename = rprof_filename,
                     chunksize = 5000,
                     memory = c("none", "both", "tseries", "stats"),
                     index = 2,
                     diff = TRUE,
                     exclude = NULL))		  
}


###################################################
cat(format(Sys.time(), "%d/%m/%y %H:%M:%S - ") ); cat("END OF ANALYSIS \n\n")

###################################################
############ END OF SCRIPT ########################
###################################################
