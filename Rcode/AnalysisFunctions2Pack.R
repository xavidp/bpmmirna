
#####################################################################################################################
### Miscel??l??nia: Funcions auxiliars
#####################################################################################################################

# Rprofmem: Initialize memory monitoring if requested

profileMemory <- function (profile_memory, eP){
  if (profile_memory)
    {
      rprof_filename = paste(eP,"_Rprof.txt", sep = "")
      Rprof(filename = rprof_filename, append = FALSE, interval = 0.02,
            memory.profiling = TRUE)
    }else{
      Rprof(NULL)
    }
}



### Instala i carrega els paquets necessaris per a Basic-Pipe


installifnot <- function (pckgName) {
  if (!(require(pckgName, character.only = TRUE))){
	biocLite(pckgName)
        require(pckgName, character.only = TRUE)
      }
}

if(!require(SortableHTMLTables)) install.packages("SortableHTMLTables")

loadPackages<- function(TRYTOINSTALL,
             chipPackAvailable, platformDesignPackAvailable,
             chipPackage,
             runMultiCore)
{
if(TRY2INSTALL)
{
  source("http://bioconductor.org/biocLite.R")
  installifnot("Biobase")
  installifnot("xtable")
  installifnot("annotate")     # Necessari per la funcio 'htmlPage'
  installifnot("AnnotationDbi")
  installifnot("affy")
  installifnot("limma")
  installifnot("VennDiagram")
  installifnot("simpleaffy")
  installifnot("affyQCReport")
  installifnot("affyPLM")
  installifnot("affycoretools")
  installifnot("rgl")
  installifnot("KEGG.db")
  installifnot("GO.db")
  installifnot("annaffy")
  installifnot("gplots")
  installifnot("GOstats")
  installifnot("marray")   # necessari pel PowerAnalysis
  installifnot("convert")  # necessari pel PowerAnalysis
  installifnot("SSPA")     # necessari pel PowerAnalysis
  #installifnot("ssize")    # alternativa a SSPA per un PowerAnalysis diferent (de moment descartat)
  installifnot("statmod")
  installifnot("vsn")
  installifnot("arrayQualityMetrics")
  installifnot("SortableHTMLTables")
  if(exonStudy) {
    installifnot("oligo")
  }
  if(chipPackAvailable){
    installifnot(old2db(chipPackage))
  }
      
}else{
# installifnot(old2db(anotPackage), lib.loc='/home/asanchez/R/x86_64-pc-linux-gnu-library/2.9')
  require("Biobase", lib.loc = mylib)
  require("limma", lib.loc = mylib)
  require("VennDiagram", lib.loc = mylib)
  require("graph", lib.loc = mylib)
  require("RBGL", lib.loc = mylib)
  require("genefilter", lib.loc = mylib)
  require("xtable", lib.loc = mylib)    
  require("AnnotationDbi", lib.loc = mylib)  
  require("annotate", lib.loc = mylib)    
  require("Category", lib.loc = mylib)
  require("GO.db", lib.loc = mylib)
  require("KEGG.db", lib.loc = mylib)  
  require("GOstats", lib.loc = mylib)
  require("biomaRt", lib.loc = mylib)
  require("affy", lib.loc = mylib)
# require("matchprobes", lib.loc = mylib)  # No es fa servir
  require("gcrma", lib.loc = mylib)
  require("annaffy", lib.loc = mylib)  
  require("affycoretools", lib.loc = mylib)  
  require("preprocessCore", lib.loc = mylib)   
  require("affyPLM", lib.loc = mylib)
  require("simpleaffy", lib.loc = mylib)
  require("RColorBrewer", lib.loc = mylib)
  require("geneplotter", lib.loc = mylib) 
  require("affyQCReport", lib.loc = mylib)
  require("Biobase", lib.loc = mylib)
  require("marray", lib.loc = mylib)   # canviat en algun punt per ggplots, pero requerit pel PowerAnalysis
  require("convert", lib.loc = mylib)  # requerit pel PowerAnalysis
  require("SSPA", lib.loc = mylib)     # requerit pel PowerAnalysis
  #require("ssize", lib.loc = mylib)   # requerit si es vol fer un PowerAnalysis diferent (de moment descartat)
  require("gtools", lib.loc = mylib)
  require("gdata", lib.loc = mylib) 
  require("gplots", lib.loc = mylib)
  require("rgl", lib.loc = mylib)
  require("statmod", lib.loc = mylib)
  require("vsn", lib.loc = mylib)
  require("arrayQualityMetrics", lib.loc = mylib)
  require("SortableHTMLTables", lib.loc = mylib)

  if(exonStudy) {
    require("oligo", lib.loc = mylib)
  }
      
  if(chipPackAvailable)
  {
    require(paste(chipPackage,"db",sep = "."), lib.loc = mylib, character.only = TRUE)
  }else{
    if(platformDesignPackAvailable)
    {
      require("statmod", lib.loc = mylib)
    }
    if(exonStudy)
    {
      require("oligo", lib.loc = mylib)
    }
    require("oligoClasses", lib.loc = mylib)
    }
}
if(runMulticore > 0 && runMulticore < 4 ){
  installifnot("multicore")
  require("multicore", lib.loc = mylib)
  if(runMulticore > 1){
    installifnot("doMC")
    require("doMC", lib.loc = mylib)
    registerDoMC()
  } 
}else if(runMulticore != 0) {
  cat(format(Sys.time(), "%d/%m/%y %H:%M:%S - ") );
  cat("Parameter runMulticore is missing or out of bounds. Revise that parameter in your procesa.XXXnnn.R file")
}
}

### Afegir un parametre a una llista

add2parsList <- function(oneList,object)
{
  pos <- length(oneList) + 1
  oneList[[pos]] <- object
  names(oneList)[pos] <- oneList[[pos]]$comparisonName
  return(oneList)
}


### Afegir a un arxiu (linksFile) el nom del fitxer, la categoria i la descripcio
### per a poder contruir el navegador ResultsFile.html

addToLinksFile <- function(linksFile, aFName, categ = "", desc = "", subcateg = "")
{
  if (!is.null(linksFile))
  {
    write(paste(aFName, categ, subcateg, desc, sep = "\t"), file = linksFile, append = TRUE)
  }
}


### Afegir extensio .db a un paquet d'anotacions

old2db <- function(anot){paste(anot, "db", sep = ".")}


### Funcio per escriure una matriu (d'expressio). Per defecte la desa en format csv2

write2csv <- function(my.data, fileName, csv = c("csv2", "csv", "txt", "xls"), outputDir)
{
  fileName<- file.path(outputDir, paste(fileName, substr(csv[1], 1, 3) , sep = "."))

  #if(require(dataframes2xls)){
  #  switch (csv,
  #    "csv" = write.csv(my.data, file=fileName, quote = F),
  #    "csv2" = write.csv2(my.data, file=fileName , quote = F),
  #    "txt" = write.table(my.data,  file=fileName, quote = F),
  #    "xls" = write.xls(my.data, file=fileName)
  #  )
  #}else{

  switch (csv[1],
          "csv" = write.csv(my.data, file = fileName, quote = F),
          "csv2" = write.csv2(my.data, file = fileName, quote = F),
          "txt" = write.table(my.data, file = fileName, quote = F))
  # }
}

# Funcio per crear una variable amb el contingut d'un arxiu binari 
# Serveix per poder donar el nom que desitgem a la variable enlloc del que tenia quan s'ha fet el "load"

loadFromFile <-function (fileName, pos=1){
  tempEnv =new("environment")
  load (fileName, tempEnv)
  varNames <-ls(tempEnv)
  myVarName <- varNames[pos]
  load (fileName)
  myVar <- eval(parse(text=myVarName))
  return(myVar)
}


### La funcio write.htmltable es va extreure del paquet arrayMagic. Aquest package el van eliminar
### de R.
### NOTA: Seria convenient estudiar les noves funcions dels paquets nous que generen taules
###       i codi HTML

write.htmltable <- function (x, filename, title = "", sortby = NULL, decreasing = TRUE, 
    open = "wt", formatNumeric = function(x) paste(signif(x, 
        3))) 
{
    if (!is.null(sortby)) {
        if (!sortby %in% colnames(x)) 
            stop(paste("Invalid argument \"sortby\": could not find a column in data frame x with name", 
                sortby))
        soby = x[, sortby]
        if (!is.numeric(soby)) 
            stop("Invalid argument \"sortby\": column is not numeric")
        x = x[order(soby, decreasing = decreasing), ]
    }
    outfile <- file(paste(filename, ".html", sep = ""), open = open)
    cat("<html>", "<STYLE>", "<!--TD { FONT-FAMILY: Helvetica,Arial; FONT-SIZE: 14px;}-->", 
        "<!--H1 { FONT-FAMILY: Helvetica,Arial; FONT-SIZE: 22px;}-->", 
        "</STYLE>", "<head>", paste("<TITLE>", title, "</TITLE>", 
            sep = ""), "</head>", "<body bgcolor=#ffffff>", file = outfile, 
        sep = "\n")
    if (title != "") 
        cat("<CENTER><H1 ALIGN=\"CENTER\">", title, " </H1></CENTER>\n", 
            file = outfile, sep = "\n")
    cat("<CENTER> \n", file = outfile)
    cat("<TABLE BORDER=0>", file = outfile, sep = "\n")
    cat("<TR>", file = outfile)
    for (j in 1:ncol(x)) cat("<TD BGCOLOR=\"", c("#e0e0ff", "#d0d0f0")[j%%2 + 
        1], "\"><B>", colnames(x)[j], "</B></TD>\n", sep = "", 
        file = outfile)
    cat("</TR>", file = outfile)
    for (i in 1:nrow(x)) {
        cat("<TR>", file = outfile)
        for (j in 1:ncol(x)) {
            txt = switch(class(x[[j]]), numeric = formatNumeric(x[i, 
                j]), as.character(x[i, j]))
            if (length(grep("^http:", txt)) > 0) {
                txt <- sub(";$", "", txt)
                s <- unlist(strsplit(txt, "[/?=]"))
                txt <- paste("<A HREF=\"", txt, "\" TARGET=\"z\">", 
                  s[length(s)], "</A>", sep = "")
            }
            cat("<TD BGCOLOR=\"", c("#e0e0ff", "#d0d0f0", "#f0f0ff", 
                "#e0e0f0")[i%%2 * 2 + j%%2 + 1], "\">", txt, 
                "</TD>\n", sep = "", file = outfile)
        }
        cat("</TR>", file = outfile)
    }
    cat("</TABLE></CENTER>", "</body>", "</html>", sep = "\n", 
        file = outfile)
    close(outfile)
}

createLinksFile <- function (linksFile)
{
  linia <- paste("FileName", "Category", "Subcategory", "Description", sep = "\t")
  write(linia, file = linksFile, append = FALSE)
}


addToLinksFile_old <- function (linksFile, aFName, categ ="", subcateg = "", desc = "")
{
  if (!is.null(linksFile))
  {
    write(paste(aFName, categ, subcateg, desc, sep = "\t"), file=linksFile, append=TRUE)
  }
}


FixLinksFile <- function(lFile)
{
  my.LnkFile <- read.table(lFile, header = TRUE, sep = "\t")

  my.newFile <- unique(my.LnkFile)
  
  write.table(my.newFile, lFile, sep = "\t", row.names = FALSE, quote = FALSE)
}


### save.sessionInfo: Desa en un arxiu de text la informacio sobre la sessio de R amb la que s'ha processat
###                   l'analisi.
###                 
###
### Parametres:
###
###     fileName : Nom del fitxer de text a on es desara la informacio. Per defecte pren com a nom "sessionInfo.txt"
###    outputDir : Ruta a on es desara el document
###
### Exemple
###
###    > save.sessionInfo(fileName = "mySessionInfo.txt", outputDir = "/mydocs/myRproject/temp")
###

save.sessionInfo <- function(fileName = "sessionInfo.txt", outputDir)
{
  sink(file = file.path(outputDir, fileName))

      cat("\n")
      cat("###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")  
      cat("### Analysis Dates\n")
      cat("###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")  
      cat("\n")
      cat("Beginning date: "); print(dataInici)
      cat("Ending date: "); print(dataFi)
      cat("Last run date: "); print(date())
      cat("\n")
      cat("\n")
      cat("###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")  
      cat("### Summary Session Info\n")
      cat("###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")  
      cat("\n")
      print(sessionInfo())
      cat("\n")
      cat("\n")
      cat("###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")  
      cat("### Extended Session Info\n")
      cat("###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")  
      cat("\n")
      cat("###------------------------------------------------------------\n")  
      cat("### R version\n")
      cat("###------------------------------------------------------------\n")  
      cat("\n")
      print(R.version)
      cat("\n")
      cat("\n")
      cat("###------------------------------------------------------------\n")  
      cat("### Installed packages\n")
      cat("###------------------------------------------------------------\n")  
      cat("\n")
      print(installed.packages())
  
  sink()
}


#####################################################################################################################
### Funcions per a l'entrada de les dades
#####################################################################################################################

readOrLoad.RawData <- function (readCELS, phenoDat, fileNames, dataFName, outputDir, exonSt=FALSE, cdf=NULL)
{
  stopifnot(require(affy))

  if (readCELS)
  {
#    print(fileNames)
    if(exonSt)
    {
      my.raw <- read.celfiles(filenames = fileNames, verbose = TRUE)
    }else{
        my.raw <- ReadAffy(filenames = fileNames, phenoData = phenoDat,
                           verbose = TRUE, cdfname=cdf)
    }
	if (is.null(platformDesignPackage))
		platformDesignPackage <<- annotation(my.raw)

	if (platformDesignPackage != annotation(my.raw))
		stop(paste("Platform entered: ", platformDesignPackage, ", Platform expected: ", annotation(my.raw),"\n"))

	switch(unlist(strsplit(platformDesignPackage, "[.]"))[2],
		"ragene" ={if(is.null(orgPackage)) orgPackage <<- "org.Rn.eg"; if(is.null(specie)) specie <<- "Rattus norvegicus"; if(is.null(organisme)) organisme <<- "rno"},
		"raex" ={if(is.null(orgPackage)) orgPackage <<- "org.Rn.eg"; if(is.null(specie)) specie <<- "Rattus norvegicus"; if(is.null(organisme)) organisme <<- "rno"},
		"aragene" = {if(is.null(orgPackage)) orgPackage <<- "org.At.tair"; if(is.null(specie)) specie <<- "Arabidopis thaliana"; if(is.null(organisme)) organisme <<- "ath"},
		"bovgene" = {if(is.null(orgPackage)) orgPackage <<- "org.Bt.eg"; if(is.null(specie)) specie <<- "Bos taurus"; if(is.null(organisme)) organisme <<- "bta"},
		"cangene" = {if(is.null(orgPackage)) orgPackage <<- "org.Cf.eg"; if(is.null(specie)) specie <<- "Canis familiaris"; if(is.null(organisme)) organisme <<- "cfa"},
		"cyngene" = {if(is.null(orgPackage)) orgPackage <<- NULL; if(is.null(specie)) specie <<- "Macaca fascicularis"; if(is.null(organisme)) organisme <<- "mcf"},
		"cyrgene" = {if(is.null(orgPackage)) orgPackage <<- NULL; if(is.null(specie)) specie <<- "Macaca fascicularis"; if(is.null(organisme)) organisme <<- "mcf"},
		"equgene" = {if(is.null(orgPackage)) orgPackage <<- NULL; if(is.null(specie)) specie <<- "Equus caballus"; if(is.null(organisme)) organisme <<- "ecb"},
		"felgene" = {if(is.null(orgPackage)) orgPackage <<- NULL; if(is.null(specie)) specie <<- "Felis catus"; if(is.null(organisme)) organisme <<- "fca"},
		"hugene" = {if(is.null(orgPackage)) orgPackage <<- "org.Hs.eg"; if(is.null(specie)) specie <<- "Homo sapiens"; if(is.null(organisme)) organisme <<- "hsa"},
		"huex" = {if(is.null(orgPackage)) orgPackage <<- "org.Hs.eg"; if(is.null(specie)) specie <<- "Homo sapiens"; if(is.null(organisme)) organisme <<- "hsa"},
		"mogene" = {if(is.null(orgPackage)) orgPackage <<- "org.Mm.eg"; if(is.null(specie)) specie <<- "Mus musculus"; if(is.null(organisme)) organisme <<- "mmu"},
		"moex" = {if(is.null(orgPackage)) orgPackage <<- "org.Mm.eg"; if(is.null(specie)) specie <<- "Mus musculus"; if(is.null(organisme)) organisme <<- "mmu"},
		"ovigene" = {if(is.null(orgPackage)) orgPackage <<- NULL; if(is.null(specie)) specie <<- "Ovis orientalis"; if(is.null(organisme)) organisme <<- NULL},
		"porgene" = {if(is.null(orgPackage)) orgPackage <<- "org.Ss.eg"; if(is.null(specie)) specie <<- "Sus scrofa"; if(is.null(organisme)) organisme <<- "ssc"},
		"rcngene" = {if(is.null(orgPackage)) orgPackage <<- NULL; if(is.null(specie)) specie <<- "Oryza sativa"; if(is.null(organisme)) organisme <<- "osa"},
		"rhegene" = {if(is.null(orgPackage)) orgPackage <<- "org.Mmu.eg"; if(is.null(specie)) specie <<- "Macaca mulatta"; if(is.null(organisme)) organisme <<- "mcc"},
		"rjpgene" = {if(is.null(orgPackage)) orgPackage <<- NULL; if(is.null(specie)) specie <<- "Oryza sativa"; if(is.null(organisme)) organisme <<- "osa"},
		"soygene" = {if(is.null(orgPackage)) orgPackage <<- NULL; if(is.null(specie)) specie <<- "Glycine max"; if(is.null(organisme)) organisme <<- "gmx"},
		"zebgene" = {if(is.null(orgPackage)) orgPackage <<- "org.Dr.eg"; if(is.null(specie)) specie <<- "Danio rerio"; if(is.null(organisme)) organisme <<- "dre"},
		{if(is.null(orgPackage)) orgPackage <<- NULL; if(is.null(specie)) specie <<- NULL; if(is.null(organisme)) organisme <<- NULL})

if (is.null(orgPackage))
	cat("a genome wide annotation package must be assigned (orgPackage)\n")

if (is.null(specie) || is.null(organisme))
	cat("unknown organism, organism and species should be entered manually\n")
  
    save(my.raw, file = file.path(outputDir, dataFName))
  }else{
    load(file = file.path(outputDir, dataFName))
  }
if (TRY2INSTALL){
        installifnot(platformDesignPackage)
  	installifnot(old2db(orgPackage))
	}
	require(platformDesignPackage, lib.loc = mylib, character.only = TRUE)
	require(old2db(orgPackage), lib.loc = mylib, character.only = TRUE)
cat("Platform: ",platformDesignPackage, "\n", "Annotation Package: ", orgPackage, "\n", "Specie: ", specie, "\n", "Specie code: ", organisme, "\n")
return(my.raw)
}


#####################################################################################################################
### Funcions per a la representacio grafica del preprocessat
#####################################################################################################################

simpleplots <- function(my.data, sampleNames, my.colors, my.method = "average", targetsinfo, my.cex = 0.7)
{

  stopifnot(require(affy))

  ### Density histogram

  hist(my.data, main = "Signal distribution", col = my.colors, lty = 1:nrow(targetsinfo))
  legend(x = "topright", legend = sampleNames, col = my.colors, lty = 1:nrow(targetsinfo), cex = my.cex)

  ### Degradation Plots

  deg <- AffyRNAdeg(my.data)
  plotAffyRNAdeg(deg, col = my.colors)
  legend(x = "bottomright", legend = sampleNames, col = my.colors, lty = 1:nrow(targetsinfo), cex = my.cex)


  ### Boxplots

  opt <- par(las = 2, cex.axis = my.cex)
  boxplot(my.data, col = my.colors, names = sampleNames, las = 2 , main = "Raw data", cex = my.cex)
  par(opt)

  ### Dendrograma
  
  opt <- par(las = 2, cex = my.cex)
  clust.euclid.average <- hclust(dist(t(exprs(my.data))), method = my.method)
  plclust(clust.euclid.average, main = "Hierarchical clustering of samples",  labels = sampleNames, hang = -1, xlab = " ", sub = " ")
  par(opt)
}


#####################################################################################################################
### Funcions per al control de qualitat
#####################################################################################################################

affyQC <- function (myRaw, grup, nomGrup, sampleNames, computeQC=TRUE, tempDir=".", outputDir, QCFile, lFile){

  if (computeQC){
      my.mas5 <- call.exprs(myRaw[,grup], "mas5")
      my.qc<- qc(myRaw[,grup],my.mas5)
      save(my.mas5,file=file.path(tempDir,paste(nomGrup,"MAS5.expres.Rda",sep="")))
      save (my.qc, file=file.path(tempDir,paste(nomGrup,"MAS5.QC.Rda",sep="")))
    }else{
      load(file.path(tempDir,paste(nomGrup,"MAS5.expres.Rda",sep="")))
      load(file.path(tempDir,paste(nomGrup,"MAS5.QC.Rda",sep="")))
    }

  scale.fac<- my.qc@scale.factors
  perc.pres<- my.qc@percent.present
  back.avg <- my.qc@average.background
  back.rang<- my.qc@maximum.background-my.qc@minimum.background
  my.qc.sumry<- cbind(scale.fac,perc.pres,back.avg,back.rang)

#  AffyQCFile <- paste(nomGrup, "QCReport.pdf",sep="")
  AffyQCFile <- paste(nomGrup, ifelse(is.null(QCFile), "QCReport.pdf", QCFile),sep="")  
  pdf(file.path(outputDir, AffyQCFile))

  plot.new()
  titlePage(myRaw[,grup])
  plot(my.qc)

  par(mfrow=c(1,4))
  boxplot(scale.fac,main="SF")
  boxplot(perc.pres,main="% Pres.")
  boxplot(back.avg,main="Mean bck")
  boxplot(back.rang,main="Bck range")
  boxplot(scale.fac,main="SF")
  text(1,scale.fac,sampleNames,cex=0.6)
  boxplot(perc.pres,main="% Pres.")
  text(1,perc.pres,sampleNames,cex=0.6)
  boxplot(back.avg,main="Mean bck")
  text(1,back.avg,sampleNames,cex=0.6)
  boxplot(back.rang,main="Bck range")
  text(1,back.rang,sampleNames,cex=0.6)
  par(mfrow=c(1,1))

  deg <- AffyRNAdeg(myRaw)
  numSamples <- nrow(pData(myRaw))
  plotAffyRNAdeg(deg,col=1:numSamples,lty=1:numSamples)
  legend(x="bottomright", legend=sampleNames, col=1:numSamples,
         lty=1:numSamples, cex=0.7)

  borderQC1(myRaw[,grup])
  borderQC2(myRaw[,grup])
  correlationPlot(myRaw[,grup])

  dev.off()

  addToLinksFile(lFile, AffyQCFile, categ = 'QC', desc = paste("Standard Affymetrix QC Report for", nomGrup, sep=" "))

  return(my.qc.sumry)
}


processaQC <- function(QCgrups, readCELS, fileNames, QCFile=NULL, lFile, phenoDat, sampleNames, computeQC, tempDir, outputDir, csvType=NULL)
{
#  stopifnot(require(Biobase))
#  stopifnot(require(simpleaffy))
#  stopifnot(require(affyQCReport))  ### Probablement no calgui aquest paquet

  my.qc.sumry <- NULL
  if (is.null(QCgrups))
    QCgrups <- list(gAll=1:nrow(pData(phenoDat)))
  for (i in 1:length(QCgrups)){
      nomGrup = names(QCgrups)[i]
      selected <- QCgrups[[i]]
      dataFileName <- paste("rawData", nomGrup, "Rda", sep = ".")

      if (readCELS){
          print(fileNames[selected])

          if (exists("my.raw")){
            rm(my.raw)
            gc()
          }

          my.raw <- ReadAffy(filenames = fileNames[selected], phenoData = phenoDat[selected,], verbose = TRUE)

          save(my.raw, file = file.path(tempDir, dataFileName))

      }else{
         load(file.path(tempDir, dataFileName))
      }

      my.qc <- affyQC(my.raw,
                      grup = 1:nrow(pData(my.raw)),
                      nomGrup = nomGrup,
                      computeQC = computeQC,
                      sampleNames = sampleNames[selected],
                      tempDir = tempDir,
                      outputDir = resultsDir,
                      QCFile = QCFile,
                       lFile = lFile)

      # QCgrups es una llista definida al fitxer "parametros.R"

      if (i==1){
         my.qc.sumry <- my.qc
      }else{
        my.qc.sumry <- rbind(my.qc.sumry, my.qc)
      }
  }

  QCSummaryFile <- "AffyQCSummary"
#  qcSummaryFName <- file.path(outputDir, QCSummaryFile)
#  write.csv2(my.qc.sumry, file = qcSummaryFName)
  csvType<-ifelse(is.null(csvType), "csv2", csvType)
  write2csv(my.qc.sumry, fileName = QCSummaryFile, csv = csvType, outputDir = outputDir)
  addToLinksFile(lFile, paste(QCSummaryFile,  substr(csvType[1], 1, 3), sep="."), categ = 'QC', desc = "Standard Affymetrix QC parameters")

}


affyPLM <- function(my.raw, fitPLModel = TRUE, sampleNames, my.colors, my.cex.axis = 0.8)
{
  if (fitPLModel){
    Pset <- fitPLM(my.raw)
    save(Pset, file = file.path(tempDir, "PLM.Rda"))
  }else{
    load(file = file.path(tempDir, "PLM.Rda"))
  }

  opt <- par(las = 2,cex.axis = my.cex.axis)

  ### RLE

  val <- RLE (Pset, type = "values")
  ylim <- range(val)
  RLE(Pset, main = "Relative Log Expression", names = sampleNames, col = my.colors, ylim = ylim)

  ### NUSE

  val <- NUSE(Pset, type = "values")
  ylim <- range(val)
  NUSE(Pset, main = "Normalized Unscaled Standard Errors", names = sampleNames, col = my.colors, ylim = ylim)

  par(opt)

  rm(Pset)  ; gc()
}


#####################################################################################################################
### Funcions per a la normalitzacio
#####################################################################################################################

### Selecciona i aplica un metode de normalitzacio (RMA, GCRMA o MAS5)

normalitza <- function (my.data, method)
{
  switch(method,
         RMA =  rma(my.data),
         GCRMA = gcrma(my.data),
         MAS5 = mas5 (my.data))
}


### Normalitzacio de les dades. Per defecte usa el metode JUSTRMA

normalization <- function(my.data = NULL,
                          method = NULL,
                          targetsinfo,
                          inputDir,
                          loadFile = FALSE,
                          normalizedFName = FALSE,
                          outputDir,
                          exonSt=FALSE)
{
  fName <- ifelse(is.null(normalizedFName),
                  paste("normalized", method, "Rda", sep = "."),
                  normalizedFName)
  if (!loadFile)
  {
    if(exonSt)
    {
      my.norm <- rma(my.data, target = ifelse(is.null(method), "core", method))
    }else{
      
#    stopifnot(require(affy))
    
      if (is.null(method) || is.null(my.data) || method == "JUSTRMA")
      {
        my.norm <- just.rma(filenames = file.path(inputDir, rownames(pData(targetsinfo))), phenoData = targetsinfo, verbose = TRUE)
      }else{
        my.norm <- normalitza(my.data, method)
      }
    }

    save(my.norm, file = file.path(outputDir, fName))
  }else{
    load(file = file.path(outputDir, fName))
  }
  
  return(my.norm)
}


#####################################################################################################################
### Funcions per a fer anotacions
#####################################################################################################################

### Anota una matriu d'expressio (normalitzada) amb Gene Symbols  (o altres IDs)

expresAndgeneSymbols <- function(expres, expresNames=colnames(expres),
                          anotPackage = NULL, SYMBOL="SYMBOL",
                          symbolsVector = NULL)
{
  if (!is.null(anotPackage))
  {
    my_SYMBOL_env <- eval(parse(text = paste(anotPackage, SYMBOL, sep="")))
    geneSymbols <- unlist(mget(rownames(expres), my_SYMBOL_env, ifnotfound=NA))
    expresWithSymbols <- data.frame(geneSymbols, expres)
    colnames(expresWithSymbols) <- c("Symbol", expresNames)
  }else{
      if (!is.null(symbolsVector)){
          geneSymbols <- symbolsVector[rownames(expres)]
          expresWithSymbols <- data.frame(geneSymbols, expres)
          colnames(expresWithSymbols) <- c("Symbol", expresNames)
      }else{
        expresWithSymbols <- expres
      }
  }
  return(expresWithSymbols)
}


### creaAnot: Construeix una taula que conte per a cada identificador (ProbeID/TranscriptID)
###           l'anotacio demanada 
###                 
###
### Parametres:
###
###        inputDir : Nom del direcori d'entrada d'on es llegira el fitxer d'anotacions 
###   annotFileName : Fitxer d'anotacions.
###           field : Nom del camp del fitxer d'anotacions que es vol extreure
###
### OBERVACIONS IMPORTANTS:
###
###   Cada linia del fitxer d'anotacions conte la informacio per a un unic gen. Les anotacions ha
###   d'estar tabulades
###

creaAnot <- function (inputDir, annotFileName, field=1)
{
  x <- read.delim(file.path(inputDir, annotFileName))
  myAnnot <- as.character(x[, field])
  names(myAnnot) <- as.character(x[, 1])  # Podria ser que dossin els rowNames, no?
  return(myAnnot)
}

### Actualitzacio 2011. Bloc provinent de 'manegaAnotacions.R'

###########################################################################################
# L'objectiu d'aquestes funcions es permetre crear les taules d'anotacions, es a dir
#   les 'EntrezTable' o 'SymbolsTable' que utilitzarem despres per anotar els resultats
#
# Hi ha tres versions segons a partir de que es creen
# 1) creaAnotFromText: Parteix d'un arxiu de text com els csv que obtenim d'affymetrix
# 2) creaAnotFromChipPackage: Parteix d'un paquet d'anotacions tipic de Bioconductor
#                      com el 'hgu133plus2.db' o els 'transcriptXXX.db'
# 3) creaAnotFromPDPackage: Parteix d'un paquet de "PlatformDesign" com els desenvolupats
#                      pel paquet 'oligo'
# Cada funcio necessita, uns parametres diferents (revisar-ho) pero retorna el mateix
# Retornen un vector de caracters que te:
# - Els noms dels identificadors (sondes/transcripts) com a 'Names' i
# - ELs noms dels entrezs o symbols o els que sigui com valors
########################################################################################
#
# Seguint la mateixa filosofia sembla raonable que el filtratge dels controls
# es pugui fer a partir d'un vector de caracters que es passi a una funcio de filtratge
# Aquest vector de caracters s'hauria de poder obtenir a partir del mateix paquet d'anotacions.
# La funcio de filtratge hauria de rebre la matriu d'expressions i el vector de controls i
# retornar la matriu "reduida" es a dir sense els controls.
#
# Probablement es podran estendre les funcions d'anotacio per a que retornin els controls
#
##############################################################################################

creaAnotFromText <- function (inputDir, annotFileName, field,
                              NAsymbol="---", cleanNAs = T,
                              multipleIDsSymbol=" /// ", removeMultipleIDs=T,
                              isControl= FALSE, ctlCode=NA, nonCtlCode=NA)
{
  x <- read.delim(file.path(inputDir, annotFileName))
  myAnotTable <- as.character(x[, field])
  names(myAnotTable) <- as.character(x[, 1])
    # Al llegir la taula amb read.delim la 1a columna conte els noms de
    # les sondes/transcripts que estaven en la 1a columna de l'arxiu de text
    # Si afegissim row.names=T la 1a columna passaria a ser els rownames...
  
  if (cleanNAs){
    if (length(grep("-", myAnotTable))>0)
       myAnotTable <- myAnotTable[-which(myAnotTable == NAsymbol)]
  }
  if (removeMultipleIDs){
    if (length(grep(multipleIDsSymbol, myAnotTable))>0) {
       if (runMulticore ==1 || runMulticore ==3) { 
         myAnotTable <- unlist(mclapply(strsplit(myAnotTable, multipleIDsSymbol),
                                    my.fun <- function(x){return(x[1])}))
       } else {
         myAnotTable <- unlist(lapply(strsplit(myAnotTable, multipleIDsSymbol),
                                    my.fun <- function(x){return(x[1])}))         
       }
    }
  }
  if (isControl){
    if ((is.na(ctlCode) && is.na(nonCtlCode))||  (!is.na(ctlCode) && !is.na(nonCtlCode)))
      stop("Either Control or Non-Control codes must be different from NA")     
    if (!is.na(ctlCode)){
      myAnotTable<- myAnotTable[myAnotTable==ctlCode]
    }else{
      myAnotTable<- myAnotTable[myAnotTable !=nonCtlCode]
    }
  }
  return(myAnotTable)
}

############################################################################
# Exemple 1: Creacio de la taula d'anotacions a partir d'un arxiu de text  
# Fins aquest canvi es el que feiem amb els arrays d'exons
############################################################################

prova1 <- function(){
  inpD <- "/media/TREKSTOR/Estudis/Microarrays/2011-01-MireiaPares-VHIR/dades"
  anotF <- "annotations.txt"
  fld <- 2 # El 2 es el Gene Symbol. El 3 es l'Entrez

  EntrezTable1 <- creaAnotFromText(inputDir =inpD, field=fld, annotFileName=anotF)
  EntrezTable13 <- creaAnotFromText(inputDir =inpD, field=fld, annotFileName=anotF, removeMultipleIDs=F)


  inpD <- "/media/TREKSTOR/Estudis/Microarrays/2011-01-MireiaPares-VHIR/dades"
  anotF <- "annotations.txt"
  fld <- 3 # El 2 es el Gene Symbol. El 3 es l'Entrez

  SymbolTable11 <- creaAnotFromText(inputDir =inpD, field=fld, annotFileName=anotF) 
  SymbolTable12 <- creaAnotFromText(inputDir =inpD, field=fld, annotFileName=anotF, cleanNAs=F)
  SymbolTable13 <- creaAnotFromText(inputDir =inpD, field=fld, annotFileName=anotF, removeMultipleIDs=F)

  inpD <- "/media/TREKSTOR/Estudis/Microarrays/2011-01-MireiaPares-VHIR/dades"
  anotF <- "annotations.txt"
  fld <- 2 # El 2 es el Gene Symbol.
         # Per fer la prova suposarem que els controls son els "---"
  isControl <- TRUE
  controlCode <- "---"
  cleanNAs<- FALSE

  controlsTable12 <- creaAnotFromText(inputDir =inpD, field=fld, annotFileName=anotF, cleanNAs=F,  isControl = T, ctlCode=controlCode)
}

#
# Per revisar : MultiplesIDs: " /// "?
# Si: El que es separa per "//" son els identificadors que anoten conceptes diferents com son l'accession, o el Gene Symbol. Cada grup d'Ids d'un mateix transcript va separat per tres "///"
# Dit aixo el que sembla logic es que un cop seleccionats els escrivim separats per dos barres, pero aixo es pot acabar de concretar.
#


#############################################################################

creaAnotFromChipPackage <- function (chipPackage, field="ENTREZ", cleanNAs=T,
                                      isControl= FALSE, ctlCode=NA, removeControls=FALSE)
{
  if (!require(old2db (chipPackage), character.only=T)){
    stop(paste("Required annotation package", chipPackage," is missing"))
  }else{
    require(old2db (chipPackage), character.only=T) # No deu caldre
  }
  if (isControl){
    cleanNAs<- FALSE
    field <- "ACCNUM"
  }
  name2extract <- paste(chipPackage, field, sep="")
  x<- eval(parse(text= name2extract))
  myAnot <- toTable(x) # Aqui hauria de venir un control d'error
                       # per si se li ha donat el nom malament
  myAnotTable <- myAnot[,2]
  names(myAnotTable) <- myAnot[,1]
  if (!cleanNAs){
     name2extract <- paste(chipPackage, "ACCNUM", sep="")
     x<- eval(parse(text= name2extract))
     myAccTable <- toTable(x)
     myNAsNames <- setdiff(myAccTable[,1], myAnot[,1])
                   # Compte! : Aquesta diferencia es assimetrica
     myNAs <- rep(NA, length(myNAsNames))
     names(myNAs) <- myNAsNames
     myAnotTable<- c(myAnotTable, myNAs)
  }
   if (isControl){
    if (is.na(ctlCode))
      stop("Control codes must be different from NA")
    controls <-sapply (myAnotTable, function (x) ifelse( regexpr(ctlCode, x) > 0, TRUE, FALSE ))
    if (!removeControls){
       myAnotTable<- myAnotTable[controls]
    }else{
      myAnotTable<- myAnotTable[!controls]
    }
  }
  return(myAnotTable)
}


##############################################################################
# Exemple 2: Creacio de la taula d'anotacions a partir d'un paquet d'anotacions  
# Fixem-nos que aqui el camp es dona coma cadena ('ENTREZ' o 'STRING')
# D'altra banda no cal un directori de treball
##############################################################################

prova2<-function(){
  source("./exemples/LPS-celltypes/RCode/AnalysisFunctions2Pack.R")
  anotP <- "hgu133plus2"
  fld1 <- "ENTREZID"
  fld2 <- "SYMBOL"

  EntrezTable21 <- creaAnotFromChipPackage (chipPackage =anotP, field=fld1)
  SymbolTable21 <- creaAnotFromChipPackage (chipPackage =anotP, field=fld2)

  controlsTable21 <- creaAnotFromChipPackage (chipPackage =anotP, field=fld2,
                                            isControl= TRUE, ctlCode='AFFX',
                                            removeControls=FALSE)
  controlsTable22 <- creaAnotFromChipPackage (chipPackage =anotP, field=fld2,
                                            isControl= TRUE, ctlCode='AFFX',
                                            removeControls=TRUE)
}

#############################################################################

creaAnotFromPDPackage <- function (dbPackage, field, fieldName=NULL, cleanNAs=T,
                                   multipleIDsSymbol=" /// ", removeMultipleIDs=T,
                                   removeControls=TRUE)
{
  #
  # dbPackage podria ser per exemple podria "pd.huex.1.0.st.v2"
  # 'field' es un nombre de 1 a 5 que descriu a
  #  quina de les parts del camp geneassigment s'extraura
  #      1. accession       - public sequence identifier for mRNA
  #      2. gene symbol     - gene name if mRNA corresponds to a known gene
  #      3. gene title      - description of gene product 
  #      4. cytoband        - cytogenetic location of gene
  #      5. entrez gene id  - Entrez Gene database identifier
  #
  
  if (!require(dbPackage, character.only=T)){
    stop(paste("Required Platfform design package", dbPackage," is missing"))
  }else{
    require(dbPackage, character.only=T) # No deu caldre
  }
  #
  # Falta implementar el cleanNAs
  # Comencem per definir la llista de transcripts, que no son controls
  #
    conn<- db(eval(parse(text=dbPackage)))
    fSetType <- dbGetQuery(conn,
            paste("SELECT DISTINCT meta_fsetid as transcript_id, type_id",
                  "FROM featureSet, core_mps, type_dict",
                  "WHERE featureSet.fsetid=core_mps.fsetid",
                  "AND featureSet.type=type_dict.type"))
      allExceptControls <- as.character(fSetType$transcript_id[fSetType$type_id=="main"])
      allControls <- as.character(fSetType$transcript_id[fSetType$type_id!="main"])
   
   # Un cop fet aixo podem mirar de recuperar els ids fent servir getNetAffx
   # getNetAffx necessita un ExpressionSet (suposa que ja ha normalitzat)
   # El fem a ma posant-li els noms dels transcripts en el camp assayData

   if  (removeControls){
     transcriptIds <- allExceptControls
   }else{
     transcriptIds <- allControls
     field <- 1 # Es l'ACCNUM que tambe existeix pels controls
     removeMultipleIDs <- FALSE
   }
  
   numTranscripts <-length(transcriptIds)
   exp <-matrix(NA, nrow=numTranscripts, ncol=1,
                   dimnames=list(transcriptIds,"Empty"))
   rownames(exp) <- transcriptIds
   aD <-assayDataNew(exprs=exp)
   nulEset <-new("ExpressionSet",  assayData =aD)
   featureNames(nulEset) <- transcriptIds
   annotation(nulEset) <-dbPackage

   # Ara invoco getNetAffx per recuperar les anotacions

   featureData(nulEset) <- getNetAffx (nulEset, "transcript")
   geneNames <-pData(featureData(nulEset))$geneassignment

   # Creem la taula d'anotacions
   
   anotTable <- data.frame(transcriptIds,
                           rep(NA, length(transcriptIds))
                           )
   if (is.null(fieldName))
       fieldName<- switch(field,
          "Accession","GeneSymbol", "Gene Title", "Cytoband","Entrez")
   colnames (anotTable) <- c("transcriptIds", fieldName)

   myAnotTable <-geneNames
   names (myAnotTable)<-transcriptIds

   stopifnot(require(gdata))
   for (i in 1:length(transcriptIds)){
     if (!is.na(geneNames[i])){
       l1<- strsplit(geneNames[i],"///")
       if (runMulticore ==1 || runMulticore ==3) { 
          l2 <- mclapply(l1, function(l) strsplit(l, "//"))
          myIds <- unique(unlist(mclapply(l2[[1]], function (x) try(trim(x[[field]])))))
         } else {
          l2 <- lapply(l1, function(l) strsplit(l, "//"))
          myIds <- unique(unlist(lapply(l2[[1]], function (x) try(trim(x[[field]])))))
       }


       anotTable[i,2] <- paste(myIds, collapse="//")
     }
   }
   if(removeMultipleIDs){
     s1<-sapply(anotTable[,2], function(s) strsplit(s, "//"))
       if (runMulticore ==1 || runMulticore ==3) { 
          s2<-mclapply (s1, function(s) return(s[1]))
         } else {
          s2<-lapply (s1, function(s) return(s[1]))         
       }

     myAnotTable <-unlist(s2)
     names(myAnotTable)<-transcriptIds
   }else{
     myAnotTable <-anotTable[,2]
     names(myAnotTable)<-transcriptIds
   }
  return(myAnotTable)
}



#######################################################################################
# Exemple 3: Creacio de la taula d'anotacions a partir d'un paquet de "Platform Design"  
# Fixem-nos que aqui el camp es dona com enter de 1 a 5 decidint que volem treure
######################################################################################

prova3<-function(){
  dbPack <- "pd.huex.1.0.st.v2"

  EntrezTable31 <-creaAnotFromPDPackage (dbPackage=dbPack, field=5, fieldName=NULL)
  EntrezTable32 <-creaAnotFromPDPackage (dbPackage=dbPack, field=5, fieldName=NULL,
                                       removeMultipleIDs=F)

  SymbolsTable31 <-creaAnotFromPDPackage (dbPackage=dbPack, field=2, fieldName=NULL)
  SymbolsTable32 <-creaAnotFromPDPackage (dbPackage=dbPack, field=2, fieldName=NULL,
                                        removeMultipleIDs=F)

  AccnumTable31 <-creaAnotFromPDPackage (dbPackage=dbPack, field=1, fieldName=NULL,
                                       removeControls=F)
  AccnumTable32 <-creaAnotFromPDPackage (dbPackage=dbPack, field=1, fieldName=NULL,
                                       removeControls=T)
}


#########################################################################


createOrLoadAnnotations <-function (loadAnnotations=FALSE,
                                    chipPackAvailable,
                                    platformDesignPackAvailable,
                                    chipPackage,
                                    platformDesignPackage,
                                    outputDir,
                                    annotationsFileName,
                                    entrezTableFileName,
                                    symbolsTableFileName,
                                    controlsTableFileName
                                    )
{
  if(!loadAnnotations){
  #
  # Primer es construeixen les anotacions
  #
    if (chipPackAvailable){
    # De moment no ho he tocat pero ho fara igual que el seguent
    entrezTable  <-creaAnotFromChipPackage (chipPackage=chipPackage,
                                            field='ENTREZID')
    symbolsTable <-creaAnotFromChipPackage (chipPackage=chipPackage,
                                            field='SYMBOL')
    controlsTable <-creaAnotFromChipPackage (chipPackage =chipPackage,
                                             field='ACCNUM',
                                             isControl= TRUE,
                                             ctlCode='AFFX',
                                             removeControls=FALSE)
  }else{
    if (platformDesignPackAvailable){
      entrezTable  <-creaAnotFromPDPackage (dbPackage=platformDesignPackage,
                                            field=5)
      symbolsTable <-creaAnotFromPDPackage (dbPackage=platformDesignPackage,
                                            field=2,
                                            removeMultipleIDs=F)
      controlsTable <-creaAnotFromPDPackage (dbPackage=platformDesignPackage,
                                             field=1, 
                                             removeControls=F)
    }else{
      symbolsTable <- NULL          # Si hi ha chipPackage ha de ser NULL
      entrezTable <- NULL           # Si hi ha chipPackage ha de ser NULL
      controlsTable <- NULL
    }
  }
  }else{
    entrezTable <-loadFromFile(file=file.path(outputDir, entrezTableFileName))
    symbolsTable <-loadFromFile(file=file.path(outputDir, symbolsTableFileName))
    controlsTable <- loadFromFile(file=file.path(outputDir, controlsTableFileName))
  }
  #
  # Un cop construides es desen en arxius binaris a la carpeta "dades"
  #
    if(!is.null(entrezTable)){
      save(entrezTable, file=file.path(outputDir, entrezTableFileName)) }
    if(!is.null(symbolsTable)){
      save(symbolsTable, file=file.path(outputDir, symbolsTableFileName))    }
    if(!is.null(controlsTable)){
      save(controlsTable, file=file.path(outputDir, controlsTableFileName))  }
    if (is.null(annotationsFileName)){
        annotationsFileName <- paste("Annotations", "txt", sep=".")
      write.table(data.frame(Entrez=entrezTable, Symbols=symbolsTable),
                  sep="\t", quote=FALSE, file=file.path(outputDir,annotationsFileName))  }
  anotacions=list(Entrez=entrezTable,Symbols=symbolsTable, controls=controlsTable)
  return(anotacions)
}

prova4<- function(){
 source("/media/TREKSTOR/microarrays/Basic-Pipe/Basic.Pipe.Package/R/AnalysisFunctions2Pack.R")
 anota11 <- createOrLoadAnnotations (loadAnnotations=FALSE,
                                    chipPackAvailable =TRUE,
                                    platformDesignPackAvailable=FALSE,
                                    chipPackage="hgu133plus2",
                                    platformDesignPackage=NULL,
                                    outputDir = "temp",
                                    annotationsFileName ="AnnotationsHGU133PLUS2",
                                    entrezTableFileName = "EntrezsHGU133.Rda",
                                    symbolsTableFileName = "SymbolsHGU133.Rda",
                                    controlsTableFileName = "controlsHGU133.Rda"
                                    )
  anota12 <- createOrLoadAnnotations (loadAnnotations=TRUE,
                                    chipPackAvailable =TRUE,
                                    platformDesignPackAvailable=FALSE,
                                    chipPackage="hgu133plus2",
                                    platformDesignPackage=NULL,
                                    outputDir = "temp",
                                    annotationsFileName ="AnnotationsHGU133PLUS2",
                                    entrezTableFileName = "EntrezsHGU133.Rda",
                                    symbolsTableFileName = "SymbolsHGU133.Rda",
                                    controlsTableFileName = "controlsHGU133.Rda"
                                    )
  anota21 <- createOrLoadAnnotations (loadAnnotations=FALSE,
                                    chipPackAvailable =FALSE,
                                    platformDesignPackAvailable=TRUE,
                                    chipPackage=NULL,
                                    platformDesignPackage="pd.huex.1.0.st.v2",
                                    outputDir = "temp",
                                    annotationsFileName ="Annotations.pd.huex",
                                    entrezTableFileName = "Entrezs.pd.huex.Rda",
                                    symbolsTableFileName = "Symbols.pd.huex.Rda",
                                    controlsTableFileName = "controls.pd.huex.Rda"
                                    )
anota21 <- createOrLoadAnnotations (loadAnnotations=TRUE,
                                    chipPackAvailable =FALSE,
                                    platformDesignPackAvailable=TRUE,
                                    chipPackage=NULL,
                                    platformDesignPackage="pd.huex.1.0.st.v2",
                                    outputDir = "temp",
                                    annotationsFileName ="Annotations.pd.huex",
                                    entrezTableFileName = "Entrezs.pd.huex.Rda",
                                    symbolsTableFileName = "Symbols.pd.huex.Rda",
                                    controlsTableFileName = "controls.pd.huex.Rda"
                                    )
}

### Bloc provinent de 'manegaAnotacions.R'

#####################################################################################################################
### Funcions per fer els grafics normalitzats
#####################################################################################################################

normplots2File <- function(my.data,
                           sampleNames,
                           my.colors,
                           my.groups,
                           my.method = "average",
                           my.cex = 0.7,
                           posText = 4,
                           dim3 = FALSE,
                           fileName,
                           outputDir,
                           PCAPlots=TRUE,
                           csv,
                           lFile = NULL)
{
if(!is.null(fileName)) pdf(file.path(outputDir,fileName))

  normplots (my.data=my.data, sampleNames= sampleNames,my.colors= my.colors, 
     my.groups= my.groups, my.method = my.method, my.cex = my.cex, posText = posText, dim3 = dim3, PCAPlots=PCAPlots, outputDir=outputDir)

  if(!is.null(fileName)) dev.off()

  addToLinksFile(linksFile=lFile, fileName, categ = 'QC', desc = "Plots of normalized data")
}


###############################################################################
### Funcions per a per fer els grafics de les dades normalitzades
###############################################################################


normplots <- function(my.data, sampleNames, my.colors, my.groups, my.method = "average", my.cex = 0.7, posText = 4, dim3 = FALSE, PCAPlots=TRUE, outputDir, csv=csv)
{

  ### Boxplots

  opt <- par(las=2, cex.axis = my.cex)
  boxplot(my.data, col = my.colors, names = sampleNames, las = 2 , main = "Normalized (RMA) data", cex.axis = my.cex)
  par(opt)

  ### Dendrograma
  
  opt <- par(las=2, cex = my.cex)
  clust.euclid.average <- hclust(dist(t(exprs(my.data))), method = my.method)
  plclust(clust.euclid.average, main = "Hierarchical clustering of samples",  labels = sampleNames, hang = -1, xlab = " ", sub = " ")
  par(opt)

  ### PCA

  doPCAplot(my.data, sampleNames, my.colors, my.groups, my.cex, PCAFFile, outputDir=outputDir, cor, comp, posText=posText,
            dim3=dim3, csv=csv, my.PCAplot=PCAPlots, x.coord = -100, y.coord = 100, pch = 19)

}


### doPCAplot: Representa graficament la PCA 2D (i 3D si es vol) el screeplot i desa en un fitxer csv
###            la informacio relativa a la variabilitat explicada per cada component
###
### Parametres:
###
###        inputDir : Nom del direcori d'entrada d'on es llegira el fitxer d'anotacions 
###   annotFileName : Fitxer d'anotacions.
###           field : Nom del camp del fitxer d'anotacions que es vol extreure
###
###         my.data : Expressions normalitzades (eset_norm)
###     sampleNames : ShortNames del targets
###       my.colors : Colores del targets
###       my.groups : Grupo del targets
###          my.cex : Valor del cex. Per defecte pren 0.8
###        PCAFFile : Nom del fitxer que contrindra la informacio relativa a la variabilitat explicada per cada component
###       outputDir : Directori a on es desara la informacio
###             cor : Per defecte es TRUE => Usara les correlacions per calcular les components principals 
###            comp : vector que contindra les components a representar. Per defecte pren la 1a i la 2a
###         posText : Posicio del text a la respresentacio
###            dim3 : Per defecte FALSE => No representar graficament la PCA 3D
###             csv : TRUE si es vol en format csv
###      my.PCAplot : Per defecte TRUE => Usem PCA plot nostre
###         x.coord = -100
###         y.coord = 100
###             pch : Tipus de punts que es representaran al grafic
###
###
### OBERVACIONS IMPORTANTS:
###

doPCAplot <- function(my.data, 
                      sampleNames ,
                      my.colors ,
                      my.groups ,
                      my.cex = 0.8,
                      PCAFFile,
                      outputDir,
                      cor = TRUE,
                      comp = c(1, 2),
                      posText = 4,
                      dim3 = FALSE,
                      csv,
                      my.PCAplot = TRUE,
                      x.coord = -100,
                      y.coord = 100,
                      pch = 19)
{
  if (my.PCAplot)
  {
   # pc.my.norm <- prcomp(t(exprs(my.data)), cor)
    pc.my.norm <- prcomp(t(exprs(my.data))) # , cor) --> Desactivat perque dona error
    pc.importance <- summary(pc.my.norm)$importance
    print(round(pc.importance, 3))
    
    # write2csv(pc.importance, fileName = PCAFile, csv, outputDir)  
        
    comp <- c(1, 2)
    plotPCA2D(pc.my.norm, pc.importance, sampleNames, my.colors, comp, posText)

    if (dim3)
    {
      comp <- c(1, 2, 3)
      plotPCA3D(pc.my.norm, pc.importance, sampleNames, my.colors, comp)
    }
    
    plotPCA(exprs(my.data), screeplot = TRUE)
  }else{
    plotPCA(my.data, groups = my.groups, groupnames = sampleNames, addtext = sampleNames,
            #x.coord = -100, y.coord = 100, 
            squarepca = TRUE, pch=pch, col = my.colors, screeplot=FALSE)
  }
}


adequancyPlot <- function(pc.importance, epsilon = 0.05)
{
  cum.prop <- 1 - pc.importance[3,]

  plot(1:dim(pc.importance)[2], cum.prop,
       main = "Adequacy of the k-dimensional representation for the PCA",
       xlab = "Principal components",
       ylab = "% of variability explained",
       type = "b", col = "blue")

  pos <- which.max(pc.importance[2,] <= epsilon)
  cutoff <- (1 - pc.importance[3, pos - 1]) - (pc.importance[2, pos] / 2)
  abline(h=cutoff, col="red")
}


plotPCA2D <- function(my.data, pc.importance, my.names, my.colors, comp = c(1, 2), posText = 4)
{
  stopifnot(require(rgl))

  scores <- my.data$x
  rownames(scores) <- my.names

  Xlims = c(min(scores[, comp[1]]-10), max(scores[, comp[1]]))
  Ylims = c(min(scores[, comp[2]]), max(scores[, comp[2]]))

  plot(scores,
         main = "Principal Components 2D Plot",
         xlab = paste("PC", comp[1], " ", round(pc.importance[2, comp[1]]*100, 1), "%", sep = ""),
         ylab = paste("PC", comp[2], " ", round(pc.importance[2, comp[2]]*100, 1), "%", sep = ""),
         xlim = Xlims, ylim = Ylims, type = "n")

  for(i in 1:dim(scores)[1])
  {
    points(scores[i, comp[1]],
           scores[i, comp[2]],
           pch = 19, col = my.colors[i])
  }

  text(scores[, comp[1]],
       scores[, comp[2]],
       my.names, cex = 0.6, pos = posText)
}


plotPCA3D <- function(my.data, pc.importance, my.names, my.colors, comp=c(1, 2, 3))
{
  stopifnot(require(rgl))

  scores <- my.data$x
  rownames(scores) <- my.names

  Xlims = c(min(scores[, comp[1]]), max(scores[, comp[1]]))
  Ylims = c(min(scores[, comp[2]]), max(scores[, comp[2]]))
  Zlims = c(min(scores[, comp[3]]), max(scores[, comp[3]]))

  plot3d(scores,
         main="Principal Components 3D Plot",
         xlab=paste("PC", comp[1], " ", round(pc.importance[2, comp[1]]*100, 1), "%", sep = ""),
         ylab=paste("PC", comp[2], " ", round(pc.importance[2, comp[2]]*100, 1), "%", sep = ""),
         zlab=paste("PC", comp[3], " ", round(pc.importance[2, comp[3]]*100, 1), "%", sep = ""),
         xlim = Xlims, ylim = Ylims, zlim = Zlims,
         type = "n", size = 4)

  for(i in 1:dim(scores)[1])
  {
    points3d(scores[i, comp[1]],
             scores[i, comp[2]],
             scores[i, comp[3]],
             size = 4, pch = 19, col = my.colors[i])
  }

  text3d(scores[, comp[1]],
         scores[, comp[2]],
         scores[, comp[3]],
         my.names, cex = 0.6, pos = 4)
}

##############################################################################################
### Funcions de filtratge
##############################################################################################

# Filtratge per prefix: treu tots els gens que comencen per una cadena

filterByString <- function(expres, filterString="AFFX")
{
  gN <- rownames(expres)
  notAffx <- sapply (gN, function (x) ifelse( regexpr(filterString, x) < 0, TRUE, FALSE ))
  expres.filtered <- expres[notAffx, ]
  return(expres.filtered)
}


filterAffy <- function(expres)
{
  gN <- rownames(expres)
  notAffx <- sapply (gN, function (x) ifelse(substr(x, 1, 4) != "AFFX", TRUE, FALSE ))
  expres.filtered <- expres[notAffx, ]

  return(expres.filtered)
}


# Filtratge basat en les intensitats

filt.by.Signal <- function(x, grups, threshold)
 {
   if( max(sapply(split(x, grups), mean), na.rm=TRUE) < threshold) return(FALSE) else return(TRUE)
 }

filterBySignal <- function(expres, groups, threshold, sigFun.Name="filt.by.signal", thr.as.perc=TRUE)
{
  if (thr.as.perc){
      signalThr <- quantile(as.vector(expres), threshold / 100)
  }else{
      signalThr <- threshold
  }
  sigFun=eval(parse(text=sigFun.Name))
  filtered.vals <- apply (expres, 1, sigFun, groups, signalThr)
  expres.filtered <- expres[filtered.vals, ]
  return(expres.filtered)
}

# Filtratge basat en la variabilitat

iqrf <- function(x, threshold) {
  #percs <- quantile(x,c(0.25,0.75))
  #if (diff(range(percs)) < threshold)
  if (IQR(x) < threshold){
    return(FALSE)
  }else{
    return(TRUE)
  }
}

sdf <- function (x, threshold) if (sd(x) < threshold) return(FALSE) else return(TRUE)

filterByVar <- function (expres, threshold, varFun.Name="sdf", thr.as.perc=TRUE, thr.fun=sd)
{
   if (thr.as.perc){
      SD <- apply(expres, 1, thr.fun)
      variabilityThr <- quantile (SD,  threshold / 100) # SD3Q
  }else{
      variabilityThr <- threshold
  }
  varFun=eval(parse(text=varFun.Name))
  filtered.vals <- apply (expres, 1, varFun, variabilityThr)
  expres.filtered <- expres[filtered.vals, ]
  return(expres.filtered)
}

check2Filter <- function (expres, toPlot=T)
{
      print("Quantiles of Intensity Values in log2 and ordinary scale")
      print(round(quants.m <- quantile(as.vector(expres), seq(0, 1, by = 0.05)), digits = 3))
      print(round(2^quants.m, digits = 3))    
      SD <- apply(expres, 1, sd)
      cat("\nQuantiles of Standard Deviations\n")
      print(round(quants <- quantile(SD, seq(0, 1, by = 0.05)), digits = 3))
      if (toPlot){
        opt<- par(mfcol=c(1,2))
        probs.1<-seq(0, 1, by = 0.001); quants.1 <-quantile(as.vector(expres), probs.1)
        plot(probs.1, quants.1, main="Quantiles of Intensities")
        probs.2<-seq(0, 1, by = 0.01); quants.2 <-quantile(SD, probs.2)
        plot(probs.2, quants.2, main="Quantiles of Standard Deviations")
        par (opt)
      }
}


### createFilteringReportFile: crea el fitxer report del filtratge, es a dir, la informacio referent al proces
###                            de filtratgesobre el filtratge
###
### Parametres:
###
###    filteringReportFName : nom de l'arxiu report. Per defecte s'anomena "FilteringReport.txt"
###               outputDir : directori a on desara l'arxiu report
###
###

createFilteringReportFile <- function(filteringReportFName, outputDir)
{
  write(paste(rep("-", 70), collapse = ""), file = file.path(outputDir, filteringReportFName), append = FALSE)
  linia <- paste("Type Of Filtering", "Threshold", "Genes Number", "Number Of Samples", sep = "\t")
  write(linia, file = file.path(outputDir, filteringReportFName), append = TRUE)
  write(paste(rep("-", 70), collapse = ""), file = file.path(outputDir, filteringReportFName), append = TRUE)
}


### filteringReport: Mira si existeix el fitxer que conte el report del filtratge i desa la informacio.
###                  Si el fitxer no existei crida a la funcio 'createFilteringReportFile' per a crear-lo
###
### Parametres:
###
###                      n.exprs : dimensions de l'objecte que conte les expressions abans de filtrar
###        n.exprs.filtered.Affy : dimensions de l'objecte que conte les expressions despres de filtrar els controls d'Affymetrix
###    n.exprs.filtered.BySignal : dimensions de l'objecte que conte les expressions despres de filtrar per senyal
###                    signalThr : threshold que s'ha utilitzat pel filtratge per senyal
###       n.exprs.filtered.ByVar : dimensions de l'objecte que conte les expressions despres de filtrar per variabilitat
###                       varThr : threshold que s'ha utilitzat pel filtratge per variabilitat
###                    outputDir : directori a on desara l'arxiu report
###         filteringReportFName : nom de l'arxiu report. Per defecte s'anomena "FilteringReport.txt"
###
###

filteringReport <- function(n.exprs = c(NA ,NA),
                            n.exprs.filtered.Affy = c(NA ,NA),
                            n.exprs.filtered.BySignal,
                            signalThr,
                            n.exprs.filtered.ByVar = c(NA ,NA),
                            varThr,
                            outputDir,
                            filteringReportFName = "FilteringReport.txt")
{
  if(!file.exists(file.path(outputDir, filteringReportFName))) createFilteringReportFile(filteringReportFName, outputDir)
  
  linia <- paste("Non-filtering", "---", n.exprs[1], n.exprs[2], sep = "\t")
  write(linia, file = file.path(outputDir, filteringReportFName), append = TRUE)
  linia <- paste("By Affymetrix Controls", "---", n.exprs.filtered.Affy[1], n.exprs.filtered.Affy[2], sep = "\t")
  write(linia, file = file.path(outputDir, filteringReportFName), append = TRUE)
  linia <- paste("By Signal", signalThr, n.exprs.filtered.BySignal[1], n.exprs.filtered.BySignal[2], sep = "\t")
  write(linia, file = file.path(outputDir, filteringReportFName), append = TRUE)
  linia <- paste("By Variability", varThr, n.exprs.filtered.ByVar[1], n.exprs.filtered.ByVar[2], sep = "\t")
  write(linia, file = file.path(outputDir, filteringReportFName), append = TRUE)
  write(paste(rep("-", 35), collapse = " "), file = file.path(outputDir, filteringReportFName), append = TRUE)    

}


filterData <- function (expres,
                        controls = NULL,
                        removeNAs = FALSE,
                        entrezs = NULL,
                        bySignal = TRUE,
                        signalThr,
                        grups = NULL,
                        sigFun.Name = "filt.by.signal",
                        sigThr.as.perc = TRUE,
                        byVar = TRUE,
                        variabilityThr,
                        varFun.Name = "sdf",
                        varThr.as.perc = TRUE,
                        pairingFun.Name = NULL,
                        targets,
                        doReport = NULL,          # Per  defecte executara la crida a la funcio que fa el report
                        outputDir = NULL,                # Directori a on desa el report del filtrat
                        filteringReportFName = NULL)     # Nom per defecte del report
                        
{
  ### Per si les mosques es un data.frame
  if (is.data.frame(expres))
    expres<-as.matrix(expres)


  ### Eliminacio dels controls i els NAs
  if(!is.null(controls)){ # si hi ha controls...
    expresNotControls <- setdiff(rownames(expres), controls) # fem la llista dels que no son controls respecte expres
    if (removeNAs &&(!(is.null(entrezs)))){ # si tambe volem eliminar els NAs i tenim entrezs...
      entrezsNotNAs <- names(entrezs[!is.na(entrezs)]) # fem la llista dels que no son NAs respecte entrezs
      # i ens quedem aquells elements de rownames(expres) que no estan a la llista de controls ni a la de NAs
      notControlsNorNAs <- intersect(expresNotControls, entrezsNotNAs)
      exprs.filtered.0 <- expres[notControlsNorNAs,] # filtrem expres deixant els que no son controls ni NAs
    }else{ # si no volem filtrar NAs o no tenim entrezs de referencia...
      exprs.filtered.0 <- expres[notControls,] # filtrem exprs deixant els que no son controls
    }
  }else { # si no hi ha controls...
    if (removeNAs &&(!(is.null(entrezs)))){ # pero volem eliminar els NAs i tenim entrezs...
      entrezsNotNAs <- names(entrezs[!is.na(entrezs)]) # fem igualment la llista dels que no son NAs respecte entrezs
      exprs.filtered.0 <- expres[entrezsNotNAs,] # filtrem expres deixant els que no son NAs
    }else{ # si tampoc volem filtrar NAs o no tenim entrezs de referencia...
      exprs.filtered.0 <- expres # tambe es crea exprs.filtered.0 pero amb tots els valors inicials
    }
  }

  
  ### Filtratge per senyal
  if(bySignal && (!is.null(grups)))
  {
    exprs.filtered.1 <- filterBySignal(exprs.filtered.0, grups, threshold= signalThr, sigFun.Name=sigFun.Name, thr.as.perc=sigThr.as.perc )
  }else{
    exprs.filtered.1 <- exprs.filtered.0
  }

  ### Aparellament. MOLTES vegades aquest pas no es fa. 
  ### EL control de si es fa o no el basem en la funcio d'aparellament que es passa com parametre
  ### Si es deixa a NULL el pas s'omet
  if (!is.null(pairingFun.Name)){
     exprs.filtered.2 <- do.call(pairingFun.Name, list(exprs.filtered.1, targets))
  }else{
    exprs.filtered.2 <- exprs.filtered.1
  }

  ### Filtratge per variabilitat
  if (byVar)
  {
    exprs.filtered.3 <- filterByVar (exprs.filtered.2, threshold= variabilityThr, varFun.Name=varFun.Name, thr.as.perc=varThr.as.perc)
  }else{
    exprs.filtered.3 <- exprs.filtered.2
  }

  ### Report del proces de filtratge
  if (!is.null(doReport))
  {
    if (doReport)
    {
      if ((!is.null(filteringReportFName)) && (!is.null(outputDir)))
      {
        filteringReport(n.exprs = dim(expres),
                        n.exprs.filtered.Affy = dim(exprs.filtered.0),
                        n.exprs.filtered.BySignal = dim(exprs.filtered.1), signalThr = signalThr, 
                        n.exprs.filtered.ByVar = dim(exprs.filtered.3), varThr = variabilityThr,
                        outputDir = outputDir,
                        filteringReportFName = filteringReportFName)
      }
    }
  }

  return(exprs.filtered.3)
}



saveData <- function (expres, expresNames=colnames(expres),
        expres.csv.FileName, csvType, description="Normalized Values",    # --> Fitxers csv
        anotPackage, SYMBOL="SYMBOL", symbolsVector = NULL,         # --> Per anotar els arxius csv 
        expres.bin.FileName,                                        # --> Arxius binaris per les dades
        linksFile, outputDir )
{
    if (!(is.null(expres.csv.FileName))){
      expres.all <- expresAndgeneSymbols(expres=expres, expresNames=expresNames, 
            anotPackage=anotPackage, SYMBOL=SYMBOL, symbolsVector=symbolsVector)
      write2csv(expres.all, fileName = expres.csv.FileName, csv = csvType, outputDir = outputDir)
      addToLinksFile(linksFile, paste(expres.csv.FileName, substr(csvType, 1, 3), sep="."), categ = 'DATA', desc= description)
    }

    if (!(is.null(expres.bin.FileName))){
#      assign(expresName, exprs.filtered)
#      save(eval(parse(text=expresName)), file=file.path(resultsDir, expresFileName))

       save(expres, file=file.path(outputDir, expres.bin.FileName))
    }
}

################ FUNCIO OBSOLETA #################
##filterDataAndSave <- function (expres,
##        bySignal = TRUE, signalThr, grups = NULL, sigFun.Name="filt.by.signal", sigThr.as.perc=TRUE,
##        byVar = TRUE, variabilityThr, varFun.Name = "sdf", varThr.as.perc=TRUE,
##        pairingFun.Name=NULL, targets,
##        is.affy=TRUE, filterString="AFFX",
##        normalizedAll.FileName, normalizedFiltered.FileName, csvType,  # --> Fitxers csv
##        anotPackage, symbolsVector = NULL, SYMBOL="SYMBOL",             # --> Per anotar els arxius csv 
##        expres.FileName, expresFiltered.FileName,                     # --> Arxius binaris per les dades
##        linksFile, outputDir )
##                              
##{   #    Aquesta versio serveix nomes amb aquests targets concrets --> Fer-la mes general (pairingFactor, etc...)
##
##    
##    ### 1. FILTRATGE
##
##    exprs.filtered <- filterData (expres = expres,
##        bySignal = bySignal, signalThr = signalThr , grups = grups, sigFun.Name=sigFun.Name, sigThr.as.perc=sigThr.as.perc,
##        byVar = byVar, variabilityThr= variabilityThr, varFun.Name = varFun.Name, varThr.as.perc = varThr.as.perc,
##        pairingFun.Name  = pairingFun.Name, targets  = targets,
##        is.affy= is.affy, filterString =filterString)
##
##    ### 2. Creacio dels .csv/.xls
##
##    # 4.1 Valors d'expression normalitzades sense filtrar i sense aparellar
##
##    saveData (expres=expres , expresNames=colnames(expres),
##        expres.csv.FileName =  normalizedFiltered.FileName, csvType=csvType, 
##          description = "Normalized Values for filtered genes" ,
##        anotPackage  = anotPackage , SYMBOL="SYMBOL", symbolsVector = symbolsVector,
##        expres.bin.FileName  = expres.FileName,                                     
##        linksFile = linksFile, outputDir = outputDir)
##    
##   saveData (expres=exprs.filtered , expresNames=colnames(exprs.filtered),
##        expres.csv.FileName =  normalizedAll.FileName, csvType=csvType, 
##          description = "Normalized Values for all genes" ,
##        anotPackage  = anotPackage , SYMBOL="SYMBOL", symbolsVector = symbolsVector,
##        expres.bin.FileName  = expresFiltered.FileName,                                     
##        linksFile = linksFile, outputDir = outputDir)
##
##    return(exprs.filtered)
##}
######################################


#################################################################
### Part 2: Funcions per trobar gens DExpressats
###     - lmAnalysis()
###     - Funcions auxiliars
#################################################################

### 2.1 Relacio de num de gens que es seleccionaran segons el cutoff


genesSelectable <- function (topTab, adj0, adj1, adj2, P1, P2)
{
  upBelowB <- sum(topTab$B > 0  & topTab$t > 0)
  downBelowB <- sum(topTab$B >0 & topTab$t < 0)

  upBelowAdj0 <- sum(topTab$adj.P.Val < adj0 & topTab$t > 0)
  downBelowAdj0 <- sum(topTab$adj.P.Val < adj0 & topTab$t < 0)
 
  upBelowAdj1 <- sum(topTab$adj.P.Val < adj1 & topTab$t > 0)
  downBelowAdj1 <- sum(topTab$adj.P.Val < adj1 & topTab$t < 0)
  
  upBelowAdj2 <- sum(topTab$adj.P.Val < adj2 & topTab$t > 0)
  downBelowAdj2 <- sum(topTab$adj.P.Val < adj2 & topTab$t < 0)
  
  upBelowP1 <- sum(topTab$P.Value < P1 & topTab$t > 0)
  downBelowP1 <- sum (topTab$P.Value < P1 & topTab$t < 0)
  
  upBelowP2 <- sum(topTab$P.Value < P2 & topTab$t > 0)
  downBelowP2 <- sum(topTab$P.Value < P2 & topTab$t < 0)
  
  return(c(upBelowB = upBelowB,downBelowB = downBelowB,
           upBelowAdj0 = upBelowAdj0, downBelowAdj0 = downBelowAdj0,
           upBelowAdj1 = upBelowAdj1, downBelowAdj1 = downBelowAdj1,
           upBelowAdj2 = upBelowAdj2, downBelowAdj2 = downBelowAdj2,
           upBelowP1 = upBelowP1, downBelowP1 = downBelowP1,
           upBelowP2 = upBelowP2, downBelowP2 = downBelowP2)) 
}

genesSelected <- function (topTab = NULL, adj0=0.01, adj1 = 0.05, adj2 = 0.25, P1 = 0.01, P2 = 0.05)
{
  if (is.null(topTab))
  {
    seleccio <- data.frame(Col1 = integer(length = 12), 
    row.names = c("upReg-B>0", "downReg-B>0",
                  paste("upReg-Adjusted-p-val", adj0, sep = " < "), 
                  paste("downReg-Adjusted-p-val", adj0, sep = " < "),      
                  paste("upReg-Adjusted-p-val", adj1, sep = " < "), 
                  paste("downReg-Adjusted-p-val", adj1, sep = " < "),
                  paste("upReg-Adjusted-p-val", adj2, sep = " < "),
                  paste("downReg-Adjusted-p-val", adj2, sep = " < "),
                  paste("upReg-P value", P1, sep = " < "),
                  paste("downReg-P value", P1, sep = " < "),
                  paste("upReg-P value", P2, sep = " < "),
                  paste("downReg-P value", P2, sep = " < ")))[, -1]
  }else{
    genesSelectable(topTab, adj0, adj1, adj2, P1, P2)
  }

  return(seleccio)
}


### 2.2 La funcio que calculara l'Expressions_And_TopTable


escriuTop_i_Express <- function(expres,
                                topTab,
                                grup1 = NULL,
                                grup2 = NULL,
                                nom1,
                                nom2, 
                                fName = NULL,
                                fit = NULL,
                                fitCoef = NULL,
                                anotPackage,
                                my.symbols = NULL,
                                my.entrezs = NULL, ### FFF: parametre afegit  2011-11-15
                                csvType,
                                outputDir)
{
###
# fit: L'objecte "fit" a partir del qual s'han generat la "topTab"
# fitCoef: El numero (o el nom entre cometes) del coeficient de la matriu de 
#     contrastos corresponent a la comparacio

#  stopifnot(require(marray)) --> No volem que en depengui

  if (!is.null(grup2)) {  # si no es null el grup2 vol dir que hi ha grup1 i grup2
    expresA <- expres[, c(grup1, grup2)]
    
    means2groups <- function(x) 
    {
      mean1 <- mean(x[1:length(grup1)]) 
      mean2 <- mean(x[(length(grup1) + 1):(length(grup1) + length(grup2))])
      foldC <- mean1 - mean2
      
      return(c(foldC, mean1, mean2))
    }
    
    meansA <- t(apply(expresA, 1, means2groups))
    colnames(meansA) <- c("logFC validation", nom1, nom2)

  }else if (!is.null(grup1)) {   # en aquest cas sols hi ha grup1
    expresA <- expres[, grup1]

    means1groups <- function(x) {
      foldC <- mean(x[1:length(grup1)]) 
      # foldC <- mean1
      return(c(foldC, foldC))
      #return(foldC)
    }        

    meansA <- t(apply(expresA, 1, means1groups))
    colnames(meansA) <- c("logFC validation", "logFC validation")
    # colnames(meansA) <- c("logFC validation") 

  ### FFF 2011-11-02: afegida aquesta part en una extensio de l'if else... 
  ### per generar CSV en tots els casos de comparacions, siguin les que siguin
  }else { # en la resta de casos (grup2 is null i grup1 is null) hem de fer-ho sense grups
    meansA <- NULL
    expresA <- expres
  }

  if (!is.null(topTab$ID)){
	topA.order <- as.character(as.integer(topTab$ID))
  }else{
  	topA.order <- rownames(topTab)
  }

  if (!is.null(meansA)) { meansA.ord <- meansA[topA.order, ] }
  expresA.ord <- expresA[topA.order, ]
  fitA.ord <- fit[topA.order, ]
  
  if (! is.null(anotPackage))
  {
    #SymbolsA <- substr(getSYMBOL(topTab$ID, old2db (anotPackage)), 1, 10)
    SymbolsA <- getSYMBOL(topA.order, old2db (anotPackage))
  }else{
    if (!is.null(my.symbols))
    {
      gNames <- topA.order
      #SymbolsA <- substr(my.symbols[gNames], 1, 10)
      SymbolsA <- my.symbols[gNames]
      EntrezsA <- my.entrezs[gNames] ### FFF: linia afegida 2011-11-15
    }
  }

  # calculem els parametres del limma
  if(!is.null(fit) && !is.null(fitCoef))
  {
    s.post <- sqrt(fitA.ord$s2.post) 
    s.unscaled <- fitA.ord$stdev.unscaled[, fitCoef]
    df.res <- fitA.ord$df.residual
 
    s.res <- fitA.ord$sigma
    s.prior <- rep(sqrt(fitA.ord$s2.prior), length(s.res))
    df.prior <- rep(fitA.ord$df.prior, length(s.res))
                          
    fit.coef <- fitA.ord$coef[, fitCoef]
    t.ord <-  fit.coef / (s.unscaled * s.res)
    p.ord <- 2 * pt (abs(t.ord), df = df.res, lower.tail = FALSE)

    limmaPars <- data.frame(s.post, s.unscaled, t.ord, p.ord, df.res, s.res, s.prior, df.prior)
   }

  # muntem la combinacio de columnes pel csv segons el cas (incloent symbols i entrezs si n'hi ha)
  if (!is.null(fit) && !is.null(fitCoef)) {  # mirem si tenim limmaPars
    if (! is.null(my.symbols)) {  # mirem si tenim symbols
      if (! is.null(my.entrezs)) {  # mirem si tenim entrezs
        if (! is.null(meansA)) {  # mirem si hi ha mitjanes (o sigui que els contrasts son simples, e.g: A - B )
          # cas 1: limmaPars+, symbols+, entrezs+, means+
	  combinedA <- cbind(SymbolsA, EntrezsA, topTab, meansA.ord, expresA.ord, limmaPars)
        }else {  # si no hi ha mitjanes (o sigui que els contrasts son simples, e.g: (A-C) - (B-C) )
          # cas 2: limmaPars+, symbols+, entrezs+, means-
          combinedA <- cbind(SymbolsA, EntrezsA, topTab, expresA.ord, limmaPars)
        }
      }else { # si no tenim entrezs
        if (! is.null(meansA)) {  # si hi ha mitjanes (o sigui que els contrasts son simples, e.g: A - B )
          # cas 3: limmaPars+, symbols+, entrezs-, means+
	  combinedA <- cbind(SymbolsA, topTab, meansA.ord, expresA.ord, limmaPars)
        }else {  # si no hi ha mitjanes (o sigui que els contrasts son simples, e.g: (A-C) - (B-C) )
          # cas 4: limmaPars+, symbols+, entrezs-, means-
          combinedA <- cbind(SymbolsA, topTab, expresA.ord, limmaPars)
        }
      }
    }else {  # si no tenim symbols
      if (! is.null(my.entrezs)) {
        if (!is.null(meansA)) {
          # cas 5: limmaPars+, symbols-, entrezs+, means+
          combinedA <- cbind(EntrezsA, topTab, meansA.ord, expresA.ord, limmaPars)
        }else {
          # cas 6: limmaPars+, symbols-, entrezs+, means-
          combinedA <- cbind(EntrezsA, topTab, expresA.ord, limmaPars)
        }
      }else {
        if (! is.null(meansA)) {  # si hi ha mitjanes (o sigui que els contrasts son simples, e.g: A - B )
          # cas 7: limmaPars+, symbols-, entrezs-, means+
	  combinedA <- cbind(topTab, meansA.ord, expresA.ord, limmaPars)
        }else {  # si no hi ha mitjanes (o sigui que els contrasts son simples, e.g: (A-C) - (B-C) )
          # cas 8: limmaPars+, symbols-, entrezs-, means-
          combinedA <- cbind(topTab, expresA.ord, limmaPars)
        }
      }
    }
  }else {  # si no tenim limmaPars
    if (!is.null(my.symbols)) {
      if (! is.null(my.entrezs)) {
        if (! is.null(meansA)) {
          # cas 9: limmaPars-, symbols+, entrezs+, means+
          combinedA <- cbind(SymbolsA, EntrezsA, topTab, meansA.ord, expresA.ord)
        }else {
          # cas 10: limmaPars-, symbols+, entrezs+, means-
          combinedA <- cbind(SymbolsA, EntrezsA, topTab, expresA.ord)
        }
      }else {
        if (! is.null(meansA)) {
          # cas 11: limmaPars-, symbols-, entrezs+, means+
          combinedA <- cbind(SymbolsA, topTab, meansA.ord, expresA.ord)
        }else {
          # cas 12: limmaPars-, symbols-, entrezs+, means-
          combinedA <- cbind(SymbolsA, topTab, expresA.ord)
        }
      }
    }else {
      if (! is.null(my.entrezs)) {
        if (!is.null(meansA)) {
          # cas 13: limmaPars-, symbols-, entrezs+, means+
          combinedA <- cbind(EntrezsA, topTab, meansA.ord, expresA.ord)
        }else {
          # cas 14: limmaPars-, symbols-, entrezs+, means-
          combinedA <- cbind(EntrezsA, topTab, expresA.ord)
        }
      }else {
        if (!is.null(meansA)) {
          # cas 15: limmaPars-, symbols-, entrezs-, means+
          combinedA <- cbind(topTab, meansA.ord, expresA.ord)
        }else {
          # cas 16: limmaPars-, symbols-, entrezs-, means-
          combinedA <- cbind(topTab, expresA.ord)
        }
      }
    }
  }


  #if (!is.null(fName)) write.csv2(combinedA, file = fName)
  if (!is.null(fName)) write2csv(combinedA, fileName = fName, csv = csvType, outputDir = outputDir)
  
  return(combinedA)
}



### 2.3 Funcio per crear la topTable anotada
### Funcions per a extreure identificadors sinonims d'una llista 

extractSinonims <- function(my.strings)
{
  my.sinonims <- list()
       if (runMulticore ==1 || runMulticore ==3) { 
            my.sinonims <- mclapply(my.strings, function(x) unlist(strsplit(x, " /// ")))

         } else {
            my.sinonims <- lapply(my.strings, function(x) unlist(strsplit(x, " /// ")))
       }

  return(my.sinonims)
}

midSinonims <- function(my.IDs)
{
  my.indexes <- grep(" /// ", my.IDs)
  my.sinonimIDs <- my.IDs[my.indexes]
    
  my.IDs[my.indexes] <- extractSinonims(my.sinonimIDs)
  return(my.IDs)
}


### knowledge2Names: Per a cada identificador dels gens contrueix un link cap al fitxer html que conte les anotacions
###
### Parametres
###
###    my.genes: vector de caracters amb els identificadors dels gens
###    knowledgeFile: nom del fitxer html que conte les anotacions dels gens
###
### Exemple
###
###    > myProbes <- c("101054_at", "101878_at", "102209_at")
###    > anotFileName <- "Annotations.html"
###    > ProbesLinked <- knowledge2gNames(myProbes, anotFileName)
###    > print(ProbesLinked)
###    [1] "<A HREF=\"Annotations.html#101054_at\" TARGET=\"_blank\">101054_at</A>"
###    [2] "<A HREF=\"Annotations.html#101878_at\" TARGET=\"_blank\">101878_at</A>"
###    [3] "<A HREF=\"Annotations.html#102209_at\" TARGET=\"_blank\">102209_at</A>"

knowledge2gNames <- function(my.genes, knowledgeFile)
{
  my.genes <- as.character(my.genes)

  my.annotation <- character()
  
  if (!is.null(knowledgeFile))
  {
    for (i in 1:length(my.genes))
    {
      my.annotation[i] <- paste("<A HREF=\"", knowledgeFile,"#", my.genes[i] ,"\" TARGET=\"_blank\">", my.genes[i],"</A>", sep = "")
    }
  }
  return(my.annotation)
}

 annotateTopTable2.Old <- function (topTab, fName, Title = "Genes selected", 
                               anotPackage, EntrezIDs = NULL, SymbolIDs = NULL, anotFileName = "Annotated.Genes.html") 
 {
    
  if (!is.null(topTab$ID)){
	gNames <- as.character(as.integer(topTab$ID))
  }else{
  	gNames <- as.character(rownames(topTab))
  }


    if (is.null(EntrezIDs)){ 
        stopifnot(require(old2db (anotPackage), character.only=T))
        myenvirENTREZID <- eval(parse(text = paste(anotPackage, "ENTREZID", sep = "")))
        EntrezIDs <- unlist(mget(gNames, env = myenvirENTREZID, ifnotfound = NA))
    }else{
        EntrezIDs <- EntrezIDs[gNames]
    }
 
    if (is.null(SymbolIDs)){ 
        stopifnot(require(old2db (anotPackage), character.only=T))
        myenvirSYMBOL <- eval(parse(text = paste(anotPackage, "SYMBOL",  sep = "")))
        SymbolIDs <- unlist(mget(gNames, env = myenvirSYMBOL, ifnotfound = NA))
    }else{
        SymbolIDs <- SymbolIDs[gNames]
    }
    linkedList <- list(en = EntrezIDs)
    topTab <- topTab[, -1]
    otherNames = cbind(affyIDs = gNames, GeneSymbols = SymbolIDs, topTab)
    htmlpage(linkedList, filename = fName, title = Title, othernames = otherNames, 
        table.head = c("EntrezID", names(otherNames)), table.center = TRUE, 
        repository = list("en"), digits=4)
}


annotateTopTable2 <- function (topTab,
                               fName,
                               Title = "Genes selected", 
                               anotPackage,
                               EntrezIDs = NULL,
                               SymbolIDs = NULL,
                               #comparison = "",
                               anotFilename = "annotations"
                               ) 
{
  #gNames <- as.character(topTab$ID)
  if (!is.null(topTab$ID)){
	gNames <- as.character(as.integer(topTab$ID))
  }else{
   	gNames <- as.character(rownames(topTab))
  }
  
  linkedGeneNanes <- knowledge2gNames(gNames, knowledgeFile =  paste(anotFilename, "html", sep="."))
  
  if (is.null(EntrezIDs))
  { 
    myenvirENTREZID <- eval(parse(text = paste(anotPackage, "ENTREZID", sep = "")))
    EntrezIDs <- unlist(mget(gNames, env = myenvirENTREZID, ifnotfound = NA))
  }else{
    EntrezIDs <- midSinonims(EntrezIDs[gNames])
  }
 
  if (is.null(SymbolIDs))
  { 
    myenvirSYMBOL <- eval(parse(text = paste(anotPackage, "SYMBOL", sep = "")))
    SymbolIDs <- unlist(mget(gNames, env = myenvirSYMBOL, ifnotfound = NA))
  }else{
    SymbolIDs <- midSinonims(SymbolIDs[gNames])
  }

  aux.SymbolIDs <- as.character(SymbolIDs)
  names(aux.SymbolIDs) <- names(SymbolIDs)
       if (runMulticore ==1 || runMulticore ==3) { 
           lS <- mclapply(aux.SymbolIDs, function(l){unlist(strsplit(l, "//"))})
         } else {
           lS <- lapply(aux.SymbolIDs, function(l){unlist(strsplit(l, "//"))})
       }
   
  linkedList <- list(en = EntrezIDs)
  topTab <- topTab[, -1]

                                        # Si existeixen sinonims "otherNames" ha de ser una llista i no un data.frame com estava
  otherNames = list(affyIDs = linkedGeneNanes, GeneSymbols = lS, topTab)  
  htmlpage(linkedList, filename = fName, title = Title, othernames = otherNames,
           table.head = c("EntrezID", "affyIDs", "GeneSymbols", names(topTab)), 
           table.center = TRUE, repository = list("en"), digits=4)
}


### Funcioneta per decidir si cal cridar l'Expressions_and_TopTable
### FFF 2011-11-02: Si la comparacio no es d'1 ni de 2 mostres tamb?? farem l'Expressions_and_TopTable

is.twoSampleContrast <-function (row.cont.matrix){
  # Rep com parametre la fila de la matriu de contrastos
  # Torna TRUE si esta formada nomes per un 1 i un -1, 
  # que es la situacio on te sentit fer un Expressions_And_Toptable
  (length(row.cont.matrix[row.cont.matrix!=0])==2) && (sum(row.cont.matrix)==0)
}

is.oneSampleContrast <-function (row.cont.matrix){
  # Rep com parametre la fila de la matriu de contrastos
  # Torna TRUE si esta formada nomes per un 1 i un -1, 
  # que es la situacio on te sentit fer un Expressions_And_Toptable
  (length(row.cont.matrix[row.cont.matrix!=0])==1) && (sum(row.cont.matrix)==1)
}

 

### 2.4 lmAnalysis(): Crida iterada al proces d'estimacio/eBayes per als
###                   contrastos seleccionats d'una matriu de contrastos donada

### La funcio lmAnalysis() fa un processat complert a partir de
### - exprs.filtered = Una matriu de dades, normalitzada i filtrada, 
### - design = Una matriu de disseny associada a la de dades
### - cont.matrix = Una matriu de contrastos
### - contrasts2Test = Un vector d'enters representant quines columnes de la matriu de contrastos (i.e. quins contrastos) volem analitzar.
### - anotPackage = Un paquet d'anotacions per la taula de resultats
### - Expressions_And_Top= Una variable per decidir si es vol fer l'Expressions and TopTable
### - showParams = Una variable per decidir si es treuen els parametres del lmFit en l'ExpressionsAndTopTable
### - Els elements per fer servir "duplicateCorrelation" (use.DupCorr, block i nDups)
### - comparison = Un nom de comparacio comu a tots els contrastos
### - outputDir= Un directori de sortida
### - ENTREZIDS = Un vector d'identificadors "ENTREZ" o be NULL si es vol que ho generi a partir del paquet d'anotacions
### - SYMBOLIDS= Un vector d'identificadors "SYMBOL" o be NULL si es vol que ho generi a partir del paquet d'anotacions
### - linksFile= Un nom d'arxiu on es posen els vincles per generar el 'Results_File'
### L'analisi es fa fent servir el limma() i retornant, de moment
### - Un objecte lmFit


lmAnalysis <- function(exprs.filtered,
                       design, cont.matrix,
                       contrasts2test = 1:ncol(cont.matrix), 
                       anotPackage,
                       Expressions_And_Top = TRUE,
                       showParams = FALSE,
                       use.dupCorr = FALSE,
                       block = NULL, nDups = 1,
                       comparison = "",
                       outputDir = ".", 
                       ENTREZIDs = NULL,
                       SYMBOLIDs = NULL,
                       linksFile,
                       fitFileName,
                       csvType,
                       rows2HTML = NULL,
                       anotFileName)
{
  categLabel ='ANALYSIS'

  ### 1. Ajust del model lineal

  if (use.dupCorr && (!is.null(block)))
  {
    stopifnot(require(statmod))

    corfit <- duplicateCorrelation(exprs.filtered, ndups = nDups, block = blocs, design = design)
    fit < -lmFit(exprs.filtered, design = design, block = blocs, cor = corfit$consensus)
  }else{
    fit <- lmFit(exprs.filtered, design)
  }

  fit.main <- contrasts.fit(fit, cont.matrix)
  fit.main <- eBayes(fit.main)               # Podriem afegir el , proportion = priori.DE.prop

  ### Taula resum dels resultats : Pendent: numGenesChanged

  numGenesChanged <- genesSelected(NULL)
  ###
  ### Creacio de les taules HTML amb els resultats
  ### Preparem les anotacions
  ###
  #########################################################################
  ### Actualitzacio 2011
  ### Aixo s'hauria de treure.
  ### De moment s'inhabilita evitant que ENTREZIDs o SYMBOLIDs siguin nuls
  ##########################################################################
  if ( (is.null(ENTREZIDs)) & (!is.null(anotPackage))){
    stopifnot(require(old2db (anotPackage), character.only = T))
    envirName <- paste(anotPackage, "ENTREZID", sep = "")
    myenvirENTREZID <- eval(parse(text = envirName))
    ENTREZIDs <- unlist(AnnotationDbi::mget(rownames(exprs.filtered), myenvirENTREZID, ifnotfound = NA))
  }

  if ( (is.null(SYMBOLIDs))  & (!is.null(anotPackage)))
  {
    stopifnot(require(old2db (anotPackage), character.only = T))
    myenvirSYMBOL <- eval(parse(text = paste(anotPackage, "SYMBOL", sep = "")))
    SYMBOLIDs <- unlist(AnnotationDbi::mget(rownames(exprs.filtered), env = myenvirSYMBOL, ifnotfound = NA))
  }

  ### Farem un bucle que recorri la matriu de contrastos
  ### No farem tots els contrastos sino els que indiquem

  for (i in contrasts2test){
    # Obte la topTable per al contrast i-essim
    top.Diff <- topTable (fit.main, coef=i, n=nrow(fit.main$t), adjust="fdr") 
    fitCoefName<- colnames(cont.matrix)[i]
    contrastTitle <- paste(comparison,colnames(cont.matrix)[i],sep=".")
    # Escriu la topTable a una taula html
    atitle<-paste("Genes analyzed in comparison",contrastTitle,sep=" ")
    aFName0<-paste(contrastTitle,"html",sep=".")
    aFName<-paste(file.path(outputDir,contrastTitle), "html",sep=".")
      
    if((!is.null(anotPackage)) |((!is.null(SYMBOLIDs)) & (!is.null(ENTREZIDs))) )
    {
      if(is.null(rows2HTML))
      {
        top.Diff2HTML <- top.Diff
      }else{
        len.rows2HTML <- length(rows2HTML)
        if (len.rows2HTML==1)
        {
          cat(paste("Only", rows2HTML, "rows selected in html gene lists... \n", sep = " "))
          top.Diff2HTML <- top.Diff[1:rows2HTML, ]
        }else{
          cat(paste("Only",length(rows2HTML), "rows selected in html gene lists... \n", sep = " "))
          top.Diff2HTML <- top.Diff[rows2HTML, ]
        }
      }

      outNUL <- annotateTopTable2(top.Diff2HTML,
                                  aFName,
                                  atitle, 
                                  anotPackage=anotPackage,
                                  EntrezIDs=ENTREZIDs,
                                  SymbolIDs=SYMBOLIDs,
                                  anotFilename=anotFileName
                                  ) 
    }
      
    addToLinksFile (linksFile, aFName0, categ=categLabel, 
                    desc = paste("Test parameters and ranked list of genes for the comparison", contrastTitle, sep=" "))
      
    # Recull els gens canviats segons diversos criteris
      
    numGenesChanged <- cbind(numGenesChanged, genesSelectable(top.Diff, 0.01, 0.05, 0.25, 0.01, 0.05))
      
    # Fa el Volcano Plot i l'escriu a un arxiu .pdf
      
    vFName0<-paste("volcano", contrastTitle, "pdf",sep=".")
    if (toTIFF == TRUE){
    	vFName<-file.path(outputDir,paste("volcano", contrastTitle, "tiff",sep="."))
	}else{
	vFName<-file.path(outputDir,paste("volcano", contrastTitle, "pdf",sep="."))
	}
    if (toTIFF == TRUE){
    	tiff(file = vFName, width = 3200, height = 3200, units = "px", res = 800)
	}else{
	pdf(vFName)
	}
    opt <- par(cex.lab = 0.7)
    #if(!is.null(anotPackage)){
    if(!is.null(SYMBOLIDs)) {
      #volcanoNames<-SYMBOLIDs 
      volcanoNames<-SYMBOLIDs[rownames(exprs.filtered)]
    }else {
      volcanoNames<-rownames(exprs.filtered)
    }

    volcanoplot(fit.main, coef=i, highlight=10, names=volcanoNames, 
                main=paste("Differentially expressed genes",contrastTitle, sep="\n"))
    abline(v=c(-1,1))
    par(opt)
    dev.off()
    addToLinksFile (linksFile, vFName0, categ = "ANALYSIS", 
                    desc = paste("Volcano Plot for the comparison", contrastTitle, sep=" "))
      
    # Calcula l'Expression_And_TopTable
    if (Expressions_And_Top) # && (is.twoSampleContrast(cont.matrix[,i]) || is.oneSampleContrast(cont.matrix[,i]))
    {
      # La versio actual de Expressions_And_TopTable necessita els grups 
      # passats com a columnes de la matriu de dades, o el que HA DE SER EL MATEIX
      # passats com files de la matriu de disseny
      if (is.twoSampleContrast(cont.matrix[,i])) { # Si la comparacio inclou 2 mostres s'ha de definir grup1 i grup2
        nom1 <- rownames(cont.matrix)[cont.matrix[,i]==1]
        nom2 <- rownames(cont.matrix)[cont.matrix[,i]==-1]
        grup1 <- which(design[, colnames(design) == nom1] == 1)
        grup2 <- which(design[, colnames(design) == nom2] == 1)
      } else if (is.oneSampleContrast(cont.matrix[,i])) { # si inclou sols 1 mostra s'ha de definir sols grup1
        nom1 <- rownames(cont.matrix)[cont.matrix[,i]==1]
        nom2 <- NULL
        grup1 <- which(design[, colnames(design) == nom1] == 1)
        grup2 <- NULL
      } else { # si la comparacio es mes complexa s'han de posar ambdos grups i els seus noms a NULL
        grup1 <- NULL
        grup2 <- NULL
        nom1 <- NULL
        nom2 <- NULL
      }

      # Una futura versio que necessites directament les files de la matriu de disseny
      # Faria servir les mostres extretes amb aquesta instruccio
      #   mostra1<-rownames(design)[design[,colnames(design)==nom1]==1]
      #   mostra2<-rownames(design)[design[,colnames(design)==nom2]==1]
#     topFileName <-paste("ExpressAndTop",contrastTitle, "csv", sep=".")

      topFileName <- paste("ExpressAndTop", contrastTitle, sep=".")
      csvType <- ifelse(is.null(csvType), "csv2", csvType)
      if (showParams) {
         combined <- escriuTop_i_Express(exprs.filtered, top.Diff, 
                                         grup1 = grup1, grup2 = grup2, nom1 = nom1, nom2 = nom2, 
                                         #fName=file.path(outputDir, topFileName),
                                         fName = topFileName,
                                         fit = fit.main, 
                                         fitCoef = fitCoefName,
                                         anotPackage = anotPackage, 
                                         my.symbols = SYMBOLIDs,
                                         my.entrezs = ENTREZIDs, ### FFF: parametre afegit  2011-11-15
                                         csvType = csvType,
                                         outputDir = outputDir)
      }else {
         combined <- escriuTop_i_Express(exprs.filtered, top.Diff, 
                                         grup1 = grup1, grup2 = grup2, nom1 = nom1, nom2 = nom2, 
                                         #fName=file.path(outputDir, topFileName),
                                         fName = topFileName,
                                         fit = NULL,
                                         fitCoef = NULL,
                                         anotPackage = anotPackage, 
                                         my.symbols = SYMBOLIDs,
                                         my.entrezs = ENTREZIDs, ### FFF: parametre afegit  2011-11-15
                                         csvType = csvType,
                                         outputDir = outputDir)
      }
      addToLinksFile (linksFile, paste(topFileName, substr(csvType, 1, 3), sep="."), categ=categLabel, 
                      desc = paste("TopTable and normalized expression values for the comparison", fitCoefName, sep=" "))
    }
  }

  # Grava la taula resum de totes les comparacions
  colnames(numGenesChanged) <- colnames(cont.matrix)[contrasts2test]
  chgdFName0 <- paste("numGenesChanged", comparison, substr(csvType, 1, 3), sep = ".")
  chgdFName  <- paste("numGenesChanged", comparison,  sep = ".")
  write2csv(numGenesChanged, fileName = chgdFName, csv = csvType, outputDir = outputDir)  

  addToLinksFile (linksFile, chgdFName0, categ = categLabel, desc = paste("Number of selectable genes in comparisons", comparison, sep = " "))

  # Grava el fit.main donant-li un nom propi de cada comparacio
  # Potser convindria que ho rebes com parametre?
  #
  # Aquestes darreres linies s'han d'acabar de polir per aconseguir que 
  # l'objecte es gravi amb el nom que li vull assignar i no com "fit.main"

  # fitMainName <- paste("fit", comparison, sep = ".")
  # Enlloc d'aixo he afegit un parametre "fitFileName"

  save(fit.main, file = file.path(outputDir, fitFileName))

  return(fit.main)
}


### 2.5 Funcio per executar el 'lmAnalysis' a partir d'un objecte parametres.

doLmAnalysis <- function (lmPar)
{
### lmPar es una llista amb la seguent composicio.
### En un futur no llunya hauria de ser una classe S4
#
#   dades = El nom de la matriu de dades
#   design = La matriu de disseny
#   contMat = La matriu de contrastos
#   whichContrasts = Quins contrastos s'inclouen en l'analisi
#   anotPackage = Quin paquet d'anotacions cal fer servir
#   outputDir = On van els resultats (
#   ExpressionsAndTop= (si esta a TRUE es calcula),
#   showLmParams = (si esta a TRUE s'afegeixen a l'Expressions and Top els parametres de l'ajust)
#   use.dupCorr =FALSE, block =NULL, nDups=1  ("Pos eso")
#   comparisonName = (El nom que es dona a la comparacio)
#   ENTREZIDs = NULL,
#   SYMBOLIDs = NULL (Si se li dona no els calcula, si son NULL els busca ell)
#   fileOfLinks (Arxiu on es desen els vincles per generar Results_File)

  p <- lmPar[[1]]

  #############################################################################
  ### Actualitzacio 2011
  ###
  ### Els ENTREZIDs i SYMBOLIDs es pasen com a text i s'avaluen internament
  ############################################################################
  
  if(!is.null(p$ENTREZIDs)){
    EntrezIDs <-  eval(parse(text = p$ENTREZIDs)) #  Seran l'entrezTable
  }

  if(!is.null(p$SYMBOLIDs)){
    SymbolIDs <-  eval(parse(text = p$SYMBOLIDs)) #  Seran l'entrezTable ### FFF: ??o la symbolsTable?
  }

  ### Aquestes variables es pasen ara a la crida a lmAnalysis
  ###############################################################################
  
  if (!is.null(p$expresFileName)){  
      # Enlloc d'aquestes dues linies que depenen de que la matriu s'hagi desat amb el nom "exprs.filtered"
      # load (file.path(p$outputDir, p$expresFileName))
      # expres <- exprs.filtered ### Aixo s'ha d'arreglar per tal que el nom es llegeixi de l'arxiu
      # Farem servir la funcio loadFromFile definida, de moment, a novesFuncions.Alex
      expres <- loadFromFile (file.path(p$outputDir, p$expresFileName))
  }else{
      if (!is.null(p$dades))
      {
        expres <- eval(parse(text = p$dades)) # Posar-hi un tryCatch per poder sortir si d error!!!
      }else{
        stop("Error, Cal definir o les dades o el nom de l'arxiu")
      }
  }
  
  if (is.null(p$whichContrasts))
  {
    contrasts2test <- 1:ncol(p$contMat)
  }else{
    contrasts2test <- p$whichContrasts
  }
       
  fitMain <- lmAnalysis(exprs.filtered = expres,
                        design = p$designMat, 
                        cont.matrix = p$contMat,
                        contrasts2test = contrasts2test,
                        anotPackage = p$anotPack,
                        outputDir = p$outputDir,
                        comparison = p$comparisonName,
                        Expressions_And_Top = p$ExpressionsAndTop,
                        showParams = p$showLmParams,
                        use.dupCorr = p$use.dupCorr,
                        block = p$block,
                        nDups = p$nDups ,
                        ENTREZIDs = EntrezIDs, # abans p$ENTREZIDs,
                        SYMBOLIDs = SymbolIDs, # abans p$SYMBOLIDs,
                        linksFile = p$fileOfLinks,
                        fitFileName = p$fitFileName,
                        csvType=p$csvType,
                        rows2HTML = p$rows2HTML,
                        anotFileName = p$anotFilename
                        )
                    
                    
  designMatrixName = paste("designMatrix",p$comparisonName, sep=".")
  contrastMatrixName = paste("contrastMatrix",p$comparisonName, sep=".")

  write2csv(p$designMat, fileName = designMatrixName, csv = p$csvType, outputDir = resultsDir)  
  write2csv(p$contMat, fileName = contrastMatrixName, csv = p$csvType, outputDir = resultsDir)

  csvType <- ifelse(is.null(p$csvType), "csv2", p$csvType)
  addToLinksFile(p$fileOfLinks, paste(designMatrixName, substr(csvType, 1, 3),  sep="."), categ = "ANALYSIS", 
                 desc = paste("Design Matrix for comparison ", p$comparisonName, sep=""))
  addToLinksFile(p$fileOfLinks, paste(contrastMatrixName, substr(csvType, 1, 3),  sep="."), categ = "ANALYSIS", 
                 desc = paste("Contrast Matrix for comparison ", p$comparisonName, sep=""))

  return (fitMain)             
}

##########################################################################33
### Funcio per crear la taula de gens anotats
##########################################################################33

### htmlids2probes: Patch que modifica un fitxer HTML desat amb la instruccio 'saveHTML' del paquet "annaffy"
###                 afegint un identificador html a cada fila que conte un gen
###
### Parametres:
###
###   gNames : vector character que conte la llista de gens per a ser anotats
###    fName : nom del fitxer desat amb la instruccio 'saveHTML'
###
### OBERVACIONS IMPORTANTS:
###
###   Aquesta funcio patch va associada amb la funcio 'annotateResults'

htmlids2probes <- function(gNames, fName)
{
  html2data.frame <- read.table(fName, sep="\t", quote="", stringsAsFactors = FALSE)
  
  my.class <- "<td class=\"aafProbe\">"
       if (runMulticore ==1 || runMulticore ==3) { 
            whichRows <- unlist(mclapply(html2data.frame, function(x){grep(my.class, x)}))
         } else {
            whichRows <- unlist(lapply(html2data.frame, function(x){grep(my.class, x)}))
       }
  new.my.class <- html2data.frame[whichRows,]

  htmlids <-   paste("<td id=\"", gNames,"\" class=\"aafProbe\">", sep="")
  for (i in 1:length(htmlids))
  {
    new.my.class[i] <- sub("<td class=\"aafProbe\">", htmlids[i], new.my.class[i])
  }

  html2data.frame[whichRows,] <-  new.my.class
  
       if (runMulticore ==1 || runMulticore ==3) { 
           whichRow.body <- unlist(mclapply(html2data.frame, function(x){grep("<body", x)}))
         } else {
           whichRow.body <- unlist(lapply(html2data.frame, function(x){grep("<body", x)}))         
       }
  html2data.frame[whichRow.body,] <- "<boby>"

  error <- write.table(html2data.frame, file = fName, sep = "\t", quote = FALSE, col.names = FALSE, row.name = FALSE)
  
  return(error)
}


### annotateResults: Crea un fitxer HTML, anomenat "Annotated.Genes.html", que conte anotacions
###                  de diferents recursos per a una llista de gens
###
### Parametres:
###
###             fit : 
###     anotPackage : Paquet d'anotacions
###  comparisonName :
###       outputDir : Nom del direcori de sortida a on es desara el fitxer amb els gens anotats 
###     fileOfLinks : Arxiu on es desen els vincles per generar Results_File
###
### OBERVACIONS IMPORTANTS:
###
###   'annotateResults' utilitza la funcio 'htmlids2probes', que es un patch, per a donar un identificador
###   html a cada fila de la taula. Aixo permet que tots els gens del fitxer "Annotated.Genes.XXX.html" siguin
###   enlla??ables des d'altres documents (p.e. des de  els documents html que contenen les llistes dels gens
###   diferencialment expressats

annotateResults <- function(gNames, anotPackage, comparisonName="", outputDir, fileOfLinks, anotationsFName=NULL)
{ 
 if (is.null(anotationsFName)){
    anotationsFName <- paste("Annotated.Genes", comparisonName, "html", sep=".")
  }else{
     anotationsFName <- paste(anotationsFName, comparisonName, "html", sep=".")
  }
 if(!is.null(anotPackage)){   
   stopifnot(require("KEGG.db"))
   stopifnot(require("GO.db"))
   stopifnot(require("annaffy"))
  
   atab <- aafTableAnn(gNames, old2db(anotPackage), aaf.handler())
 
   aFName <- file.path(outputDir, anotationsFName)

   saveHTML(atab, file = aFName)
   htmlids2probes(gNames, aFName)
  }else{
    warning("This function needs that there is a '.db' annotations package available")
    warning(paste("You must create the annotions sheet manually and name it: ",anotationsFName))
 #   aFName <- file.path(outputDir, anotationsFName)
 #   saveHTML(atab, file = aFName)
 #  htmlids2probes(gNames, aFName)    
  }
    
  addToLinksFile(fileOfLinks, anotationsFName, categ = "ANNOT", desc = "Annotations for all genes analyzed in the comparison")
}

### Funcio per executar el 'annotateResults' a partir d'un objecte 'mcParmeters'.

doAnnotateResults <- function(mcPar)
{

  p <- mcPar[[1]]

  if (!is.null(p$fitFileName)){
    fitMain <- loadFromFile(file.path(p$outputDir,p$fitFileName))
  }else{
      if (!is.null(p$fitMain)) {
	       fitMain <- eval(parse(text = p$fitMain)) # Posar-hi un tryCatch per poder sortir si dona error!!!
      }else{
        	stop("Error, Cal subministrar un nom d'arxiu o d'objecte 'fitMain'")
      }
  }
  my.genes <- fitMain[["genes"]][,1]
  annotateResults(gNames= my.genes, anotPackage = p$anotPackage, comparisonName=p$comparisonName, outputDir = p$outputDir, fileOfLinks = p$fileOfLinks)
}


##########################################################
### Part 3. Seleccio de gens despres del multiple testing
##########################################################
 
### 3.1 Selecciona les files de la matriu de decideTests on un gen surt significatiu 
###     en, com a miim una comparacio

resSelected <- function(res)
{
  sum.res.rows <- apply(abs(res), 1, sum)
  res.selected <- res[sum.res.rows != 0, ]
  if(!is.null(res.selected))
  {
    if (is.matrix(res.selected))
    {
      resSelected <- res.selected
    }else{
#      resSelected <- matrix(res.selected, nrow=1)
      if(ncol(res)>1)
        {
          resSelected <- matrix(res.selected, nrow=1)
        }else{
          resSelected <- matrix(res.selected, ncol=1)
          rownames(resSelected)<-names(res.selected)
        }

    }
  }else{
    resSelected <- NULL
  }

  return(resSelected)
}

### 3.2 Selecciona les files de la matriu de decideTests on un gen surt significatiu 
###     en, totes les comparacions

commonSelected <-function (res)
{
  sum.res.rows<-apply(abs(res),1,sum)
  common.selected<-res[sum.res.rows==ncol(res),]
  if(nrow(common.selected)>0){
    commonSelected<-common.selected
  }else{
    commonSelected <-NULL}
  return(commonSelected)
}



plotVennDiagram <- function(res.selected,
			    colsVenn,
                            whichContrasts,
                            vennColors = c("red", "yellow", "green", "blue", "pink"),
                            comparisonName,
                            titleText,
                            outputDir,
                            linksFile,
                            categLabel,
                            my.cex)
{
  if (is.null(colsVenn)) colsVenn <-1:length(whichContrasts)
  if (length(colsVenn)> 5) stop("venn.diagram cannot plot more than 5 elements at the same time")

  vennFName  <- paste("venn", comparisonName, sep = ".")
  ncomp <- as.character(length(colsVenn))
  
  switch(ncomp,
           "2" = {venn.plot <- venn.diagram(list(A = which(res.selected[, 1]!=0), B = which(res.selected[, 2]!=0)),
				category.names = colnames(res.selected)[1:2],
				fill = c(vennColors[1], vennColors[2]),
				alpha = 0.30,
				main = paste(comparisonName, titleText, sep=" "),
				resolution = 600,
				cat.cex = my.cex,
				filename = NULL)},
           "3" = {venn.plot <- venn.diagram(list(A = which(res.selected[, 1]!=0), B = which(res.selected[, 2]!=0),
                                C = which(res.selected[,3]!=0)),
				category.names = colnames(res.selected)[1:3],
				fill = c(vennColors[1], vennColors[2], vennColors[3]),
				alpha = 0.30,
				resolution = 600,
				cat.cex = my.cex,
				main = paste(comparisonName, titleText, sep=" "),
				filename = NULL)},
           "4" = {venn.plot <- venn.diagram(list(A = which(res.selected[, 1]!=0), B = which(res.selected[, 2]!=0),
				C = which(res.selected[,3]!=0), D = which(res.selected[,4]!=0)),
				category.names = colnames(res.selected)[1:4],
				fill = c(vennColors[1], vennColors[2], vennColors[3], vennColors[4]),
				alpha = 0.30,
				resolution = 600,
				cat.cex = my.cex,
				main = paste(comparisonName, titleText, sep=" "),
				filename = NULL)},
           "5" = {venn.plot <- venn.diagram(list(A = which(res.selected[, 1]!=0), B = which(res.selected[, 2]!=0),
                                C = which(res.selected[,3]!=0), D = which(res.selected[,4]!=0), W = which(res.selected[,5]!=0)),
				category.names = colnames(res.selected)[1:5],
				fill = c(vennColors[1], vennColors[2], vennColors[3], vennColors[4], vennColors[5]),
				alpha = 0.30,
				cat.cex = my.cex,
				resolution = 600,
				cat.pos = c(0,0,-140,-140,0),
				cat.default.pos = "outer",
				main = paste(comparisonName, titleText, sep=" "),
				filename = NULL)})

     pdf(file.path(outputDir,paste(vennFName, "pdf", sep = ".")))
       grid.draw(venn.plot)
     dev.off()
     addToLinksFile(linksFile, paste(vennFName, "pdf", sep = "."), categ = categLabel, desc = "Venn Diagram for genes selected from multiple comparisons")

     if(toTIFF){
          tiff(file.path(outputDir,paste(vennFName,"tiff", sep = ".")), width = 3200, height = 3200, units = "px", res = 800)
            grid.draw(venn.plot)
          dev.off()

     addToLinksFile(linksFile, paste(vennFName, "tiff", sep = "."), categ = categLabel, desc = "Venn Diagram for genes selected from multiple comparisons")
     }



}


### 3.3 Funcio principal: multipleComp

# La funcio 'multipleComp executa un 'decideTests' per als contrastos seleccionats d'un objecte 'fit'.
# Actualment retorna:
#
# 1) Una taula de decideTests amb els gens que canvien alhora entre totes les comparacions
# 2) Un diagrama de Venn en un arxiu pdf, amb el nom de la comparacio
#
# A mes a mes, cal fer:
# 
# - Que tambe retorni els gens comuns a totes les comparacions
# - Que si no pot fer el diagrama de Venn perque hi ha mes de tres contrastos
#   retorni una taula amb els comptatges de les interseccions

multipleComp <- function(fitMain,
                         whichContrasts, 
                         comparisonName,
                         titleText,
                         outputDir, 
                         anotPackage,
                         my.symbols=NULL,
                         linksFile,
                         multCompMethod = "separate",
                         adjustMethod = "none", 
                         selectionType = c("any", "all", "anyDown", "anyUp", "allDown", "allUp"),
                         P.Value.cutoff = 0.05,
                         plotVenn,
                         colsVenn = NULL,
			 vennColors,
                         cexVenn = 1,
                         geneListFName=paste("geneList", comparisonName,"Rda", sep="."),
                         csvType = NULL,
                         minLFC = 0)
{
# La informacio es treu d'un objecte 'fit' obtingut d'un analisi amb el limma.
# Si l'objecte 'fit' es nul, la matriu de contrasts que conte es nul??la 
# o nomes hi ha un contrast no es fara el decideTest i retornara NULL

  geneList <- NULL
  categLabel <- "MULTCOMP"
#  if ((!is.null(fitMain)) && (!is.null(fitMain$contrasts[,whichContrasts])) && (ncol(fitMain$contrasts[,whichContrasts])> 1))
  if ((!is.null(fitMain)) && (!is.null(fitMain$contrasts[,whichContrasts])) && (is.matrix(fitMain$contrasts[,whichContrasts])))  
  {
    res <- decideTests(fitMain[, whichContrasts], method = multCompMethod, adjust.method = adjustMethod, p.value = P.Value.cutoff, lfc = minLFC)
    res.selected <- resSelected(res)

    if (!is.null(res.selected))
    {

      if (!is.null(my.symbols))
      {
        gNames <- rownames(res.selected)

        my.indNA<-which(is.na(my.symbols[gNames]))
        my.symbols[gNames[my.indNA]]<-my.symbols[gNames[my.indNA]]

        symbols.selected <- my.symbols[gNames]
        res.selected2 <- cbind(SYMBOLS=symbols.selected, res.selected)
      }else{
        if (!is.null(anotPackage))
        {
          stopifnot(require(old2db(anotPackage), character.only = T))
          myenvirSYMBOL <- eval(parse(text = paste(anotPackage, "SYMBOL", sep = "")))
          symbols.selected <- unlist(mget(rownames(res.selected), env = myenvirSYMBOL, ifnotfound = NA))
          res.selected2 <- cbind(SYMBOLS=symbols.selected, res.selected)
        }else{
          res.selected2 <- res.selected
        }
      }
      ## if (!is.null(anotPackage)){
      ##   stopifnot(require(old2db(anotPackage), character.only = T))
      ##   myenvirSYMBOL <- eval(parse(text = paste(anotPackage, "SYMBOL", sep = "")))
      ##   symbols.selected <- unlist(mget(rownames(res.selected), env = myenvirSYMBOL, ifnotfound = NA))
      ##   res.selected2 <- cbind(SYMBOLS=symbols.selected, res.selected)
      ## }else{
      ##   res.selected2 <- res.selected
      ## }

#      selectedFName <- paste("multComp", comparisonName, "csv", sep = ".")
#      write.csv2(res.selected2, file = file.path(outputDir, selectedFName))
      selectedFName <- paste("multComp", comparisonName, sep = ".")
      
      csvType<-ifelse(is.null(csvType), "csv2", csvType)
      write2csv(res.selected2, fileName = selectedFName, csv = csvType, outputDir = outputDir)
      addToLinksFile(linksFile, paste(selectedFName, substr(csvType[1], 1, 3), sep="."), categ = categLabel, desc="Genes selected from multiple comparisons")

      if (plotVenn)
      {
        plotVennDiagram(res.selected, colsVenn, vennColors, whichContrasts = whichContrasts, comparisonName = comparisonName, titleText, outputDir = outputDir, linksFile, categLabel, my.cex = cexVenn)
      }
    }

    geneList <- rownames(res.selected)

    #   assign(genelistName, geneList, envir = .GlobalEnv)
    #   eval(parse(text=paste(genelistName, "geneList", sep="<-")))
    #   save(genelistName, file=file.path(outputDir, paste(genelistName, "Rda", sep=".")))

    save(geneList, file = file.path(outputDir, geneListFName))
  }else{
    if ((!is.null(fitMain)) && (!is.null(fitMain$contrasts[,whichContrasts])) && (!is.matrix(fitMain$contrasts[,whichContrasts])))
      {
         res <- decideTests(fitMain[, whichContrasts], method = multCompMethod, adjust.method = adjustMethod, p.value = P.Value.cutoff, lfc = minLFC)
         res.selected <- resSelected(res)
         geneList <- rownames(res.selected)
         save(geneList, file = file.path(outputDir, geneListFName))
      }
  }

  return(geneList)
}


### 2.5 Funcio per executar el 'multipleComp' a partir d'un objecte parametres.

doMultCompAnalysis <- function(mcPar)
{

  p <- mcPar[[1]]

  
  #############################################################################
  ### Actualitzacio 2011
  ###
  ### Els ENTREZIDs i SYMBOLIDs es pasen com a text i s'avaluen internament
  ############################################################################

  if(!is.null(p$my.symbols)){
    SymbolIDs <-  eval(parse(text = p$my.symbols)) #  Seran el symbolsTable
  }

  ### Aquestes variables es pasen ara a la crida a lmAnalysis
  ###############################################################################
  

  

  if (!is.null(p$fitFileName)){
    fitMain <- loadFromFile(file.path(p$outputDir,p$fitFileName))
  }else{
      if (!is.null(p$fitMain)) {
	       fitMain <- eval(parse(text = p$fitMain))
                        # Posar-hi un tryCatch per poder sortir si dona error!!!
      }else{
        	stop("Error, Cal subministrar un nom d'arxiu o d'objecte 'fitMain'")
      }
  }

  geneList <-  multipleComp(fitMain = fitMain,
                            whichContrasts = p$whichContrasts, 
                            comparisonName = p$comparisonName,
                            titleText = p$titleText, 
                            outputDir = p$outputDir, 
                            anotPackage = p$anotPackage,
                            my.symbols = SymbolIDs, # era p$my.symbols,
                            linksFile = p$fileOfLinks,
                            multCompMethod = p$multCompMethod, 
                            adjustMethod = p$adjustMethod, 
                            selectionType = p$selectionType, 
                            P.Value.cutoff = p$P.Value.cutoff, 
                            plotVenn = p$plotVenn,
                            colsVenn = p$colsVenn,
			    vennColors = p$vennColors,
                            cexVenn = p$cexVenn,
                            geneListFName=p$geneListFName,
                            csvType=p$csvType,
                            minLFC=p$minLogFC)

  return (geneList)             
}


###################################################################
### Part 4. Creacio d'un heatmap senzill amb els gens seleccionats
###################################################################

# 4.1. Funcio per seleccionar les mostres amb que farem el cluster
#       No cal fer cluster de tots els arrays sino nomes d'aquells que s'han analitzat
#       Si per exemple en una comparacio multiple nomes s'han fet servir algun contrastos
#       sembla raonable utilitzar nomes les mostres que han intervingut en aquells

selectSamples <- function(designMat, contMat, whichContrasts)
{
  selectedFactors <- which(apply(contMat[, whichContrasts], 1, sum) != 0)
  selectedSamples <- which(apply(designMat[, selectedFactors], 1, sum) != 0)

  return(selectedSamples)
}


### 4.2 Funcio per recuperar l'objecte "fit.main" de l'arxiu on es grava
###     Es una mica "parxe" perque el que voldria es haver-lo pogut gravar
###     amb un nom escollit per l'usuari

genelistFromFile <- function (directory, fileName)
{
  load(file = file.path(directory, fileName))

  return(geneList)
}


### 4.3 Funcio per fer el plot dels perfils:

plotClust <- function(mydata, ind.cut, i, scal = FALSE, ...)
{
  my.col <- ncol(mydata)
  matData <- matrix(mydata[ind.cut==i, ], ncol = my.col)

  if(scal)
  {
    matData <- (t(scale(t(matData), scale = F)) + mean(matData))
  }

  matplot(1:my.col, t(matData), type = "l", xlab = "", ylab = "log-Expression", las = 2, axes = FALSE, main = paste("Cluster ", i, sep = ""), ...)
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], border = "black")
  axis(2, las = 2, cex.axis = 0.7)
  axis(1, at = 1:my.col, lab = colnames(mydata), las = 2, cex.axis = 0.5)
}


# 4.4 Funcio 'clusterAnalysis'
#     De moment es el no va mes de la cutredat clusteristica
#     However, com te un argument "..." funcionare tan be com es vulgui que ho faci
#     a base de passar-li els parametres que li calen al heatmap

clusterAnalysis <- function(expres, 
                            genes, 
                            samples, 
                            sampleNames,
                            comparisonName,
                            anotPackage,
                            my.symbols=NULL,
                            outputDir,
                            fileOfLinks,
                            numClusters,
                            rowDistance,
                            colDistance,
                            RowVals = TRUE,
                            ColVals = TRUE,
                            colorsSet,
                            colsForGroups,
                            escala,
                            densityInfo = "none",
                            cexForColumns,
                            cexForRows,
                            plotProfiles = FALSE,
                            Title = "",
                            csvType = NULL)
{
  cat ("  - Comparison", Title, "\n") 
  categLabel <- "CLUSTER"
  exprs2Cluster <- expres[genes, samples]
  colnames(exprs2Cluster) <- sampleNames 

  ###---------------------------------------------------------------------------------
  ### Primera part: Representacio del Heatmap
  ###---------------------------------------------------------------------------------

                                        # Preparem els dendrogrames
  
  dendro <- 'both' # d'entrada esperem dibuixar tant dendro de files com de columnes
  
  if (RowVals && (!is.numeric(RowVals)))
  {
                                        # Caldria mirar si 'rowDistance' es nula, i si no ho es fer-la servir aqui
    clustRow <- hclust(as.dist(1 - cor(t(exprs2Cluster))),  method = "average")
    dendroRow <- as.dendrogram(clustRow)
  }else{
    dendroRow <- RowVals
  }

  if (ColVals && (!is.numeric(ColVals)))
  {    
                                        # Caldria mirar si 'rowDistance' es nula, i si no ho es fer-la servir aqui
    clustCol <- hclust(dist(t(exprs2Cluster)), method = "average")
    dendroCol <-  as.dendrogram(clustCol)
  }else{
    if(is.numeric(ColVals))
      {
        exprs2Cluster <- exprs2Cluster[, ColVals]
        colsForGroups <- colsForGroups[ColVals]
        ColVals <- FALSE
        dendroCol <- FALSE
      }else{
        dendroCol <- ColVals
      }
  }

  if (ColVals && RowVals)
  {
    dendro <- 'both'
  }else{
    if (ColVals && (!RowVals))
    {
      dendro <- 'column'
    }else{
      if ((!ColVals) && RowVals)
      {
        dendro <- 'row'
      }else{
        dendro <- 'none'
      }
    }
  }
  if (toTIFF == TRUE){
  	heatMapFName <- paste("HeatMap", comparisonName, "tiff", sep=".")
  }else{
	heatMapFName <- paste("HeatMap", comparisonName, "pdf", sep=".")
  }
  mainTitle <- ifelse (Title == "", comparisonName, Title )
  
  if(!is.null(outputDir))
  {
    if (toTIFF == TRUE){
    	tiff(file = file.path(outputDir, heatMapFName), width = 3200, height = 3200, units = "px", res = 800)
    	par(cex.main = 1)
    }else{
	pdf(file.path(outputDir, heatMapFName))
    }
  }

  ### FERRAN (2012/12/10)
  ### modificacions per mostrar els geneSymbols despres dels affyIds de cada fila del heat map
  ### SOLS FUNCIONA SI TENIM my.symbols !!!
  foundSymbols <- my.symbols[which(names(my.symbols) %in% rownames(exprs2Cluster))]
  newNames <- unlist(lapply(strsplit(paste(names(foundSymbols), foundSymbols ,sep="."), "//"), function(l) l[[1]][1]))
  rownames(exprs2Cluster)[rownames(exprs2Cluster) %in% names(my.symbols)] <- newNames
  #######################

  hm <- heatmap.2(exprs2Cluster, 
                  col = colorsSet,
                  ColSideColors = as.character(colsForGroups),
                  scale = escala,
                  Rowv = dendroRow,     
                  Colv = dendroCol,
                  dendrogram = dendro,
                  key = TRUE,
		  keysize=1.5,
                  symkey = FALSE,
                  density = "none",
                  trace = "none",                                       # Nomes a heamap.2
                  cexCol = cexForColumns,
		  #labRow = NA,						# Si no volen el nom dels gens
                  cexRow = cexForRows,
      cex.main =0.6,
                  main = mainTitle)
  
  if(!is.null(outputDir))
  {
    dev.off()
  }

  addToLinksFile(fileOfLinks, heatMapFName, categ = categLabel, desc = "Heatmap made from genes selected from multiple comparisons")
  
  ###---------------------------------------------------------------------------------
  ### Segona part: Si li hem passat un cert nombre de Clusters (!is.null(numClusters))
  ###              Fem dues coses: (1) El perfil mitja de cada cluster.
  ###                              (2) Escrivim un arxiu assignant els gens als clusters
  ###---------------------------------------------------------------------------------

  if  (RowVals)
  {
    cutN <- cutree(clustRow, numClusters)

    if (ColVals)
    {
      names.ord <- (clustCol$labels[clustCol$order])
    }else{
      names.ord <- 1:ncol(exprs2Cluster)
    }
    
    if(plotProfiles)
    {
      plotClustFName <- paste("ProfilePlots", comparisonName, "pdf", sep = ".")

      pdf(file.path(outputDir, plotClustFName))
         for(i in 1:numClusters)
         {
            plotClust(exprs2Cluster[, names.ord], cutN, i, scal = T)
          }
      dev.off()
      
      addToLinksFile(fileOfLinks, plotClustFName, categ = categLabel, desc = "Profiles Plots of clusterized genes selected from multiple comparisons")
    }
    
  ###---------------------------------------------------------------------------------
  ### Tercera part : Extraiem els simbols i els escrivim en format .csv.
  ###---------------------------------------------------------------------------------
  

  ### Extraiem els gens en el mateix ordre que els dona el HeatMap

  ### FERRAN (2012/12/10)
  ### correccio perque escrigui be els CSV un cop afegits els geneSymbols als rownames dels heat maps
  ### SOLS FUNCIONA SI TENIM my.symbols !!!
  raw.gNames <- rev(rownames(exprs2Cluster[hm$rowInd, ])) ### Posem rev() aqui per fer coincidir l'ordre de files del csv de sortida amb les files del heatmap
  gNames <- unlist(lapply(strsplit(raw.gNames, "[.]"), function(l) l[[1]][1]))

  if(!is.null(my.symbols))
  {
    mySymbols <- my.symbols[gNames]
    my.genes <- data.frame(symbol = mySymbols, ID = gNames, GrupID = cutN[gNames], exprs2Cluster[raw.gNames, names.ord])    
  }else{
    if (!is.null(anotPackage)){
      my_SYMBOL_env <- eval(parse(text = paste(anotPackage, "SYMBOL",sep = "")))
      mySymbols <- unlist(mget(gNames, my_SYMBOL_env, ifnotfound=NA))
      my.genes <- data.frame(symbol = mySymbols, ID = gNames, GrupID = cutN[gNames], exprs2Cluster[gNames, names.ord])
    }else{
      my.genes <- data.frame(ID = gNames, GrupID = cutN[gNames], exprs2Cluster[gNames, names.ord])
    }
  }
    
  ## if (!is.null(anotPackage)){
  ##   my_SYMBOL_env <- eval(parse(text = paste(anotPackage, "SYMBOL",sep = "")))
  ##   mySymbols <- unlist(mget(gNames, my_SYMBOL_env))
  ##   my.genes <- data.frame(symbol = mySymbols, ID = gNames, GrupID = cutN[gNames], exprs2Cluster[gNames, names.ord])
  ## }else{
  ##   my.genes <- data.frame(ID = gNames, GrupID = cutN[gNames], exprs2Cluster[gNames, names.ord])
  ## }
#  genesInClustersFName <- paste("genesInClusters", comparisonName, "csv", sep = ".")
#  write.csv2(my.genes, file.path(outputDir, genesInClustersFName))
  genesInClustersFName <- paste("genesInClusters", comparisonName, sep = ".")
  csvType<-ifelse(is.null(csvType), "csv2", csvType)
  write2csv(my.genes, fileName = genesInClustersFName, csv = csvType, outputDir = outputDir)
  addToLinksFile(fileOfLinks, paste(genesInClustersFName,  substr(csvType[1], 1, 3), sep="."),
                 categ = categLabel, desc = "Assignment of each gene to a cluster of made from genes selected from multiple comparisons")
  } ### Tot aixo es fa nomes si 'dendroRow' es cert

  return(hm)
}

### 4.3 Funcio per executar el 'clusterAnalysis' a partir d'un objecte parametres.

doClusterAnalysis <- function(clustPar)
{

  p <- clustPar[[1]]

  #############################################################################
  ### Actualitzacio 2011
  ###
  ### Els ENTREZIDs i SYMBOLIDs es pasen com a text i s'avaluen internament
  ############################################################################

  if(!is.null(p$my.symbols)){
    SymbolIDs <-  eval(parse(text = p$my.symbols)) #  Seran el symbolsTable
  }

  ### Aquestes variables es pasen ara a la crida a lmAnalysis
  ###############################################################################
  
  
  if (!is.null(p$expresFileName)){  
      expres <- loadFromFile (file.path(p$outputDir, p$expresFileName))
  }else{
      if (!is.null(p$dades)) {
	       expres <- eval(parse(text = p$dades)) # Posar-hi un tryCatch per poder sortir si d??na error!!!
      }else{
        	stop("Error, Cal definir o les dades o el nom de l'arxiu")
      }
  }

  if (!is.null(p$geneListFName)){  
      genes2cluster <- loadFromFile (file.path(p$outputDir, p$geneListFName))
  }else{
      if (is.null(p$genes2cluster)) {
        	stop("Error, Cal definir l'arxiu que conte la llista de gens o passar una variable que la contingui")
      }else{
        genes2cluster <- p$genes2cluster
      }
  }

  clust <- clusterAnalysis(expres = expres, 
                           genes = genes2cluster, 
                           samples = p$samples2cluster, 
                           sampleNames = p$sampleNames,
                           comparisonName = p$comparisonName,
                           anotPackage = p$anotPackage,
                           my.symbols = SymbolIDs, # era p$my.symbols,
                           outputDir = p$outputDir, 
                           fileOfLinks = p$fileOfLinks,                              
                           numClusters = p$numClusters,
                           rowDistance = p$rowDistance,
                           colDistance = p$colDistance,
                           RowVals = p$RowVals,
                           ColVals = p$ColVals,
                           escala = p$escala,
                           colorsSet = p$colorsSet,
                           densityInfo = p$densityInfo,
                           colsForGroups = p$colsForGroups,
                           cexForColumns = p$cexForColumns,
                           cexForRows = p$cexForRows,
                           Title = p$Title,
                           csvType = p$csvType)
  
  return(clust)
}


#####################################################################################################################
### Funcions per a fer l'analisi d'enriquiment de la GO
#####################################################################################################################

### GOTerms2Genes.sql: Per a un vector amb els EntrezID, extrau del paquet d'anotacions els Gene Symbols usant els 
###                    metodes consultes SQL
###                 
###                 
###
### Parametres:
###
###        hgResult :
###     anotPackage : Paquet d'anotacions que utilitzara
###
### OBERVACIONS IMPORTANTS:
###

GOTerms2Genes.sql <- function(hgResult, anotPackage)
{
  selectedGOTerms <- intersect(names(geneIdUniverse(hgResult)), summary(hgResult)[, 1])
  selectedGO<- geneIdUniverse(hgResult)[selectedGOTerms]
       if (runMulticore ==1 || runMulticore ==3) { 
           selectedGenes <- mclapply(selectedGO, function(x) {intersect(geneIds(hgResult),x)})
         } else {
           selectedGenes <- lapply(selectedGO, function(x) {intersect(geneIds(hgResult),x)})
       }
  
  sql.ENTREZSYMBOL <- "SELECT gene_id, symbol
                       FROM genes, gene_info
                       WHERE genes._id=gene_info._id"

  if (regexpr(".db", anotPackage) < 0)
  {
     genesAnnot <- dbGetQuery(eval(parse(text = paste(anotPackage, "_dbconn()", sep = ""))), sql.ENTREZSYMBOL)
  }else{
    genesAnnot <- dbGetQuery(eval(parse(text = paste(substr(anotPackage, 1, nchar(anotPackage) - 3), "_dbconn()", sep = ""))), sql.ENTREZSYMBOL)
  }
  
  if (runMulticore ==1 || runMulticore ==3) { 
      selectedSymb <- mclapply(selectedGenes, function(x) {genesAnnot[which(unlist(genesAnnot) %in% x), ]})
      selectedSymb <- mclapply(selectedSymb, function(x) {x <- x[, -1]})
    } else {
      selectedSymb <- lapply(selectedGenes, function(x) {genesAnnot[which(unlist(genesAnnot) %in% x), ]})
      selectedSymb <- lapply(selectedSymb, function(x) {x <- x[, -1]})         
   }

  return(selectedSymb)
}


### goLinksTest: Funcio necessaria per a afegir el link html a AmiGO d'un identificador GO
###
### Parametres:
###
###    my.GOIs : vector character que conte els identificadors GO a linkar
###
### OBERVACIONS IMPORTANTS:
###

goLinksTest <- function(my.GOIDs)
{
 if (runMulticore ==1 || runMulticore ==3) { 
    GOIDslinked <- unlist(mclapply(my.GOIDs, function(x) {paste("<a href=\"http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&query=", x, "\">",
                                                            x,
                                                            "</a>",
                                                            sep = "")}))
  } else {
    GOIDslinked <- unlist(lapply(my.GOIDs, function(x) {paste("<a href=\"http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&query=", x, "\">",
                                                            x,
                                                            "</a>",
                                                            sep = "")}))    
  }

  return(as.character(GOIDslinked))
}
######################################################################

### enrichment_Analysis: Selecciona els genes univers, els seleccionats i crida a la funcio d'enriquiment
###                 
###                 
###
### Parametres:
###
###          EntrezIDs :
###          anotFName : 
###        anotPackage : Paquet d'anotacions que utilitzara
###          outputDir :
###           universe : Per defecte prenNULL,
###               pval : Per defecte pren 0.05
###          min.count : Per defecte pren 3
###       addGeneNames : Per defecte pren TRUE
###         ontologias : Per defecte pren c("MF", "BP", "CC")
###     testDirections : Per defecte pren c("over", "under")
###      contrast.Name :
###
### OBERVACIONS IMPORTANTS:
###
###    COMPTE!!!! S'assumeix que se li pass una llista d'identificadors per a l'univers

enrichment_Analysis <- function(EntrezIDs,
                                anotFName,
                                anotPackage,
                                outputDir,
                                universe = NULL,
                                pval = 0.05,
                                min.count = 3,
                                addGeneNames = TRUE,
                                ontologias = c("MF", "BP", "CC"),
                                testDirections = c("over", "under"),
                                contrast.Name)
{
  
  universe <- unique(universe)
  EntrezIDs <- unique(EntrezIDs)

  params <- new("GOHyperGParams",
                geneIds = EntrezIDs, 
                universeGeneIds = universe,
                annotation = anotPackage,
                ontology = "MF",
                pvalueCutoff = pval,
                testDirection = "over")

  resCum <- NULL; resum <- NULL
  ontoTitle <- c(MF="Molecular Function", BP="Biological Process", CC="Cellular Component")

  for (onto in ontologias)
  {
    ontology(params) <- onto

    dir.i <- 1

    for (direction in testDirections)
    {
      testDirection(params) <- direction

      hgResult <- hyperGTest(params)

      if (addGeneNames)
      {
        EnrichedGOTerms <- as.character(summary(hgResult)[,1])
        if (length(EnrichedGOTerms) > 0)
        {
          selectedSymbols <- GOTerms2Genes.sql(hgResult, anotPackage)
          genesInHg <- sapply(selectedSymbols, function(x) paste(x, collapse=", "))
          reshgTest <- cbind(summary(hgResult), GeneNames = genesInHg)
        }else{
          reshgTest <- summary(hgResult)
        }
      }else{
        reshgTest <- summary(hgResult)
      }

      colorTest <- ifelse(direction=="over", "red;\">", "green;\">")
      reshgTest <- cbind(reshgTest,
                         OverUnder = rep(paste("<center><span style=\"color:", colorTest, direction, "</span></center>", sep = ""), dim(reshgTest)[1]))
      reshgTest[, 1] <- goLinksTest(reshgTest[, 1])
      colnames(reshgTest)[1] <- "GOID"
      
      if (dir.i==1)
      {
        reshgTest2table <- reshgTest
      }else{
        reshgTest2table <- rbind(reshgTest2table, reshgTest)
      }

      dir.i <- dir.i + 1
    }

    if (dim(reshgTest2table)[1]!=0)
    {
      reshgTest2table <- cbind(Ontology = paste("<center>", onto, "</center>", sep=""),
                               reshgTest2table)

      resum <- rbind(resum, reshgTest2table)
    }
  }

  if (!is.null(resum))
  {
    if (addGeneNames)
    {
      my.table <- resum[, c(1, 2, 8, 9, 7, 6, 5, 4, 3, 10)]
    }else{
      my.table <- resum[, c(1, 2, 8, 7, 6, 5, 4, 3, 9)]
    }
  
    write.htmltable(x = my.table, 
                    file = file.path(outputDir, anotFName), 
                    title =  paste("GO Enrichment Analysis", contrast.Name, sep = " "),
                    open = "wt")
    sortable.html.table(df = my.table,
                        output.file = paste0(anotFName, "-sortable.html"),
                        output.directory = outputDir,
                        page.title = "GO Enrichment Analysis" )
  }

}


### GOAnalysis: Selecciona els genes univers, els seleccionats i crida a la funcio d'enriquiment
###                 
###                 
###
### Parametres:
###
###            fitMain :
###     whichContrasts :
###    comparison.Name :
###          outputDir :
###        anotPackage :
###             my.IDs :
###       addGeneNames : Per defecte TRUE
###        fileOfLinks :
###       cutoffMethod : pvalor pel que se seleccionara els gens. Per defecte pren "adjusted"
###     P.Value.cutoff : Per defecte pren rep(0.05, length(whichContrasts)
###               pval : nivell de significacio pel testar si una categoria GO es significativa. Per defecta pren 0.05
###          min.count : nombre minim de vegades que ha sortit referenciada una ategoria GO per a ser reportada. Per defecte pren 3
###         ontologias : Ontologies de la GO per a les que es vol fer l'enriquiment. Per defecte pren el vector c("MF", "BP", "CC")
###     testDirections : Fara un test per les veure quines son les categories sobrerepresendades y/o infrarepresentades. Per defecte pren c("over", "under")
###
### OBSERVACIONS IMPORTANTS:
###
###   

GOAnalysis <- function(fitMain,
                       whichContrasts, 
                       comparison.Name,
                       outputDir, 
                       anotPackage,
                       my.IDs,
                       addGeneNames = TRUE,
                       fileOfLinks,
                       thrLogFC=NULL, 
                       cutoffMethod = c("adjusted", "unadjusted"), 
                       P.Value.cutoff = rep(0.05, length(whichContrasts)),
                       pval = 0.05,
                       min.count = 3,
                       ontologias = c("MF", "BP", "CC"),
                       testDirections = c("over", "under"),
                       minNumGens = 0)
{

  categLabel <- "GO"

  if ((!is.null(fitMain)) && (!is.null(fitMain$contrasts[, whichContrasts])))
  {
    if (is.null(my.IDs))
    {
      stopifnot(require(old2db (anotPackage), character.only=T))
      myenvirENTREZID <- eval(parse(text = paste(anotPackage, "ENTREZID", sep = "")))
      geneUniverse <- unlist(mget(unique( fitMain$genes[,1]), env = myenvirENTREZID, ifnotfound = NA))
    }else{
      geneUniverse <- unlist(my.IDs)
      geneUniverse <- unique(geneUniverse[geneUniverse != "---"])
    }

#   write2csv(data.frame(geneUniverse), fileName = "universe", csv = csv, outputDir = tempDir)  

    for (i in whichContrasts)
    {
      if(is.null(thrLogFC))
        {
          thrLogFC <- 0
        }else{
          thrLogFC <- abs(thrLogFC)
        }

      top.Diff <- topTable(fitMain, coef = i, n = nrow(fitMain$t), adjust = "fdr", lfc=thrLogFC)  # Seleccionem per FDR y mes
                                                                                                  # avall tallem per ajustat o no

      if (cutoffMethod=="adjusted")
      {
        top.Diff.selected.up  <- top.Diff[(top.Diff$adj.P.Val < P.Value.cutoff[i]) & (top.Diff$t >0), ]
        top.Diff.selected.down <- top.Diff[(top.Diff$adj.P.Val < P.Value.cutoff[i]) & (top.Diff$t < 0), ]
      }else{
        top.Diff.selected.up  <- top.Diff[(top.Diff$P.Value < P.Value.cutoff[i]) & (top.Diff$t >0), ]
        top.Diff.selected.down <- top.Diff[(top.Diff$P.Value < P.Value.cutoff[i]) & (top.Diff$t < 0), ]
      }
	
      if (!is.null(top.Diff$ID)){

      gNames.up <- top.Diff.selected.up$ID
      gNames.down <- top.Diff.selected.down$ID
	}else{
	      gNames.up <- rownames(top.Diff.selected.up)
      	      gNames.down <- rownames(top.Diff.selected.down)
	}
      if (!(is.null(my.IDs)))
      {
        selectedEntrezIds.up <- unlist(my.IDs[gNames.up])
        selectedEntrezIds.up <- unique(selectedEntrezIds.up[selectedEntrezIds.up != "---"])

        selectedEntrezIds.down <- unlist(my.IDs[gNames.down])
        selectedEntrezIds.down <- unique(selectedEntrezIds.down[selectedEntrezIds.down != "---"])

#       write2csv(data.frame(selectedEntrezIds.up), fileName = "up", csv = fileType, outputDir = tempDir)
#       write2csv(data.frame(selectedEntrezIds.down), fileName = "down", csv = fileType, outputDir = tempDir)

      }else{

        stopifnot(require(old2db (anotPackage), character.only=T))

        myenvirENTREZID <- eval(parse(text = paste(anotPackage, "ENTREZID", sep = "")))

        selectedEntrezIds.up <- unlist(mget(unique(gNames.up), env = myenvirENTREZID, ifnotfound = NA))
        selectedEntrezIds.down <- unlist(mget(unique(gNames.down), env = myenvirENTREZID, ifnotfound = NA))

        
#     selectedEntrezIds.up   <- unique(gNames.up)
#     selectedEntrezIds.down <- unique(gNames.down)
      }

 #    write2csv(data.frame(selectedEntrezIds.up), fileName = "up", csv = fileType, outputDir = tempDir)
 #    write2csv(data.frame(selectedEntrezIds.down), fileName = "down", csv = fileType, outputDir = tempDir)

      
      contrast.Name <- colnames(fitMain$contrasts)[i]

#      if (is.null(selectedEntrezIds.up))
      if (is.null(selectedEntrezIds.up)|(length(selectedEntrezIds.up)<=minNumGens))
      {
        warning(paste("There are not enough genes for cutoff ", P.Value.cutoff[i], " and t > 0 in comparison", contrast.Name, sep = " "))
        cat(paste("There are not enough genes for cutoff ", P.Value.cutoff[i], " and t > 0 in comparison", contrast.Name, "\n", sep = " "))        
      }else{
        eaUpFName <- paste("SignificantGO", comparison.Name, contrast.Name, "Up", sep = ".")
        eA.up <- enrichment_Analysis(EntrezIDs = selectedEntrezIds.up,
                                     anotFName = eaUpFName,
                                     anotPackage = anotPackage,
                                     outputDir = outputDir,
                                     universe = geneUniverse,
                                     pval = pval,
                                     min.count = min.count,
                                     addGeneNames = addGeneNames,
                                     ontologias = ontologias,
                                     testDirections = testDirections,
                                     contrast.Name = paste("for up-regulated genes in comparison", contrast.Name, sep = ": "))

        
        addToLinksFile(fileOfLinks,
                       paste(eaUpFName,"html", sep="."),
                       categ = categLabel,
                       desc = paste("Enrichment Analysis for up-regulated genes in comparison", contrast.Name, sep = ": "))
        
        
      }

#      if (is.null(selectedEntrezIds.down))
      if (is.null(selectedEntrezIds.down)|(length(selectedEntrezIds.down)<=minNumGens))
      {
        warning(paste("There are not enough genes for cutoff ", P.Value.cutoff[i], " and t < 0 in comparison", contrast.Name, sep = " "))
        cat(paste("There are not enough genes for cutoff ", P.Value.cutoff[i], " and t < 0 in comparison", contrast.Name, "\n", sep = " "))
      }else{
        eaDownFName <- paste("SignificantGO", comparison.Name, contrast.Name, "Down", sep = ".")
        eA.down <- enrichment_Analysis(EntrezIDs = selectedEntrezIds.down, 
                                       anotFName = eaDownFName,
                                       anotPackage = anotPackage, 
                                       outputDir = outputDir,
                                       universe = geneUniverse, 
                                       pval = pval,
                                       min.count = min.count,
                                       addGeneNames = addGeneNames,
                                       ontologias = ontologias,
                                       testDirections = testDirections,
                                       contrast.Name = paste("for down-regulated genes in comparison", contrast.Name, sep = ": "))

        addToLinksFile(fileOfLinks, paste(eaDownFName, "html", sep="."),
                       categ = categLabel,
                       desc = paste("Enrichment Analysis for down-regulated genes in comparison", contrast.Name, sep = ": "))
      }
    }
  }
  
  return(NULL)
}


### doGOAnalysis: Construeix una taula amb el resultat de fer l'analisi d'enriquiment
###                 
###                 
###
### Parametres:
###
###    GOPar: Llista amb la seguent
###
### OBERVACIONS IMPORTANTS:
###
###    En un futur no llunya hauria de ser una classe S4

doGOAnalysis <- function(GOPar)
{

  p <- GOPar[[1]]

  #############################################################################
  ### Actualitzacio 2011
  ###
  ### Els ENTREZIDs es pasen com a text i s'avaluen internament
  ############################################################################

  if(!is.null(p$my.IDs)){
    EntrezIDs <-  eval(parse(text = p$my.IDs)) #  Seran el EntrezTable
  }

  ### Aquestes variables es pasen ara a la crida a GOAnalysis
  ###############################################################################
    
 
  if (!is.null(p$fitFileName))
  {
    fitMain <- loadFromFile(file.path(p$outputDir, p$fitFileName))
  }else{
    if (!is.null(p$fitMain))
    {
      fitMain <- eval(parse(text = p$fitMain)) # Posar-hi un tryCatch per poder sortir si d??na error!!!
    }else{
      stop("Error, Cal subministrar un nom d'arxiu o d'objecte 'fitMain'")
    }
  }

  if (!is.null(p$my.IDs))
  {
    my.IDs <- midSinonims(p$my.IDs)
  }
  
  GOResult <- GOAnalysis(fitMain = fitMain,
			 whichContrasts = p$whichContrasts, 
			 comparison.Name = p$comparisonName, 
			 outputDir = p$outputDir, 
			 anotPackage = orgPackage,
			 my.IDs = EntrezIDs, # era p$my.IDs,
                         addGeneNames = p$addGeneNames,
			 fileOfLinks = p$fileOfLinks,   
                         thrLogFC = p$minLogFC,
			 cutoffMethod = p$cutoffMethod, 
			 P.Value.cutoff = p$P.Value.cutoff,
                         pval = p$pvalGOterms,
                         min.count = p$min.count,
                         ontologias = p$ontologias,
                         testDirections = p$testDirections,
                         minNumGens = p$minNumGens)

  return(GOResult)
}

#####################################


KEGGEnrAn <- function(EntrezIDs,
                      anotFName,
                      anotPackage,
                      organisme = "hsa",                      
                      outputDir,
                      universe,
                      pval,
                      min.count,
                      addGeneNames,
                      contrast.Name)
#
# Fa un test hipergeometric tant per KEGG com per GO segons els parametres que li passem
#
{
#  require("annotate")
#  require("GOstats")
#  require("hwriter")
 
  universe <- unique(universe)
  EntrezIDs <- unique(EntrezIDs)
  
  # Crea l'hiperparametre
  param <- new("KEGGHyperGParams",
               geneIds = EntrezIDs,
               universeGeneIds = universe,
               annotation = anotPackage,
               pvalueCutoff = pval)

  anotPack <- annotation(param)

  getSymbol <- function (x)
  {
    if (length(x)>0)
    {
      simbols <- getSYMBOL(x, anotPack)
    }else{
      simbols <- NULL
    }

    return(simbols)
  }
  
  hyperRes <- hyperGTest(param)
  sumari <- summary(hyperRes, p=pvalueCutoff(param))
  fName <- paste("KEGG Enrichment Analysis", contrast.Name, sep = " ") # Informe en HTML
     
  if (addGeneNames)
  {    
    EnrichedKEGGTerms <- as.character(sumari[, 1])
        
    if (length(EnrichedKEGGTerms) > 0)
    {
      selectedSymbols <- GOTerms2Genes.sql(hyperRes, anotPackage)
      genesInHg <- sapply(selectedSymbols, function(x) paste(x, collapse = ", "))
      Report <- cbind(sumari, GeneNames = genesInHg)
    }else{
      Report <- sumari
    }
  }else{
    Report <- sumari
  }

  Report[, 1] <- paste("<a href=\"http://www.genome.jp/dbget-bin/www_bget?path:", organisme, Report[, 1], "\">", Report[, 1], "</a>", sep="") 

  ReportSig <- Report[1:nrow(sumari),]

  write.htmltable(x = Report[1:nrow(sumari),],
                  file = file.path(outputDir, anotFName),
                  title =  paste("KEGG Enrichment Analysis", contrast.Name, sep = " "),)
  sortable.html.table(df = Report[1:nrow(sumari),],
                      output.file = paste0(anotFName, "-sortable.html"),
                      output.directory = outputDir,
                      page.title = "KEGG Enrichment Analysis")

}
 
### KEGGAnalysis: Selecciona els genes univers, els seleccionats i crida a la funcio d'enriquiment
###                 
###                 
###
### Parametres:
###
###            fitMain :
###     whichContrasts :
###    comparison.Name :
###          outputDir :
###        anotPackage :
###          organisme : parametre necessari per a poder crear links actius en la taula html. Per defecte pren "hsa" (Homo sapiens),
###                      pero pot pendre "mmu" (Mus musculus)
###                                      "rno" (Rattus norvegicus)
###                                      "bta" (Bos taurus)
###                                      "dre" (Danio rerio)
###                                      "ssc" (Sus scrofa)
###             my.IDs :
###       addGeneNames : Per defecte TRUE
###        fileOfLinks :
###       cutoffMethod : pvalor pel que se seleccionara els gens. Per defecte pren "adjusted"
###     P.Value.cutoff : Cutoff del pvalor (ajustat) pel que es tallara la llista de gens generats a la topTable.
###                      Per defecte pren rep(log(1), length(whichContrasts)
###           thrLogFC : Cutoff del log Fold Change pel que es tallara la llista  de gens generats a la topTable. Per defecte pren NULL
###               pval : nivell de significacio pel testar si una categoria KEGG es significativa. Per defecte pren 0.05
###
### OBSERVACIONS IMPORTANTS:
###
###   

KEGGAnalysis <- function(fitMain,
                         whichContrasts, 
                         comparison.Name,
                         outputDir, 
                         anotPackage,
                         organisme,
                         my.IDs,
                         addGeneNames = TRUE,
                         fileOfLinks,
                         cutoffMethod = c("adjusted", "unadjusted"), 
                         P.Value.cutoff = rep(0.05, length(whichContrasts)),
                         pval = 0.05,
                         thrLogFC = NULL,
                         minNumGens = 0)
{
  stopifnot(require(old2db(anotPackage), character.only = TRUE))
  
  categLabel <- "KEGG"

  if ((!is.null(fitMain)) && (!is.null(fitMain$contrasts[, whichContrasts])))
  {
    if (is.null(my.IDs))
    {
      myenvirENTREZID <- eval(parse(text = paste(anotPackage, "ENTREZID", sep = "")))
      geneUniverse <- unlist(mget(unique( fitMain$genes[,1]), env = myenvirENTREZID, ifnotfound = NA))
    }else{
      geneUniverse <- unlist(my.IDs)
      geneUniverse <- unique(geneUniverse[geneUniverse != "---"])
    }

    for (i in whichContrasts)
    {
      if(is.null(thrLogFC))
      {
        thrLogFC <- 0
      }else{
        thrLogFC <- abs(thrLogFC)
      }

      top.Diff <- topTable(fitMain, coef = i, n = nrow(fitMain$t), adjust = "fdr", lfc = thrLogFC) 

      if (cutoffMethod=="adjusted")
      {
        top.Diff.selected <- top.Diff[top.Diff$adj.P.Val < P.Value.cutoff[i], ]
      }else{
        top.Diff.selected <- top.Diff[top.Diff$P.Value < P.Value.cutoff[i], ]
      }

      if (!is.null(top.Diff.selected$ID)){
           gNames <- top.Diff.selected.up$ID
      }else{
	gNames <- rownames(top.Diff.selected)
	}


      if (!(is.null(my.IDs)))
      {
        selectedEntrezIds <- unlist(my.IDs[gNames])
        selectedEntrezIds <- unique(selectedEntrezIds[selectedEntrezIds != "---"])
      }else{
        myenvirENTREZID <- eval(parse(text = paste(anotPackage, "ENTREZID", sep = "")))
        selectedEntrezIds <- unlist(mget(unique(gNames), env = myenvirENTREZID, ifnotfound = NA))
      }
     
      contrast.Name <- colnames(fitMain$contrasts)[i]

      if (is.null(selectedEntrezIds)|(length(selectedEntrezIds) <= minNumGens))
      {
        warning(paste("There are not enough genes for cutoff ", P.Value.cutoff[i], " in comparison", contrast.Name, sep = " "))
        cat(paste("There are not enough genes for cutoff ", P.Value.cutoff[i], " in comparison", contrast.Name, "\n", sep = " "))        
      }else{
        eaFName <- paste("SignificantKEGG", comparison.Name, contrast.Name, sep = ".")
        eA <- KEGGEnrAn(EntrezIDs = selectedEntrezIds,
                        anotFName = eaFName,
                        anotPackage = anotPackage,
                        organisme = organisme,
                        outputDir = outputDir,
                        universe = geneUniverse,
                        pval = pval,
                        min.count = minNumGens,
                        addGeneNames = addGeneNames,
                        contrast.Name = paste("for regulated genes in comparison", contrast.Name, sep = ": "))
        
        addToLinksFile(fileOfLinks,
                       paste(eaFName,"html", sep="."),
                       categ = categLabel,
                       desc = paste("KEGG Enrichment Analysis regulated genes in comparison", contrast.Name, sep = ": "))
       }
    }
  }
  return(NULL)
}

### doKEGGAnalysis: Construeix una taula amb el resultat de fer l'analisi d'enriquiment
###                 
###                 
###
### Parametres:
###
###    KEGGPar: Llista amb la seguent
###
### OBERVACIONS IMPORTANTS:
###
###    En un futur no llunya hauria de ser una classe S4

doKEGGAnalysis <- function(KEGGPar)
{

  p <- KEGGPar[[1]]

  ############################################################################
  ### Actualitzacio 2011
  ###
  ### Els ENTREZIDs es pasen com a text i s'avaluen internament
  ############################################################################

  if(!is.null(p$my.IDs)){
    EntrezIDs <-  eval(parse(text = p$my.IDs)) #  Seran el EntrezTable
  }

  ### Aquestes variables es pasen ara a la crida a GOAnalysis
  ###############################################################################
  
 
  if (!is.null(p$fitFileName))
  {
    fitMain <- loadFromFile(file.path(p$outputDir, p$fitFileName))
  }else{
    if (!is.null(p$fitMain))
    {
      fitMain <- eval(parse(text = p$fitMain)) # Posar-hi un tryCatch per poder sortir si dona error!!!
    }else{
      stop("Error, cal subministrar un nom d'arxiu o d'objecte 'fitMain'")
    }
  }

#
# if (!is.null(p$my.IDs))
# {
#   my.IDs <- midSinonims(p$my.IDs)
# }
#
  
  KEGGResult <- KEGGAnalysis(fitMain = fitMain,
                             whichContrasts = p$whichContrasts, 
                             comparison.Name = p$comparisonName, 
                             outputDir = p$outputDir, 
                             anotPackage = orgPackage,
                             organisme = organisme,
                             my.IDs = EntrezIDs, # era p$my.IDs
                             addGeneNames = p$addGeneNames,
                             fileOfLinks = p$fileOfLinks,   
                             cutoffMethod = p$cutoffMethod, 
                             P.Value.cutoff = p$P.Value.cutoff,
                             pval = p$pvalKEGGterms,
                             thrLogFC = p$minLogFC,
                             minNumGens = p$minNumGens)
  
  return(KEGGResult)
}

#####################################################################################################################
### Funcions per a escriure el fitxer de Results_Files
#####################################################################################################################

### printHeader: Crea la cap??alera del fitxer html de resultats
###                 
###
### Parametres:
###
###   FileName : Nom del fitxer de resultats
###
### OBERVACIONS IMPORTANTS:
###
###   

printHeader <- function(FileName = "ResultFiles")
{
  outfile <- file(paste(FileName, ".html", sep = ""), open = "wt")
  
  txt <- "<html> \n"
  txt <- paste(txt, "<head> \n", sep = "")
  txt <- paste(txt, "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\"/>", sep = "")
  txt <- paste(txt, "<title>List of results files</Title>", sep = "")
  txt <- paste(txt, "</head> \n", sep = "")
  txt <- paste(txt, "<style> \n", sep = "")
  txt <- paste(txt, "  <!--td {font-family: Helvetica,Arial; font-size: 14px;}--> \n", sep = "")
  txt <- paste(txt, "  <!--h1 {font-family: Helvetica,Arial; font-size: 22px;}--> \n", sep = "")
  txt <- paste(txt, "</style> \n", sep = "")
  txt <- paste(txt, "<body bgcolor=\"White\" text=\"Black\"> \n", sep = "")
  
  cat(txt, file = outfile, sep = "")
  
  close(outfile)
}


### printGroupHeader: Construeix la cap??alera de presentacio que es veura a la pagina
###                 
###
### Parametres:
###
###   FileName : Nom del fitxer de resultats
###        UEB : Per defecte pren valor TRUE. Si es FALSE crea la cap??alarea adaptada a la UB
###
### OBERVACIONS IMPORTANTS:
###
###   

printGroupHeader <- function(FileName = "ResultFiles", UEB = TRUE)
{
  outfile <- file(paste(FileName, ".html", sep = ""), open = "at")

             txt <- "<table width=\"100%\"  border=\"0\"> \n"
  txt <- paste(txt, "   <tr><td width=\"24%\"> \n", sep = "")
  txt <- paste(txt, "	    <div align=\"center\"> \n", sep = "")
  
  if(UEB)
  {
    txt <- paste(txt, "		 <a href=\"http://www.vhir.org\" target=\"z\"> \n", sep = "")
    txt <- paste(txt, "			<img src=\"images/IR.jpg\" width=\"195\" height=\"73\" border=\"0\"> \n", sep= "")
    txt <- paste(txt, "		 </a> \n", sep= "")
    txt <- paste(txt, "	    </div> \n", sep= "")
    txt <- paste(txt, "	</td> \n", sep= "") 
    txt <- paste(txt, "	<td width=\"53%\"> \n", sep= "")
    txt <- paste(txt, "	    <h1 align=\"center\">Unitat d'Estad&iacute;stica i Bioinform&agrave;tica</h1> \n", sep= "")
    txt <- paste(txt, "	    <div align=\"center\">Vall d'Hebron Institut de Recerca</div> \n", sep= "")
    txt <- paste(txt, "       </td> \n", sep= "") 
    txt <- paste(txt, "	<td width=\"23%\"> \n", sep= "")
    txt <- paste(txt, "	    <div align=\"center\"> \n", sep= "")
    txt <- paste(txt, "		 <a href=\"http://ueb.vhir.org\" target=\"z\"> \n", sep= "")
    txt <- paste(txt, "			<img src=\"images/UEBblanc.jpg\" width=\"204\" height=\"48\" border=\"0\"> \n", sep= "")
  }else{

    txt <- paste(txt, "		 <a href=\"http://www.ub.edu\" target=\"z\"> \n", sep = "")
    txt <- paste(txt, "			<img src=\"images/UB.jpg\" width=\"195\" height=\"73\" border=\"0\"> \n", sep= "")
    txt <- paste(txt, "		 </a> \n", sep= "")
    txt <- paste(txt, "	    </div> \n", sep= "")
    txt <- paste(txt, "	</td> \n", sep= "") 
    txt <- paste(txt, "	<td width=\"53%\"> \n", sep= "")
    txt <- paste(txt, "	    <h1 align=\"center\">Grup de Recerca en Estad&iacute;stica i Bioinform&agrave;tica</h1> \n", sep= "")
    txt <- paste(txt, "	    <div align=\"center\">Departament d'Estad&iacute;stica - Universitat de Barcelona</div> \n", sep= "")
    txt <- paste(txt, "       </td> \n", sep= "") 
    txt <- paste(txt, "	<td width=\"23%\"> \n", sep= "")
    txt <- paste(txt, "	    <div align=\"center\"> \n", sep= "")
    txt <- paste(txt, "		 <a href=\"http://estbioinfo.stat.ub.es\" target=\"z\"> \n", sep = "")
    txt <- paste(txt, "			<img src=\"images/Logo%20EstBioinfo.jpg\" width=\"204\" height=\"48\" border=\"0\"> \n", sep= "")
  }

  txt <- paste(txt, "		 </a> \n", sep= "")
  txt <- paste(txt, "	    </div> \n", sep= "")
  txt <- paste(txt, "	</td>  \n", sep= "")
  txt <- paste(txt, "   </tr> \n", sep= "")
  txt <- paste(txt, "</table> \n", sep= "")

  cat(txt, file = outfile, sep = "")
  
  close(outfile)
}


### printAnalysisDetails: Construeix la taula que conte la informacio general de l'estudi, com dades de contacte, 
###                       titol de l'analisi
###
### Parametres
###
###    FileName : Nom del fitxer de resultats
###   Info.list : Llista amb la informacio general 
###
### Exemple
###
###    myInfolist <- list(To = "Nom i cognoms de l'usuari",
###                       Description = "Diferentially expressed genes associated with...",
###                       Analysts = "Nom del primer analista and Alex Sanchez",
###                       Contact = "Alex Sanchez (alesanchez@ir.vhebron.net)")
###    printAnalysisDetails(FileName = "ResultFiles", Info.list =  myInfoList)
###
###    # Aixo creara un fitxer que contindra una taula de l'estil
###    <table width="100%" border="0"> 
###        <tr><td height="35" colspan="2" valign="center" bgcolor="#d0d0f0"> 
###                <h1>Statistical Analysis Report</h1>
###             </td>
###        </tr> 
###        <tr height="35"><td width="15%" align="right" bgcolor="#e0e0f0"><b>To</b></td> 
###                        <td width="85%" bgcolor="#f0f0ff"><i>Nom i cognoms de l'usuari</i></td> 
###        </tr> 
###        <tr height="35"><td align="right" bgcolor="#e0e0f0"><b>Description</b></td> 
###                        <td  bgcolor="#f0f0ff"><i>Diferentially expressed genes associated with...</i></td> 
###	                   <br> 
###	   </tr> 
###	   <tr height="35"><td align="right" bgcolor="#e0e0f0"><b>Analysts</b></td> 
###	                   <td  bgcolor="#f0f0ff"><i>Nom del primer analista and Alex Sanchez</i></td>
###                        <br> 
###	   </tr> 
###	   <tr height="35"><td align="right" bgcolor="#e0e0f0"><b>Contact</b></td> 
###	                   <td bgcolor="#f0f0ff"><i>Alex Sanchez (alesanchez@ir.vhebron.net)</i></td> 
###	      <br> 
###	   </tr> 
###    </table> 
###    <hr> 

printAnalysisDetails <- function(FileName = "ResultFiles", Info.list)
{
  outfile <- file(paste(FileName, ".html", sep = ""), open = "at")

	     txt <- "<table width=\"100%\" border=\"0\"> \n"  
  txt <- paste(txt, "     <tr><td height=\"35\" colspan=\"2\" valign=\"center\" bgcolor=\"#bb66aa\"> \n", sep = "")
  txt <- paste(txt, "                <h1>Statistical Analysis Report</h1></td></tr> \n", sep = "")
  txt <- paste(txt, "     <tr height=\"35\"><td width=\"15%\" align=\"right\" bgcolor=\"#e0e0f0\"><b>To</b></td> \n", sep = "")
  txt <- paste(txt, "         <td width=\"85%\" bgcolor=\"#f0f0ff\"><i>", Info.list$To, "</i></td> \n", sep = "")
  txt <- paste(txt, "     </tr> \n", sep = "")
  txt <- paste(txt, "     <tr height=\"35\"><td align=\"right\" bgcolor=\"#e0e0f0\"><b>Description</b></td> \n", sep = "")
  txt <- paste(txt, "         <td  bgcolor=\"#f0f0ff\"><i>", Info.list$Description, "</i></td> \n", sep = "")
  txt <- paste(txt, "	      <br> \n", sep = "")
  txt <- paste(txt, "	  </tr> \n", sep = "")
  txt <- paste(txt, "	  <tr height=\"35\"><td align=\"right\" bgcolor=\"#e0e0f0\"><b>Analysts</b></td> \n", sep = "")
  txt <- paste(txt, "	      <td  bgcolor=\"#f0f0ff\"><i>", Info.list$Analysts, "</i></td> \n", sep = "")
  txt <- paste(txt, "	      <br> \n", sep = "")
  txt <- paste(txt, "	  </tr> \n", sep = "")
  txt <- paste(txt, "	  <tr height=\"35\"><td align=\"right\" bgcolor=\"#e0e0f0\"><b>Contact</b></td> \n", sep = "")
  txt <- paste(txt, "	      <td bgcolor=\"#f0f0ff\"><i>", Info.list$Contact, "</i></td> \n", sep = "")
  txt <- paste(txt, "	      <br> \n", sep = "")
  txt <- paste(txt, "	  </tr> \n", sep = "")
  txt <- paste(txt, "</table> \n", sep = "")
  txt <- paste(txt, "<hr> \n", sep = "")
  
  cat(txt, file = outfile, sep = "")

  close(outfile)
}  


### printVHIRfooter: Inserta la imatge amb logo del VHIR al peu de la pagina de resultats
###                 
###
### Parametres:
###
###   FileName : Nom del fitxer de resultats
###
### OBERVACIONS IMPORTANTS:
###   Si UEB es FALSE no hauria d'entrar a fer aixo
###   

printVHIRfooter <- function(FileName = "ResultFiles")
{
  outfile <- file(paste(FileName, ".html", sep = ""), open = "at")

  txt <- "<table width=\"100%\"  border=\"0\">\n"
  txt <- paste(txt, "   <tr height=\"50\"><td width=\"100%\"> </td></tr> \n", sep="")
  txt <- paste(txt, "   <tr><td width=\"100%\"><div align=\"center\"> \n", sep="")
  txt <- paste(txt, "      <a href=\"http://www.vhir.org\" target=\"z\"><img width=\"90%\" src=\"images/imatgelogotip.jpg\" border=\"0\"></a> \n", sep="")
  txt <- paste(txt, "   </div></td></tr> \n", sep="")

  txt <- paste(txt, "</table> \n", sep= "")

  cat(txt, file = outfile, sep = "")
  
  close(outfile)
}


### write.ResultFiles: 
###                 
###
### Parametres:
###
###    FileName :
###   Info.list :
###
### OBERVACIONS IMPORTANTS:
###
###

write.section <- function(my.info,
                          filename,
                          my.id = "",
                          sectionTitle = "Files",
                          IndexDir = "") 
{
  outfile <- file(paste(filename, ".html", sep = ""), open = "at")


  txt <- paste("<h1 id=\"", my.id, "\">", sectionTitle, "</h1> \n", sep = "")

  txt <- paste(txt, "<table border=\"0\"> \n", sep = "")

  txt <- paste(txt, "    <tr>" , sep = "")
  for (j in 1:ncol(my.info))
  {
    txt <- paste(txt, "<td bgcolor=\"", c("#cc99cc", "#cc77cc")[j%%2 + 1], "\"><b>", colnames(my.info)[j], "</b></td>\n", sep = "")
  }
  txt <- paste(txt, "    </tr> \n", sep = "")

  for (i in 1:nrow(my.info))
  {
    txt <- paste(txt, "    <tr>" , sep = "")
    for (j in 1:ncol(my.info))
    {
      my.cell <- my.info[i, j]
      
      if(j == 1)
      {
        my.cell <- paste("<a href=\"", paste(IndexDir, my.info[i, j], sep=""),"\" target=\"z\">", my.cell, "</a>", sep = "")
      }
            
      txt <- paste(txt,
                   "<td bgcolor=\"", c("#e0e0ff", "#d0d0f0", "#f0f0ff", "#e0e0f0")[i%%2 * 2 + j%%2 + 1], "\">",
                   my.cell,
                   "</td> \n",
                   sep = "")
    }
    txt <- paste(txt, "    </tr> \n" , sep = "")
  }

  txt <- paste(txt, "</table> \n", sep = "")

  cat(txt, file = outfile, sep = "")

  close(outfile)
}


### printAnalysisDetails: 
###                 
###
### Parametres:
###
###    FileName :
###
### OBERVACIONS IMPORTANTS:
###
###

closeHtml <- function(FileName = "ResultFiles")
{
  outfile <- file(paste(FileName, ".html", sep = ""), open = "at")
  
  cat("</body>", "</html>", file = outfile, sep = "\n")

  close(outfile)
}



### LinksFile2Html: 
###                 
###
### Parametres:
###
###       lFile :
###   outputDir :
###   info.list :
###    IndexDir : per defecte "ResultFiles/"
###         UEB : Mostra la cap??alera de la UEB (per defecte TRUE). En cas contrari, la de EstBioinfo
###
### OBERVACIONS IMPORTANTS:
###
###    COMPTE!!! No esta implementat el control per a tenir present les subcategories del fitxer "Links.XXXnnn.txt".
###              Per ara nomes te present el nom del fitxer (FileName), la categoria principal (Category) i
###              la descripcio del fitxer (Description). Faltaria controla la darrea columna (Subcategory)

LinksFile2Html <- function(lFile,
                           outputDir,
                           info.list,
                           IndexDir = "",
                           UEB = TRUE)
{
  FixLinksFile(lFile)
  
  my.names <- c("Result Summary",
                "Data",
                "Quality Control",
                "Pre-processing",
                "Analysis",
                "Annotations",
                "Multiple Comparisons",
                "Cluster",
                "GO Analysis",
                "KEGG Analysis",
                "Ingenuity Pathways Analysis",
                "Power Analysis")

  names(my.names) <- c("INFO", "DATA", "QC", "NORM", "ANALYSIS", "ANNOT", "MULTCOMP", "CLUSTER", "GO", "KEGG", "IPA", "POWER")
  
  my.LnkFile <- read.table(lFile, header = TRUE, sep = "\t", stringsAsFactor = FALSE)
  my.cats <- unique(my.LnkFile[, 2])
  
  resFile <- file.path(outputDir, "ResultFiles")
  
  printHeader(resFile)
  printGroupHeader(resFile, UEB = UEB)
  printAnalysisDetails(resFile, info.list)
  
  i <- 0 # Per generar el numero de categoria en el fitxer resultFiles.html de forma dinamica
  for(categ in names(my.names)[names(my.names) %in% my.cats])
  {
    i <- i+1
    my.ind <- which(my.LnkFile[, 2] == categ)
    write.section(my.info = my.LnkFile[my.ind, c(1,4)],  
                  filename = resFile,
                  my.id = categ,                  
                  sectionTitle = paste(i, ". ", my.names[categ], " Files", sep = ""),
                  IndexDir = IndexDir)
  }

  if (UEB) printVHIRfooter(resFile)

  closeHtml(resFile)
}



### trim: Elimina els espais en blanc a esquerra i/o dreta d'un vector de caracters
###
### Parametres:
###
###   x : vector de caracters
###
### Exemple:
###
###        > trim(c("  a  ", "b"))
###        [1] "a" "b"
###
###        > trim(c("  a  ", "b   "))
###        [1] "a" "b"


trim <- function(x)
{
  gsub("^[[:space:]]+|[[:space:]]+$", "", x)
}


### repofun: Donat un vector d'identificadors, el tipus d'anotacio pel que es vol crear el link
###          i l'especie si s'escau, dona una llista amb els hyperlinks a cada identificador
###          (pensant en una sortida HTML) o una concatenacio dels identificadors  (pensant en
###          una sortida txt)
###
### Parametres:
###
###      x : Vector amb els identificadors
###      y :
###      z : Nom de cientific de l'especie. Aquest parametre es necessari per a crear els links
###          en els repositoris de la KEGG i de ENSEMBL. En l'actualitat nomes es contempla les
###          seguents especies per a la KEGG: "Homo_sapiens", "Mus_musculus", "Rattus_norvegicus",
###          "Bos_taurus", "Danio_rerio", "Sus_scrofa". El parametre per defecte esta a NULL
###    ... :
###
### Exemple:
###
###

repofun <- function(x, y, z = NULL, ...)
{
  NAids <- which(is.na(x))
  blnksIDs <- which(x=="&nbsp;")
  x <- trim(x)
  
  specie <- z
  #if(!is.null(z))
  #{
  #  specie <- switch(z,
  #                   "Homo sapiens" = "hsa",
  #                   "Mus musculus" = "mmu",
  #                   "Rattus norvegicus" = "rno",
  #                   "Bos taurus" = "bta",
  #                   "Danio rerio" = "dre",
  #                   "Sus scrofa" = "ssc")
  #}
  
  if(y %in% c("SYMBOL", "GENENAME", "MAP", "ENZYME"))
  {
    out <- x
  }else{
    if(y=="PMID")
    {
      out <- paste("<a href=\"http://www.ncbi.nlm.nih.gov/pubmed/",
                   paste(x, collapse = ","), "?tool=bioconductor\" target=\"_blank\">", length(x), "</a>", sep = "")
    }else{
      out <- switch(y,
                    "ACCNUM" = paste("<a href=\"http://www.ncbi.nlm.nih.gov/protein/", x, "\" target=\"_blank\">", x, "</a>", sep = ""),
                  #  "ACCNUM" = paste("<a href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=search&db=nucleotide&term=",
                  #                   x, "[ACCN]&doptcmdl=GenBank\" target=\"_blank\">", x, "</a>", sep = ""),
                    "ENTREZ" = paste("<a href=\"http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&cmd=Retrieve&dopt=Graphics&list_uids=",
                                     x, "\" target=\"_blank\">", x, "</a>", sep = ""),
                    "REFSEQ" = paste("<a href=\"http://www.ncbi.nlm.nih.gov/nuccore/", x, "?\" target=\"_blank\">", x, "</a>", sep = ""),
                    "UNIGENE" = paste("<a href=\"http://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Hs&CID=",
                                     substr(x, 4, nchar(x)), "\" target=\"_blank\">", x, "</a>", sep = ""),
                    "ENSEMBL" = paste("<a href=\"http://www.ensembl.org/", z, "/Gene/Summary?g=", x, "\" target=\"_blank\">", x, "</a>", sep = ""),
                    "UNIPROT" = paste("<a href=\"http://www.uniprot.org/uniprot/", x, "\" target=\"_blank\">", x, "</a>", sep = ""),
                    "GO" =  paste("<a href=\"http://amigo.geneontology.org/cgi-bin/amigo/term_details?term=", x, "\" target=\"_blank\">", x, "</a>", sep = ""),
                    "PATH" = paste("<a href=\"http://www.genome.jp/dbget-bin/www_bget?", specie, x, "\" target=\"_blank\">", x, "</a>", sep = ""))
    }
  }
  
  out[NAids] <- "&nbsp;"
  out[blnksIDs] <- "&nbsp;"
  
  return(out)
}


### translateIDs: Donada un vector de de identificadors Entrez Gene torna una llista amb els identidicadors
###               o hyperlinks a un repository especific
###
###
### Parametres:
###
###        geneIDs : Vector d'identificadors d'Entrez Gene
###    anotPackage : Paquet d'anotacions (sense l'extensio .db)
###         typeID : Tipus d'anotacio que es vol extreure. En l'actualitat es pot extreure la seguent informacio
###                                   "SYMBOL" : Gene Symbol
###                                 "GENENAME" : Description
###                                   "ACCNUM" : GeneBank
###                                   "REFSEQ" : RefSeq
###                                       "GO" : Gene Ontology
###                                     "PATH" : KEGG
###                                  "UNIGENE" : UnigGene
###                                  "ENSEMBL" : ensembl!
###                                  "UNIPROT" : UniProt
###                                      "MAP" : Cytoband
###                                   "ENZYME" : Enzyme
###                                     "PMID" : Pubmed
###                  Per defecte pren "REFSEQ"
###         toHTML : Si es TRUE genera els hyperlinks al repositori sol??licitat. Aquesta es la opcio per
###                  defecte. Si es FALSE genera la concatenacio dels identificadors, seprant-los per " /// "
###         specie : Nom de cientific de l'especie. Aquest parametre es necessari per a crear els links
###                  en els repositoris de la KEGG i de ENSEMBL. En l'actualitat nomes es contempla les
###                  seguents especies per a la KEGG: "Homo_sapiens", "Mus_musculus", "Rattus_norvegicus",
###                  "Bos_taurus", "Danio_rerio", "Sus_scrofa". El parametre per defecte esta a NULL
###            ... :
###
### Exemple:
###
###

translateIDs <- function(geneIDs,
                         anotPackage,
                         typeID = "REFSEQ",
                         toHTML = TRUE,
                         specie = NULL,
                         ...)
{
  annot <- eval(parse(text = paste(anotPackage, typeID, sep = "")))
  newIDs <- mget(geneIDs, annot, ifnotfound = NA)

  if(typeID=="GO")
  {
   if (runMulticore ==1 || runMulticore ==3) { 
     goids <- mclapply(1:length(names(newIDs)), function(x, y){names(newIDs[[x]])}, y = newIDs)
    } else {
     goids <- lapply(1:length(names(newIDs)), function(x, y){names(newIDs[[x]])}, y = newIDs)         
   }
    names(goids) <- names(newIDs)
   if (runMulticore ==1 || runMulticore ==3) { 
     idxs <- which(unlist(mclapply(goids, is.null))==TRUE)
    } else {
     idxs <- which(unlist(lapply(goids, is.null))==TRUE)         
    }
    goids[idxs] <- NA
    newIDs <- goids
  }
  
  if(toHTML)
  {
   if (runMulticore ==1 || runMulticore ==3) { 
      newIDs <- mclapply(newIDs, function(x, y, z) {repofun(x, y, z)}, y = typeID, z = specie)
     } else {
      newIDs <- lapply(newIDs, function(x, y, z) {repofun(x, y, z)}, y = typeID, z = specie)         
   }

    mysep <- ", "
  }else{
    mysep <- " /// "  
  }
  
   if (runMulticore ==1 || runMulticore ==3) { 
       out <- mclapply(newIDs, FUN = function(x, y){paste(x, collapse = mysep)}, y = mysep)
     } else {
       out <- lapply(newIDs, FUN = function(x, y){paste(x, collapse = mysep)}, y = mysep)
   }

  return(out)
}



### GeneAnnotation : Donada un vector de de identificadors Entrez Gene torna una llista amb els identidicadors
###                  o hyperlinks a un repository especific
###
###
### Parametres:
###
###          egIDs : Vector d'identificadors d'Entrez Gene
###    anotPackage : Paquet d'anotacions (sense l'extensio .db)
###         toHTML : Si es TRUE genera un arxiu html amb les anotacions enlla??ades als diferents
###                  repositoris. Aquesta es la opcio per  defecte. Si es FALSE genera un arxiu de
###                  text, concatenant els sinonims amb " /// "
###
###       filename : Nom de l'arxiu de sortida sense extensio. Si toHTML = TRUE l'extensio que prendra
###                  sera .html, i si es FALSE l'extensio de sortida sera .txt
###        myTitle : Titiol de la pagina en cas de que toHTML = TRUE
###         specie :
###      info2show : Tipus d'anotacio que es vol extreure. En l'actualitat es pot extreure la seguent informacio
###
###                         "Affymetrix", "EntrezGene", "GenBank","KEGG", "GO", "RefSeq", "UniGene",
###                         "Ensembl", "UniProt", "PubMed", GeneSymbol", "GeneName", "Cytoband", "Enzyme"
###
###                  "Affymetrix", "EntrezGene" sempre es mostren i per defecte tambe es mostren els seguents
###                  repositoris 
###
###                         "Affymetrix", "EntrezGene", GeneSymbol", "GeneName", "KEGG", "GO", "PubMed"
### 
###      linksFile :
###       maxGenes : Numero maxim de gens que es mostraran per arxiu HTML d'anotacions. En cas de superar aquest
###                  nombre es crea una nova pagina amb els seguents (i.e. si maxGenes = 100 i tenim 1000 gens a
###                  anotar creara 10 pagines). Per defecte pren maxGenes = NULL i mostra tots els gens
###            ... :

GeneAnnotation <- function(egIDs,
                           anotPackage,
                           toHTML = TRUE,
                           outputDir,
                           filename = "annotations",
                           myTitle = "Annotations for all genes analyzed",
                           specie = "Homo_sapiens",
                           info2show =  c("Affymetrix", "EntrezGene"),
                           linksFile,
                           maxGenes = NULL,
                           ...)
{
  ptm <- proc.time()
  
  NAs <- egIDs[is.na(egIDs)]
  eg <- na.omit(egIDs)

  if(toHTML)
  {
    l.ai <- paste("<span id=\"", names(eg), "\">", names(eg), "</span>", sep="")
  }else{
    l.ai <- names(eg)
  }
  l.eg <- repofun(x = eg, y = "ENTREZ")
  
  a <- cbind(l.ai, l.eg)
  colnames(a) <- c("AffyID", "EntrezGene")
 
  if("GeneSymbol" %in% info2show)
  {
    colnames.a <- colnames(a)
    
    l.gs <- unlist(translateIDs(geneIDs = eg, anotPackage, typeID = "SYMBOL", toHTML))

    a  <- cbind(a, l.gs)
    colnames(a) <- c(colnames.a, "Symbol")
  }
  if("GeneName"  %in% info2show)
  {
    colnames.a <- colnames(a)
    
    l.gn <- unlist(translateIDs(geneIDs = eg, anotPackage, typeID = "GENENAME", toHTML))

    a  <- cbind(a, l.gn)  
    colnames(a) <- c(colnames.a, "Description")
  }
  if("KEGG" %in% info2show)
  {
    colnames.a <- colnames(a)
    
    l.kg <- unlist(translateIDs(geneIDs = eg, anotPackage, typeID = "PATH", toHTML, specie = specie))

    a  <- cbind(a, l.kg)
    colnames(a) <- c(colnames.a, "KEGG")
   }
  if ("GO" %in% info2show)
  {
    colnames.a <- colnames(a)
    
    l.go <- unlist(translateIDs(geneIDs = eg, anotPackage, typeID = "GO", toHTML))
    a  <- cbind(a, l.go)
    colnames(a) <- c(colnames.a, "GO")
  }
  if ("PubMed" %in% info2show)
  {
    colnames.a <- colnames(a)
    
    l.pm <- unlist(translateIDs(geneIDs = eg, anotPackage, typeID = "PMID", toHTML))

    a  <- cbind(a, l.pm)
    colnames(a) <- c(colnames.a, "Pubmed")
  }

  if("GenBank" %in% info2show)
  {
    colnames.a <- colnames(a)
    
    l.gb <- unlist(translateIDs(geneIDs = eg, anotPackage, typeID = "ACCNUM", toHTML))

    a  <- cbind(a, l.gb)
    colnames(a) <- c(colnames.a, "GenBank")
  }
  if("RefSeq" %in% info2show)
  {
    colnames.a <- colnames(a)
    
    l.rs <- unlist(translateIDs(geneIDs = eg, anotPackage, typeID = "REFSEQ", toHTML))

    a  <- cbind(a, l.rs)
    colnames(a) <- c(colnames.a, "RefSeq")
  }
  if ("UniGene" %in% info2show)
  {
    colnames.a <- colnames(a)
    
    l.ug <- unlist(translateIDs(geneIDs = eg, anotPackage, typeID = "UNIGENE", toHTML))

    a  <- cbind(a, l.ug)
    colnames(a) <- c(colnames.a, "UniGene")
  }
  if("Ensembl" %in% info2show)
  {
    colnames.a <- colnames(a)
    
    l.em <- translateIDs(geneIDs = eg, anotPackage, typeID = "ENSEMBL", toHTML, specie = specie)

    a  <- cbind(a, l.em)
    colnames(a) <- c(colnames.a, "Ensembl")
   }
  if("UniProt" %in% info2show)
  {
    colnames.a <- colnames(a)
    
    l.up <- unlist(translateIDs(geneIDs = eg, anotPackage, typeID = "UNIPROT", toHTML))

    a  <- cbind(a, l.up)
    colnames(a) <- c(colnames.a, "UniProt")
   }
  if("Cytoband" %in% info2show)
  {
    colnames.a <- colnames(a)
    
    l.cy <- unlist(translateIDs(geneIDs = eg, anotPackage, typeID = "MAP", toHTML))

    a  <- cbind(a, l.cy)
    colnames(a) <- c(colnames.a, "Cytoband")
   }
  if("Enzyme" %in% info2show)
  {
    colnames.a <- colnames(a)
    
    l.ez <- unlist(translateIDs(geneIDs = eg, anotPackage, typeID = "ENZYME", toHTML))

    a  <- cbind(a, l.ez)
    colnames(a) <- c(colnames.a, "Enzyme")
   }
  
  if(!is.null(maxGenes))
  {
    if (!toHTML) maxGenes <- nrow(a) # si volem anotacions en .txt enlloc de .html no cal tallar l'output
    k <- ceiling(nrow(a) / maxGenes) # tallarem l'output en k fitxers en funcio del maxim d'entrades per fitxer

    j <- 1 # contador dels "marges" per fer els subsets
    for(i in 1:k)
    {
      genAnnot <- a[j:ceiling((nrow(a)*(i/k))), ]
    
      outFileName <- paste(filename, paste(i, "of", k ,sep = ""), sep = ".")
    
      if(toHTML)
      {
        write.htmltable(x = genAnnot, file = file.path(outputDir, outFileName), title = myTitle)
        sortable.html.table(df = as.data.frame(genAnnot),
                            output.file = paste0(outFileName,"-sortable.html"),
                            output.directory = outputDir,
                            page.title = myTitle )

      }else{
        write.table(x = genAnnot, file = file.path(outputDir, paste(outFileName, "txt", sep = ".")),
                    row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
      }
    
      anotationsFName <- ifelse(toHTML, paste(outFileName, "html", sep="."), paste(outFileName, "txt", sep="."))
      addToLinksFile(linksFile = linksFile,
                     aFName = anotationsFName,
                     categ = 'ANNOT',
                     desc = "Gene annotations for all genes analyzed")
      j <- j + ceiling((nrow(a) / k)) # actualitzem els marges pel seguent subset
    }
  }else{
    genAnnot <- a
    
    outFileName <- filename
    
    if(toHTML)
    {
      write.htmltable(x = genAnnot, file = file.path(outputDir, outFileName), title = myTitle)
      sortable.html.table(df = as.data.frame(genAnnot),
                          output.file = paste0(outFileName,"-sortable.html"),
                          output.directory = outputDir,
                          page.title = myTitle )
    }else{
      write.table(x = genAnnot, file = file.path(outputDir, paste(outFileName, "txt", sep = ".")),
                  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    }
    
    anotationsFName <- ifelse(toHTML, paste(outFileName, "html", sep="."), paste(outFileName, "txt", sep="."))
    addToLinksFile(linksFile = linksFile,
                   aFName = anotationsFName,
                   categ = 'ANNOT',
                   desc = "Gene annotations for all genes analyzed")
  }

  return(proc.time() - ptm)
}


doGeneAnnotation <- function(AnotList)
{

  p <- AnotList

  if(!is.null(p$my.IDs))
  {
    EntrezIDs <-  eval(parse(text = p$my.IDs))
  }

  if (!is.null(p$fitFileName))
  {
    fitMain <- loadFromFile(file.path(p$outputDir, p$fitFileName))
  }else{
    if (!is.null(p$fitMain))
    {
      fitMain <- eval(parse(text = p$fitMain)) # Posar-hi un tryCatch per poder sortir si dona error!!!
    }else{
      stop("Error, cal subministrar un nom d'arxiu o d'objecte 'fitMain'")
    }
  }

 genes2annotate <- EntrezIDs[unique(rownames(fitMain$p.value))]
  
 genesAnnotated <- GeneAnnotation(egIDs = genes2annotate,
                                  anotPackage = orgPackage,
                                  toHTML = p$toHTML,
                                  outputDir = p$outputDir,
                                  filename = p$anotFilename,
                                  myTitle = p$titleAnotations,
                                  specie = organisme,
                                  info2show = p$info2show,
                                  linksFile = p$linksFile,
                                  maxGenes = p$numGenesPerPage
                                  )
  return(genesAnnotated)
}



### save.sessionInfo: Desa en un arxiu de text la informacio sobre la sessio de R amb la que s'ha processat
###                   l'analisi.
###                 
###
### Parametres:
###
###     fileName : Nom del fitxer de text a on es desara la informacio. Per defecte pren com a nom "sessionInfo.txt"
###    outputDir : Ruta a on es desara el document
###
### Exemple
###
###    > save.sessionInfo(fileName = "mySessionInfo.txt", outputDir = "/mydocs/myRproject/temp")
###

save.sessionInfo <- function(fileName = "sessionInfo.txt",
                             targetsdate,
                             beginningDate,
                             endingDate,
                             outputDir)
{
  sink(file = file.path(outputDir, fileName))

      cat("\n")
      cat("###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")  
      cat("### Analysis Dates\n")
      cat("###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")  
      cat("\n")
      cat("Creation Targets File date: "); print(targetsdate)
      cat("Beginning date: "); print(beginningDate)
      cat("Ending date: "); print(endingDate)
      cat("Last run date: "); print(date())
      cat("\n")
      cat("\n")
      cat("###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")  
      cat("### Summary Session Info\n")
      cat("###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")  
      cat("\n")
      print(sessionInfo())
      cat("\n")
      cat("\n")
      cat("###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")  
      cat("### Extended Session Info\n")
      cat("###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")  
      cat("\n")
      cat("###------------------------------------------------------------\n")  
      cat("### R version\n")
      cat("###------------------------------------------------------------\n")  
      cat("\n")
      print(R.version)
      cat("\n")
      cat("\n")
      cat("###------------------------------------------------------------\n")  
      cat("### Installed packages\n")
      cat("###------------------------------------------------------------\n")  
      cat("\n")
      print(installed.packages())
  
  sink()
}
