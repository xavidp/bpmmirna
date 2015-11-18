###################################################
### chunk number 1: processa
###################################################
#getwd()
setwd("/home/xavi/Estudis/2015-10-NuriaBarbarroja-IMIBIC-A279")
## Per anar b??, s'hauria d'especificar el working dir des d'aqu??, per poder llen??ar tot un estudi executant un sol script d'R que ho control??s tot.
# setwd("/home/ueb/Estudis/2011-05-31-foo-XXX")

## Si ??s la primera vegada que vols executar el BasicPipe, cal que posis aquest par??metre seg??ent a TRUE
TRY2INSTALL <- FALSE
require(beepr)
source("parametros.BRB279.R") # Canviar XXXnnn per les 3 lletres de l'investigador/a + l'Id numeric de l'estudi

#line 13 "ComCrearUnProjecteProcessa.rnw"
Sys.setlocale('LC_ALL','C') # Per evitar warnings deguts a l'idioma

runParam <- FALSE # From here until the next definition of runParam, sections will be set to this initial value (usally TRUE)

  createAnnotations <- runParam   # Si no les crea les ha de carregar

runParam <- TRUE # 

  readCELS <- runParam            # runParam => Llegir els fitxers .CEL
  groupColors <- FALSE #runParam         # runParam => Usar el vector de colors definits al fitxer 'targets.txt'

  preProcess <- runParam          # Decideix si es fa el QC amb arrayQualityMetrics

  doNorm <- runParam              # FALSE => Llegir el fitxer .RDa amb les dades normalitzades. Altrament, les normalitza.
    normPlots <- runParam 
    toPDF.Norm <- runParam        # NormPlots

  FILTRAR <- runParam
  processaFILTRE <- runParam      # runParam => Executa el filtrat. Altrament, llegeix el fitxer .Rda que cont?? les dades filtrades
    doReport <- runParam          # runParam => Genra el fitxer report del filtratge

  processaLM <- runParam          # Aquest proc??s sol trigar molt.
  annotateGenes <- FALSE #runParam       # S'ha fet a part amb PreparaDades.A279.R"
  processaMultComp <- runParam 
  processaCluster <- runParam     # Si nom??s hi ha una comparaci?? (e.g. ThermalStress: estr??s contra sense estr??s), aqu?? s'ha de posar a FALSE, i per tant, no far?? els diagrames de VENN.

  toTIFF <- FALSE		  # TRUE => en lloc de pdf genera un fitxer TIFF amb mes qualitat per a publicacions del Volcano, Heatmap i Venn
###--------------------------------
runParam <- TRUE  # From here onwards, sections will be set to this other value (usally FALSE)
###--------------------------------
  processaGO <- FALSE #runParam          # Alternatively, it can be hardcoded to FALSE or TRUE
  processaKEGG <- FALSE #runParam        # Alternatively, it can be hardcoded to FALSE or TRUE
  processaPower <- FALSE	  # TRUE, if power analysis is appropiate
  processaReport <- runParam      # Alternatively, it can be hardcoded to FALSE or TRUE
  processaSessionInfo <- runParam # Alternatively, it can be hardcoded to FALSE or TRUE

# Param to indicate whether the run will be held on a multicore computer and willing to parallelize the run
# Options are: 	0 : no multicore at all
#		1 : all lapply into their parallel equivalents: mclapply
#		2 : all for loops into into their parallel equivalents: foreach
#		3 : all lapply and for loops into into their parallel equivalents: mcapply and foreach
runMulticore <- 0 # 0=none | 1=mclapply | 2=foreach | 3=mcapply & foreach

###################################################
### chunk number 2: execuTar eval=FALSE
###################################################
## #line 44 "ComCrearUnProjecteProcessa.rnw"

#source("analisis.BRB279.R") # Canviar XXXnnn per les 3 lletres de l'investigador/a + l'Id numeric de l'estudi
#beep(sound = 8, expr = NULL)

