
################################################################################
####################
##### geneSelectable
####################
##### topTab: Toptable en format '.csv'
##### adj0: p.valor adjustat pel qual filtrar (de normal, 0.01)
##### adj1: segon p.valor adjustat pel qual filtrar (de normal, 0.05)
##### adj2: tercer p.valor adjustat pel qual filtrar (de normal, 0.25)
##### P1: p.valor pel qual filtrar (de normal, 0.01)
##### P2: segon p.valor pel qual filtrar (de normal, 0.05)
################################################################################


genesSelectable <- function (topTab, adj0, adj1, adj2, P1, P2,FC=1)
{
  upBelowB <- sum(topTab$B > 0  & topTab$t > 0 & abs(topTab$logFC) > FC)
  downBelowB <- sum(topTab$B > 0  & topTab$t < 0 & abs(topTab$logFC) > FC)
  
  upBelowAdj0 <- sum(topTab$adj.P.Val < adj0 & topTab$t > 0 & abs(topTab$logFC) > FC)
  downBelowAdj0 <- sum(topTab$adj.P.Val < adj0 & topTab$t < 0 & abs(topTab$logFC) > FC)
  
  upBelowAdj1 <- sum(topTab$adj.P.Val < adj1 & topTab$t > 0 & abs(topTab$logFC) > FC)
  downBelowAdj1 <- sum(topTab$adj.P.Val < adj1 & topTab$t < 0 & abs(topTab$logFC) > FC)
  
  upBelowAdj2 <- sum(topTab$adj.P.Val < adj2 & topTab$t > 0 & abs(topTab$logFC) > FC)
  downBelowAdj2 <- sum(topTab$adj.P.Val < adj2 & topTab$t < 0 & abs(topTab$logFC) > FC)
  
  upBelowP1 <- sum(topTab$P.Value < P1 & topTab$t > 0 & abs(topTab$logFC) > FC)
  downBelowP1 <- sum (topTab$P.Value < P1 & topTab$t < 0 & abs(topTab$logFC) > FC)
  
  upBelowP2 <- sum(topTab$P.Value < P2 & topTab$t > 0 & abs(topTab$logFC) > FC)
  
  table(topTab$P.Value < P2)
  sum(topTab$P.Value < P2)
  ?sum
  table(topTab$t > 0)
  table(abs(topTab$logFC) > FC)
  table(topTab$P.Value < P2 & topTab$t > 0)
  table(topTab$P.Value < P2 & topTab$t > 0 & abs(topTab$logFC) > FC)
  
  downBelowP2 <- sum(topTab$P.Value < P2 & topTab$t < 0 & abs(topTab$logFC) > FC)
  
  return(c(upBelowB = upBelowB,downBelowB = downBelowB,
           upBelowAdj0 = upBelowAdj0, downBelowAdj0 = downBelowAdj0,
           upBelowAdj1 = upBelowAdj1, downBelowAdj1 = downBelowAdj1,
           upBelowAdj2 = upBelowAdj2, downBelowAdj2 = downBelowAdj2,
           upBelowP1 = upBelowP1, downBelowP1 = downBelowP1,
           upBelowP2 = upBelowP2, downBelowP2 = downBelowP2))   
}



################################################################################
######################
##### numGeneChangedFC
######################
##### filenames: vector amb el noms dels fitxers topTable 
#####            (si estem dins de results, ex: filenames<-grep("Express",dir(),value=TRUE))
##### comparisons: vector amb els noms de les comparacions en el mateix ordre que filenames ( ex: c("NSCvsAstr","NSCvsNeur","NSCvsd6","d6vsd12"))
##### FC: valor del FoldChange pel qual es vol filtrar (si no es posa res, 0)
##### adj0: p.valor adjustat pel qual filtrar (si no es posa res, 0.01)
##### adj1: segon p.valor adjustat pel qual filtrar (si no es posa res, 0.05)
##### adj1: tercer p.valor adjustat pel qual filtrar (si no es posa res, 0.25)
##### P1: p.valor pel qual filtrar (si no es posa res, 0.01)
##### P1: segon p.valor pel qual filtrar (si no es posa res, 0.05)
################################################################################
numGeneChangedFC<- function (filenames,comparisons, FC=0, adj0=0.01, adj1=0.05, adj2=0.25, P1=0.01, P2=0.05)
{
  dat<-lapply(filenames, read.csv2,header=T)
  res<-data.frame(lapply(dat, genesSelectable,adj0=adj0, adj1=adj1, adj2=adj2, P1=P1, P2=P2,FC=FC))
  colnames(res) <- comparisons
  row.names(res) <- c("upReg-B>0", "downReg-B>0",
                      paste("upReg-Adjusted-p-val", adj0, sep = " < "), 
                      paste("downReg-Adjusted-p-val", adj0, sep = " < "),      
                      paste("upReg-Adjusted-p-val", adj1, sep = " < "), 
                      paste("downReg-Adjusted-p-val", adj1, sep = " < "),
                      paste("upReg-Adjusted-p-val", adj2, sep = " < "),
                      paste("downReg-Adjusted-p-val", adj2, sep = " < "),
                      paste("upReg-P value", P1, sep = " < "),
                      paste("downReg-P value", P1, sep = " < "),
                      paste("upReg-P value", P2, sep = " < "),
                      paste("downReg-P value", P2, sep = " < "))
  write.csv2( res, file=paste("numGenesChangedFC",FC,".csv",sep="") )
}

### EXEMPLES
#numGeneChangedFC(filenames=grep("Express",dir(),value=TRUE),comparisons= c("d6vsd12","NSCvsAstr","NSCvsd6","NSCvsNeur"),FC=0)
#numGeneChangedFC(filenames=grep("Express",dir(),value=TRUE),comparisons= c("d6vsd12","NSCvsAstr","NSCvsd6","NSCvsNeur"),FC=1)
#numGeneChangedFC(filenames=grep("Express",dir(),value=TRUE),comparisons= c("d6vsd12","NSCvsAstr","NSCvsd6","NSCvsNeur"),FC=2)






