library(VennDiagram)
basepath <- "/home/xavi/Estudis/2015-10-NuriaBarbarroja-IMIBIC-A279/"
setwd(paste0(basepath, "results"))

top1 <- read.table(file = "./results.FF.WO_19_33/Selected.Genes.in.comparison.ICvsCL.for.FF.csv",
                      header = TRUE, sep = ",", dec=".")
top2 <- read.table(file = "./results.PFA.WO_19_32/Selected.Genes.in comparison.ICvsCL.for.PFA.csv",
                   header = TRUE, sep = ",", dec=".")
top3 <- read.table(file = "./results.FFPE.WO_19_32/Selected.Genes.in.comparison.ICvsCL.for.FFPE.csv",
                   header = TRUE, sep = ",", dec=".")


dim(top1)
dim(top2)
dim(top3)
head(top1)
head(top2)
head(top3)

pass1 <- top1[which(top1$P.Value<0.01),]
tail(pass1)
pass1 <- pass1[which(pass1$logFC > 1 | pass1$logFC < -1),]
tail(pass1)

pass2 <- top2[which(top2$P.Value<0.01),]
tail(pass2)
pass2 <- pass2[which(pass2$logFC > 1 | pass2$logFC < -1),]
tail(pass2)

pass3 <- top3[which(top3$P.Value<0.01),]
tail(pass3)
pass3 <- pass3[which(pass3$logFC > 1 | pass3$logFC < -1),]
tail(pass3)

list1 <- as.character(sort(unique(pass1$X)))
length(list1)
list2 <- as.character(sort(unique(pass2$X)))
length(list2)
list3 <- as.character(sort(unique(pass3$X)))
length(list3)

list <- c(list1, list2, list3)
length(list)
list <- sort(unique(list))
length(list)

df <- data.frame(micros = list, l1 = rep(0,length(list)), l2 = rep(0,length(list)), l3 = rep(0,length(list)))
head(df)

df$l1 <- as.numeric((df$micros %in% list1)==T)
df$l2 <- as.numeric((df$micros %in% list2)==T)
df$l3 <- as.numeric((df$micros %in% list3)==T)


vennDiagram(vennCounts(df[,2:4]), names=c("ICvsCL_FF","ICvsCL_PFA","ICvsCL_FFPE"),
            circle.col=c("red", "green", "yellow"))
setwd(resultsDir)
pdf("VennDiagram.pdf")
vennDiagram(vennCounts(df[,2:4]), names=c("ICvsCL_FF","ICvsCL_PFA","ICvsCL_FFPE"),
            circle.col=c("red", "green", "yellow"))
dev.off()

)
