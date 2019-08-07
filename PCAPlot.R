MakePCA <- function(cts, coldata) {

library(DESeq2)
library(apeglm)
library(ggplot2)

#Read Files

#cts <-
  #read.table(
    #"flagCounts.txt",
    #sep = "\t",
    #header = TRUE,
    #row.names = 1,
    #check.names = FALSE
  #)

#cts <- cts[10:15]
#coldata <-
  #read.table(
    #"flagKeyFile.csv",
    #sep = ",",
    #header = TRUE,
    #row.names = 1
 # )


#No files are imported from below this point
#-------------------------------------------------------------------------------------------------------------

cts <- cts[,2:ncol(cts)]
#Construct DESseqDataObject
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Group)



#Filter low Counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Perform PCA Analysis
vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = "Group")
PCAdata <- plotPCA(vsd, intgroup = "Group", returnData = TRUE)
percentVar <- round(100 * attr(PCAdata, "percentVar"))
pngName <- "PCAPlot.png"
png(pngName)
returnPlot <- ggplot(PCAdata, aes(PC1, PC2, color = Group)) + 
  geom_point(size = 3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  coord_fixed(ylim = c(-20, 20) )

finalReturnPlot <- print(returnPlot)
dev.off()

return(pngName)

}