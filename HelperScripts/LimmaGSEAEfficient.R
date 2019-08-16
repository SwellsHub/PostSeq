

DifferentialExpression <- function(species, annotation, filter, target, minCutoff, makeVolc, mysqlUpload) {

#Import Libraries
library(edgeR)
library(limma)
library(stringi)
library(biomaRt)
library(DESeq2)
library(xml2)
library(ggplot2)
library(png)
library(data.table)


  
#No files being imported below this point
#------------------------------------------------------------------------------------------------------------


#species <- as.character(target[1, 4])
#annotation <-  as.character(target[1, 5])

#Filter out repeated values
rownames(filter) <- make.names(filter$X, unique = TRUE)

filter <- filter[2:ncol(filter)]
filter <- data.frame(filter)
filter <- filter[!(rownames(filter) == ""),]
filter <- na.omit(filter)
filter <- filter[order(rownames(filter)), , drop = FALSE]


if(stri_cmp_eq(annotation, "hgnc_symbol") || stri_cmp_eq(annotation, "mgi_symbol")) {
  
      finalsheet <- data.frame(filter, row.names = rownames(filter))
      

} else {
  #For Ensembl, trim version suffixes from IDs
  genes <- as.matrix(rownames(filter))
  colnames(genes) = annotation
  genes <- unique(sub("\\..*", "", genes))
  filter <-
    filter[unique(sub("\\..*", "", rownames(filter))), , drop = FALSE]
  rownames(filter) <- sub("\\..*", "", genes)
  filter <- filter[complete.cases(filter),]
  
  
  

ensembl <- useEnsembl(biomart = ("ensembl"), dataset = species, mirror = 'useast')
if (stri_cmp_eq(species, "hsapiens_gene_ensembl")) {
  tryCatch(
Conversion <-
  getBM(
    attributes = c(as.character(annotation), 'hgnc_symbol'),
    filters = as.character(annotation),
    values = rownames(filter),
    mart = ensembl
  ), error = function (err) {
    shinyalert("Bad News", "The Biomart servers may be down right now! If you are sure your inputs are correct, wait 30 minutes then try again. If this still fails, contact spencerwells@gmail.com", type = "error")
    return("error")
  }
  )
Conversion <-
  Conversion[order(Conversion[annotation]), , drop = FALSE]
filter <- data.frame(filter, row.names = rownames(filter))
filter$annotation <- rownames(filter)
names(filter)[ncol(filter)] <- annotation
finalsheet <- merge(filter, Conversion, by = annotation)
finalsheet <- finalsheet[!(finalsheet$hgnc_symbol == ""), ]
finalsheet <- finalsheet[complete.cases(finalsheet),]
finalsheet <-
  finalsheet[!duplicated(finalsheet$hgnc_symbol), , drop = FALSE]
rownames(finalsheet) <- finalsheet$hgnc_symbol
finalsheet <- finalsheet[2:(ncol(finalsheet) - 1)]

} else if(stri_cmp_eq(species, "mmusculus_gene_ensembl")) {
  tryCatch(
  Conversion <-
    getBM(
      attributes = c(as.character(annotation), 'mgi_symbol'),
      filters = as.character(annotation),
      values = rownames(filter),
      mart = ensembl
    ), error = function (err) {
      shinyalert("Bad News", "The Biomart servers may be down right now! If you are sure your inputs are correct, wait 30 minutes then try again. If this still fails, contact spencerwells@gmail.com", type = "error")
      return("error")
    }
    )
  Conversion <-
    Conversion[order(Conversion[annotation]), , drop = FALSE]
  filter <- data.frame(filter, row.names = rownames(filter))
  filter$annotation <- rownames(filter)
  names(filter)[ncol(filter)] <- annotation
  finalsheet <- merge(filter, Conversion, by = annotation)
  finalsheet <- finalsheet[!(finalsheet$mgi_symbol == ""), ]
  finalsheet <- finalsheet[complete.cases(finalsheet),]
  finalsheet <-
    finalsheet[!duplicated(finalsheet$mgi_symbol), , drop = FALSE]
  rownames(finalsheet) <- finalsheet$mgi_symbol
  finalsheet <- finalsheet[2:(ncol(finalsheet) - 1)]
}
} 


#Load count and phenotype label files to tables
x <- finalsheet

#Make matrix of labels for design of experiment and intended pairwise comparisons

#Factor group column of target file to construct design matrix
target <- data.frame(target)
design <- model.matrix( ~ target$Group + 0)

#Extract the relevant contrasts of the experiment from the third column of the target file
contrasts <- na.omit(contrasts)
contrasts <- as.character(unlist(target[2]))
contrasts <- noquote(contrasts)
contrasts <- contrasts[contrasts != ""]
for(y in 1:length(contrasts)) {
  contrasts[y] <- sub("^([^\r\n-]+)(-*)([^-\r\n]+)$", "\\3\\2\\1", contrasts[y])
}

#Arrange columns of design matrix in alphabetical order
target$Group <- factor(target$Group)
colnames(design) <- sort(as.character(levels(target$Group)))
contrasts <- na.omit(contrasts)


tryCatch(contrast.matrix <-
  makeContrasts(contrasts = contrasts,  levels = design), error = function(err) {
    shinyalert("Bad News", "Your experimental design is malformed. Check if there are any typos in your groups and/or comparisons, then make sure you are using syntactically valid names (cannot start with a number, no special characters besides . and _)", type = "error")
    shinyjs::toggle("downloadMessage")
    return("error")
  })






#Make DGEList from tables
#Use levels from group column of target file to feed design to DGEList
dgeFull <- DGEList(counts = x, target$Group)



#Filter out low count genes
dgeFull <-
  DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) >= minCutoff,],
          group = dgeFull$samples$group)
normCounts <- cpm(dgeFull)
eff.lib.size <-
  dgeFull$samples$lib.size * dgeFull$samples$norm.factors


#Estimate Normalization factors
dgeFull <- calcNormFactors(dgeFull, "TMM")


#The Fourth parameter in the second column of the target file (first column is rownames and does not count as column) 
#establishes whether logCPM or Voom normalization is applied to data. A 0 is logCPM, a 1 is Voom
#if (target[4, 2] == 0) {
  #Get Log counts per million
  logCPM <- cpm(dgeFull, log = TRUE, prior.count = 3)
#}
#Or apply Voom Normalization

  logCPM <- voomWithQualityWeights(dgeFull, design, plot = TRUE, normalization = "none")

#Apply and eBayes Linear trend
fit <- lmFit(logCPM, design)

##Apply to pairwise analysis
fit2 <- contrasts.fit(fit, contrast.matrix)




#Standard Method of fit calculation
# The next cell in the parameters column of the key file determines wheter eBayes or treat is applied to the linear model.
#A 0 is eBayes, a 1 is treat (treat incorporates LogFC into the calculation of adj.p.Val).
#if (target[5, 2] == 0) {
 # fit2 <- eBayes(fit2, trend = TRUE)
#}

#Method for increased importance of logFC
  fit2 <- treat(fit2, lfc = .5)


#Loop through writing for every pairwise comparison
returnList <- c()
volcReturns <- c()
for (i in seq(from = 1, to = ncol(contrast.matrix))) {
  #Write table
  name <- colnames(contrast.matrix)[i]
  
  #Standard Output
  #The cell that determines if the treat normalization is used also determines what the naming convention of the output file is.
  
  #if (target[5, 2] == 0) {
   # returnTable <- write.table(topTable(fit2, coef = i, number = 50000),
               # paste0(name, "liberal.csv"),
               # sep = ",")
    #returnList <- c(paste0(name, "liberal.csv"), returnList)
 # }
  
  #Generate Locations for each gene
  ensembl <- useEnsembl(biomart = ("ensembl"), dataset = species, mirror = 'uswest')
  
  if(species == "hsapiens_gene_ensembl") {
    tryCatch(
  LocationConversion <-
    getBM(
      attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position', 'strand'),
      filters = "hgnc_symbol",
      values = rownames(logCPM$E),
      mart = ensembl
    ), error = function (err) {
      shinyalert("Bad News", "The Biomart servers may be down right now! If you are sure your inputs are correct, wait 30 minutes then try again. If this still fails, contact spencerwells@gmail.com", type = "error")
      return("error")
    }
    )
  LocationConversion <-
    LocationConversion[order(LocationConversion['hgnc_symbol']), , drop = FALSE]
  } else {
    tryCatch(
    LocationConversion <-
      getBM(
        attributes = c('mgi_symbol', 'chromosome_name', 'start_position', 'end_position', 'strand'),
        filters = "mgi_symbol",
        values = rownames(logCPM$E),
        mart = ensembl
      ), error = function (err) {
        shinyalert("Bad News", "The Biomart servers may be down right now! If you are sure your inputs are correct, wait 30 minutes then try again. If this still fails, contact spencerwells@gmail.com", type = "error")
        return("error")
      }
    )
    LocationConversion <-
      LocationConversion[order(LocationConversion['mgi_symbol']), , drop = FALSE]
  }
  
  #Output assigning more important to LogFC
  if(species == 'hsapiens_gene_ensembl') {
    TreatTable <- topTreat(fit2, coef = i, number = 50000)
    TreatTable$hgnc_symbol <- row.names(TreatTable)
    TreatTableEff <- data.table(TreatTable, key = "hgnc_symbol")
    LocationConversionEff <- data.table(LocationConversion, key = "hgnc_symbol")
    FinalTable <- merge(TreatTableEff, LocationConversionEff, by = "hgnc_symbol")
    FinalTable <- FinalTable[!duplicated(FinalTable$hgnc_symbol),,drop = FALSE]
    FinalTable <- data.frame(FinalTable)
    row.names(FinalTable) <- FinalTable[,1]
    FinalTable <- FinalTable[,2:ncol(FinalTable)]
    returnTable <- write.table(
     FinalTable,
      paste0(name, ".csv"),
      sep = ",",
      col.names = NA,
      row.names = TRUE
    )
    
    
}
  else {
    TreatTable <- topTreat(fit2, coef = i, number = 50000)
    FinalTable <- merge(TreatTable, LocationConversion, by.TreatTable = "row.names", by.LocationConversion = "mgi_symbol")
    FinalTable <- FinalTable[!duplicated(FinalTable$hgnc_symbol),,drop = FALSE]
    FinalTable <- data.frame(FinalTable)
    row.names(FinalTable) <- FinalTable[,1]
    FinalTable <- FinalTable[,2:ncol(FinalTable)]
    returnTable <- write.table(
      FinalTable,
      paste0(name, ".csv"),
      sep = ",",
      col.names = NA,
      row.names = TRUE)
  }
  
  if(mysqlUpload == 1) {
    #Establish a connection to the database

    pool <- dbPool(drv = RMariaDB::MariaDB(), user = "root", password = "wellslab123",
                   dbname = "shinyApp", host = "localhost", port = 3306)
    
    #Get number of rows in table (number of total project)
    query <- "SELECT COUNT(*) FROM shinyApp.ProjectDetails;"
    nProjects <- dbGetQuery(pool, query)
    
    #Prepare data for upload
    uploadData <- data.frame(FinalTable[,c(1,5:8)])
    uploadData$ProjectID <- rep(nProjects[1,1], nrow(uploadData))
    uploadData$Comparison <- rep(i, nrow(uploadData))
    uploadData$ComparisonName <- rep(name, nrow(uploadData))
    uploadData$GeneID <- row.names(FinalTable)
    
    uploadData <- uploadData[,c(6:9,1:5)]
    
    #Rename necessary columns
    uploadData <- uploadData %>%
      rename(
        Chromosome_Name = chromosome_name,
        Start_Position = start_position,
        End_Position = end_position
      )
    
    #Upload table data to the database
    pool::dbWriteTable(pool, name = "DiffExp", value = uploadData, append = TRUE)
    
    poolClose(pool)
    
  }
  
    returnList <- c(paste0(name, ".csv"), returnList)
    
    
    if(makeVolc == TRUE) {
      
      res2 <- read.table(
          paste0(name, ".csv"),
          sep = ",",
          header = TRUE,
          row.names = 1
        )
      
      status <- data.frame(row.names = rownames(FinalTable))
      
      pos <- res2[(as.numeric(res2$adj.P.Val) < .05 & as.numeric(res2$logFC) > 0), ]
      if(length(pos$logFC != 0)) {
        pos$result <- 1
        pos <- pos[10]
      }
      
      neg <- res2[(as.numeric(res2$adj.P.Val) < .05 & as.numeric(res2$logFC) < 0), ]
      if(length(neg$logFC != 0)) {
        neg$result <- -1
        neg <- neg[10]
      }
      
      zeroes <- res2[as.numeric(res2$adj.P.Val) >= .05, ]
      if(length(zeroes$logFC != 0)) {
        zeroes$result <- 0
        zeroes <- zeroes[10]
      }
      
      rownames(status) <-
        append(rownames(pos), c(rownames(neg), rownames(zeroes)))
      status$stats <- append(pos$result, c(neg$result, zeroes$result))
      
      #Make sure data from all inputs are in the same order
      logCPM$E <- logCPM$E[order(rownames(logCPM$E)), , drop = FALSE]
      res2 <- res2[order(rownames(res2)), , drop = FALSE]
      status <- status[order(rownames(status)), , drop = FALSE]
      
      
      
     
      xAxis <- res2$logFC
      yAxis <- -log10(as.numeric(res2$adj.P.Val))
      
      df <- matrix(nrow = nrow(res2), ncol = 2)
      
      df <- data.frame(df)
      
      df$log2FC <- xAxis
      df$log10PFDR <- yAxis
      
      df <- df[3:4]
      
      pngName <- paste0(name,"volcPlot.png")
      
      png(pngName)
      returnVolc <- ggplot(df, mapping = aes(x = log2FC, y = log10PFDR, )) + 
        geom_point(aes(colour = factor(status$stats)),
                   size = 1) +
        
        scale_color_manual(name = "status$stats", values = c("-1" = "blue",
                                                             "0" = "black",
                                                             "1" = "red")) +
        
        geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
        geom_hline(yintercept = 1.3, linetype = "dashed", color = "black") +
        
        theme_bw() +
        theme(legend.position = "none", text = element_text(size=14,family = "sans"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.border = element_blank(), axis.line = element_line(colour = "black"), axis.ticks.length=unit(-0.25, "cm"), 
              axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )
      
      finalReturnVolc <- print(returnVolc)
      dev.off()
      
      volcReturns <- c(pngName, volcReturns)
      
    }
    
}

finalReturnList <- list()

finalReturnList$first <- returnList

if(makeVolc == TRUE) {
finalReturnList$second <- volcReturns
}

return(finalReturnList)
}
  
  #Filter data to GSEA Pre-ranked format
  #filter <-
   # read.table(
    #  paste(name, "yougottreat.csv"),
     # sep = ",",
     # header = TRUE,
     # row.names = 1
   # )
  
  
  
  #final formatting
 # filter <- filter[1]
 # colnames(filter) <- NULL
  
  #write table
 # write.table(
   # filter,
    #file = paste(name, ".rnk", sep = ""),
   # sep = "\t",
   # col.names = NA,
    #row.names = TRUE,
   # quote = FALSE
  #)
#}

#}
