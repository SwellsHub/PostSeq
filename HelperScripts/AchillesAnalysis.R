#Clustering Analysis Based on Project Achilles Data
#Spencer Wells
#2 August 2019

EssentialityClusters <- function(chosenProj, chosenComp, chosenGene, chosenCellLines) {
#Load necessary packages
library(dplyr)
library(ggpubr)

#Check if any cell lines have been selected
  browser()
if(length(chosenCellLines) == 0) {
  shinyalert("Bad News", "You didn't select any cell lines. Select some cell lines and rerun the analysis")
  return("error")
}  
  
#Load Project Achilles Data
achillesData <- read.table("ExpandedGeneZSolsCleaned.csv", sep =",", stringsAsFactors = FALSE)


#Format the Achilles Data
row.names(achillesData) <- achillesData[,1]
achillesData <- achillesData[,2:ncol(achillesData)]

tmpNames <- c()
for (i in 1:ncol(achillesData)) {
  tmpNames[i] <- as.character(achillesData[1,i])
}

colnames(achillesData) <- tmpNames
achillesData <- achillesData[2:nrow(achillesData),,drop = FALSE]

tmpFilters <- c()

#Filter to chosen cell lines
for(i in 1:length(chosenCellLines)) {
  tmpFilters[[i]] <- dplyr::select(achillesData,contains(toupper(chosenCellLines[i])))
  if(i == 1) {
    filteredCellLines <- tmpFilters[[1]]
  } else {
    filteredCellLines <- cbind(filteredCellLines, tmpFilters[[i]])
  }
}

#Obtain vector of values for chosen gene to find correlation

x <- as.numeric(as.character(filteredCellLines[row.names(filteredCellLines) == toupper(chosenGene),]),
                stringsAsFactors = FALSE)

checkNA <- is.na(x)
if(length(checkNA[checkNA == TRUE]) == length(x)) {
  shinyalert("Bad News", "That gene could not be found in the selected comparison. Check for typos, or try a different gene", type = "error")
  return("error")
}

#Calculate Pearson Correlation Coefficients for all genes relative to chosen gene
corCoefficients <- matrix(ncol = 2, nrow = nrow(filteredCellLines))
for(i in 1:nrow(filteredCellLines)) {
  corCoefficients[i,1] <- row.names(filteredCellLines)[i]
  corCoefficients[i,2] <- cor(x, as.numeric(filteredCellLines[i,]), method = c("pearson"),
                            use = "pairwise.complete.obs")
}

#Set colnames of correlation coefficient matrix
colnames(corCoefficients) <- c("GeneID", "Pearson_Correlation_Coefficient")

#Convert coefficients matrix to data frame
corCoefficients <- data.frame(corCoefficients, stringsAsFactors = FALSE)


#Connect to MySql and Query the LogFCs 
pool <- dbPool(drv = RMariaDB::MariaDB(), user = "root", password = "wellslab123",
               dbname = "shinyApp", host = "localhost", port = 3306)

query <- paste0(query <- paste0("SELECT s.GeneID, s.LogFC
                    FROM shinyApp.ProjectDetails p
                    LEFT JOIN (shinyApp.DiffExp s)
                    ON s.ProjectID = p.ProjectID
                    WHERE  p.ProjectAlias LIKE '%", chosenProj, "%'
                    AND s.ComparisonName Like '%", chosenComp, "%'
                    GROUP BY s.GeneID;"))

#Testing mode query
#query <- paste0(query <- paste0("SELECT s.GeneID, s.LogFC
                    #FROM shinyApp.ProjectDetails p
                    #LEFT JOIN (shinyApp.DiffExp s)
                    #ON s.ProjectID = p.ProjectID
                    #WHERE  p.ProjectAlias LIKE '%A%'
                    #AND s.ComparisonName Like '%sham%'
                    #AND s.`adj.P.Val` < 0.05
                    #GROUP BY s.GeneID;"))

logFCTable <- pool::dbGetQuery(pool, query)

poolClose(pool)

#Merge results in table to plot
demoTable <- merge(corCoefficients, logFCTable, by = "GeneID")
demoTable$Pearson_Correlation_Coefficient <- as.numeric(demoTable$Pearson_Correlation_Coefficient)

#Filter Table to significant results
#demoTable <- demoTable[abs(demoTable$LogFC) > 1,,drop = FALSE]

#Plot Results
png("LogFCPearsonChart.png")
visualPlot <- ggplot(demoTable, mapping = aes(x = Pearson_Correlation_Coefficient, y = LogFC)) +
  geom_point() +
  scale_x_continuous(breaks = scales::pretty_breaks())
visualPlotFinal <- print(visualPlot)
dev.off()


#Convert to Entrez Genes

ensembl <- useEnsembl(biomart = ("ensembl"), dataset = "hsapiens_gene_ensembl", mirror = "uswest")

tryCatch(
Conversion <-
  getBM(
    attributes = c("hgnc_symbol", 'entrezgene_id'),
    filters = "hgnc_symbol",
    values = demoTable$GeneID,
    mart = ensembl
  ), error = function (err) {
    shinyalert("Bad News", "The Biomart servers may be down right now! If you are sure your inputs are correct, wait 30 minutes then try again. If this still fails, contact spencerwells@gmail.com", type = "error")
    return("error")
  }
)
colnames(Conversion)[1] <- "GeneID"

demoTable <- merge(Conversion, demoTable, by = "GeneID")
demoTable <- demoTable[,2:ncol(demoTable), drop = FALSE]
colnames(demoTable)[1] <- "GeneID"

#Form weighted metric to perform pathway analysis on
demoTable$weightedMetric <- abs(demoTable$Pearson_Correlation_Coefficient)*demoTable$LogFC

statVector <- demoTable$weightedMetric
names(statVector) <- demoTable$GeneID
statVector <- na.omit(statVector)

#Retrieve MSigdb database and form index for pathway analysis
download.file(
  paste0("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata"),
  paste0("human_c2_v5p2.rdata"),
  mode = "wb"
)
pathwayData <- load(paste0("human_c2_v5p2.rdata"))
cIndex <- ids2indices(Hs.c2, demoTable$GeneID)

#Perform pathway analysis on weighted metric
cameraResults <- cameraPR(statVector, cIndex)

browser()
write.table(statVector, "rankForGSEADavid.csv", sep = ",")

#Write Pathway Analysis Table to csv (Testing mode only)
write.table(cameraResults, 
            "GeneCorrelationPathways.csv", 
            sep = ",",
            quote = FALSE,
            col.names = NA,
            row.names = TRUE)
allResults <- c("GeneCorrelationPathways.csv", "LogFCPearsonChart.png")
return(allResults)
}
