#Clustering Visualizations
#Constructs plots that visuallize differential expression clustering data
#Spencer Wells
#18 July 2019


ClusteringHistograms <- function(clusterFiles) {
  #clusterFiles <- list.files("/home/wells/Documents/A172NeighborsAllDataallFiles/ClusteringEffects")
  #Initialize list of files to return
  #setwd("/home/wells/Documents/A172NeighborsAllDataallFiles/ClusteringEffects")
  returnList <- c()
  BigGeneVec <- matrix(nrow = 1000, ncol = length(clusterFiles))
  
  #Iterate over files to generate plots
  for(i in 1:length(clusterFiles)) {
    file <- read.table(clusterFiles[i], header = TRUE, sep = ',', row.names = 1)
    
    #Format the file for the histogram
    
    #Initialize Vector of Averages and Gene Counter Vector
    avgVec <- c()
    geneCounterVec <- c()
    
    #Get average gene length for each cluster
    for(j in 1:nrow(file)) {
      #Initialize helper variables
      sum <- 0
      genePosCol <- 2
      testLoc <- 10
      geneCounter <- 0
      
      #Add lengths for each gene until end of row
      while(file[j,testLoc] != '') {
        sum <- sum + (file[j,genePosCol+3] - file[j,genePosCol+1])
        geneCounter <- geneCounter + 1
        #Don't advance counters if this row is the longest in file
        if(testLoc != ncol(file)) {
          genePosCol <- genePosCol + 3
          testLoc <- testLoc + 3
        } else {
          break
        }
      }
      
      #Compute average of lengths, add to list
      avgLength <- sum/(geneCounter)
      avgVec <- c(avgVec, avgLength)
      geneCounterVec <- c(geneCounterVec, geneCounter)
    }
    #Set the name of the histogram
    fileName <- file_path_sans_ext(clusterFiles[i])
    pngName <- paste0(fileName, "Hist.png")
    
    #Generate the histogram
    png(pngName)
    returnHist <- hist(avgVec, breaks = nrow(file))
    finalReturnHist <- print(returnHist)
    dev.off()
    
    #Add the histogram to the list of files to return
    returnList <- c(returnList, pngName)
    
    #Add the counts of genes
    BigGeneVec[1:nrow(file),i] <- geneCounterVec
  }
  #Return a list of histograms displaying clustering data
  finalReturnList <- c()
  finalReturnList[[1]] <- returnList
  finalReturnList[[2]] <- BigGeneVec
  return (finalReturnList)
}


ClusterChromosomeCharts <- function(clusterFiles, geneCounts) {
  
  #Load necessary packages
  library(purrr)
  library(chromoMap)
  library(webshot)
  browser()
  #clusterFiles <- list.files("/home/wells/Documents/A172NeighborsAllDataallFiles/ClusteringEffects")
  #Initialize list of files to return
  #setwd("/home/wells/Documents/A172NeighborsAllDataallFiles/ClusteringEffects")
  returnList <- c()
  
  #Iterate over files to generate plots
  for(i in 1:length(clusterFiles)) {
    file <- read.table(clusterFiles[i], header = TRUE, sep = ',', row.names = 1)
    
    #Appropriately format the file
    file$Chromosome_Name <- as.character(file$Chromosome_Name)
    file$Chromosome_Name <- gsub("(0)(.$)", "\\2", file$Chromosome_Name)
    file$Chromosome_Name <- as.factor(file$Chromosome_Name)
    
    #Initialize Table to construct chromosome map
    chromMap <- matrix(nrow = 1000, ncol = 3000)
    
    #Iterate over clusters in ClusterFile to generate table
    for(j in 1:nrow(file)) {
      #Get ID, Chromosome name and start position from ClusterFile
      chromMap[j,1] <- row.names(file)[j]
      chromMap[j,2] <- as.character(file[j,colnames(file) == "Chromosome_Name"])
      chromMap[j,3] <- file[j,colnames(file) == "Start_Position"]
      
      #Append blank cell to end of each row to assist in detecting index
      file[,ncol(file)+1] <- NA
      
      #Get index of cluster end position of file
      endIndex <- min(which(is.na(file[j,])))
      
      #Get end position from ClusterFile
      chromMap[j,4] <- file[j,endIndex-2]
      
      #Get number of genes in cluster
      chromMap[j,5] <- geneCounts[j,i]
      
    }
    
    #Trim data to relevant values
    chromMap <- chromMap[1:nrow(file),1:5]
    
    #Write file to be read in by ChromoMap
    chromLengths <- matrix(c(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17", "18", "19", "20", "21", "22", "X"),
                             c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                             c(247249719,242951149,199501827,191273063,180857866,170899992,158821424,146274826,140273252,135374737,134452384,132349534,114142980,106368585,100338915,88827254,78774742,76117153,63811651,62435964,46944323,49691432,154913754)), ncol = 3)
    write.table(chromLengths, "ChromSelect.txt", sep = "\t", row.names = FALSE, col.names = FALSE, na = '')
    write.table(chromMap, "mapToRead.txt", sep = "\t", row.names = FALSE, col.names = FALSE, na = '', quote = FALSE)
    
    #Set file name of Chromosome Plot
    fileName <- file_path_sans_ext(clusterFiles[i])
    pngName <- paste0(fileName, "ChromoPlot.png")
    
    #Generate Chromosomeplot
    returnChromPlot <- chromoMap("ChromSelect.txt", "mapToRead.txt", 
                                 segment_annotation = T, chr_color = "lightgrey", data_based_color_map = T,
                                 data_type = "numeric", legend = T, lg_x = 100, lg_y = 250, data_colors = list(c("blue", "red")))
    
    #Save Plot to a png
    saveWidget(returnChromPlot, "tmpWrite.html")
    webshot::webshot("tmpWrite.html", file = pngName)
    
    #Add Chromosome plot to existing list of plots
    returnList <- c(returnList, pngName)
    }
  
  #Return a list of chromosome plots
  return(returnList)
}