#Cluster Analysis
#Group genes based on proximity (neighborhoods)
#Spencer Wells
#17 July 2019

#Reads in paths to differential expression files and parameters to define a neighborhood
ClusterInfo <- function(diffExpFileNames, maxDistance, minLogFC, maxFDR, minGenes) {
  #Load necessary packages
  library(tools)
  
  #initialize list of files to return
  returnList <- c()
  
  #Loop through each differential expression file
  for(i in 1:length(diffExpFileNames)) {
    file <- read.table(diffExpFileNames[i], header = TRUE, sep = ',', row.names = 1)
    
    #Generate name for ClusterFile
    diffExpName <- file_path_sans_ext(diffExpFileNames[i])
    name <- paste0(diffExpName, "Clusters.csv")
    
    #Filter the file to geneIDs, LogFCs, Chromosome, and Start/Endpos
    file <- file[,!(colnames(file) == "AveExpr"|| colnames(file) == "t" || colnames(file) == "P.Value")]
    
    #Filter the file to only include values within the user-specified cutoffs
    file <- file[abs(file$logFC) >= minLogFC,, drop = FALSE]
    file <- file[file$adj.P.Val <= maxFDR,, drop = FALSE]
    
    #Preallocate Clustering matrix
    clustersTable <- matrix(nrow = 10000, ncol = 5000)
    rownames(clustersTable) <- 1:10000
    colnames(clustersTable) <- 1:5000
    
    #Initialize Helper Counter.
    ##Helper counter goes up if genes are within maxDistance apart.
    ##If helper counter reaches MinGenes, cluster is formed
    helperCounter <- 1
    
    #Initialize counter which increments each time new cluster row is created
    clusterRowCounter <- 0
    
    #Initialize Counter which keeps track of column to assign gene to
    geneColCounter <- 2
    
    #Prepare file for ordering
    file$chromosome_name <- as.character(file$chromosome_name)
    file[nchar(file$chromosome_name) == 1 & (file$chromosome_name != "X" & file$chromosome_name != "Y"), "chromosome_name"] <- paste0(0, file[nchar(file$chromosome_name) == 1 & (file$chromosome_name != "X" & file$chromosome_name != "Y"),
                                                "chromosome_name"])
    
    #Order DiffExp file by Chromosome name and start position
    file <- file[order(file$chromosome_name, file$start_position),, drop = FALSE]
    
    
    #Loop through DiffExp file and Generate Clusters
    for(j in 1:(nrow(file)-1)) {
      
      #Check if gene j and gene j+1 are on the same chromosome
      if(file[j,colnames(file) == "chromosome_name"] == file[j+1,colnames(file) == "chromosome_name"]) {
        
        #Check if gene start positions are less than minDistance base pairs apart
        if((file[j+1,colnames(file) == "start_position"] - file[j,colnames(file) == "end_position"]) <= maxDistance) {
          helperCounter = helperCounter + 1
          
          #Only initiazlize/add to cluster if helperCounter is >= minGenes
          if(helperCounter >= minGenes) {
            
            
            #Initialize new cluster if helperCounter == MinGenes
            if(helperCounter == minGenes) {
              clusterRowCounter = clusterRowCounter + 1
              
              #Create counter that indicates how many genes back to go to start creating cluster
              initHelper <- minGenes - 2
              
              #Initialize chromosome name, startpos, and endpos of cluster
              rownames(clustersTable)[clusterRowCounter] <- paste0("Cluster", as.character(clusterRowCounter))
              clustersTable[clusterRowCounter,1] <- file[j-initHelper,colnames(file) == "chromosome_name"]
              
              #Add all genes that should be in cluster
              for(l in seq((initHelper),0,by = -1)) {
                clustersTable[clusterRowCounter,geneColCounter] <- file[j-l,colnames(file) == "start_position"]
                clustersTable[clusterRowCounter,geneColCounter+1] <- file[j-l,colnames(file) == "end_position"]
                clustersTable[clusterRowCounter,geneColCounter+2] <- row.names(file)[j-l]
                geneColCounter = geneColCounter + 3
              }
            } 
            
            #Add gene j+1 to the cluster
            clustersTable[clusterRowCounter,geneColCounter] <- file[j+1,colnames(file) == "start_position"]
            clustersTable[clusterRowCounter,geneColCounter+1] <- file[j+1,colnames(file) == "end_position"]
            clustersTable[clusterRowCounter,geneColCounter+2] <- row.names(file)[j+1]
            geneColCounter = geneColCounter + 3
          }
        } else {
          #Check if cluster is currently being made
          if(helperCounter >= minGenes) {
            
            #Reset the Gene Column Counter
            geneColCounter <- 2
          }
          #Reset the helperCounter
          helperCounter <- 1
          
        } 
        
      }  else {
        #Check if cluster is currently being made
        if(helperCounter >= minGenes) {
          
          #Reset the Gene Column Counter
          geneColCounter <- 2
        }
        #Reset the helperCounter
        helperCounter <- 1
      }
    }

    #Add end position of final cluster if necessary
    if(helperCounter >= minGenes) {
      clustersTable[clusterRowCounter,3] <- file[j,colnames(file) == "end_position"]
    }
    
    #Trim the table to relevant values
    clustersTable <- clustersTable[rowSums(is.na(clustersTable)) != ncol(clustersTable),]
    clustersTable <- clustersTable[,colSums(is.na(clustersTable)) != nrow(clustersTable)]
    
   #Add relevant colnames
   colnames(clustersTable)[1:4] <- c("Chromosome_Name", "Start_Position", "End_Position", "Gene_Name")
   
    ##Label the Table correctly for Chromosome X - Deprecated as of 
    ##converting Chromosome names to character vector - July 18 2019
    ##clustersTable[clustersTable[, "Chromosome_Name"] == "68", colnames(clustersTable) == "Chromosome_Name"] <- "X"
   ##Label the table correctly for special case of chromosome 6
   clustersTable[clustersTable[, "Chromosome_Name"] == "CHR_HSCHR6_MHC_COX_CTG1", colnames(clustersTable) == "Chromosome_Name"] <- "6"
   
    #Write the CSV to be returned to the app
    write.table(clustersTable, name,
                sep = ",",
                col.names = NA,
                row.names = TRUE,
                na = "")
    
    #Add CSV to list of files to return
    returnList <- c(name, returnList)
  }
  
  #Return list of CSVs giving clusters for each comparison
  return(returnList)
}