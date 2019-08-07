#Pathway Analysis with CAMERA

PathwayAnalysis <- function(species, annotation, filter, target, database, minCutoff, makeBars, makePDFs) {
  
  
  #import libraries
  library(edgeR)
  library(limma)
  library(stringi)
  library(biomaRt)
  library(CAMERA)
  library(UpSetR)
  library(textclean)
  
  #Set Working directory
  
  
  #Load count and phenotype label/parameter files to tables
  #x <-
  #read.table(
  #"hEAT RNAseq 2018.txt",
  # sep = "\t",
  #header = TRUE,
  #row.names = 1
  #)
  #target <-
  #read.table(
  # "hEAT RNAseq key file.csv",
  #header = TRUE,
  # sep = ",",
  # row.names = 1
  # )
  
  
  #No files being imported below this point
  #-----------------------------------------------------------------------------------------------------
  
  rownames(filter) <- make.names(filter$X, unique = TRUE)
  
  filter <- filter[2:ncol(filter)]
  filter <- data.frame(filter)
  filter <- filter[!(rownames(filter) == ""),]
  filter <- na.omit(filter)
  filter <- filter[order(rownames(filter)), , drop = FALSE]
  
  
  if(stri_cmp_eq(annotation, "entrezgene")) {
    
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
    
    
    
    
    ensembl <- useEnsembl(biomart = ("ensembl"), dataset = species, mirror = "uswest")
    if (stri_cmp_eq(species, "hsapiens_gene_ensembl")) {
      Conversion <-
        getBM(
          attributes = c(as.character(annotation), 'entrezgene_id'),
          filters = as.character(annotation),
          values = rownames(filter),
          mart = ensembl
        )
      Conversion <-
        Conversion[order(Conversion[annotation]), , drop = FALSE]
      filter$annotation <- rownames(filter)
      names(filter)[ncol(filter)] <- annotation
      finalsheet <- merge(filter, Conversion, by = annotation)
      finalsheet <- finalsheet[!(finalsheet$entrezgene == ""), ]
      finalsheet <- finalsheet[complete.cases(finalsheet),]
      finalsheet <-
        finalsheet[!duplicated(finalsheet$entrezgene_id), , drop = FALSE]
      rownames(finalsheet) <- finalsheet$entrezgene_id
      finalsheet <- finalsheet[2:(ncol(finalsheet) - 1)]
      
    } 
  } 
  
  
  
  x <- finalsheet
  
  #Make matrix of labels for design of experiment and intended pairwise comparisons
  
  #Factor group column of target file to construct design matrix
  
  design <- model.matrix(~target$Group + 0)
  
  #Extract the relevant contrasts of the experiment from the third column of the target file
  contrasts <- as.character(unlist(target[2]))
  contrasts <- noquote(contrasts)
  contrasts <- contrasts[contrasts != ""]
  contrastsSwap <- contrasts
  for(y in 1:length(contrastsSwap)) {
    contrastsSwap[y] <- sub("^([^\r\n-]+)(-*)([^-\r\n]+)$", "\\3\\2\\1", contrastsSwap[y])
  }
  #Arrange columns of design matrix in alphabetical order
  target$Group <- factor(target$Group)
  colnames(design) <- sort(as.character(levels(target$Group)))
  contrasts <- na.omit(contrasts)
  contrast.matrix <-
    makeContrasts(contrasts = contrastsSwap,  levels = design)
  
  #Make DGEList from tables
  dgeFull <- DGEList(counts = x, target$Group)
  
  #logCPM out low count genes
  dgeFull <-
    DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) >= minCutoff, ],
            group = dgeFull$samples$group)
  normCounts <- cpm(dgeFull)
  eff.lib.size <-
    dgeFull$samples$lib.size * dgeFull$samples$norm.factors
  
  #Estimate Normalization factors
  dgeFull <- calcNormFactors(dgeFull, "TMM")
  
  #Normalize Data with voom
  logCPM <-
    voom(
      dgeFull,
      design,
      plot = TRUE,
      header = FALSE,
      row.names = FALSE
    )
  
  #Download and load desired gene sets
  download.file(
    paste0("http://bioinf.wehi.edu.au/software/MSigDB/human_",database,"_v5p2.rdata"),
    paste0("human_",database,"_v5p2.rdata"),
    mode = "wb"
  )
  pathwayData <- load(paste0("human_",database,"_v5p2.rdata"))
  browser()
  #prep camera test with gene set and index of probe identifiers
  
  if(stri_cmp_eq(database, "c2")) {
    cIndex <- ids2indices(Hs.c2, rownames(logCPM))
  }
  
  else if(stri_cmp_eq(database, "c1")) {
    cIndex <- ids2indices(Hs.c1, rownames(logCPM))
  }
  
  else if(stri_cmp_eq(database, "H")) {
    cIndex <- ids2indices(Hs.H, rownames(logCPM))
  }
  
  
  #Prepare for Bar Graph Analysis
  
  xUp <- vector("list", ncol(contrast.matrix))
  xDown <- vector("list", ncol(contrast.matrix))
  LUp <- vector("list", ncol(contrast.matrix))
  LDown <- vector("list", ncol(contrast.matrix))
  
  
  #Conduct Camera test for each comparison
  
  finalReturnList <- c()
  returnListCamera = c()
  for (i in seq(from = 1, to = ncol(contrast.matrix))) {
    name <- paste0(contrasts[i],database)
    cameraTest <-
      camera(
        logCPM,
        index = cIndex,
        design = design,
        contrast = contrast.matrix[, i],
        inter.gene.cor = 0.01
      )
    
    cameraTest$names <- row.names(cameraTest)
    #Write Camera results to excel
    
    returnTable <- write.table(
      cameraTest,
      file = paste0(
        name,
        ".csv"
      ),
      sep = ",",
      quote = FALSE,
      col.names = NA,
      row.names = TRUE
    )
    
    
    returnListCamera <- c(returnListCamera, paste0(name, ".csv"))
    
    if(makeBars != "") {
      Pu <-
        cameraTest[(cameraTest$Direction == "Up" &
                      cameraTest$FDR < .05), ]
      xUp[[i]] <- Pu
      xNamesU <- as.matrix(xUp[[i]]$names)
      LUp[[i]] <- xNamesU
      names(LUp)[[i]] <- paste(contrasts[i], "Up")
      
      #Add downregulated sets to down array
      Pd <-
        cameraTest[(cameraTest$Direction == "Down" &
                      cameraTest$FDR < .05), ]
      xDown[[i]] <- Pd
      xNamesD <- as.matrix(xDown[[i]]$names)
      LDown[[i]] <- xNamesD
      names(LDown)[[i]] <- paste(contrasts[i], "Down")
    }
    
  }
  returnListCSV <- c()
  returnListChart <- c()
  if(makeBars != "") {
    userInput <- makeBars
    x <- 0
    returnListCSV <- c()
    returnListChart <- c()
    
    for (l in seq(from = 1, to = ncol(userInput))) {
      tmpName <- paste0("col", toString(l));
      
      checkSup <- makeBars[, l]
      checkSup <- checkSup[checkSup != ""]
      InputTable <- vector("list", length(checkSup))
      k <- 1
      for (h in seq(from = 1, to = length(checkSup))) {
        for (j in seq(from = 1, to = ncol(contrast.matrix))) {
          if (identical(names(LUp[j]), toString(userInput[h, l]))) {
            InputTable[k] <- LUp[j]
            names(InputTable)[[k]] <- names(LUp[j])
            k <- k + 1
          }
          if (identical(names(LDown[j]), toString(userInput[h, l]))) {
            InputTable[k] <- LDown[j]
            names(InputTable)[[k]] <- names(LDown[j])
            k <- k + 1
          }
        }
      }
      #Construct graph to represent overlaps
      sets <- names(InputTable)
      setsd <- data.frame("names" = sets)
      
      #Subset to last two characters of each element
      setsdSub <-
        substr(setsd$names, nchar(as.character(setsd$names)) - 1, nchar(as.character(setsd$names)))
      
      
      #Coerce Data frames to greppable metadata for UpsetR (distinguish based on Up/Down)
      setlevels <- factor(setsdSub)
      setlevelsd <- data.frame("names" = setlevels)
      metadata <- as.data.frame(cbind(sets, setlevelsd))
      
      
      
      #Extract lists of Up and Down Regulated Sets to prepare calculation of overlaps
      grepgdown <- InputTable[grep("Down", setsd$names, value = TRUE)]
      grepgup <-  InputTable[grep("Up", setsd$names, value = TRUE)]
      n <- length(grepgdown)
      picName = paste0(
        toString(x), database,
        "Overlaps.png"
      )
      browser()
      Overlaps <-
        list(
          Down_Down = intersect(grepgdown[[1]], grepgdown[[2]]),
          Up_Up = intersect(grepgup[[2]], grepgup[[1]]),
          Down_Up = intersect(grepgdown[[1]], grepgup[[2]]),
          Up_Down = intersect(grepgdown[[2]], grepgup[[1]])
        )
      names(Overlaps) <-
        list(
          paste(names(InputTable[3]), "_", names(InputTable[4])),
          paste(names(InputTable[1]), "_", names(InputTable[2])),
          paste(names(InputTable[3]), "_", names(InputTable[2])),
          paste(names(InputTable[1]), "_", names(InputTable[4]))
        )
      noOverlaps <-
        list(
          setdiff(grepgup[[1]], c(grepgdown[[2]], grepgup[[2]])),
          setdiff(grepgup[[2]], c(grepgdown[[1]], grepgup[[1]])),
          setdiff(grepgdown[[1]], c(grepgdown[[2]], grepgup[[2]])),
          setdiff(grepgdown[[2]], c(grepgdown[[1]], grepgup[[1]]))
        )
      if (length(names(Overlaps)) == length(names(InputTable))) {
        names(noOverlaps) <- names(InputTable)
        breakdown <- append(Overlaps, noOverlaps)
        Overlapdf <-
          data.frame(lapply(breakdown, 'length<-', max(lengths(Overlaps))))
        
        #Write an Excel specifying which pathways correspond to which overlap
        
        
        if(makePDFs == TRUE) {
          pdfName <-  paste0(
            tmpName, database,
            toString(x),
            "InclusionExclusion.csv"
          )
          returnCSV <- write.csv(
            Overlapdf,
            pdfName,
            na = ''
          )
          
          returnListCSV <- c(pdfName, returnListCSV)
        }
      }
      
      
      #Determine how to generate the plot based on what comparisons exist in data
      if ((length(intersect(grepgdown[[1]], grepgdown[[2]])) != 0) &
          (length(intersect(grepgup[[1]], grepgup[[2]])) != 0 &
           length(grepgup) <= 2)) {
        png(file = picName, width = 862, height = 622)
        print({
          upset(
            fromList((InputTable)),
            sets = names(InputTable),
            keep.order = T,
            set.metadata = list(data = metadata, plots = list(
              list(
                type = "matrix_rows",
                column = "names",
                colors = c(Up = "navy", wn = "black"),
                alpha = .4
              )
            )),
            queries = list(
              list(
                query = intersects,
                params = grep("Up", setsd$names, value = TRUE),
                color = "blue",
                active = T
              ),
              list(
                query = intersects,
                params = grep("Down", setsd$names, value = TRUE),
                color = "blue",
                active = T
              )
            ),
            nsets = length(InputTable),
            text.scale = 2,
            main.bar.color = "black",
            mainbar.y.label = "number of sets"
          )
        })
      } else if (length(intersect(grepgup[[1]], grepgup[[2]])) != 0 &
                 length(intersect(grepgdown[[1]], grepgdown[[2]])) == 0 &
                 length(grepgup) <= 2) {
        png(file = picName, width = 862, height = 622)
        print({
          upset(
            fromList((InputTable)),
            sets = names(InputTable),
            keep.order = T,
            set.metadata = list(data = metadata, plots = list(
              list(
                type = "matrix_rows",
                column = "names",
                colors = c(Up = "navy", wn = "black"),
                alpha = .4
              )
            )),
            queries = list(
              list(
                query = intersects,
                params = grep("Up", setsd$names, value = TRUE),
                color = "blue",
                active = T
              )
            ),
            nsets = length(InputTable),
            text.scale = 2,
            main.bar.color = "black",
            mainbar.y.label = "number of sets"
          )
        })
      } else if (length(intersect(grepgdown[[1]], grepgdown[[2]])) != 0 &
                 length(intersect(grepgup[[1]], grepgup[[2]])) == 0 &
                 length(grepgup) <= 2) {
        png(file = picName, width = 862, height = 622)
        print({
          upset(
            fromList((InputTable)),
            sets = names(InputTable),
            keep.order = T,
            set.metadata = list(data = metadata, plots = list(
              list(
                type = "matrix_rows",
                column = "names",
                colors = c(Up = "navy", wn = "black"),
                alpha = .4
              )
            )),
            queries = list(
              list(
                query = intersects,
                params = grep("Down", setsd$names, value = TRUE),
                color = "blue",
                active = T
              )
            ),
            nsets = length(InputTable),
            text.scale = 2,
            main.bar.color = "black",
            mainbar.y.label = "number of sets"
          )
        })
        
      } else if (length(grepgup) > 2) {
        png(file = picName, width = 862, height = 622)
        print({
          upset(
            fromList((InputTable)),
            sets = names(InputTable),
            keep.order = T,
            set.metadata = list(data = metadata, plots = list(
              list(
                type = "matrix_rows",
                column = "names",
                colors = c(Up = "navy", wn = "black"),
                alpha = .4
              )
            )),
            nsets = length(InputTable),
            text.scale = 2,
            main.bar.color = "black",
            mainbar.y.label = "number of sets"
          )
        })
        
      }
      
      
      #Ouput graphs
      
      
      
      returnListChart <- c(picName, returnListChart)
      while (!is.null(dev.list()))  dev.off()
      
      x <- x + 1
      
      
    }
    
    
    
  }
  
  browser()
  finalReturnList <- list()
  
  finalReturnList$first <- returnListCamera
  finalReturnList$second <- returnListCSV
  finalReturnList$third <- returnListChart
  
  return(finalReturnList)
}