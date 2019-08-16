#Chromosome Map Generator from Direct Positional Data
#Spencer Wells
#1 August 2019

ChromPlot <- function(chromMap) {
browser()
#Import Necessary Package
library(chromoMap)

chromLengths <- matrix(c(c("01","02","03","04","05","06","07","08","09","10","11","12","13",
                           "14","15","16","17", "18", "19", "20", "21", "22","M", "X","Y"),
                           c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                           c(247249719,242951149,199501827,191273063,180857866,170899992,
                             158821424,146274826,140273252,135374737,134452384,132349534,
                             114142980,106368585,100338915,88827254,78774742,76117153,
                             63811651,62435964,46944323,49691432,16569,154913754,57227415)),
                       ncol = 3)

#chromMap <- read.table("/home/wells/Documents/MillieChromSelect.csv", sep = ",")
chromMap$placeHolder <- 1:nrow(chromMap)
chromMap <- chromMap[,c(ncol(chromMap),1:(ncol(chromMap)-1))]
chromMap$placeHolder <- as.character(chromMap$placeHolder)
chromMap$V1 <- as.character(chromMap$V1)
chromMap$V2 <- as.character(chromMap$V2)
chromMap$V3 <- as.character(chromMap$V3)
#chromMap$V1 <- sub("^...", "", chromMap$V1)
colnames(chromMap)[2:3] <- c("chromosome_name", "start_position")


chromMap$chromosome_name <- as.character(chromMap$chromosome_name)
chromMap[nchar(chromMap$chromosome_name) == 1 & (chromMap$chromosome_name != "X" & chromMap$chromosome_name != "Y" & chromMap$chromosome_name != "M"), "chromosome_name"] <- paste0(0, chromMap[nchar(chromMap$chromosome_name) == 1 & (chromMap$chromosome_name != "X" & chromMap$chromosome_name != "Y" & chromMap$chromosome_name != "M"),
                                                                                                                                          "chromosome_name"])

#Order DiffExp file by Chromosome name and start position
chromMap <- chromMap[order(chromMap$chromosome_name),, drop = FALSE]

write.table(chromLengths, "ChromSelect.txt", sep = "\t", row.names = FALSE, col.names = FALSE, na = '')
write.table(chromMap, "mapToRead.txt", sep = "\t", row.names = FALSE, col.names = FALSE, na = '', quote = FALSE)


returnChromPlot <- chromoMap("ChromSelect.txt", "mapToRead.txt", 
                   segment_annotation = T, chr_color = "lightgrey",
                   lg_x = 100, lg_y = 250, anno_col = c("Darkblue"))

#Define name for plot to export
pngName <- paste0("ChromoPlot.png")

#Save Plot to a png
saveWidget(returnChromPlot, "tmpWrite.html")
webshot::webshot("tmpWrite.html", file = pngName)

return(pngName)
}