<!--
*** Thanks for checking out this README Template. If you have a suggestion that would
*** make this better, please fork the repo and create a pull request or simply open
*** an issue with the tag "enhancement".
*** Thanks again! Now go create something AMAZING! :D
-->





<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->

[![MIT License][license-shield]][license-url]




<!-- PROJECT LOGO -->
<br />
<p align="center">



  <h3 align="center">PostSeq</h3>

  <p align="center">
    Intuitive Post-Processing for RNAseq!
    <br />
  </p>
</p>



<!-- TABLE OF CONTENTS -->
## Table of Contents

* [About the Project](#about-the-project)
  * [Built With](#built-with)
* [Getting Started](#getting-started)
  * [Manual Install Prerequisites](#manual-installation-prerequisites)
  * [Manual Installation](#manual-installation)
* [Usage](#usage)
* [Contributing](#contributing)
* [License](#license)
* [Contact](#contact)



<!-- ABOUT THE PROJECT -->
## About The Project

PostSeq is a webapp hosted at 134.68.250.10:3838/shinyApp which increases the effciency and approachability of RNAseq analysis.
PostSeq allows you to complete differential expression, PCA, pathway analysis and more through a streamlined graphical interface that requires no coding.

Additionally, PostSeq allows you to upload your sequencing data to a MySQL database where it can be retrieved or further analyzed at any time.

### Built With
PostSeq is built in Shiny, a popular package used to build webapps with R code. PostSeq is also linked to a MySQL database to host and retrieve sequencing data
* [RShiny](https://shiny.rstudio.com/)
* [MySQL](https://www.mysql.com/)




<!-- GETTING STARTED -->
## Getting Started

PostSeq can be accessed with no installation at 134.68.250.10:3838/shinyApp . This is the recommended way to use the app and allows you to access all functionalities of the program. 

If you want to host a local instance of PostSeq, the setup process consists of launching a shiny server on your own linux server and populating it with the files in this Repo. Detailed instructions on how to install and maintain a shiny server can be found at [https://docs.rstudio.com/shiny-server/](https://docs.rstudio.com/shiny-server/)

Further instructions for the process of hosting a local instance of PostSeq can be found [here](#manual-installation-prerequisites).


<!-- USAGE EXAMPLES -->
## Usage

**Overview**

PostSeq's main workflow consists of uploading your sequencing data (comma delimited counts file or BAM file), generating the experimental design from this data, choosing the desired forms of analysis, and submitting.

Once you have uploaded a data set to MySQL additional forms of analysis are opened up in the MySQL and Essentiality tabs of the app.

**Sequence Data Upload**

If you are uploading a personal counts file, it needs to be a csv with the first column consisting of gene ids and the first row consisting of sample names. Each subsequent row should contain the raw counts for the correct gene and sample.

If you are uploading BAM files, choose all bam files to be processed after selecting the uploading BAM files checkbox and the 
"choose all BAM files to upload" file input. Please note that due to the size of BAM files, this upload process will take some time. Once Bam files have been uploaded, the generated counts file can be downloaded so the BAM files will not have to be uploaded in future analysis.

**Experimental Design Generation**

After uploading your sequence data, it is time to generate your experimental design (design matrix).
Selecting "Let me autofill sample data" will automatically populate the sample names boxes on the right side of the app with the column names from your counts file. This is recommended as it avoids any naming inconsistencies between the counts and experimental design file. If you wish to input samples manually, modify the "Number of Samples" numeric input and then enter each sample name in the sample names boxes. Note that the sample names in these boxes and your counts file must match **exactly**!

**Selecting Classifiers**

Now that your main inputs are sorted, you need to choose the classifiers that match your input data. Choose the species and gene annotation that match your sequencing data. If your counts are from a BAM file, your annotation format will be "HGNC_Symbol".
Currently supported species are human and mouse, and currently supported annotation formats are Ensembl_gene_id, hgnc_symbol, entrezid, and MGI symbol.

**Choosing Whether to Upload to MySQL**

If you want to rerun your analysis in the future or use any of the extended analysis functions in the MySQL and Essentiality tabs of the app, you should upload your data to the MysQL database. After choosing to upload your data, you will need to input a ProjectID, description, and date to be stored in the database. When the analysis is run, your counts file, experimental design file, and differential expression file will be stored in the database. You can then retrieve and analyze the updated data as described [here](#extended-analysis). Important Note: do not upload the same analysis data multiple times. This will cause the database to be filled with redundant information.

**Choosing Outputs to Include in Analysis**

There are a variety of outputs you can include with your analysis to explore your data. The four primary outputs are differential expression using limma, pathway analysis using either CAMERA or GO, principal component analysis, and volcano plots. Some of these primary outputs will have specific options and parameters. In addition to these primary outputs, there are a number of secondary outputs which are generated from the results of the primary analysis.
These options and secondary analysis are described here and grouped by their respective primary analysis:

 1. Differential Expression
  * Options:
    * Gene Filtering: Choose to filter out genes where the counts across all samples do not meet a certain threshold.
    
  * Secondary Analysis:
    * Clustering Effects: The differential expression data can be clustered into groups of 3 or more genes sufficiently close to
    each other on a chromosome. This analysis includes a csv of all the detected clusters, a histogram with the average distance
    between genes in each cluster, and a chromosomal map that will visualize the positions of the clusters. The clustering
    analysis has the following options:
      * MinGenes: Minimum number of genes to consider a cluster
      * MaxDistance: Maximum distance between two genes in a cluster (base pairs) 
      * MaxFDR: Only consider genes with an FDR at or below this value
      * MinLogFC: Only consider genes with a LogFC magnitude (absolute value) at or above this value
      
 2. Pathway Analysis - CAMERA
  * Options:
    * Database Select: Choose which databases to perform pathway analysis on. You have the following three options:
      * c1: Positional data
      * c2: Large set of many pathway from MSigdb database
      * H: Hallmark - curated set of just a few pathways
  * Secondary Analysis
    * Wordclouds: The program will rank order the words in the pathway analysis csv and create a wordcloud of the most common ones (after some filtering of extraneous words)
    * Inclusion/Exclusion Analysis: You will compare the overlaps in the pathways affected in multiple comparisons. You determine what pathways you want compared to each other by dragging the generated boxes on the right side of the screen to containers. Each container represents a set of comparisons that will be tested for overlaps similar to a Venn diagram. The inclusion/exclusion analysis produces the following outputs:
      * Bar graphs: Bar graphs are made with the UpSetR package and function like venn diagrams. There is a legend below the bars to represent what overlap each bar represents
      * csv of the pathways: This will show a matrix of 0s and 1s where pathways affected in a given comparison will be labelled with a 1 and unaffected pathways will be labelled with a 0.
      * Pairwise Comparisons: The program will take all the pairwise comparisons between treatments and calculate which comparisons have the most similar set of the pathways. A csv will be returned that lists these comparisons in order of similarity
      * Pairwise Rankings: The program will return a csv with a quantitative metric of the measured similiarity for all pairwise comparisons
      
 3. Pathway Analysis - GO
  * Options
    * Subontology Select: You can select what subontology go performs pathway analysis on. You can choose from the following three options:
      * Biological Process
      * Cellular Component
      * Molecular Function
      
 4. Principal Component Analysis
  * Returns a png that groups your samples based on PCA and colors them based on their group in the experimental design file
  
 5. Volcano Plot
  * Returns a plot of -Log10FDR vs LogFC and color codes genes either red, blue, or gray depending on the gene being up, down, or insignificant
  


### Manual Installation Prerequisites


### Manual Installation


<!-- CONTRIBUTING -->
## Contributing

Public contributions are the backbone of the open source community. If you want to add something, that's **awesome**.
Here are the instructions on how to do that:

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact

Spencer Wells - spencerwells@gmail.com

Project Link: [https://github.com/SwellsHub/PostSeq](https://github.com/SwellsHub/PostSeq)



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[license-shield]: https://img.shields.io/github/license/othneildrew/Best-README-Template.svg?style=flat-square
[license-url]: https://github.com/othneildrew/Best-README-Template/blob/master/LICENSE.txt
