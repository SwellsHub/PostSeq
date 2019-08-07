# GUI Interface for RNA Sequencing Applications
# Spencer Wells
# 21/5/19

# Load Packages
library("shiny")
library("limma")
library("biomaRt")
library("mzR")
library("CAMERA")
library("Glimma")
library("quantmod")
library("shinythemes")
library("shinyjs")
library("htmlwidgets")
library("readr")
library("readxl")
library("shinyjqui")
library("DBI")
library("pool")
library("data.table")
library("tidyverse")
library("DT")
library("grid")
library("Rsubread")
library("ggplot2")
library("chromoMap")
library("filesstrings")



# Source Relevant R Scripts
source("LimmaGSEAEfficient.R")
source("CameraPathwayAnalysis.R")
source("PCAPlot.R")
source("wordClouds.R")
source("GeneNeighbors.R")
source("ClusterPlots.R")
source("ChromChart.R")
source("AchillesAnalysis.R")

#Read in supplementary files for wordclouds
blackList <-
  read.table("data/blacklist.csv", sep = ",")



Modifiers <-
  read.csv(
    "data/modifiers.csv",
    stringsAsFactors = FALSE
  )

tweaks <- 
  list(tags$head(tags$style(HTML("
                                 .multicol { 
                                   height: 200px;
                                   -webkit-column-count: 3; /* Chrome, Safari, Opera */ 
                                   -moz-column-count: 3;    /* Firefox */ 
                                   column-count: 3; 
                                   -moz-column-fill: auto;
                                   -column-fill: auto;
                                 } 
                                 ")) 
  ))

#Specify UI of application
ui <- fluidPage(theme = shinytheme("simplex"),
                
                #Use JS    
                useShinyjs(),
                
                #Include JS Script
                includeScript("www/check.js"),
                
                #Theme Selector
                
                #MulticolCheckboxGroup
                tweaks,
                
                #CSS Styles
                tags$head(
                  tags$style(
                    HTML(
                      ".checkbox-inline { 
                    margin-left: 0px;
                    margin-right: 10px;
          }
         .checkbox-inline+.checkbox-inline {
                    margin-left: 0px;
                    margin-right: 10px;
          }
        "
                    )
                  ) 
                ),
                
                tags$style(HTML("hr {
   display: block;
    margin-before: 0.5em;
    margin-after: 0.5em;
    margin-start: auto;
    margin-end: auto;
    overflow: hidden;
    border-style: inset;
    border-width: 1px;
    }")),
                tags$head(tags$style(
                  type="text/css",
                  "img {width: 100%; min-width: 200px; min-height: 200px; max-width: 550px; max-height: 500px; height: 100%;}"
                )),
                tags$head(tags$style(HTML("#helpLink {color: blue}"))),
                tags$head(tags$style(HTML("a {color: gray}"))),
                
                
                
                #JS Script
                tags$head(tags$script('
                        var width = 0;
                        $(document).on("shiny:connected", function(e) {
                          width = window.innerWidth;
                          Shiny.onInputChange("width", width);
                        });
                        $(window).resize(function(e) {
                          width = window.innerWidth;
                          Shiny.onInputChange("width", width);
                        });
                        ')),
                
                
                
                
                titlePanel("RNA Sequencing Post Processing"),
                
                tabsetPanel(id = "allTabs", type = "tabs",
                            
                            tabPanel(title = "Input", value = "introInput",                     
                                     
                                     sidebarLayout(
                                       
                                       sidebarPanel(
                                         conditionalPanel(condition = "input.check.data[1] == true",
                                                          div(style = "color: red;","Please switch to Chrome for full functionality")),
                                         strong("Provide the experimental inputs and the desired mode of analysis"),
                                         
                                         br(),
                                         br(),
                                         
                                         textInput("projectName", label = "Project Name"),
                                         
                                         checkboxInput("mysqlUpload", "Check to upload differential expression results to MySQL database"),
                                         
                                         conditionalPanel(condition = "input.mysqlUpload == 1",
                                         textInput("projectID", label = "Unique Project ID"),
                                         textAreaInput("projectDescr", label = "Project Description", resize = "vertical"),
                                         dateInput("projectDate", label = "Project Date")
                                         ),
                                         
                                         selectInput("species", label = "Species of Samples", choices = 
                                                       list(Human = "hsapiens_gene_ensembl", Mouse = "mmusculus_gene_ensembl"), selected = "Human"),
                                         
                                         
                                         selectInput("nomenclature", label = "Input Format", choices =
                                                       list(Ensembl = "ensembl_gene_id",  HGNC_Symbol = "hgnc_symbol", Entrez = "entrezgene_id", MGI = 'mgi_symbol')),
                                         
                                         conditionalPanel(condition = "input.convertBam == 0 && input.uploadCore == 0",
                                         checkboxInput("pickMysql", label = "Check to retrieve information from previous run")
                                         ),
                                         
                                         conditionalPanel(condition = "input.pickMysql == 0 && input.uploadCore == 0",
                                         checkboxInput("convertBam", "Check if uploading BAM files (Long wait times)", value = FALSE)
                                         ),
                                         
                                         conditionalPanel(condition = "input.convertBam == 0 && input.pickMysql == 0",
                                                          checkboxInput("uploadCore", "Check if uploading counts file from the core")
                                         ),
                                         
                                         
                                         br(),
                                         
                                         conditionalPanel(condition = "input.convertBam == 0", 
                                         div(style = "height: 75px;",
                                             fileInput("fileSelectExpr", label = "Gene Counts File")
                                         )),
                                         
                                         conditionalPanel(condition = "input.convertBam == 1", 
                                                          div(style = "height: 75px;",
                                                              fileInput("bamInputs", label = "Choose all BAM files to upload", multiple = TRUE)
                                                          )),
                                         
                                         
                                         actionLink("helpLink", label = "Need Help?"),
                                         
                                         conditionalPanel(condition = "input.convertBam == 1", 
                                         downloadButton("bamCounts", label = "Download counts from BAM")
                                         ),
                                        
                                          
                                         fluidRow(column(width = 5,
                                                         div(id = "compCheckBox",
                                                             checkboxInput("expEntryMode", "Input Experiment Design Manually", value = TRUE), style = "font-style:bold;")
                                         ), column(width = 2, helpText("or")),
                                         column(width = 4,
                                                div(id = "cvUploadBox", checkboxInput("dispCVInput", "Upload Design CSV", value = FALSE), style = "font-style:bold;"))
                                         ),
                                         
                                         
                                         conditionalPanel(condition = "input.expEntryMode == 0",
                                                          fileInput("fileSelectComp", label = "Experimental Design File")
                                         ),
                                         conditionalPanel(condition = "input.expEntryMode == 1 ",
                                                          
                                                          conditionalPanel(condition = "input.convertBam == 0", 
                                                          checkboxInput("sampleGuess", "Let me autofill sample data from counts file")
                                                          ),
                                                          div(id = "insertGuessLocation"),
                                                          
                                                          numericInput("nSamples", "Number of Samples:", value = 1, min = 1),
                                                          
                                                          numericInput("nComparisons", "Number of Comparisons", value = 1, min = 1)
                                         ),
                                         
                                         
                                         div(id = "compFileInput"),
                                         
                                         
                                         checkboxInput("advancedToggle", "Show Advanced Parameters", value = FALSE),
                                         
                                         conditionalPanel(condition = "input.advancedToggle == 1",
                                                          sliderInput("minCutoff", label = "Filter genes with counts below:", min = 0, max = 200, value = 1)
                                                          
                                                          
                                                          
                                         ),
                                         
                                         
                                         checkboxGroupInput("outputParameters", "Select files to include in output:", choiceNames = 
                                                              c("Differential Expression - Limma", "Pathway Analysis - CAMERA", "Principal Component Analysis", "Volcano Plot"), choiceValues =
                                                              c("DE", "PA", "PCA", "VP")),
                                         
                                         
                                         conditionalPanel(condition = "input.outputParameters.includes('DE')",
                                                          checkboxInput("clusterCheck", "Include CSV of gene clustering effects (Execution Time >= 30 min)", value = FALSE),
                                                          conditionalPanel(condition = "input.clusterCheck == 1",
                                                                           numericInput("maxClusterDistance", "Maximum distance between two genes in cluster (bp)", min = 0, value = 1000000),
                                                                           sliderInput("minClusterLogFC", "Minimum logFC magnitude", min = 0, max = 10, value = 1),
                                                                           sliderInput("maxClusterFDR", "FDR cutoff", min = 0, max = 1, value = .05),
                                                                           sliderInput("minClusterGenes", "Minimum number of genes to consider a cluster", min = 3, max = 25, value = 3)
                                                                           )
                                                          ),
                                         
                                         conditionalPanel(condition = "input.outputParameters.includes('PA')",
                                                          checkboxGroupInput("databaseSelect", "Database(s) for pathway analysis:", choiceNames = c("c1","c2 (recommended)", "H"), choiceValues = c("c1","c2","H")),
                                                          strong("Additional outputs to include"),
                                                          checkboxInput("wordClouds", "Wordclouds"),
                                                          checkboxInput("barGraphs", "Inclusion/Exclusion Bar Graphs")
                                         ),
                                         
                                         conditionalPanel(condition = "input.barGraphs == 1", 
                                                          checkboxInput("autoBars", "Check to upload CSV of desired inclusion/exclusion charts"),
                                                          checkboxInput("barPDF", "Include CSV of inclusion/exclusion data in output"),
                                                          checkboxInput("overlapCloud", "Include wordclouds of overlap data (Experimental feature. Use at your own risk.)")
                                                          ),
                                                          
                                         
                                         conditionalPanel(condition = "input.autoBars == 1",
                                                          fileInput("barUpload", "Desired inclusion/exclusion charts")
                                         ),
                                         
                                         
                                         downloadButton("submitButton", label = "Run and Save Analysis"),
                                         
                                         div(style = "color: red;", "Analysis started. May take several minutes.", id = "downloadMessage"),
                                         div(style = "color: green;", "Analysis complete! Files Downloaded.", id = "downloadComplete"),
                                         tags$script(JS("$('#downloadMessage').toggle()")),
                                         tags$script(JS("$('#downloadComplete').toggle()"))
                                         
                                         
                                         
                                       ),
                                       
                                       
                                       
                                       mainPanel(conditionalPanel(condition = "input.expEntryMode == 1", wellPanel(id = "mainWell", 
                                                                                                                   h1("Experimental Design"),
                                                                                                                   
                                                                                                                   br(),
                                                                                                                   
                                                                                                                   fluidRow( column(width = 4, h4(strong("Sample Names"))),
                                                                                                                             column(width = 4, h4(strong("Sample Groups"))),
                                                                                                                             column(width = 4, h4(strong("Comparisons")),
                                                                                                                                    strong("(Control-Treatment)"))
                                                                                                                             
                                                                                                                   ),
                                                                                                                   
                                                                                                                   fluidRow( column(width = 4,
                                                                                                                                    uiOutput("sampleInputs")),
                                                                                                                             column(width = 4,
                                                                                                                                    uiOutput("groupInputs")),
                                                                                                                             column(width = 4,
                                                                                                                                    uiOutput("compInputs"),
                                                                                                                                    downloadButton("designDownload", label = "Save Configuration", align = "left"))
                                                                                                                             
                                                                                                                   )
                                                                                                                   
                                                                                                                   
                                       )
                                       ),
                                       
                                       
                                       conditionalPanel(condition = "input.barGraphs == 1", wellPanel(id = "barPanel",
                                                                                                      h1("Bar Graph Comparisons"),
                                                                                                      
                                                                                                      br(),
                                                                                                      
                                                                                                      fluidRow(
                                                                                                        div( class = "col-12",
                                                                                                             uiOutput("compSource")),
                                                                                                        div( class = "col-12",
                                                                                                             uiOutput("compSamples"),
                                                                                                             actionButton("clearBars", "Clear")),
                                                                                                        downloadButton("makeBarComps", "Save Configuration")
                                                                                                        
                                                                                                        
                                                                                                      )
                                                                                                      
                                                                                                      
                                       )
                                       
                                       )
                                       
                                       )
                                     )
                            ),
                            
                            tabPanel(title = "Mysql", value = 'mysql',
                                     
                                     wellPanel(
                                       
                                       selectInput("searchMode", label = "Results to display", choices = 
                                                     list(Project_Details = "ProjectDetails", Project_Results = "DiffExp", Gene_Tracks = "GeneTracks")
                                                     ),
                                       conditionalPanel(condition = "input.searchMode != 'GeneTracks'",
                                       selectInput("searchType", label = "Parameter to search by", choices = 
                                                     list(Project_ID = "p.ProjectAlias", Comparison = "s.ComparisonName", GeneID = "s.GeneID" ), selected = "Human")
                                       ),
                                       conditionalPanel(condition = "input.searchMode == 'GeneTracks'",
                                                        textInput("searchTypeProject", label = "Selected Project"),
                                                        dataTableOutput("searchTypeProjTable")
                                       ),
                                       conditionalPanel(condition = "input.searchMode == 'GeneTracks'",
                                                        textInput("searchTypeComparison", label = "Selected Comparison"),
                                                        dataTableOutput("searchTypeCompTable")
                                       ),
                                       fluidRow(
                                         column(width = 10,
                                                textInput("dataSearch", label = "Enter search term")
                                                ),
                                         column(width = 2,
                                                br(),
                                                actionButton("submitQuery", label = "Search")
                                                )
                                         
                                       ),
                                       
                                       dataTableOutput("searchData"),
                                       
                                       uiOutput("searchDownload"),
                                      
                                       plotOutput("geneTracks", width = "100%"),
                                       
                                       uiOutput("trackDownload")
                                       
                                       
                                     )
                                   
                                     
                              
                            ),
                            
                            tabPanel(title = "ChromoPlot", value = "chromoplot",
                                     
                                     wellPanel(
                                       fileInput("chromPlotInput", "Upload positional data"),
                                       downloadButton("chromPlotDownload"),
                                       div(style = "color: red;", "Analysis started. May take several minutes.", id = "chromMessage"),
                                       div(style = "color: green;", "Analysis complete! Files Downloaded.", id = "chromComplete"),
                                       tags$script(JS("$('#chromMessage').toggle()")),
                                       tags$script(JS("$('#chromComplete').toggle()"))
                                       
                                       
                                     )
                                     ),
                            
                            tabPanel(title = "Essentiality", value = "essentiality", id="essentiality",
                                     
                                     wellPanel(
                                       
                                       textInput("eProjInput", "Selected project"),
                                       dataTableOutput("eProjInputTable"),
                                       textInput("eCompInput", "Selected comparison"),
                                       dataTableOutput("eCompInputTable"),
                                       textInput("eGeneInput", strong("Enter GeneID to compare against")),
                                       dataTableOutput("eGeneInputTable"),
                                       tags$div(align = "left", class = "multicol",
                                       checkboxGroupInput("eCellLineInput",
                                                          "Select cell lines to include",
                                                          choices = c("Bone", "Breast",
                                                                          "Central_Nervous_System",
                                                                          "Cervix", "Endometrium",
                                                                          "Haematopoietic_and_lymphoid_tissue", "Kidney",
                                                                          "Large_intestine", "Liver", "Lung",
                                                                          "Oesophagus", "Ovary", "Pancreas",
                                                                          "Pleura", "Prostate", "Skin",
                                                                          "Small_intestine", "Soft_tissue",
                                                                          "Stomach", "Upper_aerodigestive_tract",
                                                                          "Urinary_tract"), inline = FALSE)),
                                       downloadButton("corDownload", "Submit and Download"),
                                       div(style = "color: red;", "Analysis started. May take several minutes.", id = "essMessage"),
                                       div(style = "color: green;", "Analysis complete! Files Downloaded.", id = "essComplete"),
                                       tags$script(JS("$('#essMessage').toggle()")),
                                       tags$script(JS("$('#essComplete').toggle()"))
                                       
                                     )
                                     
                                     ),
                            
                            tabPanel(title = "Instructions", value = "instructions",
                                     wellPanel(
                                       
                                       strong(h1("Instructions and Reference")),
                                       hr(),
                                       fluidRow(
                                         div( class = "col-12 col-lg-6",
                                              h3("Layout of the Gene Counts File"),
                                              fluidRow(
                                                column(width = 4, 
                                                       p("Upload the Gene Counts file as a .csv with the following specifications:"),
                                                       
                                                       p("     1.First column contains signature names"),
                                                       p("     2.Each row contains raw counts for a signature"),
                                                       p("     3.The First row contains the names of each sample"),
                                                       p("     4.Each column contains a unique sample"),
                                                       p("     5.Do not include dashes (-) in the names of the samples")
                                                ),
                                                
                                                column(width = 8,
                                                       fluidRow(
                                                         tags$img(src = 'sampleCountsFile.png', height = 250, width = 350, align = "left")),
                                                       fluidRow(
                                                         helpText("Sample gene counts file")
                                                       )
                                                )
                                                
                                                
                                              )
                                         ),
                                         div( class = "col-12 col-lg-6",
                                              conditionalPanel(condition = "input.width < 1200",
                                                               hr()
                                              ),
                                              fluidRow(
                                                h3("Layout of the Experimental Design File"),
                                                p(paste0("The experimental data can be entered manually into the provided boxes to generate a design file for future usage.",
                                                         " Alternatively, a design.csv file can be created and uploaded according to the following specifications:")),
                                                
                                                
                                                
                                                column(width = 4,
                                                       
                                                       p("     1.First column contains sample names (must be identical to column names in counts file)"),
                                                       p("     2.Second column contains sample groups (each set of replicates should form a group)"),
                                                       p("     3.Third column contains desired comparisons in the form 'SampleNameofTreatment-SampleNameofControl'"),
                                                       p("     4.Do not include dashes (-) in the group names. Only use dashes in the comparisons")
                                                ),
                                                column(width = 8,
                                                       fluidRow(
                                                         tags$img(src = 'sampleDesignFile.png', height = 250, width = 250, align = "left")
                                                       ),
                                                       fluidRow(
                                                         helpText("sample experimental design file")
                                                       )
                                                )
                                              )
                                         )
                                       )
                                     )
                            )
                )
)







server <- function(input,output, session) {
  
  shinyjs::disable("submitButton")
  
  shouldParseVal <- reactiveValues(shouldParse = FALSE)
  
  globalCounter <<- 1
  globalCompCounter <<- 1
  globalCompCounterTrack <<- 1
  
  observe({
    #Processing for Essentiality Tab
    if(input$allTabs == "essentiality" | input$allTabs == "mysql") {
    
      
      
      pool <- dbPool(drv = RMariaDB::MariaDB(), user = "root", password = "wellslab123",
                     dbname = "shinyApp", host = "localhost", port = 3306)
      
      query <- paste0("SELECT p.ProjectAlias
                    FROM shinyApp.ProjectDetails p;")
      
      
      displayData <- dbGetQuery(pool, query)
      
      if(globalCounter == 1) {
      #DT::renderDataTable(relevantResults, selection = list(mode = "single", selected = 1)
      output$eProjInputTable <- DT::renderDataTable(displayData,
                                                    selection = list(mode = "single",
                                                                     selected = 1)
      )
      output$searchTypeProjTable <- DT::renderDataTable(displayData,
                                                        selection = list(mode = "single",
                                                                         selected = 1)
      )
      
      globalCounter <<- 2
      }
      
      
      if(input$allTabs == "essentiality") {
        if(is.null(input$eProjInputTable_row_last_clicked)) {
          inputSelect <- 1
        } else {
          inputSelect <- input$eProjInputTable_row_last_clicked
        }
      }
      
      if(input$allTabs == "mysql") {
        if(is.null(input$searchTypeProjTable_row_last_clicked)) {
          inputSelect <- 1
        } else {
          inputSelect <- input$searchTypeProjTable_row_last_clicked
        }
      }
      
      query <- paste0("SELECT s.ComparisonName
                    FROM shinyApp.DiffExp s
                    LEFT JOIN (shinyApp.ProjectDetails p)
                    ON s.ProjectID = p.ProjectID
                    WHERE p.ProjectID = ",inputSelect,"
                    GROUP BY s.ComparisonName;")
      
      displayDataComp <- dbGetQuery(pool, query)
      
      if(globalCompCounter == 1) {
        if(input$allTabs == "essentiality") {
      output$eCompInputTable <- DT::renderDataTable(displayDataComp, selection = list(mode = "single", selected = 1))
      globalCompCounter <<- 2
        }
      }
      if(globalCompCounterTrack == 1) {
        if(input$allTabs == "mysql") {
        output$searchTypeCompTable <- DT::renderDataTable(displayDataComp, selection = list(mode = "single", selected = 1))
        globalCompCounterTrack <<- 2
        }
      }
      
      if(input$allTabs == "essentiality") {
        if(is.null(input$eCompInputTable_row_last_clicked)) {
          inputSelectComp <- 1
        } else {
          inputSelectComp <- input$eCompInputTable_row_last_clicked
        }
      }
      
      if(input$allTabs == "mysql") {
        if(is.null(input$searchTypeCompTable_row_last_clicked)) {
          inputSelectComp <- 1
        } else {
          inputSelectComp <- input$searchTypeCompTable_row_last_clicked
        }
      }
      
      query <- paste0("SELECT p.ProjectAlias, s.ComparisonName
                      FROM shinyApp.ProjectDetails p
                      LEFT JOIN(shinyApp.DiffExp s)
                      ON s.ProjectID = p.ProjectID
                      WHERE p.ProjectID = ",inputSelect,"
                      AND s.Comparison = ",inputSelectComp,"
                      GROUP BY s.ComparisonName;")
      
      textUpdateFields <- dbGetQuery(pool, query)
      
      #Update text fields based on selected row in table
      if(input$allTabs == "essentiality") {
        updateTextInput(session, "eProjInput", value = as.character(textUpdateFields$ProjectAlias))
        updateTextInput(session, "eCompInput", value = as.character(textUpdateFields$ComparisonName))
      }
      if(input$allTabs == "mysql") {
        updateTextInput(session, "searchTypeProject", value = as.character(textUpdateFields$ProjectAlias))
        updateTextInput(session, "searchTypeComparison", value = as.character(textUpdateFields$ComparisonName))
      }
      
      poolClose(pool)
    }
    
    
    if((!(is.null(input$fileSelectExpr) & is.null(input$bamInputs)) && !is.null(input$outputParameters)) && !("PA" %in% input$outputParameters && is.null(input$databaseSelect)) && !((is.null(input$fileSelectComp)) && (input$expEntryMode == 0))) {
      
      shinyjs::enable("submitButton")
      
    } else {
      shinyjs::disable("submitButton")
      validate(need(!is.null(input$fileSelectExpr), "Please upload a gene counts file!"))
      validate(need(!is.null(input$fileSelectComp), "Please upload an experimental design file"))
    }
    
  })
  
  observeEvent({input$eProjInputTable_row_last_clicked}, {
    
    globalCompCounter <<- 1
    
  })
  
  observeEvent({input$searchTypeProjTable_row_last_clicked}, {
    
    globalCompCounterTrack <<- 1
    
  })
  
  options(shiny.maxRequestSize=300000*1024^2) 
  
  compUIExist = FALSE
  
  rv <- reactiveValues(download_flag = 0)
  
  
  
  
  
  observeEvent({input$dispCVInput}, {
    if(!compUIExist && input$dispCVInput == 1) {
      updateCheckboxInput(session, "expEntryMode", value = !input$dispCVInput)
    }
    if(compUIExist && input$dispCVInput == 0) {
      updateCheckboxInput(session, "expEntryMode", value = !input$dispCVInput)
    }
  })
  
  
  observeEvent({input$expEntryMode}, { 
    
    
    if(!compUIExist && input$expEntryMode == 0) {
      updateCheckboxInput(session, "dispCVInput", value = !input$expEntryMode)
      
      
      compUIExist <<- TRUE
    }
    
    if(compUIExist && input$expEntryMode == 1) {
      updateCheckboxInput(session, "dispCVInput", value = !input$expEntryMode)
      
      
      compUIExist <<- FALSE
    }
    
  })  
  
  observeEvent(input$helpLink, {
    
    updateTabsetPanel(session, "allTabs", selected = "instructions")
    
  })
  
  observeEvent(input$sampleGuess, {
    
    if(is.null(input$fileSelectExpr$name) && input$sampleGuess == TRUE) {
      
      insertUI(selector = "#insertGuessLocation", where = "afterEnd", ui =
                 div(id = "helpGuess", helpText("Upload counts file and recheck box")))
      
      shouldRemove = TRUE
      
    } else {
      shouldRemove = FALSE
    }
    if((input$sampleGuess == 1) && !(is.null(input$fileSelectExpr) && is.null(input$bamInputs))) {
      shouldParseVal$shouldParse = TRUE
      parseSamples()
      
      if(exists("shouldRemove")) {
        removeUI(selector = "#helpGuess")
      }
      
    } else if(input$sampleGuess == 0){
      shouldParseVal$shouldParse == FALSE
      clearSamples()
      
      if(exists("shouldRemove")) {
        removeUI(selector = "#helpGuess")
      }
      
    } 
      
    
    
   
    
  })
  
  observeEvent(input$submitQuery, {
    pool <- dbPool(drv = RMariaDB::MariaDB(), user = "root", password = "wellslab123",
                   dbname = "shinyApp", host = "localhost", port = 3306)
    
    
    if(input$searchMode == "ProjectDetails") {
            query <- paste0("SELECT DISTINCT p.ProjectID, p.ProjectAlias, p.ProjectDescr, p.ProjectDate 
                    FROM shinyApp.ProjectDetails p
                    LEFT JOIN (shinyApp.DiffExp s)
                    ON s.ProjectID = p.ProjectID
                    WHERE ", input$searchType,  " LIKE '%", input$dataSearch, "%';")
    
            relevantResults <- dbGetQuery(pool, query)
            poolClose(pool)
            #Render search results in interactive data table
            output$searchData <- DT::renderDataTable(relevantResults, selection = list(mode = "single", selected = 1),
            )
    } else if(input$searchMode == "DiffExp"){
            query <- paste0("SELECT GeneID, LogFC, ProjectAlias 
                    FROM shinyApp.DiffExp s
                    LEFT JOIN (shinyApp.ProjectDetails p)
                    ON s.ProjectID = p.ProjectID
                    WHERE ", input$searchType,  " LIKE '%", input$dataSearch, "%';")
    
            relevantResults <- dbGetQuery(pool, query)
            poolClose(pool)
            #Render search results in interactive data table
            output$searchData <- DT::renderDataTable(relevantResults, selection = list(mode = "single", selected = 1),
            )
    } else {
      query <- paste0("SELECT GeneID, LogFC, Chromosome_Name, Start_Position, End_Position, `adj.P.Val`
                      FROM shinyApp.DiffExp s
                      LEFT JOIN (shinyApp.ProjectDetails p)
                      on s.ProjectID = p.ProjectID
                      WHERE ComparisonName LIKE '%",input$searchTypeComparison, "%' 
                      AND p.ProjectAlias LIKE '%", input$searchTypeProject, "%'
                      GROUP BY GeneID;")
      relevantResults <- dbGetQuery(pool, query)
      
      poolClose(pool)
      
      relevantResults$Chromosome_Name <- as.character(relevantResults$Chromosome_Name)
      relevantResults[nchar(relevantResults$Chromosome_Name) == 1 & (relevantResults$Chromosome_Name != "X" & relevantResults$Chromosome_Name != "Y"), "Chromosome_Name"] <- paste0(0, relevantResults[nchar(relevantResults$Chromosome_Name) == 1 & (relevantResults$Chromosome_Name != "X" & relevantResults$Chromosome_Name != "Y"),
                                                                                                                                                "Chromosome_Name"])
      
      #Order DiffExp file by Chromosome name and start position
      relevantResults <- relevantResults[order(relevantResults$Chromosome_Name, relevantResults$Start_Position),, drop = FALSE]
      
      relevantResults$Significance <- rep("FDR>.05", nrow(relevantResults))
      relevantResults[relevantResults$adj.P.Val < .05, colnames(relevantResults) == "Significance"] <- "FDR<.05"
      
      index <- which(relevantResults$GeneID == toupper(input$dataSearch))
      
      relevantResults <- relevantResults[c((index-3):(index+3)),c(1,2,4,5, ncol(relevantResults))]
      
      relevantResults$GeneID <- factor(relevantResults$GeneID, levels = relevantResults$GeneID)
      
      
      Distances <- c()
      annotDistances <- c()
      for (i in 1:(nrow(relevantResults)-1)) {
        if((relevantResults$Start_Position[i+1] - relevantResults$End_Position[i]) >= 0) {
          Distances[[i]] <- relevantResults$Start_Position[i+1] - relevantResults$End_Position[i]
        } else {
          Distances[[i]] <- "Overlap"
        }
        annotDistances[[i]] <- grid.text(Distances[i], x = ((i*1.4)/10), y = .02, gp=gpar(fontsize = 11, fontface = "bold"))
      }
      
      svgName <- paste0(input$dataSearch, "Track.png")
      browser()
        png(svgName)
      plotTokeep <- ggplot(data = relevantResults,
                           aes(x = GeneID, y = LogFC, names.arg)) + 
        geom_bar(aes(fill = Significance), stat="identity") +
        scale_fill_manual(values = c("black", "grey")) +
        theme_bw() +
        theme( text = element_text(size=14,family = "sans"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.border = element_blank(), axis.line = element_line(colour = "black"))+
        annotation_custom(annotDistances[[1]])+
        annotation_custom(annotDistances[[2]])+
        annotation_custom(annotDistances[[3]])+
        annotation_custom(annotDistances[[4]])+
        annotation_custom(annotDistances[[5]])+
        annotation_custom(annotDistances[[6]])
      
      plotToKeepFinal <- print(plotTokeep)
      dev.off()
            
      output$geneTracks <- renderPlot({ggplot(data = relevantResults,
                                               aes(x = GeneID, y = LogFC, names.arg)) + 
                                               geom_bar(aes(fill = Significance), stat="identity") +
                                                scale_fill_manual(values = c("black", "grey")) +
                                                theme_bw() +
                                                
                                                theme( text = element_text(size=20,family = "sans"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                panel.border = element_blank(), axis.line = element_line(colour = "black")) +
                                                annotation_custom(annotDistances[[1]])+
                                                annotation_custom(annotDistances[[2]])+
                                                annotation_custom(annotDistances[[3]])+
                                                annotation_custom(annotDistances[[4]])+
                                                annotation_custom(annotDistances[[5]])+
                                                annotation_custom(annotDistances[[6]])
    }, width = 707, height = 400)
      
      
      
    }
    if(input$searchMode == "ProjectDetails") {
      output$searchDownload <- renderUI(
        downloadButton("searchInputsDownload", "Download Project Inputs (~3 min)")
      )
    } else if(input$searchMode == "GeneTracks"){
      output$trackDownload <- renderUI(
        downloadButton("trackDownloadButton", "Download picture of track")
      )
    } else {
      removeUI(selector = "#searchInputsDownload")
      removeUI(selector = "#trackDownloadButton")
    }
  })
  
  
  
  output$sampleInputs <- renderUI({
    if(isNumeric(input$nSamples)) {
      if(input$nSamples >= 1) {
        textInputs <- lapply(1:input$nSamples, function(number){
          
          textId <- paste0("sample.", number)
          textInputId <- paste0("sampleInput.", number)
          div(
            textInput(inputId = textInputId, label = NULL),
            style = "font-size:14px;"
          )
          
        })
      }  else {
        helpText("Invalid Entry for Number of Samples")
      }
    } else {
      helpText("Invalid Entry for Number of Samples")
    }
    
    
    
  })  
  
  observeEvent(input[[paste0("sampleInput.", input$nSamples)]], {
    
    if(shouldParseVal$shouldParse == TRUE) {
      
      geneCounts <- getCounts()
      
      for(i in 2:ncol(geneCounts)) {
        
        tmpId <- paste0("sampleInput.", i-1)
        updateTextInput(session, tmpId, value = names(geneCounts)[i])
        
      }
      
      shouldParseVal$shouldParse = FALSE
    }
    
  })
  
  output$groupInputs <- renderUI({
    if(isNumeric(input$nSamples)) {
      if(input$nSamples >= 1) {
        groupInputs <- lapply(1:input$nSamples, function(number){
          
          groupId <- paste0("group.", number)
          groupInputId <- paste0("groupInput.", number)
          div(
            textInput(inputId = groupInputId, label = NULL),
            style = "font-size:14px;"
          )
          
        })
      } 
    }
    
    
    
  })
  
  output$compInputs <- renderUI({
    
    if(isNumeric(input$nComparisons)) {
      if(input$nComparisons >= 1) {
        compInputs <- lapply(1:input$nComparisons, function(number){
          
          compId <- paste0("comp.", number)
          compInputId <- paste0("compInput.", number)
          div(
            textInput(inputId = compInputId, label = NULL),
            style = "font-size:14px;"
          )
          
        })
      } else {
        helpText("Invalid entry for number of comparisons")
      }
      
      
    } else {
      helpText("Invalid entry for number of comparisons")
    }
  })
  
  
  
  
  
  listCompsRender <- c()
  
  nCharts <- reactiveValues(nChartWells = 1)
  checkN <- reactiveValues(checkNumber = 0)
  checkPanel <- reactiveValues(checkPanelId = "comparisons0_order")
  checkTheOld <- reactiveValues(checkOld = 1)
  barId <- reactiveValues(barIds = c("comparisons0"))
  preventionMeasures <- reactiveValues(prevention = 1)
  secondPrevention <- reactiveValues(second = 1)
  thirdPrevention <- reactiveValues(third = 1)
  
  output$compSource <- renderUI({
    if(((input$expEntryMode == 0) && (!is.null(input$fileSelectComp))) || ((input$expEntryMode == 1) && (input$compInput.1 != ""))) {
      if(input$expEntryMode == 0) {
        expDesign <-
          read.table(
            input$fileSelectComp$datapath,
            header = TRUE,
            sep = ",",
            row.names = 1
          )
        listComps <- expDesign[2]
        
        
      } else {
        listComps <- matrix(nrow = input$nComparisons, ncol = 1)
        listComps <- as.data.frame(listComps)
        
        for(x in 1:input$nComparisons) {
          compId <- paste0("compInput.", x)
          
          listComps[x,1] <- input[[compId]]
        }
      }
      colnames(listComps) <- "Comparison"
      listCompsUp <- listComps
      listCompsUp$Comparison <- paste0(listCompsUp$Comparison, " Up")
      listCompsDown <- listComps
      listCompsDown$Comparison <- paste0(listCompsDown$Comparison, " Down")
      
      listCompsUp <- listCompsUp[listCompsUp$Comparison != " Up",]
      listCompsUp <- as.list(listCompsUp)
      listCompsDown <- listCompsDown[listCompsDown$Comparison != " Down",]
      listCompsDown <- as.list(listCompsDown)
      
      
      listCompsFinal <- c(listCompsUp, listCompsDown)
      listCompsRender <- listCompsFinal
      for(j in 1:length(listCompsFinal)) {
        names(listCompsFinal)[j] <- listCompsFinal[j]
      }
      
      wellPanel(id = "sourcePanel",
                orderInput('source', div(strong('Source'), style="font-size: 200%;"), items = listCompsRender, as_source = TRUE, connect = barId$barIds)
      )
    } else {
      helpText("You need to input comparison data")
    }
  })
  
  output$compSamples <- renderUI({
    if(((input$expEntryMode == 0) && (!is.null(input$fileSelectComp))) || ((input$expEntryMode == 1) && (input$compInput.1 != ""))) {
      div(id = "compsPanel0", style = "border-style: solid; border-width: 2px; border-color: lightgrey;",
          orderInput('comparisons0', div(strong('Comparison 1'), style="font-size: 125%;"), items = NULL, placeholder = "Drag comparisons here", connect = "source")
      )
      
    }
  })
  
  
  necComps <- reactiveValues(necComp = list())
  
  
  
  observeEvent({input[[checkPanel$checkPanelId]]}, {
    
    
    
    if(0 == 1) {
      nCharts$nChartWells <- 1
      checkN$checkNumber <- 0
      checkPanel$checkPanelId <- "comparisons0_order"
      barId$barIds <- c("comparisons0")
    }
    
    if(checkN$checkNumber == 0) {
      necComps$necComp[[1]] <- list(input[[("comparisons0_order")]])
    }
    
    
    if(thirdPrevention$third == 1) {
      if(preventionMeasures$prevention == 1) {
        
        
        if(length(necComps$necComp[[max(1, checkN$checkNumber)]]) != 0) {
          
          for(x in 0:checkN$checkNumber) {
            necComps$necComp[[x+1]] <- list(input[[(paste0("comparisons", x, "_order"))]])
          }
          
          if(checkN$checkNumber == 0) {
            necComps$necComp <- necComps$necComp[[1]]
          }
          
          if(length(necComps$necComp[[checkN$checkNumber+1]]) != 0) {
            
            
            
            panelId <- paste0("compsPanel", checkN$checkNumber + 1)
            ordInput <- paste0("comparisons", checkN$checkNumber + 1)
            selectorId = paste0("#compsPanel", (checkN$checkNumber))
            
            
            barId$barIds[nCharts$nChartWells+1] <- paste0("comparisons", nCharts$nChartWells)
            
            nCharts$nChartWells <- nCharts$nChartWells + 1
            
            
            
            
            insertUI(selector = selectorId, where = "afterEnd",
                     ui = div(id = panelId, style = "border-style: solid; border-width: 2px; border-color: lightgrey;",
                              orderInput(ordInput, label = div(strong(paste0("Comparison ", checkN$checkNumber + 2)), style="font-size: 125%;"), items = NULL, placeholder = "Drag comparisons here", connect = "source")
                     )
            )
            
            
            
            
            #necComps$necComp <- c(input[["comparisons0_order"]], necComps$necComp)
            
            checkN$checkNumber <- length(barId$barIds) - 1
            
            updateCheck = paste0("comparisons", checkN$checkNumber, "_order")
            checkPanel$checkPanelId = updateCheck
            if (secondPrevention$second == 1) {
              checkTheOld$checkOld <- checkN$checkNumber
            }
            
            
            if(secondPrevention$second == 0 && checkTheOld$checkOld == checkN$checkNumber) {
              finalDestruction()
            }
          } 
        } 
      } else {
        preventionMeasures$prevention <- 1
      }
    }  else {
      
      updateCheckboxInput(session, "expEntryMode", value = 0)
      updateCheckboxInput(session, "expEntryMode", value = 1)
      
      thirdPrevention$third <- 1
      
      
    } 
    
  })
  
  finalDestruction <- function() {
    for(x in 1:checkN$checkNumber) {
      paneToRemove = paste0("#compsPanel", x)
      slotToRemove = paste0("#comparisons", x)
      idsToRemove = paste0("#comparisons", x, "_order")
      removeUI(selector = slotToRemove, multiple = TRUE, immediate = FALSE)
      removeUI(selector = paneToRemove, multiple = TRUE, immediate = FALSE)
      
    }
    secondPrevention$second <- 1
    thirdPrevention$third <- 0
    nCharts$nChartWells <- 1
    checkN$checkNumber <- 0
    barId$barIds <- "comparisons0"
    checkPanel$checkPanelId <- "comparisons0_order"
    
    
    
  }
  
  
  
  observeEvent({input$clearBars}, {
    
    necComps$necComp <- list()
    
    for(x in 1:checkN$checkNumber) {
      paneToRemove = paste0("#compsPanel", x)
      slotToRemove = paste0("#comparisons", x)
      idsToRemove = paste0("#comparisons", x, "_order")
      removeUI(selector = slotToRemove, multiple = TRUE, immediate = TRUE)
      removeUI(selector = paneToRemove, multiple = TRUE, immediate = TRUE)
    }
    
    
    removeUI(selector = "#comparisons0", multiple = TRUE, immediate = TRUE)
    removeUI(selector = "#compsPanel0", multiple = TRUE, immediate = TRUE)
    
    insertUI(selector = "#sourcePanel", where = "afterEnd", ui =
               div(id = "compsPanel0", style = "border-style: solid; border-width: 2px; border-color: lightgrey;",
                   orderInput('comparisons0', div(strong('Comparison 1'), style="font-size: 125%;"), items = NULL, placeholder = "Drag comparisons here", connect = "source")))
    
    
    
    
    
    nCharts$nChartWells <- 1
    checkN$checkNumber <- 0
    barId$barIds <- "comparisons0"
    checkPanel$checkPanelId <- "comparisons0_order"
    
    secondPrevention$second <- 0
    
  })
  
  
  output$submitButton <- downloadHandler(
    
    
    
    filename = function() {
      paste0(input$projectName, "Analysis.zip.zip")
      
    },
    
    content = function(file) {
      
      rv$download_flag <- rv$download_flag + 1
      
      if(rv$download_flag == 1) {
        shinyjs::toggle("downloadMessage")
      } else if(rv$download_flag > 1) {
        shinyjs::toggle("downloadMessage")
        shinyjs::toggle("downloadComplete")
      }
      
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      files <- NULL;
      
      #Read in counts file
      geneCounts <- getCounts()

      
      names(geneCounts)[1] <- "X"
      
      #Read in comparisons matrix
      if(input$expEntryMode == 0) {
        expDesign <-
          read.table(
            input$fileSelectComp$datapath,
            header = TRUE,
            sep = ",",
            row.names = 1
          )
      } else if(input$expEntryMode == 1) {
        expDesign <- makeDesign()
        
      }
      
      if(input$barPDF == 1) {
          makePDFs = TRUE
        } else {
          makePDFs = FALSE
        }
      
      if(input$barGraphs == 1) {
        
        if(input$autoBars == 1) {
          barDesign <- read.table(
            input$barUpload$datapath,
            header = FALSE,
            sep = ",",
            row.names = NULL
          )
          barDesign <- barDesign[2:nrow(barDesign),]
          
        } else if(input$autoBars == 0) {
          barDesign <- makeBarCSV()
          barDesign <- as.data.frame(barDesign)
          barDesign <- sapply(barDesign, as.character)
          barDesign[is.na(barDesign)] <- ""
          barDesign <- noquote(barDesign)
        }
        
        expDesign <- as.data.frame(expDesign)
        
        minCutoff <- input$minCutoff
        
        if(input$barPDF == 1) {
          makePDFs = TRUE
        } else {
          makePDFs = FALSE
        }
        
      } 
        
      
      if(input$mysqlUpload == 1) {
        
        
        pool <- dbPool(drv = RMariaDB::MariaDB(), user = "root", password = "wellslab123",
                       dbname = "shinyApp", host = "localhost", port = 3306)
        
        #Get total number of projects
        query <- "SELECT COUNT(*) FROM shinyApp.ProjectDetails;"
        nProjects <- dbGetQuery(pool, query)
        #Flatten counts file and experimental design file
        countStore <- c()
        countStore[[1]] <- colnames(geneCounts)
        countStore[[1]] <- paste(unlist(countStore[[1]]), collapse = ";")
        
        countStore[[2]] <- geneCounts[,1]
        countStore[[2]] <- paste(unlist(countStore[[2]]), collapse = ";")
        
        countStore[[3]] <- geneCounts[,2:ncol(geneCounts), drop = TRUE]
        countStore[[3]] <- paste(unlist(countStore[[3]]), collapse = ";")
        
        countStoreFlat <- paste(countStore[[1]], countStore[[2]], countStore[[3]], collapse = ",")
        
        #countStoreFlat <- c()
        #countStoreFlat[[1]] <- countStore
        
        expStore <- c()
        expStore[[1]] <- row.names(expDesign)
        expStore[[1]] <- paste(unlist(expStore[[1]]), collapse = ";")
        
        expStore[[2]] <- expDesign[,1]
        expStore[[2]] <- paste(unlist(expStore[[2]]), collapse = ";")
        
        expStore[[3]] <- expDesign[,2]
        expStore[[3]] <- paste(unlist(expStore[[3]]), collapse = ";")
        
        expStoreFlat <- paste(expStore[[1]], expStore[[2]], expStore[[3]], collapse = ",")
        
        #Prepare project details to be uploaded into mysql table
        projectDetails <- c((as.character(nProjects[1,1]+1)),input$projectID,
                            input$projectDescr, as.character(input$projectDate[[1]]), countStoreFlat, expStoreFlat)
        projectDetails <- data.frame(projectDetails)
        projectDetails <- data.table::transpose(projectDetails)
        colnames(projectDetails) <- c("ProjectID", "ProjectAlias", "ProjectDescr", "ProjectDate", "Counts", "ExpDesign")
        
        #Add project details to database
        pool::dbWriteTable(pool, name = "ProjectDetails", value = projectDetails, append = TRUE)
        
        poolClose(pool)
      }
      
      
      files = c()
      allFilesDir <- paste0(input$projectName, "allFiles")
      unlink(allFilesDir, recursive = TRUE)
      
      #Create output directory
      dir.create(allFilesDir)
      
      #Check which outputs to produce
      if(("DE" %in% input$outputParameters) || ("VP" %in% input$outputParameters)) {
        
        
        unlink("tmpDir1", recursive = TRUE) 
        
        
        if("VP" %in% input$outputParameters) {
          makeVolc <- TRUE
        } else {
          makeVolc <- FALSE
        }
        
        
        diffExpFiles = c(DifferentialExpression(input$species,input$nomenclature,geneCounts, expDesign, input$minCutoff, makeVolc, input$mysqlUpload), files)
        
        if("DE" %in% input$outputParameters) {
          dir.create(paste0(allFilesDir, "/DifferentialExpression"))  
          file.move(diffExpFiles[[1]], paste0(allFilesDir, "/DifferentialExpression"))
        }
        if("VP" %in% input$outputParameters) {
          dir.create(paste0(allFilesDir, "/VolcanoPlots"))
          file.move(diffExpFiles[[2]], paste0(allFilesDir, "/VolcanoPlots"))
        }
        
        if(input$clusterCheck == 1) {
          #Generate CSVs of Clusters
          clusterCSVs <- c()
          dir.create(paste0(allFilesDir, "/ClusteringEffects"))
          clusterCSVs <- ClusterInfo(diffExpFiles[[1]], input$maxClusterDistance, input$minClusterLogFC, input$maxClusterFDR, input$minClusterGenes)
          file.move(clusterCSVs, paste0(allFilesDir, "/ClusteringEffects"))
          
          
          #Generate Histograms of Clusters
          clusterHists <- c()
          dir.create(paste0(allFilesDir, "/ClusterHistograms"))
          clusterData <- ClusteringHistograms(clusterCSVs)
          
          file.move(clusterData[[1]], paste0(allFilesDir, "/ClusterHistograms"))
          
          #Generate Chromosome Plots
          clusterChroms <- c()
          dir.create(paste0(allFilesDir, "/ClusterChromosomePlots"))
          clusterChroms <- ClusterChromosomeCharts(clusterCSVs, clusterData[[2]])
          file.move(clusterChroms, paste0(allFilesDir, "/ClusterChromosomePlots"))
        }
      }
        
      if("PA" %in% input$outputParameters) {
        
        if(input$barGraphs == 0) {
          barDesign = ""
        } else {
          barCharts <- c()
          ieCSVs = c()
        }
        
        
        for(x in 1:length(input$databaseSelect)) {
          pathwayFiles = c()
          
          paAllFiles <- PathwayAnalysis(input$species,input$nomenclature,geneCounts, expDesign, input$databaseSelect[x], input$minCutoff, barDesign, makePDFs)
          
          pathwayFiles = c(paAllFiles[[1]], pathwayFiles)
          
          dir.create(paste0(allFilesDir, "/PathwayAnalysis"))
          
          file.move(pathwayFiles, paste0(allFilesDir, "/PathwayAnalysis"))
          
          if(input$barGraphs == 1) {
            browser()
            dir.create(paste0(allFilesDir, "/InclusionExclusion"))
            if(length(paAllFiles) == 3) {
            barCharts <- c(paAllFiles[[3]], barCharts)
            ieCSVs <- c(paAllFiles[[2]], ieCSVs)
            } else {
              barCharts <- c(paAllFiles[[2]], barCharts)
            }
            file.move(barCharts, paste0(allFilesDir, "/InclusionExclusion"))
          }
          
          if(makePDFs == TRUE) {
            if(length(paAllFiles) == 3) {
            dir.create(paste0(allFilesDir, "/IECSVs"))
            file.move(ieCSVs, paste0(allFilesDir, "/IECSVs"))
            }
          }
          
        }
        
      }
      if("PCA" %in% input$outputParameters) {
        pcaPlot <- MakePCA(geneCounts, expDesign)
        dir.create(paste0(allFilesDir, "/PCA"))
        file.move(pcaPlot, paste0(allFilesDir, "/PCA"))
      }
      
      
      if(input$wordClouds == 1 || input$overlapCloud == 1) {
        
        
        if(input$overlapCloud == 1) {
          overlapClouds <- c()
          overlapClouds <- MakeWordClouds(blackList, Modifiers, ieCSVs, TRUE)
          dir.create(paste0(allFilesDir, "/InclusionExclusionClouds"))
          file.move(overlapClouds, paste0(allFilesDir, "/InclusionExclusionClouds"))
        }
        
        if(input$wordClouds == 1) {
          wordClouds <- c()
          wordClouds <- MakeWordClouds(blackList, Modifiers, pathwayFiles, FALSE)
          dir.create(paste0(allFilesDir, "/Wordclouds"))
          file.move(wordClouds, paste0(allFilesDir, "/Wordclouds"))
        }
        
      }
      
      zip(file, allFilesDir)
      shinyjs::toggle("downloadMessage")
      shinyjs::toggle("downloadComplete")
      
      
    },
    
    
    
    contentType = "application/zip"
    
    
  )
  
  output$designDownload <- downloadHandler(
    
    filename <- function() {
      
      paste0("ExpDesign.csv")
    },
    
    content <- function(file) {
      expArray <- makeDesign()
      browser()
      write.csv(expArray, file, row.names = TRUE, na = "")
    }
  )
  
  
  output$searchInputsDownload <- downloadHandler(
    
    filename <- function() {
      
      paste0("ProjectInputs.zip.zip")
      
    },
    
    content <- function(file) {
      
      if(is.null(input$searchData_row_last_clicked)) {
        inputSelect <- 1
      } else {
        inputSelect <- input$searchData_row_last_clicked
      }
      
      pool <- dbPool(drv = RMariaDB::MariaDB(), user = "root", password = "wellslab123",
                     dbname = "shinyApp", host = "localhost", port = 3306)
      
      query <- paste0("SELECT Counts, ExpDesign 
                  FROM shinyApp.ProjectDetails s
                  WHERE s.ProjectID = ", inputSelect, ";")
      
      projectInputs <- dbGetQuery(pool, query)
      
      query <- paste0(paste0("SELECT ProjectAlias 
                  FROM shinyApp.ProjectDetails s
                  WHERE s.ProjectID = ", inputSelect, ";"))
      
      projectName <- dbGetQuery(pool, query)
      
      df <- read.table(text = projectInputs[,1], sep = " ")
      tmpString <- as.character(df[1,3])
      countsGen <- read.table(text = tmpString, sep = ";")
      countsGen <- as.data.frame(matrix(countsGen,
                                        ncol = length(strsplit(as.character(df[1,1]), ";")[[1]])-1, byrow = FALSE))
      
      tmpString <- as.character(df[1,1])
      countsGen$X <- 1:nrow(countsGen)
      countsGen <- countsGen[,c(ncol(countsGen),1:(ncol(countsGen)-1))]
      colnames(countsGen) <- as.character(unlist(read.table(text = tmpString, sep = ";")))
      countsGen <- countsGen[,2:ncol(countsGen)]
      
      tmpString <- as.character(df[1,2])
      countsGen$X <- as.character(unlist(read.table(text = tmpString, sep = ";")))
      countsGen <- countsGen[,c(ncol(countsGen),1:(ncol(countsGen)-1))]
      
      df <- read.table(text = projectInputs[,2], sep = " ")
      tmpString <- as.character(df[1,1])
      expGen <- as.character(unlist(read.table(text = tmpString, sep = ";")))
      expGen <- as.data.frame(matrix(expGen, ncol = 1))
      colnames(expGen) <- "Sample"
      tmpString <- as.character(df[1,2])
      expGen$Group <- as.character(unlist(read.table(text = tmpString, sep = ";")))
      tmpString <- as.character(df[1,3])
      expGen$Comparison <- as.character(unlist(read.table(text = tmpString, sep = ";",
                                                          stringsAsFactors = FALSE)))
      
      projInputFiles <- paste0("Inputs")
      if(exists(projInputFiles)) {
        unlink(projInputFiles, recursive = TRUE)
      }
      
      #Create output directory
      dir.create(projInputFiles)
      
      #Set output file names
      countsFileName <- paste0(projectName, "Counts.csv")
      expDesignName <- paste0(projectName, "ExpDesign.csv")
      
      #Create output files
      fwrite(countsGen, countsFileName)
      fwrite(expGen, expDesignName)
      
      file.move(c(countsFileName, expDesignName), projInputFiles)
      
      zip(file, projInputFiles)
      
      poolClose(pool)
      
    },
    
    contentType = "application/zip"
    
  )
  
  
  makeDesign <- function() {
    browser()
    expArray <- matrix(ncol = 3, nrow = input$nSamples)
    for(i in 1:2) {  
      for(x in 1:input$nSamples) {
        if(i == 1) {
          sampleId <- paste0("sampleInput." ,x)
          expArray[x,i] <- input[[sampleId]]
          
        } else if(i == 2) {
            groupId <- paste0("groupInput.", x)
          expArray[x,i] <- input[[groupId]]
        } 
        
      }
    }
    
    for(x in 1:input$nComparisons) {
      compId <- paste0("compInput.", x)
      expArray[x,3] <- input[[compId]]
    }
    
    colnames(expArray) <- c("Sample", "Group", "Comparison")
    
    row.names(expArray) <- expArray[,1]
    expArray <- expArray[,2:ncol(expArray)]
    return(expArray)
    
    
  }
  
  output$makeBarComps <- downloadHandler(
    
    filename <- function() {
      paste0("InclusionExclusionCharts.csv")
    },
    
    content <- function(file) {
      barArray <- makeBarCSV()
      
      write.csv(barArray, file, row.names = FALSE, na = "")
      
    }
    
  )
  
  output$bamCounts <- downloadHandler(
    
    filename <- function() {
      
      paste0("countsFromBam.csv")
    },
    
    content <- function(file) {
    geneCounts <- getCounts()
    write.csv(geneCounts, file, row.names = TRUE, na = "")
      
      
    }
    
    
  )
  
  
  makeBarCSV <- function() {
    
    for(x in 0:checkN$checkNumber) {
      necComps$necComp[[x+1]] <- list(input[[(paste0("comparisons", x, "_order"))]])
    }
    
    lengthChecks <- list()
    for(x in 1:length(necComps$necComp)) {
      lengthChecks[[x]] <- length(unlist(necComps$necComp[[x]]))
    }
    
    barArray <- matrix(ncol = length(necComps$necComp)-1, nrow = max(unlist(lengthChecks)));
    
    for(i in 1:ncol(barArray)) {
      
      for(j in 1:nrow(barArray)) {
        barArray[j,i] <- unlist(necComps$necComp[[i]][[1]])[j];
    }
    }
    return(barArray)
    
    
  }
  
  
  parseSamples <- function() {
    
    
    geneCounts <- getCounts()
    
    updateNumericInput(session, "nSamples", value = ncol(geneCounts)-1)
    
   
  }
  
  
  
  clearSamples <- function() {
    
    updateNumericInput(session, 'nSamples', value = 1)
    updateTextInput(session, "sampleInput.1", value = "")
    
    
  }
  
  alreadyDone = 0
  tmpCounts <- c()
  getCounts <- function() {
    
    if(input$uploadCore == 0 && input$convertBam == 0) {
      geneCounts <-
        read.table(
          input$fileSelectExpr$datapath,
          sep = ",",
          header = TRUE,
          row.names = NULL
        )
    } else if(input$convertBam == 1) {
        #Convert Bam files to count files
      

        if(alreadyDone == 1) {
          geneCounts <- tmpCounts$counts
          geneCounts <- data.frame(geneCounts)
          geneCounts$X <- row.names(geneCounts)
          geneCounts <- geneCounts[,c(ncol(geneCounts),(1:ncol(geneCounts)-1))]
        } else {
        if(input$species == "hsapiens_gene_ensembl") {
        tmpCounts <<- featureCounts(input$bamInputs$datapath, annot.inbuilt = "hg38",
                                   isPairedEnd = TRUE,
                                   annot.ext = "/srv/shiny-server/guiWork/shinyApp/data/gencode.v31.annotation.gtf",
                                   isGTFAnnotationFile = TRUE)
        } else {
          tmpCounts <<- featureCounts(input$bamInputs$datapath, annot.inbuilt = "mm10",
                                     isPairedEnd = TRUE,
                                     annot.ext = "/srv/shiny-server/guiWork/shinyApp/data/gencode.vM22.annotation.gtf",
                                     isGTFAnnotationFile = TRUE)
        }
        colnames(tmpCounts$counts) <<- file_path_sans_ext(input$bamInputs$name)
        geneCounts <- tmpCounts$counts
        geneCounts <- data.frame(geneCounts)
        geneCounts$X <- row.names(geneCounts)
        geneCounts <- geneCounts[,c(ncol(geneCounts),(1:ncol(geneCounts)-1))]
        alreadyDone <<- 1
      }
    } else {
      rawSheet <- read_excel(input$fileSelectExpr$datapath, sheet = "DE_analysis")
      
      firstTestCol <- which(colnames(rawSheet) == "Raw counts")
      colsToKeep <- c(1,firstTestCol:ncol(rawSheet))
      rawSheet <- rawSheet[,colsToKeep]
      colnames(rawSheet) <- rawSheet[1,]
      rawSheet <- rawSheet[2:nrow(rawSheet),]
      colnames(rawSheet)[1] <- "X"
      rawSheet <- as.data.frame(rawSheet)
      geneCounts <- rawSheet
      
      }
    return(geneCounts)
  }
  
  output$trackDownloadButton <- downloadHandler(
   
    filename <- function() {
      paste0(input$dataSearch, "Track.png")
    }, 
    
    content <- function(file) {
      
      svgName <- paste0(input$dataSearch, "Track.png")
      file.copy(svgName, file)
      file.remove(svgName)
    }
  )
  
  output$chromPlotDownload <- downloadHandler(
    
    filename <- function() {
      paste0("ChromoPlot.png")
    },
    
    content <- function(file) {
      shinyjs::toggle("chromMessage")
      
      chromData <- read.table(input$chromPlotInput$datapath, sep = ",")
      chromoPlot <- ChromPlot(chromData)
      shinyjs::toggle("chromComplete")
      file.copy(chromoPlot, file)
      file.remove(chromoPlot)
    }
  )
  
  output$corDownload <- downloadHandler(
    
    filename <- function() {
      paste0("CorAnalysis.zip")
    },
    
    content <- function(file) {
      shinyjs::toggle("essMessage")
      unlink("EssentialityAnalysis", recursive = TRUE)
      corData <- EssentialityClusters(input$eProjInput, input$eCompInput, input$eGeneInput, input$eCellLineInput)
      corCloud <- MakeWordClouds(blackList, Modifiers, corData[[1]], FALSE)
      #Process complete, generate necessary files
      shinyjs::toggle("essComplete")
      shinyjs::toggle("essMessage")
      dir.create("EssentialityAnalysis")
      file.move(corData, "EssentialityAnalysis")
      file.move(corCloud, "EssentialityAnalysis")
      zip(file, "EssentialityAnalysis")
      
    },
    
    contentType = "application/zip"
    
    
    
  )
  

}

shinyApp(ui, server)
