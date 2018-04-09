############################################################################
# Copyright CNRS 2017
# Contributor : Marie Locard-Paulet (21/11/2017) [marie.locard@ipbs.fr]
# This software is a computer program whose purpose is to visualize and inspect deconvoluted MS 2D maps.
# This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software. You can  use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info". 
# As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive licensors have only limited liability. In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or reproducing the software by the user in light of its specific status of free software,that may mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data to be ensured and, more generally, to use and operate it in the  same conditions as regards security. 
# The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
############################################################################



# Packages:
############################################################################

library(shiny)
# devtools::install_github('hadley/ggplot2')
library(ggplot2)
library(plotly)
library(dplyr)
library(RColorBrewer)

library(shinyBS)


############################################################################
# Functions:
############################################################################

RBindList <- function(l) {
  # combine tables of a list l using rbind
  if (length(l)==1) {tab <- l[[1]]} else {
    if (class(l[[1]])!="matrix" | class(l[[1]])!="data.frame") {
      tab <- rbind(l[[1]], l[[2]])
      if (length(l)>2) {
        for (i in 3:length(l)) {
          tab <- rbind(tab, l[[i]])
        }
      } 
    } else {
      l <- l[2:length(l)]
      RBindList(l)
    }
  }
  return(tab)
}

RenameRoWinPro <- function(tab) {
  # Filter and rename the tables RoWinPro style
  names(tab)[3] <- "intensity"
  names(tab)[2] <- "Mass"
  names(tab)[1] <- "RT"
  tab
}

RenameBioPharma <- function(tab) {
  # Filter and rename the tables BioPharma style
  vec <- names(tab)
  vec[3] <- "intensity"
  vec[2] <- "Mass"
  vec[1] <- "RT"
  vec[4] <- "PeakStart"
  vec[5] <- "PeakStop"
  names(tab) <- vec
  return(tab)
}

ThresholdCleaning <- function(l, threshold) {
  # removes the points with lowest intensity according to the threshold chosen by the user:
  threshold <- threshold/100
  l1 <- list()
  for (x in seq_along(l)) {
    # Filter the intensities according to threshold:
    temp <- cbind(l[[x]], "File" = rep(names(l)[x], nrow(l[[x]])))
    temp <- temp[order(temp[,3], decreasing = T),]
    temp <- temp[!is.na(temp[,3]),]
    thresh <- floor(threshold * nrow(temp))
    temp <- temp[c(1:thresh),]
    l1[[x]] <- temp
  }
  names(l1) <- names(l)
  return(l1)
}
############################################################################

# App:
############################################################################
## UI:
############################################################################

ui <- fluidPage(
  titlePanel("VisioProt-MS"),
  #checkboxInput("MSModeCheck", "MS data only", TRUE), # to switch from MS data to MS2 mode
  radioButtons("MSModeCheck", "MS mode:",
               c("MS data only" = 'MS',
                 "MS2 overlay" = 'MS2'),
               selected = 'MS',
               inline = TRUE
  ), # to switch from MS data to MS2 mode
  bsTooltip("MSModeCheck", 
            "Choose between plotting MS data only or overlaying results of Top-Down searches",
            "right"),
  # Control side bar:
  sidebarLayout( 
    sidebarPanel( 
      # Conditional panels:
      # Part 1:
      conditionalPanel(condition="input.MSModeCheck== 'MS'",
                       # File selection:
                       fluidRow(
                         column(11,
                                fileInput("fileMS", "Select input file(s):",  
                                          accept = c(
                                            "text/csv",
                                            "text/comma-separated-values,text/plain",
                                            ".txt",
                                            ".ms1ft",
                                            ".csv"),
                                          multiple = T,
                                          width = "100%"
                                )),
                         column(1, style = "margin-top: 30px",
                                tipify(bsButton("bsfileMS", "?", style = "default", size = "extra-small"),
                                       "Upload a file with deconvoluted MS data for plotting",
                                       placement = "right")
                         )),
                       checkboxInput("TestModeCheck", "Using test mode", FALSE),
                       bsTooltip("TestModeCheck", 
                                 "Check to test the application without uploading any file. Then click on a button to upload a single test file or several for overlay",
                                 "right"), # to switch from user data to test mode
                       # Modifying output when passing in test mode:
                       conditionalPanel(condition="input.TestModeCheck==true",
                                        actionButton("TestFile1", "Single test file"),
                                        actionButton("TestFile2", "Multiple test files")
                       )
      ),
      # Part 2:
      conditionalPanel(condition="input.MSModeCheck== 'MS2'", 
                       checkboxInput("MS2TestModeCheck", "Using test mode", FALSE),
                       bsTooltip("MS2TestModeCheck", 
                                 "Check to test the application without uploading any file",
                                 "right"), # to switch from user data to test mode
                       # Modifying output when passing in test mode:
                       conditionalPanel(condition = "input.MS2TestModeCheck==true",
                                        tags$span(style="color:red", "You are in test mode. Uncheck to exit."),
                                        br()
                       ),
                       conditionalPanel(condition = "input.MS2TestModeCheck==false",
                                        fluidRow(
                                          column(11,
                                                 fileInput("fileMS2", "Select input file for MS:",  
                                                           accept = c(
                                                             "text/csv",
                                                             "text/comma-separated-values,text/plain",
                                                             ".csv",
                                                             ".txt",
                                                             ".ms1ft"),
                                                           multiple = F,
                                                           width = "100%"
                                                 )),
                                          column(1, style = "margin-top: 30px", 
                                                 tipify(bsButton(inputId = "fileMS2", label = "?", style = "default", size = "extra-small"),
                                                        title = "Upload a file with deconvoluted MS data for plotting with the corresponding results of a top-down analysis",
                                                        options = NULL)
                                          )),
                                        radioButtons("PDPFModeCheck", "Origin of the MS2 files:",
                                                     c("Proteome Discoverer" = 'PD',
                                                       "MSPathFinder" = 'PF'),
                                                     selected = 'PD',
                                                     inline = TRUE
                                        ),
                                        bsTooltip("PDPFModeCheck", 
                                                  "Choose the software utilized for analysing of the top-down data",
                                                  "right"),
                                        conditionalPanel(condition = "input.PDPFModeCheck== 'PD'", 
                                                         fluidRow(
                                                           column(11,
                                                                  fileInput("MS2file", "Choose MS2 File:", 
                                                                            accept = c(
                                                                              "text/csv",
                                                                              "text/comma-separated-values,text/plain",
                                                                              ".txt")
                                                                  )),
                                                           column(1, style = "margin-top: 30px",
                                                                  tipify(bsButton("MS2file", "?", style = "default", size = "extra-small"),
                                                                         "Upload the corresponding MSMS file from Proteome Discoverer",
                                                                         placement = "right")
                                                           )),
                                                         fluidRow(
                                                           column(11,
                                                                  fileInput("PSMfile", "Choose PSM File:", 
                                                                            accept = c(
                                                                              "text/csv",
                                                                              "text/comma-separated-values,text/plain",
                                                                              ".txt")
                                                                  )),
                                                           column(1, style = "margin-top: 30px",
                                                                  tipify(bsButton("PSMfile", "?", style = "default", size = "extra-small"),
                                                                         "Upload the corresponding PSMs file from Proteome Discoverer",
                                                                         placement = "right")
                                                           ))),
                                        conditionalPanel(condition = "input.PDPFModeCheck== 'PF'", 
                                                         fluidRow(
                                                           column(11,
                                                                  fileInput("MS2filePF", "Choose IcTarget or IcTda File from MSPathFinder:", 
                                                                            accept = c(
                                                                              "text/csv",
                                                                              "text/comma-separated-values,text/plain",
                                                                              ".tsv")
                                                                  )),
                                                           column(1, style = "margin-top: 30px",
                                                                  tipify(bsButton("MS2filePF", "?", style = "default", size = "extra-small"),
                                                                         "Upload the corresponding IcTarget or IcTda file from MSPathFinder",
                                                                         placement = "right")
                                                           ))
                                        ),
                                        selectInput("SelectProt", "Select the ID to highlight:", 
                                                    NULL,
                                                    multiple = TRUE),
                                        bsTooltip("SelectProt", 
                                                  "Select among the identified proteins which one(s) to highlight on the plot",
                                                  "right"),
                                        checkboxInput("HideMSMS", "Hide MSMS withouth ID", FALSE),
                                        bsTooltip("HideMSMS", 
                                                  "Removes the MSMS spectra from the top-down analysis that were not matched to a protein.",
                                                  "right"),
                                        checkboxInput("MSTrace", "Display the MS trace", TRUE),
                                        bsTooltip("MSTrace", 
                                                  "Adds the MS trace to the plot.",
                                                  "right")
                       )
      ),
      checkboxInput("DataPoints", "Show data labels (slower)", FALSE), # To switch between ggplot and plotly.
      bsTooltip("DataPoints", 
                "Switch to \"data\" mode: data appears on hovering",
                "right"),
      # Parameters for the plot:
      fluidRow(
        column(5,
               # Selection of the colour scales. This depends on the number of input files:
               # With updateSelectInput:
               selectInput("colourscale", "Colour scale:", # for continuous scales
                           c("Spectral" = "Spectral",
                             "Red/yellow/blue" = "RdYlBu", 
                             "Red/yellow/green" = "RdYlGn",
                             "yellow to red" = "YlOrRd"
                           )),
               bsTooltip("colourscale", 
                         "Select the colour scale for the MS data.",
                         "right")),
        column(3,
               numericInput("pch", label = "Point size:", value = 1, min = 0.1, step = 0.1, max = 10),
               bsTooltip("pch", 
                         "Define the size of the point (from 0.1 to 10).",
                         "right")),
        column(4,
               numericInput("IntensityThresh", label = "Threshold:", value = 20, min = 0, max = 100, step = 1),
               bsTooltip("IntensityThresh", 
                         "Define the percentage of highest intensity features of the MS data to display.",
                         "right"))
      ),
      # Information regarding how to zoom (depends on the plotting type):
      htmlOutput("ZoomParam"),
      br(),
      actionButton("DeZoom", "Unzoom one step", style='padding:8px; font-size:150%'),
      bsTooltip("DeZoom", 
                "Unzoom to previous window (only once).",
                "right"),
      actionButton("TotalDeZoom", "Total unzoom", style='padding:8px; font-size:150%'),
      bsTooltip("TotalDeZoom", 
                "Total unzoom.",
                "right"),
      br(),
      br(),
      # Buttons for download:
      downloadButton("Download", "Download .pdf"),
      downloadButton("Download1", "Download .png"),
      downloadButton("Download2", "Download .svg"),
      br(),
      a("Help", href="Help/VisioProtHelp.html", target="blank") # Access to help
    ),
    # Main panel for plotting (output different in function of the checkbox DataPoints):
    mainPanel(
      uiOutput("plotUI")
    )
  ),
  # Footer
  tabsetPanel(
    tabPanel(
      HTML('<footer><font size="0.8">copyright 2017 - CNRS - All rights reserved - VisioProt-MS V2.0</font></footer>')
    )
  )
)
############################################################################
## Server:
############################################################################

server <- function(input, output, clientData, session) {
  
  # UI modifications:
  ###################
  # Text output to describe how to zoom (dependent on the checkbox DataPoints):
  output$ZoomParam <- renderUI({
    if (input$DataPoints) {
      HTML("<h4>Zoom in: select the ranges of interest.<br/>Zoom out: Click on the unzoom button.</h4>")
    } else {
      HTML("<h4>Zoom in: select the ranges of interest and double click.<br/>Zoom out: double click.</h4>")
    }
  })
  
  # Change colour scale in function of the number of scales:
  observe({
    if (!is.null(linput())) {
      x <- linput()
      if (x > 1) {
        updateSelectInput(session, "colourscale",
                          "Colour scale:",
                          c("Set1" = "Set1",
                            "Set2" = "Set2",
                            "Set3" = "Set3",
                            "Dark2" = "Dark2", 
                            "Paired" = "Paired",
                            "Accent" = "Accent"
                          ))
      } else {
        updateSelectInput(session, "colourscale",
                          "Colour scale:",
                          c("Spectral" = "Spectral",
                            "Red/yellow/blue" = "RdYlBu", 
                            "Red/yellow/green" = "RdYlGn",
                            "yellow to red" = "YlOrRd"
                          ))
      }
    }
  })
  colval <- reactiveVal()
  observe({
    x <- input$colourscale
    colval(x)
  })
  
  # Add the protein IDs to select to highlight them in the plot:
  observe({
    if (!is.null(filedataMS2())) {
      if (length(filedataMS2()$PSMfile$Master.Protein.Descriptions[!is.na(filedataMS2()$PSMfile$Master.Protein.Descriptions)]) > 0) {
        updateSelectInput(session, "SelectProt",
                          "Select the ID to highlight:",
                          sort(unique(filedataMS2()$PSMfile$Master.Protein.Descriptions[!is.na(filedataMS2()$PSMfile$Master.Protein.Descriptions)]))
        )
      }
    }
    if (!is.null(filedataMS2PF)) { # PathFinder
      if (length(filedataMS2PF()$Protein.Descriptions[!is.na(filedataMS2PF()$Protein.Descriptions)]) > 0) {
        updateSelectInput(session, "SelectProt",
                          "Select the ID to highlight:",
                          sort(unique(filedataMS2PF()$Protein.Descriptions[!is.na(filedataMS2PF()$Protein.Descriptions)]))
        )
      }
    }
  })
  
  ###################
  # When input MS file:
  #####################
  InputFileMS <- reactiveVal(NULL)
  
  observeEvent(input$fileMS, {
    ranges$x <- NULL
    ranges$y <- NULL
    if (input$MSModeCheck == "MS" & !is.null(input$fileMS) & input$TestModeCheck == FALSE & input$MS2TestModeCheck == FALSE) {
      InputFileMS(input$fileMS)
    } else {
      InputFileMS(NULL)
    }
  })
  
  observeEvent(input$fileMS2, {
    if (input$MSModeCheck == "MS2" & !is.null(input$fileMS2)) {
      InputFileMS(input$fileMS2)
    } else {
      InputFileMS(NULL)
    }
  })
  
  #####################
  # When input MS2 file:
  #####################
  InputFilesMS2 <- reactiveVal(NULL) # For BioPharma input
  observeEvent(c(input$MS2file, input$PSMfile), {
    if (input$MSModeCheck == "MS2" & !is.null(input$MS2file) & !is.null(input$PSMfile)) {
      InputFilesMS2(list("MS2file" = input$MS2file, "PSMfile" = input$PSMfile))
    } else {
      InputFilesMS2(NULL)
    }
  })
  InputFilesMS2PF <- reactiveVal(NULL) # For PathFinder input
  observeEvent(input$MS2filePF, {
    if (input$MSModeCheck == "MS2" & !is.null(input$MS2filePF)) {
      InputFilesMS2PF(input$MS2filePF)
    } else {
      InputFilesMS2PF(NULL)
    }
  })
  
  #####################
  # Plotting MS trace:
  #####################
  
  # Number of input file(s) from the same type:
  linput <- reactiveVal()
  
  # Test files input:
  testfileinput <- reactiveVal(0) # 0: no test file; 1: single file; 2: multiple file; 3: MS2 mode test.
  
  filetype <- reactiveValues(RoWinPro = 0, BioPharma = 0, ProMex = 0) # Number of files of each type. Bruker files fall into the "RoWinPro" category once recognised and opened properly.
  
  # test mode / test files:
  observeEvent(input$TestFile1, {
    ranges$x <- NULL
    ranges$y <- NULL
    linput(1)
    testfileinput(1)
    filetype$RoWinPro <- 1
    filetype$BioPharma <- 0
    filetype$ProMex <- 0
    colval("Spectral")
    InputFileMS(NULL)
    InputFilesMS2(NULL)
    InputFilesMS2PF(NULL)
  })
  
  observeEvent(input$TestFile2, {
    ranges$x <- NULL
    ranges$y <- NULL
    linput(4)
    testfileinput(2)
    filetype$RoWinPro <- 4
    filetype$BioPharma <- 0
    filetype$ProMex <- 0
    colval("Set1")
    InputFileMS(NULL)
    InputFilesMS2(NULL)
    InputFilesMS2PF(NULL)
  })
  
  observeEvent(input$MS2TestModeCheck, {
    ranges$x <- NULL
    ranges$y <- NULL
    if (input$MS2TestModeCheck == TRUE) {
      linput(1)
      testfileinput(3)
      filetype$RoWinPro <- 1
      filetype$BioPharma <- 0
      filetype$ProMex <- 0
      colval("Spectral")
    }
    InputFileMS(NULL)
    InputFilesMS2(NULL)
    InputFilesMS2PF(NULL)
  })
  
  # When uploading an MS file in MS mode:
  observeEvent(c(input$fileMS, input$fileMS2), { 
    InputFileMS <- InputFileMS()
    testfileinput(0)
    if (!is.null(InputFileMS)) {
      l <- list()
      l2 <- list()
      l3 <- list()
      for(i in 1:nrow(InputFileMS)){
        l[[i]] <- grepl("Monoisotopic Mass", readLines(InputFileMS[i, 'datapath'])[1]) & grepl("Apex RT", readLines(InputFileMS[i, 'datapath'])[1]) & grepl("Sum Intensity", readLines(InputFileMS[i, 'datapath'])[1]) & grepl("Start Time (min)", readLines(InputFileMS[i, 'datapath'])[1], fixed = T) & grepl("Stop Time (min)", readLines(InputFileMS[i, 'datapath'])[1], fixed = T)  # TRUE if Biopharma
        l2[[i]] <- substr(readLines(InputFileMS[i, 'datapath'])[2], 0, 13) == "Compound Name" # TRUE if Bruker
        l3[[i]] <- grepl(".ms1ft", InputFileMS$name[i], fixed = T) # TRUE if ProMex
      }
      l <- unlist(l)
      l2 <- unlist(l2)
      l3 <- unlist(l3)
      filetype$RoWinPro <- sum(l==F & l3==F) # Bruker files too
      filetype$BioPharma <- sum(l==T & l2==F & l3==F)
      filetype$ProMex <- sum(l==F & l2==F & l3==T)
      #linput(max(as.numeric(table(l))))
      linput(max(as.numeric(c(filetype$RoWinPro, filetype$BioPharma, filetype$ProMex))))
      if (linput() == 1 & length(c(filetype$RoWinPro, filetype$BioPharma, filetype$ProMex)[c(filetype$RoWinPro, filetype$BioPharma, filetype$ProMex)!=0])>1) {
        linput(sum(as.numeric((c(filetype$RoWinPro, filetype$BioPharma, filetype$ProMex)))))
      }
      if (linput() > 1) {
        colval("Set1")
      } else {
        colval("Spectral")
      }
    } else {
      filetype$RoWinPro <- 0
      filetype$BioPharma <- 0
      filetype$ProMex <- 0
    }
  })
  #####################
  # Prevent plotting when modifying modes:
  ########################################
  # Remove plot when getting out of test mode.
  observeEvent(input$TestModeCheck, {
    InputFileMS(NULL)
    testfileinput(0)
    ranges$x <- NULL
    ranges$y <- NULL
    filetype$RoWinPro <- 0
    filetype$BioPharma <- 0
    filetype$ProMex <- 0
  })
  
  # Remove plot when getting out of MS2 test mode.
  observeEvent(input$MS2TestModeCheck, {
    if (input$MS2TestModeCheck==FALSE) {
      InputFileMS(NULL)
      InputFilesMS2(NULL)
      InputFilesMS2PF(NULL)
      testfileinput(0)
      updateSelectInput(session, "SelectProt",
                        "Select the ID to highlight:",
                        c(""))
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
  # Remove plot when getting out of MS or MS2 mode.
  observeEvent(input$MSModeCheck, {
    InputFileMS(NULL)
    testfileinput(0)
    ranges$x <- NULL
    ranges$y <- NULL
  })
  
  # Remove plot when getting out of PD or PF mode.
  observeEvent(input$PDPFModeCheck, {
    InputFilesMS2(NULL)
    InputFilesMS2PF(NULL)
    testfileinput(0)
    ranges$x <- NULL
    ranges$y <- NULL
    updateSelectInput(session, "SelectProt",
                      "Select the ID to highlight:",
                      c(""))
  })
  ########################################
  
  ftype <- reactive({
    if (is.null(InputFileMS()) & testfileinput() == 0) {
      return(NULL)
    } else {
      InputFileMS <- InputFileMS()
      l <- list()
      for(i in 1:nrow(InputFileMS)){
        val <- grepl("Monoisotopic Mass", readLines(InputFileMS[i, 'datapath'])[1]) & grepl("Apex RT", readLines(InputFileMS[i, 'datapath'])[1]) & grepl("Sum Intensity", readLines(InputFileMS[i, 'datapath'])[1]) & grepl("Start Time (min)", readLines(InputFileMS[i, 'datapath'])[1], fixed = T) & grepl("Stop Time (min)", readLines(InputFileMS[i, 'datapath'])[1], fixed = T) # T for BioPharma, F for RoWinPro
        val2 <- substr(readLines(InputFileMS[i, 'datapath'])[2], 0, 13) == "Compound Name" # TRUE if Bruker
        val3 <- grepl(".ms1ft", InputFileMS$name[i]) # TRUE if ProMex
        val <- ifelse(val,  "BioPharma", "RoWinPro")
        val[val2] <- "Bruker"
        val[val3] <- "ProMex"
        # Check it is the correct input format:
        if (val == "BioPharma") { # Find "Apex RT" in BioPharma files
          if (!grepl("Apex RT", readLines(InputFileMS[i, 'datapath'])[1])) {
            val <- "DoNotApply"
            validate(
              need(val!="DoNotApply", "Incorrect input format.\nVisioProt-MS accepts the following input files:\n- outputs from RoWinPro (Gersch et al. 2015).\n- outputs from DataAnalysis 4.2 (Bruker).\n- BioPharma Finder 3.0 (Thermo Fisher Scientific) tables that have been exported at \"Component Level Only\" before being converted in tab-separated files.\n- ProMex exports in \".ms1ft\".")
            )
          }
        }
        if (val == "RoWinPro") { # Four columns in RoWinPro files
          if (length(as.numeric(gregexpr("\t", readLines(InputFileMS[i, 'datapath'])[1])[[1]]))!=3) {
            val <- "DoNotApply"
            validate(
              need(val!="DoNotApply", "Incorrect input format.\nVisioProt-MS accepts the following input files:\n- outputs from RoWinPro (Gersch et al. 2015).\n- outputs from DataAnalysis 4.2 (Bruker).\n- BioPharma Finder 3.0 (Thermo Fisher Scientific) tables that have been exported at \"Component Level Only\" before being converted in tab-separated files.\n- ProMex exports in \".ms1ft\".")
            )
          }
        }
        if (val == "Bruker") { # First line has no "," and second line contains the column header " RT / min"
          if (!(!grepl(",", readLines(InputFileMS[i, 'datapath'])[1]) & grepl(" RT", readLines(InputFileMS[i, 'datapath'])[2]))) {
            val <- "DoNotApply"
            validate(
              need(val!="DoNotApply", "Incorrect input format.\nVisioProt-MS accepts the following input files:\n- outputs from RoWinPro (Gersch et al. 2015).\n- outputs from DataAnalysis 4.2 (Bruker).\n- BioPharma Finder 3.0 (Thermo Fisher Scientific) tables that have been exported at \"Component Level Only\" before being converted in tab-separated files.\n- ProMex exports in \".ms1ft\".")
            )
          }
        }
        if (val == "ProMex") { # Contains the columns MonoMass, ApexIntensity and Min/MaxElutionTime
          if (!grepl("MinElutionTime", readLines(InputFileMS[i, 'datapath'])[1])) {
            val <- "DoNotApply"
            validate(
              need(val!="DoNotApply", "Incorrect input format.\nVisioProt-MS accepts the following input files:\n- outputs from RoWinPro (Gersch et al. 2015).\n- outputs from DataAnalysis 4.2 (Bruker).\n- BioPharma Finder 3.0 (Thermo Fisher Scientific) tables that have been exported at \"Component Level Only\" before being converted in tab-separated files.\n- ProMex exports in \".ms1ft\".")
            )
          }
        }
        l[[i]] <- val
      }
      vec <- unlist(l)
      return(vec)
    }
  })
  
  # Input the data table:
  filedata0 <- reactive({
    #This function is repsonsible for loading in the selected file
    if (testfileinput() == 0) { # no input test file
      validate(
        need((input$TestModeCheck==FALSE & input$MSModeCheck == "MS") | (input$MS2TestModeCheck == FALSE & input$MSModeCheck == "MS2"), "Uncheck test mode before loading a MS file.")
      )
      if (is.null(InputFileMS())) {
        # User has not uploaded a file yet
        return(NULL)
      } else {
        # Warning if trying to plot several types of data AND several files:
        validate(
          need(!(max(table(ftype())) > 1 & length(unique(ftype())) > 1), "Visioprot-MS can compare several files from the same deconvolution tool (RoWinPro, BioPharma Finder, MSPathFinder or Bruker). If you want to compare files from different input types, please only input one file per tool.")
        )
        InputFileMS <- InputFileMS()
        lfiles <- list()
        for(i in 1:nrow(InputFileMS)){
          if (ftype()[i] == "BioPharma") { # If the file is from Thermo BioPharma
            lfiles[[i]] <- read.table(InputFileMS[i, 'datapath'], sep = "\t", header = T)
            lfiles[[i]] <- lfiles[[i]][,c("Apex.RT", "Monoisotopic.Mass", "Sum.Intensity", "Start.Time..min.", "Stop.Time..min.")] # Map the columns as in RoWinPro format, but with apex RT, start and stop instead of all the points of the peak.
          } else if (ftype()[i] == "RoWinPro")  { # RoWinPro output
            lfiles[[i]] <- read.table(InputFileMS[i, 'datapath'], sep = "\t", header = F)
            lfiles[[i]] <- cbind(lfiles[[i]][,1:3], " Temp1" = rep(NA, nrow(lfiles[[i]])), "Temp2" = rep(NA, nrow(lfiles[[i]]))) # add one more column to allow row binding later on 
          } else if (ftype()[i] == "Bruker")  { # Bruker output
            lfiles[[i]] <- read.table(InputFileMS[i, 'datapath'], sep = ",", header = F, skip = 2)
            lfiles[[i]] <- cbind(lfiles[[i]][,2:4], " Temp1" = rep(NA, nrow(lfiles[[i]])), "Temp2" = rep(NA, nrow(lfiles[[i]]))) # add one more column to allow row binding later on 
          } else if (ftype()[i] == "ProMex")  { # ProMex output
            lfiles[[i]] <- read.table(InputFileMS[i, 'datapath'], sep = "\t", header = T)
            lfiles[[i]] <- cbind("RT" = (lfiles[[i]]$MinElutionTime + ((lfiles[[i]]$MaxElutionTime - lfiles[[i]]$MinElutionTime)/2)), lfiles[[i]][,c("MonoMass", "ApexIntensity", "MinElutionTime", "MaxElutionTime")]) # Map the columns as in RoWinPro format, but with start and stop instead of all the points of the peak. I add a first column with the middle of the peak for zooming (plotly_select needs points, not ranges).
          }
        }
        names(lfiles) <- InputFileMS()$name
        return(lfiles)
      }
    } else if (testfileinput() == 1) { # single file test
      infile <- list.files("files/Unique/", pattern = ".csv", full.names = T)
      lfiles <- list()
      for(i in 1){
        lfiles[[i]] <- read.table(infile[i], sep = "\t", header = F)
      }
      names(lfiles) <- c("test data")
      return(lfiles)
      testfileinput(0)
    } else if (testfileinput() == 2) { # Multiple file tests
      infile <- list.files("files/Multiple/", pattern = ".csv", full.names = T)
      lfiles <- list()
      for(i in 1:length(infile)){
        lfiles[[i]] <- read.table(infile[i], sep = "\t", header = F)
        lfiles[[i]] <- cbind(lfiles[[i]][,1:3], " Temp1" = rep(NA, nrow(lfiles[[i]])), "Temp2" = rep(NA, nrow(lfiles[[i]]))) # add one more column to allow row binding later on 
      }
      names(lfiles) <- c("test data 1", "test data 2", "test data 3", "test data 4")
      return(lfiles)   
      testfileinput(0)
    } else if (testfileinput() == 3) { # Test file mode MS2
      infile <- list.files("files/MS2/", pattern = ".csv", full.names = T)
      print(infile)
      lfiles <- list()
      for(i in 1){
        lfiles[[i]] <- read.table(infile[i], sep = "\t", header = F)
      }
      names(lfiles) <- c("test data")
      return(lfiles)
    }
  })
  
  filedataMS2 <- reactive({
    if (is.null(InputFilesMS2()) & testfileinput() != 3) {
      # User has not uploaded a file yet
      return(NULL)
    } else {
      if (testfileinput() == 3) {
        infileMS2 <- list.files("files/MS2/", pattern = "MSMS", full.names = T)[[1]]
        infilePSM <- list.files("files/MS2/", pattern = "PSM", full.names = T)[[1]]
        PSM <- read.table(infilePSM, sep = "\t", header = T)
        MS2 <- read.table(infileMS2, sep = "\t", header = T)
      } else {
        PSM <- read.table(InputFilesMS2()$PSMfile$datapath, sep = "\t", header = T)
        MS2 <- read.table(InputFilesMS2()$MS2file$datapath, sep = "\t", header = T)
        validate(
          need(sum(grepl("Master.Protein.Descriptions", names(PSM))) == 1 & sum(grepl("RT.in.min", names(MS2))) == 1, "Error in file format for plotting MS2 data.\nYou have to upload the following files:\n- A MSMSSpectrumInfo.txt file from BioPharma Finder (in the \"MS2 File\" field.\n- The corresponding PSMs.txt file (in the \"PSM File\" field).")
        )
      }
    }
    return(list("MS2file" = MS2, "PSMfile" = PSM))
  })
  
  #filedataMS2PF <- reactiveVal(NULL)
  #observeEvent(c(input$MS2filePF, input$fileMS2), {
  filedataMS2PF <- reactive({
    if (is.null(input$MS2filePF) | is.null(input$fileMS2)) {
      return(NULL)
    } else if (input$MSModeCheck == "MS2" & input$PDPFModeCheck == "PF")  {
      validate(
        need(!is.null(InputFileMS()), "You need to upload the MS2 with the associated MS file to plot MS2 results from MSPathFinder."
        )
      )
      MS2PF <- read.table(InputFilesMS2PF()$datapath, sep = "\t", header = F, skip = 1)
      vec <- lapply(unique(MS2PF[,15]), function(x) {
        MS2PF[MS2PF[,15]==x,8]
      })
      vec <- lapply(vec, function(x) {
        length(unique(x))
      })
      vec <- unlist(vec)
      val <- max(vec)
      validate(
        need(grepl("IcT", InputFilesMS2PF()$name, fixed = T), "Error in file format for plotting MS2 data.\nYou have to upload the \"IcTarget\" or \"IcTda\" output file from MSPathFinder associated with the deconvoluted MS masses uploaded as \"input file for MS\".")
      )
      validate(
        need(grepl(".ms1ft", InputFileMS()$name, fixed = T), "Error in file format for plotting MS data.\nYou have to upload the \"IcTarget\" or \"IcTda\" output file from MSPathFinder associated with the deconvoluted MS masses uploaded as \"input file for MS\".")
      )
      validate(
        need(val==1, "Several IDs have been attributed to the same MS feature.")
      )
      #names(MS2PF)[8] <- "Master.Protein.Descriptions"
      names(MS2PF)[8] <- "Protein.Descriptions"
      names(MS2PF)[15] <- "FeatureID"
      MS <- read.table(InputFileMS()$datapath, sep = "\t", header = T)
      MS2PF <- merge(MS, MS2PF, by = "FeatureID", all = T)
      MS2PF <- MS2PF[!is.na(MS2PF[,3]),]
      MS2PF <- cbind("RT" = (MS2PF$MinElutionTime + ((MS2PF$MaxElutionTime - MS2PF$MinElutionTime)/2)), MS2PF[,c("MonoMass", "ApexIntensity", "MinElutionTime", "MaxElutionTime", "Protein.Descriptions")]) 
      names(MS2PF)[2] <- "Mass"
      names(MS2PF)[3] <- "intensity"
      names(MS2PF)[4] <- "PeakStart"
      names(MS2PF)[5] <- "PeakStop"
      return(MS2PF) 
    }
  })
  
  # Create table for plotting:
  filedata <- function() {
    validate(
      need(input$IntensityThresh <= 100, "Threshold value too high")
    )
    if (is.null(filedata0())) {
      return(NULL)
    } else {
      lfiles <- filedata0()
      lfiles <- ThresholdCleaning(lfiles, input$IntensityThresh)
      if (filetype$BioPharma == 0 & filetype$ProMex == 0) { # Only RoWinPro files
        l <- list()
        for (i in seq_along(lfiles)) {
          l[[i]] <- RenameBioPharma(lfiles[[i]])
        }
        names(l) <- names(lfiles)
        lfiles <- l
        return(RBindList(lfiles))
      } else if (filetype$RoWinPro == 0) { # No RoWinPro, so BioPharma or ProMex
        lfiles <- lapply(lfiles, function(x) {
          RenameBioPharma(x)
        })
        return(RBindList(lfiles))
      } else { # More than one type of files
        lfiles <- lapply(lfiles, function(x) {
          RenameBioPharma(x) # The files from ProMex will have an empty column named "RT".
        })
        return(lfiles)
      }
    }
  }
  
  
  
  #####################
  # For zoomable plot:
  #####################
  ranges <- reactiveValues(x = NULL, y = NULL)
  observeEvent(input$DeZoom, {
    ranges$x <- oldranges$x
    ranges$y <- oldranges$y
  })
  observeEvent(input$TotalDeZoom, {
    if (is.null(filedataMS2())) { # Only MS trace, no MS2
      if (class(filedata()) != "list") { # one table
        if (filetype$ProMex > 0 | filetype$BioPharma > 0) {
          ranges$x <- c(0, range(filedata()[,5])[2])
          ranges$y <- range(filedata()[,2])
        } else {
          ranges$x <- c(0, range(filedata()[,1])[2])
          ranges$y <- range(filedata()[,2])
        }
      } else { # two tables because two types of files
        x1 <- sapply(filedata(), function(z) {
          z$RT
        })
        x2 <- sapply(filedata(), function(z) {
          z$PeakStart
        })
        x3 <- sapply(filedata(), function(z) {
          z$PeakStop
        })
        x <- c(unlist(x1), unlist(x2), unlist(x3))
        x <- x[!is.na(x)]
        y <- sapply(filedata(), function(z) {
          z$Mass
        })
        ranges$x <- c(0, range(x)[2])
        ranges$y <- range(y)
      }
    } else if (is.null(filedata())) { # Only MS2 data, no MS trace
      ranges$x <- c(0, range(filedataMS2()$MS2file$RT.in.min)[2])
      ranges$y <- range(filedataMS2()$MS2file$Precursor.MHplus.in.Da)
    } else { # Overlay of MS2 and Ms data
      if (filetype$ProMex > 0 | filetype$BioPharma > 0) {
        ranges$x <- c(0, range(c(filedataMS2()$MS2file$RT.in.min, filedata()[,5]))[2])
      } else {
        ranges$x <- c(0, range(c(filedataMS2()$MS2file$RT.in.min, filedata()[,1]))[2])
      }
      ranges$y <- range(c(filedataMS2()$MS2file$Precursor.MHplus.in.Da, filedata()[,2]))
    }
  })
  
  oldranges <- reactiveValues(x = NULL, y = NULL)
  observeEvent(event_data("plotly_selected"), {
    oldranges$x <- ranges$x
    oldranges$y <- ranges$y
    if (input$DataPoints) {
      newdata <- event_data("plotly_selected")
      if (!is.null(newdata) & class(newdata)=="data.frame") {
        if (class(filedata()) != "list" & (filetype$ProMex > 0 | filetype$BioPharma > 0)) {
          ranges$x <- c(min(filedata()[filedata()[,5] >= min(newdata$x),4]), max(filedata()[filedata()[,4] <= max(newdata$x),5]))
          ranges$y <- range(newdata$y)  
        } else if (class(filedata()) == "list") { # multiple file types
          tab <- filedata()
          if (sum(ftype()=="ProMex" | ftype()=="BioPharma") > 1) {
            tab2 <- RBindList(tab[ftype()=="ProMex" | ftype()=="BioPharma"])
          } else {
            tab2 <- tab[ftype()=="ProMex" | ftype()=="BioPharma"][[1]]
          }
          if (sum(!(ftype()=="ProMex" | ftype()=="BioPharma")) > 1) {
            tab <- RBindList(tab[!(ftype()=="ProMex" | ftype()=="BioPharma")])
          } else {
            tab <- tab[!(ftype()=="ProMex" | ftype()=="BioPharma")][[1]]
          }
          minx <- min(tab2[tab2[,5] >= min(newdata$x),4])
          minx <- min(minx, min(tab$RT))
          maxx <- max(tab2[tab2[,4] <= max(newdata$x),5])
          maxx <- max(maxx, max(tab$RT))
          ranges$x <- c(minx, maxx)
          ranges$y <- range(newdata$y)  
        } else {
          ranges$x <- range(newdata$x)
          ranges$y <- range(newdata$y)      
        }
      } else {
        newdata <- NULL
      }
    }
  })
  #####################
  
  #####################
  # Plot:
  #####################
  
  defineranges <- function() {
    if (!is.null(ranges$x) & !is.null(ranges$y)) {
      rangesx <- ranges$x
      rangesy <- ranges$y
    } else if (is.null(filedataMS2())) { # Only MS trace, no MS2
      if (class(filedata()) != "list") { # one table
        rangesx <- range(filedata()[,1])
        rangesy <- range(filedata()[,2])
      } else { # two tables because two types of files
        x1 <- sapply(filedata(), function(z) {
          z$RT
        })
        x2 <- sapply(filedata(), function(z) {
          z$PeakStart
        })
        x3 <- sapply(filedata(), function(z) {
          z$PeakStop
        })
        x <- c(unlist(x1), unlist(x2), unlist(x3))
        x <- x[!is.na(x)]
        y <- sapply(filedata(), function(z) {
          z$Mass
        })
        rangesx <- range(x)
        rangesy <- range(y)
      }
    } else if (is.null(filedata())) { # Only MS2 data, no MS trace
      rangesx <- range(filedataMS2()$MS2file$RT.in.min)
      rangesy <- range(filedataMS2()$MS2file$Precursor.MHplus.in.Da)
    } else { # Overlay of MS2 and MS data
      if (filetype$ProMex > 0 | filetype$BioPharma > 0) {
        rangesx <- range(c(filedataMS2()$MS2file$RT.in.min, filedata()[,5]))
      } else {
        rangesx <- range(c(filedataMS2()$MS2file$RT.in.min, filedata()[,1]))
      }
      rangesy <- range(c(filedataMS2()$MS2file$Precursor.MHplus.in.Da, filedata()[,2]))
    }
    return(list(rangesx, rangesy))
  }
  
  plotInput1 <- function(){
    validate(
      need(input$pch <= 10, "Please define a smaller size of points (max. 10).")
    )
    #if (!is.null(linput()) | input$MSModeCheck == "MS2") {
    if (is.null(filedata()) & is.null(filedataMS2()) & is.null(filedataMS2PF())) {
      return(NULL)
    } else {
      if (!is.null(linput())) {
        rangesx <- defineranges()[[1]]
        rangesy <- defineranges()[[2]]
        if (filetype$BioPharma == 0 & filetype$ProMex == 0) { # Only RoWinPro
          gtab <- filedata()
          if (linput() >= 2) { # if comparing several plots
            g <- ggplot() + 
              geom_point(data = gtab, aes(x = RT, y = Mass, col = File, text = paste(RT, "min\n", Mass, "Da\nSignal:", intensity, "\n", File)), alpha = 0.7, size = input$pch) +
              geom_text(parse = TRUE) +
              coord_cartesian(xlim = rangesx, ylim = rangesy, expand = TRUE) +
              theme_bw() + 
              scale_colour_brewer(palette = colval()) + 
              ylab("Protein mass (Da)") + 
              xlab("Retention time (min)")+
              theme(legend.title=element_blank())
          } else { # only one plot, no overlay
            g <- ggplot() + 
              geom_point(data = gtab, aes(x = RT, y = Mass, col = log10(intensity), text = paste(RT, "min\n", Mass, "Da\nSignal:", intensity)), alpha = 0.7, size = input$pch) +
              coord_cartesian(xlim = rangesx, ylim = rangesy, expand = TRUE) +
              theme_bw() + 
              scale_colour_distiller(palette = colval()) + 
              ylab("Protein mass (Da)") + 
              xlab("Retention time (min)")
          }
        } else if (filetype$RoWinPro == 0) { # Only type BioPharma/Promex
          gtab <- filedata()
          gtab <- gtab[gtab$PeakStart >= (rangesx[1]-(rangesx[1]*0.01)) & gtab$PeakStop <= (rangesx[2]+(rangesx[2]*0.01)),]
          if (linput() >= 2) { # if comparing several plots
            g <- ggplot() + 
              geom_segment(data = gtab, aes(y = Mass, x = PeakStart, col = File, yend = Mass, xend = PeakStop, text = paste("Start:", PeakStart, "min\n", "Stop:", PeakStop, "min\n", Mass, "Da\nSignal:", intensity, "\n", File)), alpha = 0.7, size = input$pch) + 
              geom_point(data = gtab, aes(x = RT, y = Mass, text = paste("Start:", PeakStart, "min\n", "Stop:", PeakStop, "min\n", Mass, "Da\nSignal:", intensity, "\n", File), col = File), alpha = 0) +
              coord_cartesian(xlim = rangesx, ylim = rangesy, expand = TRUE) + 
              theme_bw() + 
              scale_colour_brewer(palette = colval()) + 
              ylab("Protein mass (Da)") + 
              xlab("Retention time (min)")
          } else {
            g <- ggplot() + 
              geom_segment(data = gtab, aes(x = PeakStart, y = Mass, xend = PeakStop, yend = Mass, col = log10(intensity), text = paste("Start:", PeakStart, "min\n", "Stop:", PeakStop, "min\n", Mass, "Da\nSignal:", intensity)), alpha = 0.7, size = input$pch) +
              geom_point(data = gtab, aes(x = RT, y = Mass, col = log10(intensity), text = paste("Start:", PeakStart, "min\n", "Stop:", PeakStop, "min\n", Mass, "Da\nSignal:", intensity)), alpha = 0) +
              coord_cartesian(xlim = rangesx, ylim = rangesy, expand = TRUE) + 
              theme_bw() + 
              scale_colour_distiller(palette = colval()) + 
              ylab("Protein mass (Da)") + 
              xlab("Retention time (min)")
          }
        } else { # several types of input format
          gtabRWP <- RBindList(filedata()[ftype()=="RoWinPro" | ftype()=="Bruker"])
          gtabBP <- RBindList(filedata()[ftype()=="BioPharma" | ftype()=="ProMex"])
          
          gtabBP <- gtabBP[gtabBP$PeakStart >= (rangesx[1]-(rangesx[1]*0.01)) & gtabBP$PeakStop <= (rangesx[2]+(rangesx[2]*0.01)),]
          
          # Define the ranges for margins in the plot:
          rangesyB <- c(min(gtabBP$PeakStart, na.rm = T), max(gtabBP$PeakStop, na.rm = T))
          rangesyB <- c(min(rangesyB[1], rangesx[1]), max(rangesyB[2], rangesx[2]))
          
          g <- ggplot() + 
            geom_segment(data = gtabBP, aes(x = PeakStart, y = Mass, col = File, xend = PeakStop, yend = Mass, text = paste("Start:", PeakStart, "min\n", "Stop:", PeakStop, "min\n", Mass, "Da\nSignal:", intensity, "\n", File)), size = input$pch, alpha = 0.7) + 
            geom_point(data = gtabBP, aes(x = RT, y = Mass, text = paste("Start:", PeakStart, "min\n", "Stop:", PeakStop, "min\n", Mass, "Da\nSignal:", intensity, "\n", File), col = File), alpha = 0) +
            coord_cartesian(xlim = rangesyB, ylim = rangesy, expand = TRUE) + 
            theme_bw() + 
            scale_colour_brewer(palette = colval()) + 
            geom_point(data = gtabRWP, aes(y = Mass, x = RT, col = File, text = paste(RT, "min\n", Mass, "Da\nSignal:", intensity, "\n", File))) + 
            ylab("Protein mass (Da)") + 
            xlab("Retention time (min)")
        }
      }
      if (is.null(filedataMS2()) & is.null(filedataMS2PF())) {
        return(g)
      } else { # MS2 files loaded
        # When in MS2 mode: overlay of the MS2 values:
        if (input$MSModeCheck == "MS2" & (!is.null(filedataMS2()) | !is.null(filedataMS2PF()))) {
          if (input$PDPFModeCheck == "PD" | testfileinput() == 3) {
            PSM <- filedataMS2()$PSM
            MS2 <- filedataMS2()$MS2
            PSM$ID <- paste0(PSM$Spectrum.File, "|", PSM$First.Scan)
            MS2$ID <- paste0(MS2$Spectrum.File, "|", MS2$First.Scan)
            # Retrieve protein IDs in the MS2 table:
            MS2$Master.Protein.Descriptions <- PSM$Master.Protein.Descriptions[match(MS2$ID, PSM$ID)]
            # Plot:
            gtabMS2 <- MS2[,c("RT.in.min", "Precursor.MHplus.in.Da", "Precursor.Intensity", "Master.Protein.Descriptions")]
            gtabMS2$Identification <- ifelse(!is.na(gtabMS2$Master.Protein.Descriptions), "IDed", "Not IDed")
            
            # Action button:
            if (input$HideMSMS == TRUE) {
              gtabMS2 <- gtabMS2[gtabMS2$Identification == "IDed",]
            }
            
            names(gtabMS2)[3] <- "intensity"
            gtabMS2 <- gtabMS2[order(gtabMS2$Identification, decreasing = T),]
            
            if (!is.null(input$SelectProt) & input$PDPFModeCheck == "PD") {
              vec <- unique(gtabMS2$Master.Protein.Descriptions[gtabMS2$Master.Protein.Descriptions %in% input$SelectProt])
              vec <- vec[!is.na(vec)]
              getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
            }
            names(gtabMS2)[names(gtabMS2)=="Master.Protein.Descriptions"] <- "Protein.Descriptions"
            
            if (is.null(filedata0()) | input$MSTrace == FALSE) { # No MS trace
              g <- ggplot() + 
                geom_point(data = gtabMS2, aes(x = RT.in.min, y = Precursor.MHplus.in.Da, shape = Identification, text = paste(RT.in.min, "min\n", Precursor.MHplus.in.Da, "Da\nSignal:", intensity, "\n", Protein.Descriptions)), alpha = 0.8, size = input$pch, col = "grey30", show.legend = FALSE) + 
                coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)  + 
                theme_bw() + 
                scale_shape_manual(values = c(16, 1)) + 
                ylab("Protein mass (Da)") + 
                xlab("Retention time (min)")
              if (!is.null(input$SelectProt)) {
                g <- g + 
                  geom_point(data = gtabMS2[gtabMS2$Protein.Descriptions %in% input$SelectProt[!is.na(input$SelectProt)],], aes(x = RT.in.min, y = Precursor.MHplus.in.Da, fill = Protein.Descriptions, text = paste(RT.in.min, "min\n", Precursor.MHplus.in.Da, "Da\nSignal:", intensity, "\n", Protein.Descriptions)), shape = 21, size = input$pch, alpha = 0.8, stroke = 0, col = "white") +
                  scale_fill_manual(values = getPalette(length(vec))) 
              }
            } else if (input$MSTrace == TRUE) { # Overlay on MS trace
              g <- g +
                geom_point(data = gtabMS2, aes(x = RT.in.min, y = Precursor.MHplus.in.Da, shape = Identification, text = paste(RT.in.min, "min\n", Precursor.MHplus.in.Da, "Da\nSignal:", intensity, "\n", Protein.Descriptions)), alpha = 0.8, size = input$pch, col = "grey30", show.legend = FALSE) + 
                theme_bw() + 
                scale_shape_manual(values = c(16, 1)) + 
                ylab("Protein mass (Da)") + 
                xlab("Retention time (min)")
              if (!is.null(input$SelectProt)) {
                g <- g + 
                  geom_point(data = gtabMS2[gtabMS2$Protein.Descriptions %in% input$SelectProt[!is.na(input$SelectProt)],], aes(x = RT.in.min, y = Precursor.MHplus.in.Da, fill = Protein.Descriptions, text = paste(RT.in.min, "min\n", Precursor.MHplus.in.Da, "Da\nSignal:", intensity, "\n", Protein.Descriptions)), shape = 21, size = input$pch, alpha = 0.8, stroke = 0, col = "white") +
                  scale_fill_manual(values = getPalette(length(vec)))
              }
            }
          } else { # MSPathFinder
            gtabMS2 <- filedataMS2PF()
            gtabMS2$Identification <- ifelse(!is.na(as.character(gtabMS2$Protein.Descriptions)), "IDed", "Not IDed")
            # Action button:
            if (input$HideMSMS == TRUE) {
              gtabMS2 <- gtabMS2[gtabMS2$Identification == "IDed",]
            }
            gtabMS2 <- gtabMS2[order(gtabMS2$Identification, decreasing = T),]    
            if (!is.null(input$SelectProt)) {
              vec <- unique(gtabMS2$Protein.Descriptions[gtabMS2$Protein.Descriptions %in% input$SelectProt])
              vec <- vec[!is.na(vec)]
              getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
            }
            
            if (input$MSTrace == TRUE) {
              g <- g +
                geom_point(data = gtabMS2, aes(x = RT, y = Mass, shape = Identification, text = paste("Start:", PeakStart, "min\n", "Stop:", PeakStop, "min\n", Mass, "Da\nSignal:", intensity, "\n", Protein.Descriptions)), alpha = 0.8, size = input$pch, col = "grey30", show.legend = FALSE) + 
                theme_bw() + 
                scale_shape_manual(values = c(16, 1)) + 
                ylab("Protein mass (Da)") + 
                xlab("Retention time (min)")
              if (!is.null(input$SelectProt)) {
                g <- g + 
                  geom_point(data = gtabMS2[gtabMS2$Protein.Descriptions %in% input$SelectProt[!is.na(input$SelectProt)],], aes(x = RT, y = Mass, fill = Protein.Descriptions, text = paste("Start:", PeakStart, "min\n", "Stop:", PeakStop, "min\n", Mass, "Da\nSignal:", intensity, "\n", Protein.Descriptions)), shape = 21, size = input$pch, alpha = 0.8, stroke = 0, col = "white") +
                  scale_fill_manual(values = getPalette(length(vec)))
              }
            } else {
              g <- ggplot() +
                geom_point(data = gtabMS2, aes(x = RT, y = Mass, shape = Identification, text = paste("Start:", PeakStart, "min\n", "Stop:", PeakStop, "min\n", Mass, "Da\nSignal:", intensity, "\n", Protein.Descriptions)), alpha = 0.8, size = input$pch, col = "grey30", show.legend = FALSE) + 
                theme_bw() + 
                scale_shape_manual(values = c(16, 1)) + 
                ylab("Protein mass (Da)") + 
                xlab("Retention time (min)")
              if (!is.null(input$SelectProt)) {
                g <- g + 
                  geom_point(data = gtabMS2[gtabMS2$Protein.Descriptions %in% input$SelectProt[!is.na(input$SelectProt)],], aes(x = RT, y = Mass, fill = Protein.Descriptions, text = paste("Start:", PeakStart, "min\n", "Stop:", PeakStop, "min\n", Mass, "Da\nSignal:", intensity, "\n", Protein.Descriptions)), shape = 21, size = input$pch, alpha = 0.8, stroke = 0, col = "white") +
                  scale_fill_manual(values = getPalette(length(vec)))
              }
            }
          }
          return(g)
        }
      }
      
    }
    #}
  }
  
  
  # Plotly output if DataPoints == T
  output$plot1 <- renderPlotly({
    validate(
      need(!is.null(plotInput1()), '')
    )
    if (input$DataPoints == TRUE) {
      g <- plotInput1() + 
        theme(legend.title = element_blank())
      p <- ggplotly(g, tooltip = "text", height = 800) %>%
        layout(dragmode = "select") %>%
        config(displayModeBar = F) %>%
        layout(xaxis=list(fixedrange=TRUE)) %>%
        layout(yaxis=list(fixedrange=TRUE)) %>%
        layout(margin = list(l = 110, b = 40, r = 10, t = 10, pad = -2))  %>% 
        layout(title = "")
      # Remove IDed and Not IDed from the legend:
      if (input$MSModeCheck == "MS2") {
        if (input$MSTrace & !is.null(filedata()) & filetype$RoWinPro > 0) {
          if (input$HideMSMS) {
            p <- style(p, showlegend = FALSE, traces = 2)
          } else {
            p <- style(p, showlegend = FALSE, traces = 2:3)
          }
        } else {
          print(plotly_build(g)$data)
          if (input$HideMSMS) {
            p <- style(p, showlegend = FALSE, traces = 1)
          } else {
            p <- style(p, showlegend = FALSE, traces = 1:2)
          }
        }
      }
      p
    } else {
      plotly_empty()
    }
  })
  
  # Plotly output if DataPoints == F
  output$plot2 <- renderPlot({
    validate(
      need(!is.null(plotInput1()), '')
    )
    if (input$DataPoints == F) {
      plotInput1()  
    }
  }, height = 800)
  
  # For rendering the plot in the UI, in function of DataPoints:
  output$plotUI <- renderUI({
    if (input$DataPoints) {
      plotlyOutput("plot1")
    } else {
      plotOutput("plot2",
                 click = "plot_click",
                 dblclick = "plot_dblclick",
                 brush = brushOpts(
                   id = "plot_brush",
                   resetOnNew = TRUE))
    }
  })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$plot_dblclick, {
    oldranges$x <- ranges$x
    oldranges$y <- ranges$y
    if (input$DataPoints == FALSE) {
      brush <- input$plot_brush
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)       
      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    }
  })
  
  # For export:
  ############
  # pdf output:
  output$Download <- downloadHandler(
    filename = function(){
      if (is.null(InputFilesMS2())) {
        paste0("VisioProt-MS_", substring(InputFileMS()$name, first = 1, last = (nchar(InputFileMS()$name)-4)), "_", Sys.Date(), ".pdf")
      } else {
        paste0("VisioProt-MS_", substring(InputFilesMS2()[[1]]$name, first = 1, last = (nchar(InputFilesMS2()[[1]]$name)-4)), "_", Sys.Date(), ".pdf")
      }
    },
    content = function(file) {
      device <- function(..., width, height) {
        grDevices::pdf(..., width = 10, height = 8)
      }
      ggsave(file, plot = plotInput1(), device = device)
    })
  # png output:
  output$Download1 <- downloadHandler(
    filename = function(){
      if (is.null(InputFilesMS2())) {
        paste0("VisioProt-MS_", substring(InputFileMS()$name, first = 1, last = (nchar(InputFileMS()$name)-4)), "_", Sys.Date(), ".png")
      } else {
        paste0("VisioProt-MS_", substring(InputFilesMS2()[[1]]$name, first = 1, last = (nchar(InputFilesMS2()[[1]]$name)-4)), "_", Sys.Date(), ".png")
      }
    },
    content = function(file) {
      device <- function(..., width, height) {
        grDevices::png(..., width = 1000, height = 800, res = 120)
      }
      ggsave(file, plot = plotInput1(), device = device)
    })
  # svg output:
  output$Download2 <- downloadHandler(
    filename = function(){
      if (is.null(InputFilesMS2())) {
        paste0("VisioProt-MS_", substring(InputFileMS()$name, first = 1, last = (nchar(InputFileMS()$name)-4)), "_", Sys.Date(), ".svg")
      } else {
        paste0("VisioProt-MS_", substring(InputFilesMS2()[[1]]$name, first = 1, last = (nchar(InputFilesMS2()[[1]]$name)-4)), "_", Sys.Date(), ".svg")
      }
    },
    content = function(file) {
      device <- function(..., width, height) {
        grDevices::svg(..., width = 10, height = 8)
      }
      ggsave(file, plot = plotInput1(), device = device)
    })
  ############
  
  ############
  
}

############################################################################

shinyApp(ui = ui, server = server)

############################################################################

# rsconnect::deployApp("T:/RRelatedWork/VisioProt-MS")

