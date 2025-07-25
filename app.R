############################################################################
# VisioProt-MS: Interactive 2D maps from intact protein mass spectrometry
# Copyright CNRS 2017
# Contributor : Marie Locard-Paulet (21/11/2017) [marie.locard@ipbs.fr]
# This software is a computer program whose purpose is to visualize and inspect deconvoluted MS 2D maps.
# This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software. You can  use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info". 
# As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive licensors have only limited liability. In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or reproducing the software by the user in light of its specific status of free software,that may mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data to be ensured and, more generally, to use and operate it in the  same conditions as regards security. 
# The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
############################################################################


#=============================================================================
# STEP 1: LOAD REQUIRED PACKAGES AND SET GLOBAL OPTIONS
#=============================================================================

# Load required R packages for the Shiny application
library(shiny)         # Core Shiny framework for web applications
# devtools::install_github('hadley/ggplot2')
library(ggplot2)       # Grammar of graphics plotting system
library(plotly)        # Interactive web-based data visualization
library(dplyr)         # Data manipulation and transformation
library(RColorBrewer)  # Color palettes for data visualization
library(shinyBS)       # Bootstrap components for Shiny (tooltips, modals, etc.)
library(data.table)    # High-performance data manipulation

# Set global options for the Shiny application
options(shiny.maxRequestSize=90*1024^2) # Set maximum file upload size to 90MB

#=============================================================================
# STEP 2: SOURCE EXTERNAL HELPER FUNCTIONS
#=============================================================================

# Load custom R functions from external files
# These functions handle specific data processing tasks for different MS formats
source("files/Rfiles/Keep_highest_signal.R", local = TRUE)        # Filter data to keep highest intensity features
source("files/Rfiles/custom_rbindlist.R", local = TRUE)           # Custom function to combine data frames
source("files/Rfiles/Parse_TopPIC_input.R", local = TRUE)         # Parse TopPIC software output files
source("files/Rfiles/rename_input_coulumns.R", local = TRUE)      # Standardize column names across different formats

#=============================================================================
# STEP 3: DEFINE USER INTERFACE (UI)
#=============================================================================

# Create the main user interface layout using Shiny's fluidPage
ui <- fluidPage(
  #---------------------------------------------------------------------------
  # Header Section: Application title and help button
  #---------------------------------------------------------------------------
  fluidRow(
    # Application title in a 3-column span
    column(3, titlePanel("VisioProt-MS")
           ),
    # Help button in a 1-column span that opens documentation in new tab
    column(1, actionButton(inputId='ab1', label="?", 
                           onclick ="window.open('https://masstools.ipbs.fr/visioprothelp.html', '_blank')",
                           style="color: #fff; background-color: #673a49; border-color: #000000")
           )
  ),
  # Custom CSS styling for the help button
  tags$style(type='text/css', "#ab1 { width:80%; margin-top: 25px; font-family : Cursive; font-weight: 900; font-size: 160%;}"),
  
  #---------------------------------------------------------------------------
  # Mode Selection: Choose between MS-only or MS/MS analysis
  #---------------------------------------------------------------------------
  #checkboxInput("MSModeCheck", "MS data only", TRUE), # Legacy checkbox - now replaced by radio buttons
  radioButtons("MSModeCheck", "MS mode:",
               c("MS" = 'MS',           # MS data only mode
                 "MS/MS" = 'MS2'),      # MS/MS overlay mode
               selected = 'MS',         # Default to MS mode
               inline = TRUE            # Display options horizontally
  ),
  # Add tooltip explaining the mode selection
  bsTooltip("MSModeCheck", 
            "Choose between plotting MS data only or overlaying results of Top-Down searches",
            "right"),
  
  #---------------------------------------------------------------------------
  # Main Layout: Sidebar for controls and main panel for plots
  #---------------------------------------------------------------------------
  sidebarLayout( 
    #-------------------------------------------------------------------------
    # Sidebar Panel: Contains all user controls and file inputs
    #-------------------------------------------------------------------------
    sidebarPanel(width = 4,
                 #---------------------------------------------------------------------
                 # Conditional Panels: Different UI elements based on selected mode
                 #---------------------------------------------------------------------
                 
                 #===================================================================
                 # MS Mode Panel: For uploading and processing MS data only
                 #===================================================================
                 conditionalPanel(condition="input.MSModeCheck== 'MS'",
                                  #---------------------------------------------------------
                                  # File Input Section for MS Mode
                                  #---------------------------------------------------------
                                  fileInput("fileMS", "Select input file(s):",  
                                            accept = c(
                                              "text/csv",                              # CSV files
                                              "text/comma-separated-values,text/plain", # Text files
                                              ".txt",                                  # Text files
                                              ".ms1ft",                               # ProMex format
                                              ".csv",                                 # CSV files
                                              ".msalign"),                            # TopPIC format
                                            multiple = T,                             # Allow multiple file selection
                                            width = "100%"
                                  ),
                                  #---------------------------------------------------------
                                  # Test Mode Controls for MS Mode
                                  #---------------------------------------------------------
                                  checkboxInput("TestModeCheck", "Using test mode", FALSE),
                                  bsTooltip("TestModeCheck", 
                                            "Check to test the application without uploading any file. Then click on a button to upload a single test file or several for overlay",
                                            "right"),
                                  # Show test file buttons only when test mode is enabled
                                  conditionalPanel(condition="input.TestModeCheck==true",
                                                   actionButton("TestFile1", "Single test file"),    # Load one test file
                                                   actionButton("TestFile2", "Multiple test files")  # Load multiple test files for comparison
                                  )
                 ),
                 
                 #===================================================================
                 # MS/MS Mode Panel: For uploading MS data with MS/MS identification overlay
                 #===================================================================
                 conditionalPanel(condition="input.MSModeCheck== 'MS2'", 
                                  #---------------------------------------------------------
                                  # Test Mode Controls for MS/MS Mode
                                  #---------------------------------------------------------
                                  checkboxInput("MS2TestModeCheck", "Using test mode", FALSE),
                                  bsTooltip("MS2TestModeCheck", 
                                            "Check to test the application without uploading any file",
                                            "right"),
                                  # Display test mode message when enabled
                                  conditionalPanel(condition = "input.MS2TestModeCheck==true",
                                                   tags$span(style="color:red", "You are in test mode. Click on a button to select a single test file or multiple test files.\nUncheck to exit and upload your own data."),
                                                   br()
                                  ),
                                  
                                  #---------------------------------------------------------
                                  # File Input Section for MS/MS Mode (when not in test mode)
                                  #---------------------------------------------------------
                                  conditionalPanel(condition = "input.MS2TestModeCheck==false",
                                                   # MS data file input
                                                   fileInput("fileMS2", "Select input file for MS:",  
                                                             accept = c(
                                                               "text/csv",
                                                               "text/comma-separated-values,text/plain",
                                                               ".csv",
                                                               ".txt",
                                                               ".ms1ft",                    # ProMex format
                                                               ".msalign"),                 # TopPIC format
                                                             multiple = F,                  # Single file only for MS background
                                                             width = "100%"
                                                   ),
                                                   
                                                   #-------------------------------------------------
                                                   # Software Selection: Choose MS/MS analysis software
                                                   #-------------------------------------------------
                                                   radioButtons("PDPFModeCheck", "Origin of the MS/MS files:",
                                                                c("Proteome Discoverer" = 'PD',    # Thermo Fisher software
                                                                  "MSPathFinder" = 'PF',           # PNNL software
                                                                  "TopPIC" = 'TP'),                # TopPIC software
                                                                selected = 'PD',                   # Default to Proteome Discoverer
                                                                inline = TRUE
                                                   ),
                                                   bsTooltip("PDPFModeCheck", 
                                                             "Choose the software utilized for analysing of the top-down data",
                                                             "right"),
                                                   
                                                   #=================================================
                                                   # Proteome Discoverer File Inputs
                                                   #=================================================
                                                   conditionalPanel(condition = "input.PDPFModeCheck== 'PD'", 
                                                                    # MS/MS spectrum information file
                                                                    fileInput("MS2file", "Choose MSMSSpectrumInfo File:", 
                                                                              accept = c(
                                                                                "text/csv",
                                                                                "text/comma-separated-values,text/plain",
                                                                                ".txt")
                                                                    ),
                                                                    # Peptide/protein identification file
                                                                    fileInput("PSMfile", "Choose PSM File:", 
                                                                              accept = c(
                                                                                "text/csv",
                                                                                "text/comma-separated-values,text/plain",
                                                                                ".txt")
                                                                    )),
                                                   
                                                   #=================================================
                                                   # MSPathFinder File Inputs
                                                   #=================================================
                                                   conditionalPanel(condition = "input.PDPFModeCheck== 'PF'", 
                                                                    # MSPathFinder output file
                                                                    fileInput("MS2filePF", "Choose IcTarget or IcTda File from MSPathFinder:", 
                                                                              accept = c(
                                                                                "text/csv",
                                                                                "text/comma-separated-values,text/plain",
                                                                                ".tsv")                           # Tab-separated values
                                                                    )
                                                   ),
                                                   
                                                   #=================================================
                                                   # TopPIC File Inputs
                                                   #=================================================
                                                   conditionalPanel(condition = "input.PDPFModeCheck== 'TP'", 
                                                                    # TopPIC MS/MS file
                                                                    fileInput("MS2fileTP", "Choose MS/MS File:", 
                                                                              accept = c(
                                                                                "text",
                                                                                "text/comma-separated-values,text/plain",
                                                                                ".msalign")                       # TopPIC MS/MS format
                                                                    ),
                                                                    # TopPIC identification results file
                                                                    fileInput("IDfileTP", "Choose the OUTPUT_TABLE/_ms2_toppic (saved at tab-delimited .txt) from TopPIC:", 
                                                                              accept = c(
                                                                                "text",
                                                                                "text/comma-separated-values,text/plain",
                                                                                ".OUTPUT_TABLE")                  # TopPIC output format
                                                                    )
                                                   )
                                  ),
                                  
                                  #---------------------------------------------------------
                                  # Protein Selection and Display Options for MS/MS Mode
                                  #---------------------------------------------------------
                                  # Dropdown to select which identified proteins to highlight
                                  selectInput("SelectProt", "Select the ID to highlight:", 
                                              NULL,                                           # Options populated dynamically
                                              multiple = TRUE),                               # Allow multiple protein selection
                                  bsTooltip("SelectProt", 
                                            "Select among the identified proteins which one(s) to highlight on the plot",
                                            "right"),
                                  
                                  # Option to hide unidentified MS/MS spectra
                                  checkboxInput("HideMSMS", "Hide MS/MS withouth ID", FALSE),
                                  bsTooltip("HideMSMS", 
                                            "Removes the MS/MS spectra from the top-down analysis that were not matched to a protein.",
                                            "right"),
                                  
                                  # Option to show/hide the underlying MS trace
                                  checkboxInput("MSTrace", "Display the MS trace", TRUE),
                                  bsTooltip("MSTrace", 
                                            "Adds the MS trace to the plot.",
                                            "right")
                 ),
                 
                 #---------------------------------------------------------------------
                 # Common Controls: Available in both MS and MS/MS modes
                 #---------------------------------------------------------------------
                 
                 # Toggle between static plot and interactive plot with data labels
                 checkboxInput("DataPoints", "Show data labels (slower)", FALSE),
                 bsTooltip("DataPoints", 
                           "Switch to \"data\" mode: data appears on hovering",
                           "right"),
                 
                 #---------------------------------------------------------------------
                 # Plot Customization Parameters
                 #---------------------------------------------------------------------
                 fluidRow(
                   #-------------------------------------------------------------------
                   # Color Scale Selection
                   #-------------------------------------------------------------------
                   column(5,
                          # Color palette dropdown (options change based on number of files)
                          selectInput("colourscale", "Colour scale:",
                                      c("Spectral" = "Spectral",                # Spectral color scheme
                                        "Red/yellow/blue" = "RdYlBu",           # Diverging color scheme
                                        "Red/yellow/green" = "RdYlGn",          # Diverging color scheme
                                        "yellow to red" = "YlOrRd"              # Sequential color scheme
                                      )),
                          bsTooltip("colourscale", 
                                    "Select the colour scale for the MS data.",
                                    "right")),
                   
                   #-------------------------------------------------------------------
                   # Point Size Control
                   #-------------------------------------------------------------------
                   column(3,
                          numericInput("pch", label = "Point size:", value = 1, min = 0.1, step = 0.1, max = 10),
                          bsTooltip("pch", 
                                    "Define the size of the point (from 0.1 to 10).",
                                    "right")),
                   
                   #-------------------------------------------------------------------
                   # Intensity Threshold Control
                   #-------------------------------------------------------------------
                   column(4,
                          numericInput("IntensityThresh", label = "Threshold:", value = 20, min = 0, max = 100, step = 1),
                          bsTooltip("IntensityThresh", 
                                    "Define the percentage of highest intensity features of the MS data to display.",
                                    "right"))
                 ),
                 
                 #---------------------------------------------------------------------
                 # Zoom and Navigation Controls
                 #---------------------------------------------------------------------
                 # Dynamic instructions that change based on plot type (static vs interactive)
                 htmlOutput("ZoomParam"),
                 br(),
                 
                 # Zoom control buttons
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
                 
                 #---------------------------------------------------------------------
                 # Export Controls: Download plot in different formats
                 #---------------------------------------------------------------------
                 downloadButton("Download", "Download .pdf"),      # PDF format for publications
                 downloadButton("Download1", "Download .png"),     # PNG format for presentations
                 downloadButton("Download2", "Download .svg"),     # SVG format for vector graphics
                 
                 #---------------------------------------------------------------------
                 # Citation Information
                 #---------------------------------------------------------------------
                 HTML(paste("<br/><br/>",
                   h4("If you use VisioProt-MS for your research please cite:"),
                   "<br/>Marie Locard-Paulet, Julien Parra, Renaud Albigot, Emmanuelle Mouton-Barbosa, Laurent Bardi, Odile Burlet-Schiltz, Julien Marcoux; VisioProt-MS: interactive 2D maps from intact protein mass spectrometry, Bioinformatics, bty680, https://doi.org/10.1093/bioinformatics/bty680")
                 )
    ),
    
    #-------------------------------------------------------------------------
    # Main Panel: Contains the plot output area
    #-------------------------------------------------------------------------
    mainPanel(width = 8,
              # Dynamic UI that switches between static and interactive plots
              uiOutput("plotUI")
    )
  ),
  
  #---------------------------------------------------------------------------
  # Footer Section
  #---------------------------------------------------------------------------
  tabsetPanel(
    tabPanel(
      HTML('<footer><font size="0.8">copyright 2017 - CNRS - All rights reserved - VisioProt-MS V2.2</font></footer>')
    )
  )
)
#=============================================================================
# STEP 4: DEFINE SERVER LOGIC
#=============================================================================

server <- function(input, output, clientData, session) {
  
  #---------------------------------------------------------------------------
  # SECTION 4.1: DYNAMIC UI MODIFICATIONS
  #---------------------------------------------------------------------------
  
  #---------------------------------------------------------------------------
  # 4.1.1: Update zoom instructions based on plot type
  #---------------------------------------------------------------------------
  # Display different zoom instructions for static vs interactive plots
  output$ZoomParam <- renderUI({
    if (input$DataPoints) {
      # Instructions for interactive plots (plotly)
      HTML("<h4>Zoom in: select the ranges of interest.<br/>Zoom out: Click on the unzoom button.</h4>")
    } else {
      # Instructions for static plots (ggplot2)
      HTML("<h4>Zoom in: select the ranges of interest and double click.<br/>Zoom out: double click.</h4>")
    }
  })
  
  #---------------------------------------------------------------------------
  # 4.1.2: Dynamic color scale updates based on number of files
  #---------------------------------------------------------------------------
  # Change available color palettes based on whether single or multiple files are loaded
  observe({
    if (!is.null(linput())) {
      x <- linput()  # Get number of input files
      if (x > 1) {
        # Multiple files: use qualitative color schemes for file comparison
        updateSelectInput(session, "colourscale",
                          "Colour scale:",
                          c("Paired" = "Paired",         # Qualitative palette
                            "Set1" = "Set1",             # Qualitative palette
                            "Set2" = "Set2",             # Qualitative palette
                            "Set3" = "Set3",             # Qualitative palette
                            "Set4" = "Dark2",            # Qualitative palette
                            "Set5" = "Accent"            # Qualitative palette
                          ))
      } else {
        # Single file: use sequential/diverging color schemes for intensity visualization
        updateSelectInput(session, "colourscale",
                          "Colour scale:",
                          c("Spectral" = "Spectral",                # Diverging palette
                            "Red/yellow/blue" = "RdYlBu",           # Diverging palette
                            "Red/yellow/green" = "RdYlGn",          # Diverging palette
                            "yellow to red" = "YlOrRd"              # Sequential palette
                          ))
      }
    }
  })
  
  # Store the selected color value in a reactive variable
  colval <- reactiveVal()
  observe({
    x <- input$colourscale
    colval(x)
  })
  
  #---------------------------------------------------------------------------
  # 4.1.3: Dynamic protein ID selection for MS/MS mode
  #---------------------------------------------------------------------------
  # Populate protein selection dropdown based on loaded identification data
  observe({
    # For Proteome Discoverer data
    if (!is.null(filedataMS2())) {
      if (length(filedataMS2()$PSMfile$Majority.Protein.Accessions[!is.na(filedataMS2()$PSMfile$Majority.Protein.Accessions)]) > 0) {
        updateSelectInput(session, "SelectProt",
                          "Select the ID to highlight:",
                          sort(unique(filedataMS2()$PSMfile$Majority.Protein.Accessions[!is.na(filedataMS2()$PSMfile$Majority.Protein.Accessions)]))
        )
      }
    }
    # For MSPathFinder data
    if (!is.null(filedataMS2PF())) {
      if (length(filedataMS2PF()$Protein.Descriptions[!is.na(filedataMS2PF()$Protein.Descriptions)]) > 0) {
        updateSelectInput(session, "SelectProt",
                          "Select the ID to highlight:",
                          sort(unique(filedataMS2PF()$Protein.Descriptions[!is.na(filedataMS2PF()$Protein.Descriptions)]))
        )
      }
    }
    # For TopPIC data
    if (!is.null(filedataMS2TP())) {
      if (length(filedataMS2TP()$Protein.Descriptions[!is.na(filedataMS2TP()$Protein.Descriptions)]) > 0) {
        updateSelectInput(session, "SelectProt",
                          "Select the ID to highlight:",
                          sort(unique(filedataMS2TP()$Protein.Descriptions[!is.na(filedataMS2TP()$Protein.Descriptions)]))
        )
      }
    }
  })
  
  #---------------------------------------------------------------------------
  # 4.1.4: Track number of selected proteins for plot layout
  #---------------------------------------------------------------------------
  # Count selected proteins to adjust plot dimensions for export
  nProtSelection <- reactiveVal(0)
  observeEvent(c(plotInput1, input$SelectProt), {
    nProtSelection(length(input$SelectProt))
  })
  
  #---------------------------------------------------------------------------
  # SECTION 4.2: FILE INPUT HANDLING
  #---------------------------------------------------------------------------
  
  #---------------------------------------------------------------------------
  # 4.2.1: MS file input processing
  #---------------------------------------------------------------------------
  # Reactive variable to store MS file inputs
  InputFileMS <- reactiveVal(NULL)
  
  # Handle MS file uploads in MS mode
  observeEvent(input$fileMS, {
    ranges$x <- NULL  # Reset zoom ranges when new file uploaded
    ranges$y <- NULL
    if (input$MSModeCheck == "MS" & !is.null(input$fileMS) & input$TestModeCheck == FALSE & input$MS2TestModeCheck == FALSE) {
      InputFileMS(input$fileMS)
    } else {
      InputFileMS(NULL)
    }
  })
  
  # Handle MS file uploads in MS/MS mode (background MS data)
  observeEvent(input$fileMS2, {
    if (input$MSModeCheck == "MS2" & !is.null(input$fileMS2)) {
      InputFileMS(input$fileMS2)
    } else {
      InputFileMS(NULL)
    }
  })
  
  #---------------------------------------------------------------------------
  # 4.2.2: MS/MS file input processing for different software platforms
  #---------------------------------------------------------------------------
  
  # Proteome Discoverer files (MSMSSpectrumInfo + PSM files)
  InputFilesMS2 <- reactiveVal(NULL)
  observeEvent(c(input$MS2file, input$PSMfile), {
    if (input$MSModeCheck == "MS2" & !is.null(input$MS2file) & !is.null(input$PSMfile)) {
      InputFilesMS2(list("MS2file" = input$MS2file, "PSMfile" = input$PSMfile))
    } else {
      InputFilesMS2(NULL)
    }
  })
  
  # MSPathFinder files
  InputFilesMS2PF <- reactiveVal(NULL)
  observeEvent(input$MS2filePF, {
    if (input$MSModeCheck == "MS2" & !is.null(input$MS2filePF)) {
      InputFilesMS2PF(input$MS2filePF)
    } else {
      InputFilesMS2PF(NULL)
    }
  })
  
  # TopPIC files (MS/MS + identification files)
  InputFilesMS2TP <- reactiveVal(NULL)
  observeEvent(c(input$MS2fileTP, input$IDfileTP), {
    if (input$MSModeCheck == "MS2" & !is.null(input$MS2fileTP) & !is.null(input$IDfileTP)) {
      InputFilesMS2TP(list("MS2file" = input$MS2fileTP, "IDfile" = input$IDfileTP))
    } else {
      InputFilesMS2TP(NULL)
    }
  })
  
  #---------------------------------------------------------------------------
  # SECTION 4.3: TEST MODE AND FILE TYPE DETECTION
  #---------------------------------------------------------------------------
  
  #---------------------------------------------------------------------------
  # 4.3.1: Test mode variables and file type tracking
  #---------------------------------------------------------------------------
  # Track number of input files
  linput <- reactiveVal()
  
  # Test mode states: 0=no test, 1=single file, 2=multiple files, 3=MS2 test
  testfileinput <- reactiveVal(0)
  
  # Count files by deconvolution software type
  filetype <- reactiveValues(
    RoWinPro = 0,     # RoWinPro, Bruker, TopPIC files (point-based data)
    BioPharma = 0,    # BioPharma Finder files (peak range data)
    ProMex = 0        # ProMex files (peak range data)
  )
  
  #---------------------------------------------------------------------------
  # 4.3.2: Test mode file loading
  #---------------------------------------------------------------------------
  
  # Load single test file
  observeEvent(input$TestFile1, {
    ranges$x <- NULL  # Reset zoom
    ranges$y <- NULL
    linput(1)
    testfileinput(1)
    filetype$RoWinPro <- 1
    filetype$BioPharma <- 0
    filetype$ProMex <- 0
    colval("Spectral")  # Use spectral color scheme for single file
    InputFileMS(NULL)
    InputFilesMS2(NULL)
    InputFilesMS2PF(NULL)
  })
  
  # Load multiple test files for comparison
  observeEvent(input$TestFile2, {
    ranges$x <- NULL  # Reset zoom
    ranges$y <- NULL
    linput(4)
    testfileinput(2)
    filetype$RoWinPro <- 4
    filetype$BioPharma <- 0
    filetype$ProMex <- 0
    colval("Paired")  # Use qualitative color scheme for multiple files
    InputFileMS(NULL)
    InputFilesMS2(NULL)
    InputFilesMS2PF(NULL)
  })
  
  # Enable MS/MS test mode
  observeEvent(input$MS2TestModeCheck, {
    ranges$x <- NULL  # Reset zoom
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
  
  #---------------------------------------------------------------------------
  # 4.3.3: Automatic file type detection and counting
  #---------------------------------------------------------------------------
  # Detect file formats when MS files are uploaded
  observeEvent(c(input$fileMS, input$fileMS2), { 
    InputFileMS <- InputFileMS()
    testfileinput(0)  # Exit test mode
    
    if (!is.null(InputFileMS)) {
      # Initialize detection arrays
      l <- list()   # BioPharma format detection
      l2 <- list()  # Bruker format detection  
      l3 <- list()  # ProMex format detection
      l4 <- list()  # TopPIC format detection
      
      # Check each uploaded file for format markers
      for(i in 1:nrow(InputFileMS)){
        # BioPharma detection: check for specific column headers
        l[[i]] <- grepl("Mass", readLines(InputFileMS[i, 'datapath'])[1]) & 
                  grepl("Apex RT", readLines(InputFileMS[i, 'datapath'])[1]) & 
                  grepl("Sum Intensity", readLines(InputFileMS[i, 'datapath'])[1]) & 
                  grepl("Start Time (min)", readLines(InputFileMS[i, 'datapath'])[1], fixed = T) & 
                  grepl("Stop Time (min)", readLines(InputFileMS[i, 'datapath'])[1], fixed = T)
        
        # Bruker detection: check for specific header pattern
        l2[[i]] <- substr(readLines(InputFileMS[i, 'datapath'])[2], 0, 13) == "Compound Name"
        
        # ProMex detection: check file extension
        l3[[i]] <- grepl(".ms1ft", InputFileMS$name[i], fixed = T)
        
        # TopPIC detection: check file extension  
        l4[[i]] <- grepl("ms1.msalign", InputFileMS$name[i], fixed = T)
      }
      
      # Convert to logical vectors
      l <- unlist(l)
      l2 <- unlist(l2)
      l3 <- unlist(l3)
      l4 <- unlist(l4)
      
      # Count files by type
      filetype$RoWinPro <- sum(l==F & l3==F)      # RoWinPro + Bruker + TopPIC files
      filetype$BioPharma <- sum(l==T & l2==F & l3==F)  # BioPharma files only
      filetype$ProMex <- sum(l==F & l2==F & l3==T)     # ProMex files only
      
      # Determine appropriate number of files for color scaling
      linput(max(as.numeric(c(filetype$RoWinPro, filetype$BioPharma, filetype$ProMex))))
      
      # Handle mixed file types
      if (linput() == 1 & length(c(filetype$RoWinPro, filetype$BioPharma, filetype$ProMex)[c(filetype$RoWinPro, filetype$BioPharma, filetype$ProMex)!=0])>1) {
        linput(sum(as.numeric((c(filetype$RoWinPro, filetype$BioPharma, filetype$ProMex)))))
      }
      
      # Set appropriate color scheme based on file count
      if (linput() > 1) {
        colval("Paired")    # Multiple files: qualitative colors
      } else {
        colval("Spectral")  # Single file: continuous colors
      }
    } else {
      # Reset file type counts when no files are loaded
      filetype$RoWinPro <- 0
      filetype$BioPharma <- 0
      filetype$ProMex <- 0
    }
  })
  #---------------------------------------------------------------------------
  # SECTION 4.4: STATE MANAGEMENT AND RESET FUNCTIONS
  #---------------------------------------------------------------------------
  
  #---------------------------------------------------------------------------
  # 4.4.1: Reset plot when switching modes or test settings
  #---------------------------------------------------------------------------
  
  # Reset when exiting test mode in MS mode
  observeEvent(input$TestModeCheck, {
    InputFileMS(NULL)
    testfileinput(0)
    ranges$x <- NULL  # Clear zoom settings
    ranges$y <- NULL
    # Reset file type counters
    filetype$RoWinPro <- 0
    filetype$BioPharma <- 0
    filetype$ProMex <- 0
  })
  
  # Reset when exiting test mode in MS/MS mode
  observeEvent(input$MS2TestModeCheck, {
    if (input$MS2TestModeCheck==FALSE) {
      InputFileMS(NULL)
      InputFilesMS2(NULL)
      InputFilesMS2PF(NULL)
      testfileinput(0)
      # Clear protein selection dropdown
      updateSelectInput(session, "SelectProt",
                        "Select the ID to highlight:",
                        c(""))
      # Reset to default software selection
      updateRadioButtons(session, "PDPFModeCheck",
                        selected = 'PD')
      ranges$x <- NULL  # Clear zoom settings
      ranges$y <- NULL
    }
  })
  
  # Reset when switching between MS and MS/MS modes
  observeEvent(input$MSModeCheck, {
    InputFileMS(NULL)
    testfileinput(0)
    ranges$x <- NULL  # Clear zoom settings
    ranges$y <- NULL
  })
  
  # Reset when switching between different MS/MS software types
  observeEvent(input$PDPFModeCheck, {
    InputFilesMS2(NULL)
    InputFilesMS2PF(NULL)
    testfileinput(0)
    ranges$x <- NULL  # Clear zoom settings
    ranges$y <- NULL
    # Clear protein selection dropdown
    updateSelectInput(session, "SelectProt",
                      "Select the ID to highlight:",
                      c(""))
  })
  
  # Reset zoom when uploading new files
  observeEvent(c(input$fileMS, input$fileMS2), {
    ranges$x <- NULL
    ranges$y <- NULL
  })
  
  #---------------------------------------------------------------------------
  # SECTION 4.5: FILE TYPE VALIDATION AND FORMAT DETECTION
  #---------------------------------------------------------------------------
  
  #---------------------------------------------------------------------------
  # 4.5.1: Comprehensive file format validation
  #---------------------------------------------------------------------------
  # Reactive function to validate and classify input file formats
  ftype <- reactive({
    if (is.null(InputFileMS()) & testfileinput() == 0) {
      return(NULL)
    } else {
      InputFileMS <- InputFileMS()
      l <- list()
      
      # Error message for invalid file formats
      validationText <- "Incorrect input format.\nVisioProt-MS accepts the following input files:\n- outputs from RoWinPro (Gersch et al. 2015).\n- outputs from DataAnalysis 4.2 (Bruker).\n- BioPharma Finder 3.0 (Thermo Fisher Scientific) tables that have been exported at \"Component Level Only\" before being converted in tab-separated files.\n- ProMex exports in \".ms1ft\".\n- TopPic export of the deconvoluted MS data: \"_ms1.msalign\" files."
      
      # Check each uploaded file
      for(i in 1:nrow(InputFileMS)){
        # Initial format detection
        val <- grepl("Mass", readLines(InputFileMS[i, 'datapath'])[1]) & 
               grepl("Apex RT", readLines(InputFileMS[i, 'datapath'])[1]) & 
               grepl("Sum Intensity", readLines(InputFileMS[i, 'datapath'])[1]) & 
               grepl("Start Time (min)", readLines(InputFileMS[i, 'datapath'])[1], fixed = T) & 
               grepl("Stop Time (min)", readLines(InputFileMS[i, 'datapath'])[1], fixed = T)
        
        val2 <- substr(readLines(InputFileMS[i, 'datapath'])[2], 0, 13) == "Compound Name"  # Bruker format
        val3 <- grepl(".ms1ft", InputFileMS$name[i])                                       # ProMex format
        val4 <- grepl("ms1.msalign", InputFileMS$name[i])                                  # TopPIC format
        
        # Classify file type based on detection results
        val <- ifelse(val,  "BioPharma", "RoWinPro")
        val[val2] <- "Bruker"
        val[val3] <- "ProMex"
        val[val4] <- "TopPic"
        
        #=======================================================================
        # Format-specific validation checks
        #=======================================================================
        
        # BioPharma format validation
        if (val == "BioPharma") {
          if (!grepl("Apex RT", readLines(InputFileMS[i, 'datapath'])[1])) {
            val <- "DoNotApply"
            validate(need(val!="DoNotApply", validationText))
          }
        }
        
        # RoWinPro format validation (should have exactly 4 columns)
        if (val == "RoWinPro") {
          if (length(as.numeric(gregexpr("\t", readLines(InputFileMS[i, 'datapath'])[1])[[1]]))!=3) {
            val <- "DoNotApply"
            validate(need(val!="DoNotApply", validationText))
          }
        }
        
        # Bruker format validation
        if (val == "Bruker") {
          if (!(!grepl(",", readLines(InputFileMS[i, 'datapath'])[1]) & grepl(" RT", readLines(InputFileMS[i, 'datapath'])[2]))) {
            val <- "DoNotApply"
            validate(need(val!="DoNotApply", validationText))
          }
        }
        
        # ProMex format validation
        if (val == "ProMex") {
          if (!grepl("MinElutionTime", readLines(InputFileMS[i, 'datapath'])[1])) {
            val <- "DoNotApply"
            validate(need(val!="DoNotApply", validationText))
          }
        }
        
        # TopPIC format validation
        if (val == "TopPic") {
          if (!grepl("#TopFD", readLines(InputFileMS[i, 'datapath'])[1], fixed = T)) {
            val <- "DoNotApply"
            validate(need(val!="DoNotApply", validationText))
          }
        }
        
        l[[i]] <- val
      }
      vec <- unlist(l)
      return(vec)
    }
  })
  #---------------------------------------------------------------------------
  # SECTION 4.6: DATA LOADING AND PROCESSING FUNCTIONS
  #---------------------------------------------------------------------------
  
  #---------------------------------------------------------------------------
  # 4.6.1: Main MS data loading function
  #---------------------------------------------------------------------------
  # Load and process MS data files based on format type and test mode
  filedata0 <- reactive({
    if (testfileinput() == 0) { 
      # Regular file upload mode (not test mode)
      validate(
        need((input$TestModeCheck==FALSE & input$MSModeCheck == "MS") | (input$MS2TestModeCheck == FALSE & input$MSModeCheck == "MS2"), 
             "You are in test mode. Click on a button to select a single test file or multiple test files.\nUncheck to exit and upload your own data.")
      )
      
      if (is.null(InputFileMS())) {
        return(NULL)  # No file uploaded yet
      } else {
        # Validation: prevent mixing different file types with multiple files
        validate(
          need(!(max(table(ftype())) > 1 & length(unique(ftype())) > 1), 
               "Visioprot-MS can compare several files from the same deconvolution tool (RoWinPro, BioPharma Finder, MSPathFinder or Bruker). If you want to compare files from different input types, please only input one file per tool.")
        )
        
        InputFileMS <- InputFileMS()
        lfiles <- list()
        
        # Process each uploaded file according to its detected format
        for(i in 1:nrow(InputFileMS)){
          if (ftype()[i] == "BioPharma") {
            #===================================================================
            # BioPharma format processing
            #===================================================================
            lfiles[[i]] <- read.table(InputFileMS[i, 'datapath'], sep = "\t", header = T)
            # Standardize mass column name
            names(lfiles[[i]])[names(lfiles[[i]])=="Monoisotopic.Mass"|names(lfiles[[i]])=="Average.Mass"] <- "Mass"
            # Select relevant columns and map to standard format
            lfiles[[i]] <- lfiles[[i]][,c("Apex.RT", "Mass", "Sum.Intensity", "Start.Time..min.", "Stop.Time..min.")]
            
          } else if (ftype()[i] == "RoWinPro")  {
            #===================================================================
            # RoWinPro format processing
            #===================================================================
            lfiles[[i]] <- read.table(InputFileMS[i, 'datapath'], sep = "\t", header = F)
            # Add placeholder columns for peak start/stop times (not available in RoWinPro format)
            lfiles[[i]] <- cbind(lfiles[[i]][,1:3], " Temp1" = rep(NA, nrow(lfiles[[i]])), "Temp2" = rep(NA, nrow(lfiles[[i]])))
            
          } else if (ftype()[i] == "Bruker")  {
            #===================================================================
            # Bruker format processing
            #===================================================================
            lfiles[[i]] <- read.table(InputFileMS[i, 'datapath'], sep = ",", header = F, skip = 2)
            # Select relevant columns (RT, Mass, Intensity) and add placeholder columns
            lfiles[[i]] <- cbind(lfiles[[i]][,2:4], " Temp1" = rep(NA, nrow(lfiles[[i]])), "Temp2" = rep(NA, nrow(lfiles[[i]])))
            
          } else if (ftype()[i] == "ProMex")  {
            #===================================================================
            # ProMex format processing
            #===================================================================
            lfiles[[i]] <- read.table(InputFileMS[i, 'datapath'], sep = "\t", header = T)
            # Calculate peak center time and map columns to standard format
            lfiles[[i]] <- cbind("RT" = (lfiles[[i]]$MinElutionTime + ((lfiles[[i]]$MaxElutionTime - lfiles[[i]]$MinElutionTime)/2)), 
                                lfiles[[i]][,c("MonoMass", "ApexIntensity", "MinElutionTime", "MaxElutionTime")])
            
          } else if (ftype()[i] == "TopPic")  {
            #===================================================================
            # TopPIC format processing
            #===================================================================
            lfiles[[i]] <- TopPicMS1Parsing(InputFileMS[i,'datapath'])  # Use custom parsing function
          }
        }
        names(lfiles) <- InputFileMS()$name  # Set file names as list names
        return(lfiles)
      }
      
    } else if (testfileinput() == 1) {
      #=========================================================================
      # Test mode: Single file
      #=========================================================================
      infile <- list.files("files/Unique/", pattern = ".csv", full.names = T)
      lfiles <- list()
      for(i in 1){
        lfiles[[i]] <- read.table(infile[i], sep = "\t", header = F)
      }
      names(lfiles) <- c("test data")
      return(lfiles)
      testfileinput(0)
      
    } else if (testfileinput() == 2) {
      #=========================================================================
      # Test mode: Multiple files
      #=========================================================================
      infile <- list.files("files/Multiple/", pattern = ".csv", full.names = T)
      lfiles <- list()
      for(i in 1:length(infile)){
        lfiles[[i]] <- read.table(infile[i], sep = "\t", header = F)
        # Add placeholder columns for consistency
        lfiles[[i]] <- cbind(lfiles[[i]][,1:3], " Temp1" = rep(NA, nrow(lfiles[[i]])), "Temp2" = rep(NA, nrow(lfiles[[i]])))
      }
      names(lfiles) <- c("test data 1", "test data 2", "test data 3", "test data 4")
      return(lfiles)   
      testfileinput(0)
      
    } else if (testfileinput() == 3) {
      #=========================================================================
      # Test mode: MS/MS mode
      #=========================================================================
      infile <- list.files("files/MS2/", pattern = ".csv", full.names = T)
      lfiles <- list()
      for(i in 1){
        lfiles[[i]] <- read.table(infile[i], sep = "\t", header = F)
      }
      names(lfiles) <- c("test data")
      return(lfiles)
    }
  })
  
  #---------------------------------------------------------------------------
  # 4.6.2: MS/MS data loading for Proteome Discoverer
  #---------------------------------------------------------------------------
  # Load and process Proteome Discoverer MS/MS identification data
  filedataMS2 <- reactive({
    if (is.null(InputFilesMS2()) & testfileinput() != 3) {
      return(NULL)  # No MS/MS files uploaded
    } else {
      if (testfileinput() == 3) {
        #=====================================================================
        # Test mode: Load test MS/MS files
        #=====================================================================
        infileMS2 <- list.files("files/MS2/", pattern = "MSMS", full.names = T)[[1]]
        infilePSM <- list.files("files/MS2/", pattern = "SMs.txt", full.names = T)[[1]]
        PSM <- read.table(infilePSM, sep = "\t", header = T, comment.char = "#")
        MS2 <- read.table(infileMS2, sep = "\t", header = T)
      } else {
        #=====================================================================
        # Regular mode: Load user-uploaded files
        #=====================================================================
        PSM <- read.table(InputFilesMS2()$PSMfile$datapath, sep = "\t", header = T)
        MS2 <- read.table(InputFilesMS2()$MS2file$datapath, sep = "\t", header = T, comment.char = "#")
        
        # Validate file format for Proteome Discoverer
        validate(
          need((sum(grepl("Master.Protein.Descriptions", names(PSM))) == 1 & sum(grepl("RT.in.min", names(MS2))) == 1) | 
               (sum(grepl("Protein.Accessions", names(PSM))) == 1 & sum(grepl("RT..min.", names(MS2))) == 1), 
               "Error in file format for plotting MS2 data.\nYou have to upload the following files:\n- A MSMSSpectrumInfo.txt file from BioPharma Finder (in the \"MS/MS File\" field).\n- The corresponding PSMs.txt or PrSMs.txt file (in the \"PSM File\" field).")
        )
        
        # Standardize column names for consistency
        PSM <- RenamePDPrSM(PSM)
        MS2 <- RenamePDMSMS(MS2)
      }
    }
    return(list("MS2file" = MS2, "PSMfile" = PSM))
  })
  
  #---------------------------------------------------------------------------
  # 4.6.3: MS/MS data loading for MSPathFinder
  #---------------------------------------------------------------------------
  # Load and process MSPathFinder identification data
  filedataMS2PF <- reactive({
    if (is.null(input$MS2filePF) | is.null(input$fileMS2)) {
      return(NULL)
    } else if (input$MSModeCheck == "MS2" & input$PDPFModeCheck == "PF")  {
      # Validate that both MS and MS/MS files are uploaded
      validate(
        need(!is.null(InputFileMS()), "You need to upload the MS2 with the associated MS file to plot MS2 results from MSPathFinder.")
      )
      
      # Load MSPathFinder results
      MS2PF <- read.table(InputFilesMS2PF()$datapath, sep = "\t", header = F, skip = 1)
      
      # Validate that each MS feature has only one identification
      vec <- lapply(unique(MS2PF[,15]), function(x) {
        MS2PF[MS2PF[,15]==x,8]
      })
      vec <- lapply(vec, function(x) {
        length(unique(x))
      })
      vec <- unlist(vec)
      val <- max(vec)
      
      # File format validation
      validate(
        need(grepl("IcT", InputFilesMS2PF()$name, fixed = T), 
             "Error in file format for plotting MS2 data.\nYou have to upload the \"IcTarget\" or \"IcTda\" output file from MSPathFinder associated with the deconvoluted MS weights uploaded as \"input file for MS\".")
      )
      validate(
        need(grepl(".ms1ft", InputFileMS()$name, fixed = T), 
             "Error in file format for plotting MS data.\nYou have to upload the \"IcTarget\" or \"IcTda\" output file from MSPathFinder associated with the deconvoluted MS weights uploaded as \"input file for MS\".")
      )
      validate(
        need(val==1, "Several IDs have been attributed to the same MS feature.")
      )
      
      # Standardize column names
      names(MS2PF)[8] <- "Protein.Descriptions"
      names(MS2PF)[15] <- "FeatureID"
      
      # Merge with MS data and process
      MS <- read.table(InputFileMS()$datapath, sep = "\t", header = T)
      MS2PF <- merge(MS, MS2PF, by = "FeatureID", all = T)
      MS2PF <- MS2PF[!is.na(MS2PF[,3]),]  # Remove entries without identification
      
      # Calculate retention time and standardize column names
      MS2PF <- cbind("RT" = (MS2PF$MinElutionTime + ((MS2PF$MaxElutionTime - MS2PF$MinElutionTime)/2)), 
                     MS2PF[,c("MonoMass", "ApexIntensity", "MinElutionTime", "MaxElutionTime", "Protein.Descriptions")]) 
      names(MS2PF)[2] <- "Mass"
      names(MS2PF)[3] <- "intensity"
      names(MS2PF)[4] <- "PeakStart"
      names(MS2PF)[5] <- "PeakStop"
      return(MS2PF) 
    }
  })
  
  #---------------------------------------------------------------------------
  # 4.6.4: MS/MS data loading for TopPIC
  #---------------------------------------------------------------------------
  # Load and process TopPIC identification data
  filedataMS2TP <- reactive({
    if (is.null(InputFilesMS2TP())) {
      return(NULL)
    } else if (input$MSModeCheck == "MS2" & input$PDPFModeCheck == "TP")  {
      # Validate TopPIC file formats
      validate(
        need(grepl("_ms2.msalign", InputFilesMS2TP()$MS2file$name, fixed = T), 
             "Error in file format for plotting MS2 data.\nYou have to upload the \"_ms2.msalign\" output file from TopPic associated with the \"_ms2.OUTPUT_TABLE\".")
      )
      validate(
        need(grepl("_ms2.OUTPUT_TABLE", InputFilesMS2TP()$IDfile$name, fixed = T) | 
             grepl("_ms2_toppic", InputFilesMS2TP()$IDfile$name, fixed = T), 
             "Error in file format for plotting ID data.\nYou have to upload the \"_ms2.OUTPUT_TABLE\", or \"_ms2_toppic\" output file from TopPic associated with the deconvoluted MS2 weights uploaded as \"input file for MS2\".")
      )
      
      # Parse TopPIC files using custom functions
      IDTP <- TopPicIDParsing(InputFilesMS2TP()$IDfile$datapath)
      MS2TP <- TopPicMS2Parsing(InputFilesMS2TP()$MS2file$datapath)
      
      # Merge identification and MS/MS data
      names(IDTP)[names(IDTP) == "Spectrum ID"] <- "Scan"
      dat <- merge(MS2TP, IDTP, by = "Scan", all = T)
      names(dat)[names(dat)=="Protein accession"] <- "Protein.Descriptions"
      dat$Mass <- as.numeric(dat$Mass)
      dat$Identification <- ifelse(!is.na(as.character(dat$Protein.Descriptions)), "IDed", "Not IDed")
      return(dat)
    }
  })
  #---------------------------------------------------------------------------
  # SECTION 4.7: DATA PROCESSING AND PREPARATION FOR PLOTTING
  #---------------------------------------------------------------------------
  
  #---------------------------------------------------------------------------
  # 4.7.1: Final data processing and formatting for visualization
  #---------------------------------------------------------------------------
  # Process loaded data for plotting, including intensity filtering and column standardization
  filedata <- function() {
    # Validate intensity threshold parameter
    validate(
      need(input$IntensityThresh <= 100, "Threshold value too high")
    )
    
    if (is.null(filedata0())) {
      return(NULL)
    } else {
      lfiles <- filedata0()
      
      # Apply intensity threshold filtering to keep only high-intensity features
      lfiles <- ThresholdCleaning(lfiles, input$IntensityThresh)
      
      if (filetype$BioPharma == 0 & filetype$ProMex == 0) {
        #=====================================================================
        # Only RoWinPro-type files (point-based data)
        #=====================================================================
        l <- list()
        for (i in seq_along(lfiles)) {
          l[[i]] <- RenameBioPharma(lfiles[[i]])  # Standardize column names
        }
        names(l) <- names(lfiles)
        lfiles <- l
        return(RBindList(lfiles))  # Combine into single data frame
        
      } else if (filetype$RoWinPro == 0) {
        #=====================================================================
        # Only BioPharma/ProMex files (peak range data)
        #=====================================================================
        lfiles <- lapply(lfiles, function(x) {
          RenameBioPharma(x)  # Standardize column names
        })
        return(RBindList(lfiles))  # Combine into single data frame
        
      } else {
        #=====================================================================
        # Mixed file types - keep as separate list for different processing
        #=====================================================================
        lfiles <- lapply(lfiles, function(x) {
          RenameBioPharma(x)  # Standardize column names
        })
        return(lfiles)  # Return as list to handle different data types
      }
    }
  }
  
  #---------------------------------------------------------------------------
  # SECTION 4.8: ZOOM AND INTERACTION HANDLING
  #---------------------------------------------------------------------------
  
  #---------------------------------------------------------------------------
  # 4.8.1: Zoom range management
  #---------------------------------------------------------------------------
  # Reactive values to store current plot zoom ranges
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  # Store previous zoom state for "unzoom one step" functionality
  oldranges <- reactiveValues(x = NULL, y = NULL)
  
  # Handle "Unzoom one step" button - restore previous zoom level
  observeEvent(input$DeZoom, {
    ranges$x <- oldranges$x
    ranges$y <- oldranges$y
  })
  
  # Handle "Total unzoom" button - reset to full data range
  observeEvent(input$TotalDeZoom, {
    if (is.null(filedataMS2())) {
      #=========================================================================
      # MS-only mode: Calculate ranges from MS data
      #=========================================================================
      if (class(filedata()) != "list") {
        # Single data table
        if (filetype$ProMex > 0 | filetype$BioPharma > 0) {
          # Peak range data: use peak stop times for x-axis
          ranges$x <- c(0, range(filedata()[,5])[2])
          ranges$y <- range(filedata()[,2])
        } else {
          # Point data: use retention times for x-axis  
          ranges$x <- c(0, range(filedata()[,1])[2])
          ranges$y <- range(filedata()[,2])
        }
      } else {
        # Multiple data tables (mixed file types)
        x1 <- sapply(filedata(), function(z) { z$RT })
        x2 <- sapply(filedata(), function(z) { z$PeakStart })
        x3 <- sapply(filedata(), function(z) { z$PeakStop })
        x <- c(unlist(x1), unlist(x2), unlist(x3))
        x <- x[!is.na(x)]
        y <- sapply(filedata(), function(z) { z$Mass })
        ranges$x <- c(0, range(x)[2])
        ranges$y <- range(y)
      }
    } else if (is.null(filedata())) {
      #=========================================================================
      # MS/MS-only mode: Calculate ranges from MS/MS data
      #=========================================================================
      if (input$PDPFModeCheck == "PD") {
        ranges$x <- c(0, range(filedataMS2()$MS2file$RT.in.min)[2])
        ranges$y <- range(filedataMS2()$MS2file$Precursor.MHplus.in.Da)
      } else {
        # TopPIC data
        ranges$x <- c(0, range(filedataMS2TP()$MS2file$RT)[2])
        ranges$y <- range(filedataMS2TP()$MS2file$Mass)
      }
    } else {
      #=========================================================================
      # Combined MS + MS/MS mode: Calculate ranges from both datasets
      #=========================================================================
      if (filetype$ProMex > 0 | filetype$BioPharma > 0) {
        ranges$x <- c(0, range(c(filedataMS2()$MS2file$RT.in.min, filedata()[,5]))[2])
        ranges$y <- range(c(filedataMS2()$MS2file$Precursor.MHplus.in.Da, filedata()[,2]))
      } else if (input$PDPFModeCheck == "PD") {
        ranges$x <- c(0, range(c(filedataMS2()$MS2file$RT.in.min, filedata()[,1]))[2])
        ranges$y <- range(c(filedataMS2()$MS2file$Precursor.MHplus.in.Da, filedata()[,2]))
      } else if (input$PDPFModeCheck == "TP") {
        ranges$x <- c(0, range(c(filedataMS2TP()$MS2file$RT, filedata()[,1]))[2])
        ranges$y <- range(c(filedataMS2()$MS2file$Mass, filedata()[,2]))
      }
    }
  })
  
  #---------------------------------------------------------------------------
  # 4.8.2: Interactive plot selection handling (Plotly mode)
  #---------------------------------------------------------------------------
  # Handle zoom selection in interactive plots
  observeEvent(event_data("plotly_selected"), {
    oldranges$x <- ranges$x  # Store current ranges before changing
    oldranges$y <- ranges$y
    
    if (input$DataPoints) {
      newdata <- event_data("plotly_selected")
      if (!is.null(newdata) & class(newdata)=="data.frame") {
        if (class(filedata()) != "list" & (filetype$ProMex > 0 | filetype$BioPharma > 0)) {
          #=================================================================
          # Peak range data: Convert point selection to peak ranges
          #=================================================================
          ranges$x <- c(min(filedata()[filedata()[,5] >= min(newdata$x),4]), 
                       max(filedata()[filedata()[,4] <= max(newdata$x),5]))
          ranges$y <- range(newdata$y)  
          
        } else if (class(filedata()) == "list") {
          #=================================================================
          # Mixed file types: Handle both point and range data
          #=================================================================
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
          
          # Calculate combined ranges from both data types
          minx <- min(tab2[tab2[,5] >= min(newdata$x),4])
          minx <- min(minx, min(tab$RT))
          maxx <- max(tab2[tab2[,4] <= max(newdata$x),5])
          maxx <- max(maxx, max(tab$RT))
          ranges$x <- c(minx, maxx)
          ranges$y <- range(newdata$y)  
          
        } else {
          #=================================================================
          # Point data: Direct range assignment
          #=================================================================
          ranges$x <- range(newdata$x)
          ranges$y <- range(newdata$y)      
        }
      } else {
        newdata <- NULL
      }
    }
  })
  
  #---------------------------------------------------------------------------
  # SECTION 4.9: PLOTTING FUNCTIONS
  #---------------------------------------------------------------------------
  
  #---------------------------------------------------------------------------
  # 4.9.1: Helper function to define plot ranges
  #---------------------------------------------------------------------------
  # Calculate appropriate plot ranges based on current data and zoom state
  defineranges <- function() {
    if (!is.null(ranges$x) & !is.null(ranges$y)) {
      # Use current zoom ranges if they exist
      rangesx <- ranges$x
      rangesy <- ranges$y
      
    } else if (is.null(filedataMS2())) {
      #=========================================================================
      # MS-only mode: Calculate ranges from MS data
      #=========================================================================
      if (class(filedata()) != "list") {
        # Single data table
        rangesx <- range(filedata()[,1])  # Retention time range
        rangesy <- range(filedata()[,2])  # Mass range
      } else {
        # Multiple data tables (mixed file types)
        x1 <- sapply(filedata(), function(z) { as.numeric(z$RT) })
        x2 <- sapply(filedata(), function(z) { as.numeric(z$PeakStart) })
        x3 <- sapply(filedata(), function(z) { as.numeric(z$PeakStop) })
        x <- c(unlist(x1), unlist(x2), unlist(x3))
        x <- x[!is.na(x)]
        y <- sapply(filedata(), function(z) { as.numeric(z$Mass) })
        rangesx <- range(x)
        rangesy <- range(y)
      }
      
    } else if (is.null(filedata()) | input$MSTrace == FALSE) {
      #=========================================================================
      # MS/MS-only mode: Calculate ranges from MS/MS data
      #=========================================================================
      rangesx <- range(filedataMS2()$MS2file$RT.in.min)
      rangesy <- range(filedataMS2()$MS2file$Precursor.MHplus.in.Da)
      
    } else {
      #=========================================================================
      # Combined MS + MS/MS mode: Calculate ranges from both datasets
      #=========================================================================
      if (filetype$ProMex > 0 | filetype$BioPharma > 0) {
        rangesx <- range(c(filedataMS2()$MS2file$RT.in.min, filedata()[,5]))
        rangesy <- range(c(filedataMS2()$MS2file$Precursor.MHplus.in.Da, filedata()[,2]))
      }  else if (input$PDPFModeCheck == "PD") {
        ranges$x <- range(c(filedataMS2()$MS2file$RT.in.min, filedata()$RT))
        ranges$y <- range(c(filedataMS2()$MS2file$Precursor.MHplus.in.Da, filedata()$Mass))
      } else if (input$PDPFModeCheck == "TP") {
        ranges$x <- range(c(filedataMS2TP()$RT, filedata()[,1]))[2]
        ranges$y <- range(c(filedataMS2TP()$Mass, filedata()[,2]))
      }
    }
    return(list(rangesx, rangesy))
  }
  
  #---------------------------------------------------------------------------
  # 4.9.2: Main plotting function
  #---------------------------------------------------------------------------
  # Generate the main plot based on data type and user settings
  plotInput1 <- function(){
    # Validate user input parameters
    validate(
      need(input$pch <= 10 & input$pch >= 0.1, "Please define a size of points between 0.1 and 10.")
    )
    validate(
      need(input$IntensityThresh <= 100 & input$IntensityThresh >= 0.1, 
           "Please define a threshold between 0.1 and 100. This value corresponds to the proportion of the points in the data set that you want to upload (highest intensities).")
    )
    
    # Check if any data is available for plotting
    if (is.null(filedata()) & is.null(filedataMS2()) & is.null(filedataMS2PF()) & is.null(filedataMS2TP())) {
      return(NULL)
    } else {
      if (!is.null(linput())) {
        # Get plot ranges
        rangesx <- defineranges()[[1]]
        rangesy <- defineranges()[[2]]
        
        #=====================================================================
        # MS DATA PLOTTING SECTION
        #=====================================================================
        
        if (filetype$BioPharma == 0 & filetype$ProMex == 0) {
          #===================================================================
          # RoWinPro-type data (point-based): Plot as points
          #===================================================================
          gtab <- filedata()
          
          if (linput() >= 2) {
            #=================================================================
            # Multiple files: Color by file name
            #=================================================================
            g <- ggplot() + 
              geom_point(data = gtab, 
                        aes(x = RT, y = Mass, col = File, 
                            text = paste(RT, "min\n", Mass, "Da\nSignal:", intensity, "\n", File)), 
                        alpha = 0.7, size = input$pch) +
              geom_text(parse = TRUE) +
              coord_cartesian(xlim = rangesx, ylim = rangesy, expand = TRUE) +
              theme_bw() + 
              scale_colour_brewer(palette = colval()) + 
              ylab("Molecular Weight (Da)") + 
              xlab("Retention time (min)")
          } else {
            #=================================================================
            # Single file: Color by intensity
            #=================================================================
            g <- ggplot() + 
              geom_point(data = gtab, 
                        aes(x = RT, y = Mass, col = log10(intensity), 
                            text = paste(RT, "min\n", Mass, "Da\nSignal:", intensity)), 
                        alpha = 0.7, size = input$pch) +
              coord_cartesian(xlim = rangesx, ylim = rangesy, expand = TRUE) +
              theme_bw() + 
              scale_colour_distiller(palette = colval()) + 
              ylab("Molecular Weight (Da)") + 
              xlab("Retention time (min)")
          }
          
        } else if (filetype$RoWinPro == 0) {
          #===================================================================
          # BioPharma/ProMex data (peak ranges): Plot as line segments
          #===================================================================
          gtab <- filedata()
          # Filter data to current view range (with small margin)
          gtab <- gtab[gtab$PeakStart >= (rangesx[1]-(rangesx[1]*0.01)) & 
                      gtab$PeakStop <= (rangesx[2]+(rangesx[2]*0.01)),]
          
          if (linput() >= 2) {
            #=================================================================
            # Multiple files: Color by file name
            #=================================================================
            g <- ggplot() + 
              geom_segment(data = gtab, 
                          aes(y = Mass, x = PeakStart, col = File, yend = Mass, xend = PeakStop, 
                              text = paste("Start:", PeakStart, "min\n", "Stop:", PeakStop, "min\n", 
                                          Mass, "Da\nSignal:", intensity, "\n", File)), 
                          alpha = 0.7, size = input$pch) + 
              # Invisible points for plotly interaction
              geom_point(data = gtab, 
                        aes(x = RT, y = Mass, 
                            text = paste("Start:", PeakStart, "min\n", "Stop:", PeakStop, "min\n", 
                                        Mass, "Da\nSignal:", intensity, "\n", File), col = File), 
                        alpha = 0) +
              coord_cartesian(xlim = rangesx, ylim = rangesy, expand = TRUE) + 
              theme_bw() + 
              scale_colour_brewer(palette = colval()) + 
              ylab("Molecular Weight (Da)") + 
              xlab("Retention time (min)")
          } else {
            #=================================================================
            # Single file: Color by intensity
            #=================================================================
            g <- ggplot() + 
              geom_segment(data = gtab, 
                          aes(x = PeakStart, y = Mass, xend = PeakStop, yend = Mass, col = log10(intensity), 
                              text = paste("Start:", PeakStart, "min\n", "Stop:", PeakStop, "min\n", 
                                          Mass, "Da\nSignal:", intensity)), 
                          alpha = 0.7, size = input$pch) +
              # Invisible points for plotly interaction
              geom_point(data = gtab, 
                        aes(x = RT, y = Mass, col = log10(intensity), 
                            text = paste("Start:", PeakStart, "min\n", "Stop:", PeakStop, "min\n", 
                                        Mass, "Da\nSignal:", intensity)), 
                        alpha = 0) +
              coord_cartesian(xlim = rangesx, ylim = rangesy, expand = TRUE) + 
              theme_bw() + 
              scale_colour_distiller(palette = colval()) + 
              ylab("Molecular Weight (Da)") + 
              xlab("Retention time (min)")
          }
          
        } else {
          #===================================================================
          # Mixed file types: Combine point and range data visualization
          #===================================================================
          gtabRWP <- RBindList(filedata()[ftype()=="RoWinPro" | ftype()=="Bruker" | ftype()=="TopPic"])
          gtabBP <- RBindList(filedata()[ftype()=="BioPharma" | ftype()=="ProMex"])
          
          # Filter peak range data to current view
          gtabBP <- gtabBP[gtabBP$PeakStart >= (rangesx[1]-(rangesx[1]*0.01)) & 
                          gtabBP$PeakStop <= (rangesx[2]+(rangesx[2]*0.01)),]
          
          # Define plot ranges including both data types
          rangesyB <- c(min(gtabBP$PeakStart, na.rm = T), max(gtabBP$PeakStop, na.rm = T))
          rangesyB <- c(min(rangesyB[1], rangesx[1]), max(rangesyB[2], rangesx[2]))
          
          # Create combined plot with segments and points
          g <- ggplot() + 
            geom_segment(data = gtabBP, 
                        aes(x = PeakStart, y = Mass, col = File, xend = PeakStop, yend = Mass, 
                            text = paste("Start:", PeakStart, "min\n", "Stop:", PeakStop, "min\n", 
                                        Mass, "Da\nSignal:", intensity, "\n", File)), 
                        size = input$pch, alpha = 0.7) + 
            geom_point(data = gtabBP, 
                      aes(x = RT, y = Mass, 
                          text = paste("Start:", PeakStart, "min\n", "Stop:", PeakStop, "min\n", 
                                      Mass, "Da\nSignal:", intensity, "\n", File), col = File), 
                      alpha = 0) +
            coord_cartesian(xlim = rangesyB, ylim = rangesy, expand = TRUE) + 
            theme_bw() + 
            scale_colour_brewer(palette = colval()) + 
            geom_point(data = gtabRWP, 
                      aes(y = Mass, x = RT, col = File, 
                          text = paste(RT, "min\n", Mass, "Da\nSignal:", intensity, "\n", File)), 
                      size = input$pch) + 
            ylab("Molecular Weight (Da)") + 
            xlab("Retention time (min)")
        }
      }
      
      #=========================================================================
      # MS/MS DATA OVERLAY SECTION
      #=========================================================================
      
      if (is.null(filedataMS2()) & is.null(filedataMS2PF()) & is.null(filedataMS2TP())) {
        # No MS/MS data: return the MS plot only
        return(g)
      } else {
        # MS/MS files loaded: add MS/MS overlay to the plot
        if (input$MSModeCheck == "MS2" & (!is.null(filedataMS2()) | !is.null(filedataMS2PF()) | !is.null(filedataMS2TP()))) {
          #===================================================================
          # PROTEOME DISCOVERER MS/MS DATA PROCESSING AND PLOTTING
          #===================================================================
          if (input$PDPFModeCheck == "PD" | testfileinput() == 3) {
            PSM <- filedataMS2()$PSM  # Peptide/protein identification data
            MS2 <- filedataMS2()$MS2  # MS/MS spectrum information
            
            # Create unique identifiers for matching PSM and MS/MS data
            if (sum(grepl("First.Scan", names(PSM))) == 1 & sum(grepl("First.Scan", names(MS2))) == 1) {
              # Use spectrum file + scan number as identifier
              PSM$ID <- paste0(PSM$Spectrum.File, "|", PSM$First.Scan)
              MS2$ID <- paste0(MS2$Spectrum.File, "|", MS2$First.Scan)
            } else {
              # Use spectrum file + m/z value as identifier
              PSM$ID <- paste0(PSM$Spectrum.File, "|", PSM$m.z..Da.)
              MS2$ID <- paste0(MS2$Spectrum.File, "|", MS2$Precursor.m.z..Da.)
              MS2$Master.Protein.Descriptions <- PSM$Master.Protein.Descriptions[match(MS2$ID, PSM$ID)]
            }
            
            # Map protein identifications to MS/MS spectra
            MS2$Master.Protein.Descriptions <- PSM$Master.Protein.Descriptions[match(MS2$ID, PSM$ID)]
            
            # Prepare data for plotting
            gtabMS2 <- MS2[,c("RT.in.min", "Precursor.MHplus.in.Da", "Master.Protein.Descriptions")]
            gtabMS2$Identification <- ifelse(!is.na(gtabMS2$Master.Protein.Descriptions), "IDed", "Not IDed")
            
            # Apply filter to hide unidentified MS/MS if requested
            if (input$HideMSMS == TRUE) {
              gtabMS2 <- gtabMS2[gtabMS2$Identification == "IDed",]
            }
            
            # Order data to show identified spectra on top
            gtabMS2 <- gtabMS2[order(gtabMS2$Identification, decreasing = T),]
            
            # Prepare color palette for selected proteins
            if (!is.null(input$SelectProt) & input$PDPFModeCheck == "PD") {
              vec <- unique(gtabMS2$Master.Protein.Descriptions[gtabMS2$Master.Protein.Descriptions %in% input$SelectProt])
              vec <- vec[!is.na(vec)]
              getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
            }
            names(gtabMS2)[names(gtabMS2)=="Master.Protein.Descriptions"] <- "Protein.Descriptions"
            
            if (is.null(filedata0()) | input$MSTrace == FALSE) {
              g <- ggplot() + 
                geom_point(data = gtabMS2, aes(x = RT.in.min, y = Precursor.MHplus.in.Da, shape = Identification, text = paste(RT.in.min, "min\n", Precursor.MHplus.in.Da, "\n", Protein.Descriptions)), alpha = 0.8, size = input$pch, col = "grey30", show.legend = FALSE) + 
                coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = TRUE)  + 
                theme_bw() + 
                scale_shape_manual(values = c(16, 1)) + 
                ylab("Molecular Weight (Da)") + 
                xlab("Retention time (min)")
              if (!is.null(input$SelectProt)) {
                g <- g + 
                  geom_point(data = gtabMS2[gtabMS2$Protein.Descriptions %in% input$SelectProt[!is.na(input$SelectProt)],], aes(x = RT.in.min, y = Precursor.MHplus.in.Da, fill = Protein.Descriptions, text = paste(RT.in.min, "min\n", Precursor.MHplus.in.Da, "Da\n", Protein.Descriptions)), shape = 21, size = input$pch+1, alpha = 0.8, stroke = 0, col = alpha("black", 1)) +
                  scale_fill_manual(values = getPalette(length(vec))) 
              }
            } else if (input$MSTrace == TRUE) { # Overlay on MS trace
              g <- g +
                geom_point(data = gtabMS2, aes(x = RT.in.min, y = Precursor.MHplus.in.Da, shape = Identification, text = paste(RT.in.min, "min\n", Precursor.MHplus.in.Da, "Da\n", Protein.Descriptions)), alpha = 0.8, size = input$pch, col = "grey30", show.legend = FALSE) + 
                theme_bw() + 
                scale_shape_manual(values = c(16, 1)) + 
                ylab("Molecular Weight (Da)") + 
                xlab("Retention time (min)")
              if (!is.null(input$SelectProt)) {
                g <- g + 
                  geom_point(data = gtabMS2[gtabMS2$Protein.Descriptions %in% input$SelectProt[!is.na(input$SelectProt)],], aes(x = RT.in.min, y = Precursor.MHplus.in.Da, fill = Protein.Descriptions, text = paste(RT.in.min, "min\n", Precursor.MHplus.in.Da, "Da\n", Protein.Descriptions)), shape = 21, size = input$pch+1, alpha = 0.8, stroke = 0, col = alpha("black", 1)) +
                  scale_fill_manual(values = getPalette(length(vec)))
              }
            }
          } else if (input$PDPFModeCheck == "PF") { # MSPathFinder
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
                ylab("Molecular Weight (Da)") + 
                xlab("Retention time (min)")
              if (!is.null(input$SelectProt)) {
                g <- g + 
                  geom_point(data = gtabMS2[gtabMS2$Protein.Descriptions %in% input$SelectProt[!is.na(input$SelectProt)],], aes(x = RT, y = Mass, fill = Protein.Descriptions, text = paste("Start:", PeakStart, "min\n", "Stop:", PeakStop, "min\n", Mass, "Da\nSignal:", intensity, "\n", Protein.Descriptions)), shape = 21, size = input$pch+1, alpha = 0.8, stroke = 0, col = alpha("black", 1)) +
                  scale_fill_manual(values = getPalette(length(vec)))
              }
            } else {
              g <- ggplot() +
                geom_point(data = gtabMS2, aes(x = RT, y = Mass, shape = Identification, text = paste("Start:", PeakStart, "min\n", "Stop:", PeakStop, "min\n", Mass, "Da\nSignal:", intensity, "\n", Protein.Descriptions)), alpha = 0.8, size = input$pch, col = "grey30", show.legend = FALSE) + 
                theme_bw() + 
                scale_shape_manual(values = c(16, 1)) + 
                ylab("Molecular Weight (Da)") + 
                xlab("Retention time (min)")
              if (!is.null(input$SelectProt)) {
                g <- g + 
                  geom_point(data = gtabMS2[gtabMS2$Protein.Descriptions %in% input$SelectProt[!is.na(input$SelectProt)],], aes(x = RT, y = Mass, fill = Protein.Descriptions, text = paste("Start:", PeakStart, "min\n", "Stop:", PeakStop, "min\n", Mass, "Da\nSignal:", intensity, "\n", Protein.Descriptions)), shape = 21, size = input$pch+1, alpha = 0.8, stroke = 0, col = alpha("black", 1)) +
                  scale_fill_manual(values = getPalette(length(vec)))
              }
            }
          } else if (input$PDPFModeCheck == "TP") { # TopPic
            gtabMS2 <- filedataMS2TP()
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
            if (is.null(filedata0()) | input$MSTrace == FALSE) { # No MS trace
              g <- ggplot() + 
                geom_point(data = gtabMS2, aes(x = RT, y = Mass, shape = Identification, text = paste(RT, "min\n", Mass, "Da\nSignal:", intensity, "\n", Protein.Descriptions)), alpha = 0.8, size = input$pch, col = "grey30", show.legend = FALSE) + 
                coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = TRUE)  + 
                theme_bw() + 
                scale_shape_manual(values = c(16, 1)) + 
                ylab("Molecular Weight (Da)") + 
                xlab("Retention time (min)")
              if (!is.null(input$SelectProt)) {
                g <- g + 
                  geom_point(data = gtabMS2[gtabMS2$Protein.Descriptions %in% input$SelectProt[!is.na(input$SelectProt)],], aes(x = RT, y = Mass, fill = Protein.Descriptions, text = paste(RT, "min\n", Mass, "Da\nSignal:", intensity, "\n", Protein.Descriptions)), shape = 21, size = input$pch+1, alpha = 0.8, stroke = 0, col = alpha("black", 1)) +
                  scale_fill_manual(values = getPalette(length(vec))) 
              }
            } else if (input$MSTrace == TRUE) { # Overlay on MS trace
              g <- g +
                geom_point(data = gtabMS2, aes(x = RT, y = Mass, shape = Identification, text = paste(RT, "min\n", Mass, "Da\nSignal:", intensity, "\n", Protein.Descriptions)), alpha = 0.8, size = input$pch, col = "grey30", show.legend = FALSE) + 
                theme_bw() + 
                scale_shape_manual(values = c(16, 1)) + 
                ylab("Molecular Weight (Da)") + 
                xlab("Retention time (min)")
              if (!is.null(input$SelectProt)) {
                g <- g + 
                  geom_point(data = gtabMS2[gtabMS2$Protein.Descriptions %in% input$SelectProt[!is.na(input$SelectProt)],], aes(x = RT, y = Mass, fill = Protein.Descriptions, text = paste(RT, "min\n", Mass, "Da\nSignal:", intensity, "\n", Protein.Descriptions)), shape = 21, size = input$pch+1, alpha = 0.8, stroke = 0, col = alpha("black", 1)) +
                  scale_fill_manual(values = getPalette(length(vec)))
              }
            }
          }
          return(g)
        }
      }
    }
  }
  
  #---------------------------------------------------------------------------
  # SECTION 4.10: PLOT OUTPUT RENDERING
  #---------------------------------------------------------------------------
  
  #---------------------------------------------------------------------------
  # 4.10.1: Interactive plot output (Plotly)
  #---------------------------------------------------------------------------
  # Render interactive plot when DataPoints mode is enabled
  output$plot1 <- renderPlotly({
    validate(
      need(!is.null(plotInput1()), '')  # Only render if plot data exists
    )
    if (input$DataPoints == TRUE) {
      g <- plotInput1() 
      g <- g + theme(legend.title = element_blank())  # Remove legend title
      
      # Convert ggplot to interactive plotly object
      p <- ggplotly(g, tooltip = "text", height = 800) %>%
        layout(dragmode = "select",           # Enable selection mode for zooming
               xaxis=list(fixedrange=TRUE),   # Disable plotly's default zoom
               yaxis=list(fixedrange=TRUE),   # Disable plotly's default zoom
               title = "") %>%                # Remove default title
        config(displayModeBar = F) %>%       # Hide plotly toolbar
        layout(margin = list(l = 110, b = 40, r = 3, t = 10, pad = -2))  # Adjust margins
      
      # Remove legend entries for identification status in MS/MS mode
      if (input$MSModeCheck == "MS2") {
        # Hide legends for all traces to avoid clutter
        p <- style(p, showlegend = FALSE, traces = seq_len(length(p$x$data)))
      }
      p
    } else {
      plotly_empty()  # Return empty plot when in static mode
    }
  })
  
  #---------------------------------------------------------------------------
  # 4.10.2: Static plot output (ggplot2)
  #---------------------------------------------------------------------------
  # Render static plot when DataPoints mode is disabled
  output$plot2 <- renderPlot({
    validate(
      need(!is.null(plotInput1()), '')  # Only render if plot data exists
    )
    if (input$DataPoints == F) {
      if (!is.null(linput())) {
        if ((linput() > 1 & input$MSModeCheck == "MS")) {
          # Multiple MS files: place legend at bottom with multiple columns
          plotInput1() + 
            theme(legend.direction ="vertical", legend.position="bottom") +
            guides(fill=guide_legend(ncol=2))
        } else if (nProtSelection() > 0 & input$MSModeCheck == "MS2") {
          # MS/MS mode with selected proteins: place legend at bottom
          plotInput1() + 
            theme(legend.direction ="vertical", legend.position="bottom") +
            guides(fill=guide_legend(ncol=2))
        } else {
          # Single file or no protein selection: place legend on right
          plotInput1()+ 
            theme(legend.direction ="vertical", legend.position="right") 
        }
      } else if (nProtSelection() > 0 & input$MSModeCheck == "MS2") {
          # MS/MS mode with selected proteins: place legend at bottom
          plotInput1() + 
            theme(legend.direction ="vertical", legend.position="bottom") +
            guides(fill=guide_legend(ncol=2))
      } else {
        # Default: place legend on right
        plotInput1()+ 
          theme(legend.direction ="vertical", legend.position="right") 
      }
    }
  }, height = 800)  # Set plot height to 800 pixels
  
  #---------------------------------------------------------------------------
  # 4.10.3: Dynamic plot UI selection
  #---------------------------------------------------------------------------
  # Dynamically switch between interactive and static plot outputs
  output$plotUI <- renderUI({
    if (input$DataPoints) {
      # Interactive mode: use plotly with click/brush interactions disabled
      plotlyOutput("plot1")
    } else {
      # Static mode: use ggplot with click, double-click, and brush interactions
      plotOutput("plot2",
                 click = "plot_click",        # Single click detection
                 dblclick = "plot_dblclick",  # Double click for zoom
                 brush = brushOpts(           # Brush selection for zoom area
                   id = "plot_brush",
                   resetOnNew = TRUE))        # Clear brush on new selection
    }
  })
  
  #---------------------------------------------------------------------------
  # 4.10.4: Static plot zoom interaction handling
  #---------------------------------------------------------------------------
  # Handle double-click zoom functionality in static plots
  observeEvent(input$plot_dblclick, {
    oldranges$x <- ranges$x  # Store current ranges before changing
    oldranges$y <- ranges$y
    
    if (input$DataPoints == FALSE) {
      brush <- input$plot_brush
      if (!is.null(brush)) {
        # Zoom to brushed area
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)       
      } else {
        # No brush: reset zoom
        ranges$x <- NULL
        ranges$y <- NULL
      }
    }
  })
  
  #---------------------------------------------------------------------------
  # SECTION 4.11: PLOT EXPORT FUNCTIONALITY
  #---------------------------------------------------------------------------
  
  #---------------------------------------------------------------------------
  # 4.11.1: PDF export handler
  #---------------------------------------------------------------------------
  # Generate PDF download with automatic filename and size adjustment
  output$Download <- downloadHandler(
    filename = function(){
      # Generate filename based on input file name and current date
      if (is.null(InputFilesMS2())) {
        paste0("VisioProt-MS_", substring(InputFileMS()$name, first = 1, last = (nchar(InputFileMS()$name)-4)), "_", Sys.Date(), ".pdf")
      } else {
        paste0("VisioProt-MS_", substring(InputFilesMS2()[[1]]$name, first = 1, last = (nchar(InputFilesMS2()[[1]]$name)-4)), "_", Sys.Date(), ".pdf")
      }
    },
    content = function(file) {
      # Calculate plot height based on content
      h <- 8  # Base height in inches
      if (nProtSelection() > 0) {
        if (!is.null(filedata0()) & input$MSTrace) {
          h <- h + 1.6  # Extra height for MS trace overlay
        }
        if (nProtSelection() > 5) {
          h <- h+(nProtSelection()*0.152)  # Extra height for protein legend
        }
      } else if (linput() > 1) {
        h <- h+(linput()*0.152)  # Extra height for file comparison legend
      }
      
      # Configure plot layout for export
      if ((linput() > 1 & input$MSModeCheck == "MS")|(nProtSelection() > 0 & input$MSModeCheck == "MS2")) {
        g <- plotInput1() + 
          theme(legend.direction ="vertical", legend.position="bottom")
      } else {
        g <- plotInput1()+ 
          theme(legend.direction ="vertical", legend.position="right") 
      }
      # Save plot to PDF with calculated dimensions
      ggsave(file, plot = g, device = "pdf", width = 10, height = h)
    })
  
  #---------------------------------------------------------------------------
  # 4.11.2: PNG export handler
  #---------------------------------------------------------------------------
  # Generate PNG download with automatic filename and size adjustment
  output$Download1 <- downloadHandler(
    filename = function(){
      # Generate filename based on input file name and current date
      if (is.null(InputFilesMS2())) {
        paste0("VisioProt-MS_", substring(InputFileMS()$name, first = 1, last = (nchar(InputFileMS()$name)-4)), "_", Sys.Date(), ".png")
      } else {
        paste0("VisioProt-MS_", substring(InputFilesMS2()[[1]]$name, first = 1, last = (nchar(InputFilesMS2()[[1]]$name)-4)), "_", Sys.Date(), ".png")
      }
    },
    content = function(file) {
      # Calculate plot height in pixels
      h <- 800  # Base height in pixels
      if (nProtSelection() > 0) {
        if (!is.null(filedata0()) & input$MSTrace) {
          h <- h + 160  # Extra height for MS trace overlay
        }
        if (nProtSelection() > 5) {
          h <- h+(nProtSelection()*15.2)  # Extra height for protein legend
        }
      } else if (linput() > 1) {
        h <- h+(linput()*15.2)  # Extra height for file comparison legend
      }
      
      # Define PNG device with specific resolution
      device <- function(..., width, height) {
        grDevices::png(..., width = 1000, height = h, res = 120)
      }
      
      # Configure plot layout for export
      if ((linput() > 1 & input$MSModeCheck == "MS")|(nProtSelection() > 0 & input$MSModeCheck == "MS2")) {
        g <- plotInput1() + 
          theme(legend.direction ="vertical", legend.position="bottom")
      } else {
        g <- plotInput1()+ 
          theme(legend.direction ="vertical", legend.position="right") 
      }
      # Save plot to PNG with calculated dimensions
      ggsave(file, plot = g, device = device)
    })
  
  #---------------------------------------------------------------------------
  # 4.11.3: SVG export handler
  #---------------------------------------------------------------------------
  # Generate SVG download with automatic filename and size adjustment
  output$Download2 <- downloadHandler(
    filename = function(){
      # Generate filename based on input file name and current date
      if (is.null(InputFilesMS2())) {
        paste0("VisioProt-MS_", substring(InputFileMS()$name, first = 1, last = (nchar(InputFileMS()$name)-4)), "_", Sys.Date(), ".svg")
      } else {
        paste0("VisioProt-MS_", substring(InputFilesMS2()[[1]]$name, first = 1, last = (nchar(InputFilesMS2()[[1]]$name)-4)), "_", Sys.Date(), ".svg")
      }
    },
    content = function(file) {
      # Calculate plot height based on content
      h <- 8  # Base height in inches
      if (nProtSelection() > 0) {
        if (!is.null(filedata0()) & input$MSTrace) {
          h <- h + 1.6  # Extra height for MS trace overlay
        }
        if (nProtSelection() > 5) {
          h <- h+(nProtSelection()*0.152)  # Extra height for protein legend
        }
      } else if (is.null(linput())) {
        if (linput() > 1) {
          h <- h+(linput()*0.152)  # Extra height for file comparison legend
        }
      }
      
      # Define SVG device
      device <- function(..., width, height) {
        grDevices::svg(..., width = 10, height = h)
      }
      
      # Configure plot layout for export
      if (is.null(linput())) {
        g <- plotInput1()+ 
          theme(legend.direction ="vertical", legend.position="right") 
      } else {
        if ((linput() > 1 & input$MSModeCheck == "MS")|(nProtSelection() > 0 & input$MSModeCheck == "MS2")) {
          g <- plotInput1() + 
            theme(legend.direction ="vertical", legend.position="bottom")
        } else {
          g <- plotInput1()+ 
            theme(legend.direction ="vertical", legend.position="right") 
        }
      }
      # Save plot to SVG with calculated dimensions
      ggsave(file, plot = g, device = device)
    })
  
} # End of server function

#=============================================================================
# STEP 5: LAUNCH THE SHINY APPLICATION
#=============================================================================

# Initialize and run the VisioProt-MS Shiny application
# This combines the UI and server components to create the complete web application
shinyApp(ui = ui, server = server)

#=============================================================================
# OPTIONAL: DEPLOYMENT COMMANDS (COMMENTED)
#=============================================================================

# Uncomment and modify the following line for deployment to shinyapps.io or other hosting platforms
# rsconnect::deployApp("T:/RRelatedWork/VisioProt-MS")

