############################################################################
# USER INTERFACE (UI) DEFINITION
############################################################################
# Define the web application's user interface using Shiny's fluidPage layout
# The UI consists of a header, mode selection, sidebar controls, and main plot panel

# Core Shiny framework for building interactive web applications
library(shiny)

# Visualization and UI enhancement packages
library(shinyBS)       # Bootstrap components for Shiny (tooltips, modals, etc.)

ui <- fluidPage(
  #--------------------------------------------------------------------------
  # APPLICATION HEADER
  #--------------------------------------------------------------------------
  # Create header row with application title and help button
  fluidRow(
    # Application title in a responsive column
    column(3, titlePanel("VisioProt-MS")),
    
    # Help button that opens documentation in a new browser tab
    column(1, actionButton(inputId='ab1', label="",
                           icon = icon("fa-solid fa-circle-question"),
                           onclick ="window.open('https://masstools.ipbs.fr/visioprothelp.html', '_blank')",
                           style="color: #fff; background-color: #673a49; border-color: #000000"),
           bsTooltip("ab1", 
                     "go to help page",
                     placement = "bottom")
    )
  ),
  
  # Custom CSS styling for the help button
  tags$style(type='text/css', "#ab1 { width:80%; margin-top: 25px; font-family : Cursive; font-weight: 900; font-size: 160%;}"),
  
  #--------------------------------------------------------------------------
  # MODE SELECTION
  #--------------------------------------------------------------------------
  # Allow users to choose between MS-only visualization or MS/MS overlay mode
  #checkboxInput("MSModeCheck", "MS data only", TRUE), # Legacy checkbox implementation
  
  # Radio buttons for mode selection (replaces checkbox for better UX)
  radioButtons("MSModeCheck", "MS mode:",
               c("MS" = 'MS',        # MS data visualization only
                 "MS/MS" = 'MS2'),   # MS data with MS/MS identification overlay
               selected = 'MS',      # Default to MS mode
               inline = TRUE         # Display options horizontally
  ),
  
  # Tooltip explaining the mode selection options
  bsTooltip("MSModeCheck", 
            "Choose between plotting MS data only or overlaying results of Top-Down searches",
            "right"),
  
  #--------------------------------------------------------------------------
  # MAIN LAYOUT: SIDEBAR AND PLOT PANEL
  #--------------------------------------------------------------------------
  # Create two-panel layout with sidebar controls and main plotting area
  sidebarLayout( 
    #------------------------------------------------------------------------
    # SIDEBAR PANEL: USER CONTROLS AND FILE INPUTS
    #------------------------------------------------------------------------
    sidebarPanel(width = 4,
                 
                 #====================================================================
                 # CONDITIONAL PANELS: DIFFERENT CONTROLS FOR EACH MODE
                 #====================================================================
                 
                 #--------------------------------------------------------------------
                 # MS MODE PANEL: For MS data visualization only
                 #--------------------------------------------------------------------
                 conditionalPanel(condition="input.MSModeCheck== 'MS'",
                                  
                                  #--------------------------------------------------------
                                  # File Input Section for MS Mode
                                  #--------------------------------------------------------
                                  # Allow users to upload deconvoluted MS data files
                                  fileInput("fileMS", "Select input file(s):",  
                                            accept = c(
                                              "text/csv",                              # CSV files
                                              "text/comma-separated-values,text/plain", # Text files  
                                              ".txt",                                  # Text files
                                              ".ms1ft",                               # ProMex format
                                              ".csv",                                 # CSV files
                                              ".msalign"),                            # TopPIC format
                                            multiple = T,      # Allow multiple files for comparison
                                            width = "100%"
                                  ),
                                  
                                  #--------------------------------------------------------
                                  # Test Mode Controls for MS Mode
                                  #--------------------------------------------------------
                                  # Enable test mode to try the application with sample data
                                  checkboxInput("TestModeCheck", "Using test mode", FALSE),
                                  bsTooltip("TestModeCheck", 
                                            "Check to test the application without uploading any file. Then click on a button to upload a single test file or several for overlay",
                                            "right"),
                                  
                                  # Show test file buttons only when test mode is enabled
                                  conditionalPanel(condition="input.TestModeCheck==true",
                                                   actionButton("TestFile1", "Single test file"),    # Load one sample file
                                                   actionButton("TestFile2", "Multiple test files")  # Load multiple files for comparison
                                  )
                                  
                 ),
                 
                 #--------------------------------------------------------------------
                 # MS/MS MODE PANEL: For MS data with identification overlay
                 #--------------------------------------------------------------------
                 conditionalPanel(condition="input.MSModeCheck== 'MS2'", 
                                  
                                  #--------------------------------------------------------
                                  # Test Mode Controls for MS/MS Mode
                                  #--------------------------------------------------------
                                  checkboxInput("MS2TestModeCheck", "Using test mode", FALSE),
                                  bsTooltip("MS2TestModeCheck", 
                                            "Check to test the application without uploading any file",
                                            "right"),
                                  
                                  # Display warning message when in test mode
                                  conditionalPanel(condition = "input.MS2TestModeCheck==true",
                                                   tags$span(style="color:red", "You are in test mode. Click on a button to select a single test file or multiple test files.\nUncheck to exit and upload your own data."),
                                                   br()
                                  ),
                                  
                                  #--------------------------------------------------------
                                  # File Input Section for MS/MS Mode (User Data)
                                  #--------------------------------------------------------
                                  conditionalPanel(condition = "input.MS2TestModeCheck==false",
                                                   
                                                   # MS data file input (background for MS/MS overlay)
                                                   fileInput("fileMS2", "Select input file for MS:",  
                                                             accept = c(
                                                               "text/csv",
                                                               "text/comma-separated-values,text/plain",
                                                               ".csv",
                                                               ".txt",
                                                               ".ms1ft",    # ProMex format
                                                               ".msalign"), # TopPIC format
                                                             multiple = F,  # Single file only for background MS
                                                             width = "100%"
                                                   ),
                                                   
                                                   #------------------------------------------------
                                                   # Software Selection for MS/MS Analysis
                                                   #------------------------------------------------
                                                   # Choose which top-down software was used for analysis
                                                   radioButtons("PDPFModeCheck", "Origin of the MS/MS files:",
                                                                c("Proteome Discoverer" = 'PD',  # Thermo Fisher
                                                                  "MSPathFinder" = 'PF',         # PNNL
                                                                  "TopPIC" = 'TP'),             # UCSD
                                                                selected = 'PD',
                                                                inline = TRUE
                                                   ),
                                                   bsTooltip("PDPFModeCheck", 
                                                             "Choose the software utilized for analysing of the top-down data",
                                                             "right"),
                                                   
                                                   #============================================
                                                   # PROTEOME DISCOVERER FILE INPUTS
                                                   #============================================
                                                   # Proteome Discoverer requires two files:
                                                   # 1. MSMSSpectrumInfo: MS/MS spectrum metadata
                                                   # 2. PSM: Peptide spectrum matches with protein IDs
                                                   conditionalPanel(condition = "input.PDPFModeCheck== 'PD'", 
                                                                    fileInput("MS2file", "Choose MSMSSpectrumInfo File:", 
                                                                              accept = c(
                                                                                "text/csv",
                                                                                "text/comma-separated-values,text/plain",
                                                                                ".txt")
                                                                    ),
                                                                    fileInput("PSMfile", "Choose PSM File:", 
                                                                              accept = c(
                                                                                "text/csv",
                                                                                "text/comma-separated-values,text/plain",
                                                                                ".txt")
                                                                    )),
                                                   
                                                   #============================================
                                                   # MSPATHFINDER FILE INPUTS
                                                   #============================================
                                                   # MSPathFinder requires IcTarget or IcTda results file
                                                   conditionalPanel(condition = "input.PDPFModeCheck== 'PF'", 
                                                                    fileInput("MS2filePF", "Choose IcTarget or IcTda File from MSPathFinder:", 
                                                                              accept = c(
                                                                                "text/csv",
                                                                                "text/comma-separated-values,text/plain",
                                                                                ".tsv")  # Tab-separated values
                                                                    )
                                                   ),
                                                   
                                                   #============================================
                                                   # TOPPIC FILE INPUTS
                                                   #============================================
                                                   # TopPIC requires two files:
                                                   # 1. MS/MS data file (.msalign format)
                                                   # 2. Identification results (OUTPUT_TABLE)
                                                   conditionalPanel(condition = "input.PDPFModeCheck== 'TP'", 
                                                                    fileInput("MS2fileTP", "Choose MS/MS File:", 
                                                                              accept = c(
                                                                                "text",
                                                                                "text/comma-separated-values,text/plain",
                                                                                ".msalign")  # TopPIC MS/MS format
                                                                    ),
                                                                    fileInput("IDfileTP", "Choose the OUTPUT_TABLE/_ms2_toppic (saved at tab-delimited .txt) from TopPIC:", 
                                                                              accept = c(
                                                                                "text",
                                                                                "text/comma-separated-values,text/plain",
                                                                                ".OUTPUT_TABLE")  # TopPIC results format
                                                                    )
                                                   )
                                  ),
                                  
                                  #--------------------------------------------------------
                                  # Protein Selection and Display Options for MS/MS Mode
                                  #--------------------------------------------------------
                                  # Dropdown to select which identified proteins to highlight in the plot
                                  selectInput("SelectProt", "Select the ID to highlight:", 
                                              NULL,          # Options will be populated dynamically based on loaded data
                                              multiple = TRUE), # Allow multiple protein selection
                                  bsTooltip("SelectProt", 
                                            "Select among the identified proteins which one(s) to highlight on the plot",
                                            "right"),
                                  
                                  # Option to hide MS/MS spectra without protein identification
                                  checkboxInput("HideMSMS", "Hide MS/MS withouth ID", FALSE),
                                  bsTooltip("HideMSMS", 
                                            "Removes the MS/MS spectra from the top-down analysis that were not matched to a protein.",
                                            "right"),
                                  
                                  # Option to show/hide the underlying MS1 trace
                                  checkboxInput("MSTrace", "Display the MS trace", TRUE),
                                  bsTooltip("MSTrace", 
                                            "Adds the MS trace to the plot.",
                                            "right")
                 ),
                 
                 #====================================================================
                 # COMMON CONTROLS: Available in both MS and MS/MS modes
                 #====================================================================
                 
                 # Toggle between static ggplot2 and interactive plotly visualization
                 checkboxInput("DataPoints", "Show data labels (slower)", FALSE),
                 bsTooltip("DataPoints", 
                           "Switch to \"data\" mode: data appears on hovering",
                           "right"),
                 
                 #--------------------------------------------------------------------
                 # PLOT CUSTOMIZATION PARAMETERS
                 #--------------------------------------------------------------------
                 # Arrange plot parameters in a responsive row layout
                 fluidRow(
                   #------------------------------------------------------------------
                   # Color Scale Selection (Column 1)
                   #------------------------------------------------------------------
                   column(5,
                          # Color palette selection - options change based on number of files loaded
                          selectInput("colourscale", "Colour scale:",
                                      c("Spectral" = "Spectral",           # For single files: continuous scales
                                        "Red/yellow/blue" = "RdYlBu", 
                                        "Red/yellow/green" = "RdYlGn",
                                        "yellow to red" = "YlOrRd"
                                      )),
                          bsTooltip("colourscale", 
                                    "Select the colour scale for the MS data.",
                                    "right")),
                   
                   #------------------------------------------------------------------
                   # Point Size Control (Column 2)
                   #------------------------------------------------------------------
                   column(3,
                          numericInput("pch", label = "Point size:", value = 1, min = 0.1, step = 0.1, max = 10),
                          bsTooltip("pch", 
                                    "Define the size of the point (from 0.1 to 10).",
                                    "right")),
                   
                   #------------------------------------------------------------------
                   # Intensity Threshold Control (Column 3)
                   #------------------------------------------------------------------
                   column(4,
                          numericInput("IntensityThresh", label = "Threshold:", value = 20, min = 0, max = 100, step = 1),
                          bsTooltip("IntensityThresh", 
                                    "Define the percentage of highest intensity features of the MS data to display.",
                                    "right"))
                 ),
                 
                 #--------------------------------------------------------
                 # Signal sum for MS Mode
                 #--------------------------------------------------------
                 # Click to calculate the sum of signal in plot window
                 # conditionalPanel(condition="input.MSModeCheck== 'MS' && output.filetype.RoWinPro == 0 && !(output.filetype.BioPharma == 0 && output.filetype.ProMex == 0)",
                 conditionalPanel(condition="input.MSModeCheck== 'MS' && output.show_calcsum",
                                  actionButton("CalcSum", "Sum signal in plot",
                                               style="color: #fff; background-color: #3c3c3c; border-color:  #3c3c3c"),
                                  bsTooltip("CalcSum", 
                                            "Clic to calculate sum of signal of plotted points",
                                            "right"),
                                  actionButton("CalcSumReset", "Reset",
                                               style="color: #fff; background-color: #3c3c3c; border-color:  #3c3c3c"),
                                  bsTooltip("CalcSumReset", 
                                            "Clic to reset values",
                                            "right"),
                                  tableOutput('SumSignal')
                 ),
                 
                 #--------------------------------------------------------------------
                 # ZOOM AND NAVIGATION CONTROLS
                 #--------------------------------------------------------------------
                 # Dynamic instructions that change based on plot type (static vs interactive)
                 hr(),
                 htmlOutput("ZoomParam"),
                 br(),
                 
                 # Zoom control buttons
                 actionButton("DeZoom", "Unzoom one step", style='padding:8px; font-size:100%'),
                 bsTooltip("DeZoom", 
                           "Unzoom to previous window (only once).",
                           "right"),
                 actionButton("TotalDeZoom", "Total unzoom", style='padding:8px; font-size:100%'),
                 bsTooltip("TotalDeZoom", 
                           "Total unzoom.",
                           "right"),
                 br(),
                 br(),
                 
                 #--------------------------------------------------------------------
                 # EXPORT CONTROLS: Download plot in different formats
                 #--------------------------------------------------------------------
                 downloadButton("Download", "Download .pdf"),   # PDF format for publications
                 downloadButton("Download1", "Download .png"),  # PNG format for presentations  
                 downloadButton("Download2", "Download .svg"),  # SVG format for vector graphics
                 downloadButton("Download3", "Download table"),  # tab-delimited file with values
                 
                 #--------------------------------------------------------------------
                 # CITATION INFORMATION
                 #--------------------------------------------------------------------
                 HTML(paste("<br/><br/>",
                            h5("If you use VisioProt-MS for your research please cite:"),
                            "Marie Locard-Paulet, Julien Parra, Renaud Albigot, Emmanuelle Mouton-Barbosa, Laurent Bardi, Odile Burlet-Schiltz, Julien Marcoux; VisioProt-MS: interactive 2D maps from intact protein mass spectrometry, Bioinformatics, bty680, https://doi.org/10.1093/bioinformatics/bty680")
                 )
    ),
    
    #--------------------------------------------------------------------------
    # MAIN PANEL: Plot display area
    #--------------------------------------------------------------------------
    # Main plotting area that adapts based on user's visualization choice
    # (static ggplot2 vs interactive plotly)
    mainPanel(width = 8,
              uiOutput("plotUI"),  # Dynamic UI element for plot output
              ## if MS1 panel, plot the sum of intensities in the plot window
              # if (print_sum_results) { #TODO: conditional: only if MS trace panel 
              uiOutput("sum_results")  # Dynamic UI element for printing intensity summarization results
              # }
    )
    
  ),
  
  #----------------------------------------------------------------------------
  # FOOTER
  #----------------------------------------------------------------------------
  # Application footer with copyright information
  tabsetPanel(
    tabPanel(
      HTML('<footer><font size="0.8">copyright 2017 - CNRS - All rights reserved - VisioProt-MS V2.4</font></footer>')
    )
  )
)