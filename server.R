# Advanced plotting libraries
# devtools::install_github('hadley/ggplot2')  # Use development version if needed
library(ggplot2)    # Grammar of graphics for static plots
library(plotly)     # Interactive web-based data visualization
library(dplyr)      # Data manipulation and transformation

# Visualization and UI enhancement packages
library(RColorBrewer)  # Color palettes for data visualization
library(shinyBS)       # Bootstrap components for Shiny (tooltips, modals, etc.)
library(data.table)    # High-performance data manipulation and reading

############################################################################
# EXTERNAL HELPER FUNCTIONS
############################################################################
# Load custom R functions for data processing and format handling
# These functions are stored in separate files for modularity and reusability

# Filter data to keep only the highest intensity features for each mass/RT pair
source("files/Rfiles/Keep_highest_signal.R", local = TRUE)

# Custom function to efficiently combine multiple data frames
source("files/Rfiles/custom_rbindlist.R", local = TRUE)

# Parse and process TopPIC software output files
source("files/Rfiles/Parse_TopPIC_input.R", local = TRUE)

# Standardize column names across different deconvolution software formats
source("files/Rfiles/rename_input_coulumns.R", local = TRUE)

# Check input files with adapted error messages
source("files/Rfiles/check_input_tables.R", local = TRUE)


############################################################################
# SERVER LOGIC
############################################################################
# Define the server-side logic that handles user inputs, processes data,
# and generates reactive outputs for the user interface

server <- function(input, output, clientData, session) {
  
  #--------------------------------------------------------------------------
  # SECTION 1: DYNAMIC UI MODIFICATIONS
  #--------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------
  # 1.1: Update zoom instructions based on plot type
  #--------------------------------------------------------------------------
  # Provide different zoom instructions for static vs interactive plots
  output$ZoomParam <- renderUI({
    if (input$DataPoints) {
      # Instructions for interactive plotly plots
      HTML("<h5>Zoom in: select the ranges of interest.<br/>Zoom out: Click on the unzoom button.</h5>")
    } else {
      # Instructions for static ggplot2 plots
      HTML("<h5>Zoom in: select the ranges of interest and double click.<br/>Zoom out: double click.</h5>")
    }
  })
  
  #--------------------------------------------------------------------------
  # 1.2: Dynamic color scale updates based on number of files
  #--------------------------------------------------------------------------
  # Automatically switch between continuous and qualitative color palettes
  # based on whether single or multiple files are loaded
  observe({
    if (!is.null(linput())) {
      x <- linput()  # Get number of input files
      if (x > 1) {
        # Multiple files: use qualitative color schemes for file comparison
        updateSelectInput(session, "colourscale",
                          "Colour scale:",
                          c("Paired" = "Paired",     # Qualitative palettes for distinguishing files
                            "Set1" = "Set1",
                            "Set2" = "Set2",
                            "Set3" = "Set3",
                            "Set4" = "Dark2", 
                            "Set5" = "Accent"
                          ))
      } else {
        # Single file: use continuous color schemes for intensity visualization
        updateSelectInput(session, "colourscale",
                          "Colour scale:",
                          c("Spectral" = "Spectral",        # Continuous palettes for intensity
                            "Red/yellow/blue" = "RdYlBu", 
                            "Red/yellow/green" = "RdYlGn",
                            "yellow to red" = "YlOrRd"
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
  
  #--------------------------------------------------------------------------
  # 1.3: Dynamic protein ID selection for MS/MS mode
  #--------------------------------------------------------------------------
  # Populate protein selection dropdown based on loaded identification data
  observe({
    # For Proteome Discoverer data: extract protein descriptions from PSM file
    if (!is.null(filedataMS2())) {
      if (length(filedataMS2()$PSMfile$Master.Protein.Descriptions[!is.na(filedataMS2()$PSMfile$Master.Protein.Descriptions)]) > 0) {
        updateSelectInput(session, "SelectProt",
                          "Select the ID to highlight:",
                          sort(unique(filedataMS2()$PSMfile$Master.Protein.Descriptions[!is.na(filedataMS2()$PSMfile$Master.Protein.Descriptions)]))
        )
      }
    }
    
    # For MSPathFinder data: extract protein descriptions
    if (!is.null(filedataMS2PF())) {
      if (length(filedataMS2PF()$Protein.Descriptions[!is.na(filedataMS2PF()$Protein.Descriptions)]) > 0) {
        updateSelectInput(session, "SelectProt",
                          "Select the ID to highlight:",
                          sort(unique(filedataMS2PF()$Protein.Descriptions[!is.na(filedataMS2PF()$Protein.Descriptions)]))
        )
      }
    }
    
    # For TopPIC data: extract protein descriptions
    if (!is.null(filedataMS2TP())) {
      if (length(filedataMS2TP()$Protein.Descriptions[!is.na(filedataMS2TP()$Protein.Descriptions)]) > 0) {
        updateSelectInput(session, "SelectProt",
                          "Select the ID to highlight:",
                          sort(unique(filedataMS2TP()$Protein.Descriptions[!is.na(filedataMS2TP()$Protein.Descriptions)]))
        )
      }
    }
  })
  
  #--------------------------------------------------------------------------
  # 1.4: Track number of selected proteins for plot layout
  #--------------------------------------------------------------------------
  # Count selected proteins to adjust plot dimensions for export
  nProtSelection <- reactiveVal(0)
  observeEvent(c(plotInput1, input$SelectProt), {
    nProtSelection(length(input$SelectProt))
  })
  
  #--------------------------------------------------------------------------
  # SECTION 2: FILE INPUT HANDLING
  #--------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------
  # 2.1: MS file input processing
  #--------------------------------------------------------------------------
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
  
  #--------------------------------------------------------------------------
  # 2.2: MS/MS file input processing for different software platforms
  #--------------------------------------------------------------------------
  
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
  
  #--------------------------------------------------------------------------
  # SECTION 3: TEST MODE AND FILE TYPE DETECTION
  #--------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------
  # 3.1: Test mode variables and file type tracking
  #--------------------------------------------------------------------------
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
  
  #--------------------------------------------------------------------------
  # 3.2: Test mode file loading
  #--------------------------------------------------------------------------
  
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
  
  #--------------------------------------------------------------------------
  # 3.3: Automatic file type detection and counting
  #--------------------------------------------------------------------------
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
        l[[i]] <- is_biopharama_MS_input(input = InputFileMS[i, 'datapath'])
        
        # Bruker detection: check for specific header pattern
        l2[[i]] <- is_Bruker_MS_input(input = InputFileMS[i, 'datapath'])
        
        # ProMex detection: check file extension
        l3[[i]] <- is_ProMex_MS_input(input = InputFileMS$name[i])
        
        # TopPIC detection: check file extension  
        l4[[i]] <- is_TopFD_input(input = InputFileMS$name[i])
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
  
  #--------------------------------------------------------------------------
  # SECTION 4: STATE MANAGEMENT AND RESET FUNCTIONS
  #--------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------
  # 4.1: Reset plot when switching modes or test settings
  #--------------------------------------------------------------------------
  
  # Reset when exiting test mode in MS mode
  observeEvent(input$TestModeCheck, {
    InputFileMS(NULL)
    testfileinput(0)
    ranges$x <- NULL  # Clear zoom settings
    ranges$y <- NULL
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
  
  #--------------------------------------------------------------------------
  # SECTION 5: FILE TYPE VALIDATION AND FORMAT DETECTION
  #--------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------
  # 5.1: Comprehensive file format validation
  #--------------------------------------------------------------------------
  # Reactive function to validate and classify input file formats
  ftype <- reactive({
    if (is.null(InputFileMS()) & testfileinput() == 0) {
      return(NULL)
    } else {
      InputFileMS <- InputFileMS()
      l <- list()
      
      # Check each uploaded file
      for(i in 1:nrow(InputFileMS)){
        
        # Initial format detection
        val <- is_biopharama_MS_input(input = InputFileMS[i, 'datapath'])
        val2 <- is_Bruker_MS_input(input = InputFileMS[i, 'datapath'])  # Bruker format
        val3 <- is_ProMex_MS_input(input = InputFileMS$name[i]) # ProMex format
        val4 <- is_TopFD_input(InputFileMS$name[i]) # TopPIC format
        
        # Classify file type based on detection results
        val <- ifelse(val,  "BioPharma", "RoWinPro")
        val[val2] <- "Bruker"
        val[val3] <- "ProMex"
        val[val4] <- "TopPic"
        
        #=======================================================================
        # Format-specific validation checks
        #=======================================================================
        
        # # BioPharma format validation 
        # if (val == "BioPharma") { 
        #
        # }
        
        # RoWinPro format validation (should have exactly 4 columns)
        if (val == "RoWinPro") {
          RoWinPro_MS_check(input = InputFileMS[i, 'datapath'])
        }
        
        # Bruker format validation
        if (val == "Bruker") {
          Bruker_MS_check(input = InputFileMS[i, 'datapath'])
        }
        
        # ProMex format validation
        if (val == "ProMex") {
          ProMex_MS_check(input = InputFileMS[i, 'datapath'])
        }
        
        # TopPIC format validation
        if (val == "TopPic") {
          TopFD_MS_check(input = InputFileMS[i, 'datapath'])
        }
        
        l[[i]] <- val
      }
      vec <- unlist(l)
      # print(vec)
      return(vec)
    }
  })
  
  #--------------------------------------------------------------------------
  # SECTION 6: DATA LOADING AND PROCESSING FUNCTIONS
  #--------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------
  # 6.1: Main MS data loading function
  #--------------------------------------------------------------------------
  # Load and process MS data files based on format type and test mode
  filedata0 <- reactive({
    if (testfileinput() == 0) { 
      # Regular file upload mode (not test mode)
      validate(
        need((input$TestModeCheck==FALSE & input$MSModeCheck == "MS") | (input$MS2TestModeCheck == FALSE & input$MSModeCheck == "MS2"), 
             "You are in test mode. Click on a button to select a single test file or multiple test files.\nUncheck to exit and upload your own data."
        )
      )
      if (is.null(InputFileMS())) {
        # User has not uploaded a file yet
        return(NULL)
      } else {
        # Validation: prevent mixing different file types with multiple files
        validate(
          need(!(max(table(ftype())) > 1 & length(unique(ftype())) > 1), 
               "With multiple files, all files need to be from the same deconvolution software."
          )
        )
        InputFileMS <- InputFileMS()
        lfiles <- list()
        
        # Process each uploaded file according to its detected format
        for(i in 1:nrow(InputFileMS)){
          if (ftype()[i] == "BioPharma") {
            # BioPharma Finder format: read tab-delimited with headers
            lfiles[[i]] <- read.table(InputFileMS[i, 'datapath'], sep = "\t", header = T)
            # Standardize mass column name (could be Monoisotopic.Mass or Average.Mass)
            names(lfiles[[i]])[names(lfiles[[i]])=="Monoisotopic.Mass"|names(lfiles[[i]])=="Average.Mass"] <- "Mass"
            # Map to standard column format: RT, Mass, Intensity, Start Time, Stop Time
            lfiles[[i]] <- lfiles[[i]][,c("Apex.RT", "Mass", "Sum.Intensity", "Start.Time..min.", "Stop.Time..min.")]
            
          } else if (ftype()[i] == "RoWinPro")  {
            # RoWinPro format: simple tab-delimited without headers
            lfiles[[i]] <- read.table(InputFileMS[i, 'datapath'], sep = "\t", header = F)
            # Add placeholder columns for consistency with BioPharma format
            lfiles[[i]] <- cbind(lfiles[[i]][,1:3], " Temp1" = rep(NA, nrow(lfiles[[i]])), "Temp2" = rep(NA, nrow(lfiles[[i]])))
            
          } else if (ftype()[i] == "Bruker")  {
            # Bruker DataAnalysis format: comma-separated, skip header lines
            lfiles[[i]] <- read.table(InputFileMS[i, 'datapath'], sep = ",", header = F, skip = 2)
            lfiles[[i]] <- cbind(lfiles[[i]][,2:4], " Temp1" = rep(NA, nrow(lfiles[[i]])), "Temp2" = rep(NA, nrow(lfiles[[i]]))) # add one more column to allow row binding later on 
          } else if (ftype()[i] == "ProMex")  { # ProMex output
            lfiles[[i]] <- read.table(InputFileMS[i, 'datapath'], sep = "\t", header = T)
            lfiles[[i]] <- cbind("RT" = (lfiles[[i]]$MinElutionTime + ((lfiles[[i]]$MaxElutionTime - lfiles[[i]]$MinElutionTime)/2)), lfiles[[i]][,c("MonoMass", "ApexIntensity", "MinElutionTime", "MaxElutionTime")]) # Map the columns as in RoWinPro format, but with start and stop instead of all the points of the peak. I add a first column with the middle of the peak for zooming (plotly_select needs points, not ranges).
          } else if (ftype()[i] == "TopPic")  { # TopPic output
            lfiles[[i]] <- TopPicMS1Parsing(InputFileMS[i,'datapath'])
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
      lfiles <- list()
      for(i in 1){
        lfiles[[i]] <- read.table(infile[i], sep = "\t", header = F)
      }
      names(lfiles) <- c("test data")
      return(lfiles)
    }
  })
  
  #--------------------------------------------------------------------------
  # 6.2: MS/MS data loading for Proteome Discoverer
  #--------------------------------------------------------------------------
  # Load and process Proteome Discoverer MS/MS identification data
  filedataMS2 <- reactive({
    if (is.null(InputFilesMS2()) & testfileinput() != 3) {
      return(NULL)  # No MS/MS files uploaded
    } else {
      if (testfileinput() == 3) {
        # Test mode: load sample MS/MS files
        infileMS2 <- list.files("files/MS2/", pattern = "MSMS", full.names = T)[[1]]
        infilePSM <- list.files("files/MS2/", pattern = "SMs.txt", full.names = T)[[1]]
        PSM <- read.table(infilePSM, sep = "\t", header = T, comment.char = "#")
        MS2 <- read.table(infileMS2, sep = "\t", header = T)
      } else {
        # User uploaded files: read PSM and MS/MS spectrum files
        PSM <- read.table(InputFilesMS2()$PSMfile$datapath, sep = "\t", header = T)
        MS2 <- read.table(InputFilesMS2()$MS2file$datapath, sep = "\t", header = T, comment.char = "#")
        
        # Validate file formats for Proteome Discoverer compatibility
        PD_MS2_check(PSM_tab = PSM, MSMS_tab = MS2)
        
        PSM <- RenamePDPrSM(PSM)
        MS2 <- RenamePDMSMS(MS2)
      }
    }
    return(list("MS2file" = MS2, "PSMfile" = PSM))
  })
  
  #--------------------------------------------------------------------------
  # 6.3: MS/MS data loading for MSPathFinder
  #--------------------------------------------------------------------------
  # Load and process MSPathFinder identification data
  filedataMS2PF <- reactive({
    if (is.null(input$MS2filePF) | is.null(input$fileMS2)) {
      return(NULL)
    } else if (input$MSModeCheck == "MS2" & input$PDPFModeCheck == "PF")  {
      # Validate that both MS and MS/MS files are uploaded
      validate(
        need(!is.null(InputFileMS()), 
             "You need to upload the MS2 with the associated MS file to plot MS2 results from MSPathFinder.")
      )
      
      # File format validation
      MSPathFinder_MS2_check(inputMSMS = InputFilesMS2PF(), inputMS = InputFileMS())
      
      # Load MSPathFinder results
      MS2PF <- read.table(InputFilesMS2PF()$datapath, sep = "\t", header = F, skip = 1)
      
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
  
  #--------------------------------------------------------------------------
  # 6.4: MS/MS data loading for TopPIC
  #--------------------------------------------------------------------------
  # Load and process TopPIC identification data
  filedataMS2TP <- reactive({
    if (is.null(InputFilesMS2TP())) {
      return(NULL)
    } else if (input$MSModeCheck == "MS2" & input$PDPFModeCheck == "TP")  {
      # Validate TopPIC file formats
      TopFD_MS2_check(input = InputFilesMS2TP())
      
      # Parse TopPIC files using custom functions
      IDTP <- TopPicIDParsing(InputFilesMS2TP()$IDfile$datapath)
      MS2TP <- TopPicMS2Parsing(InputFilesMS2TP()$MS2file$datapath)
      
      # Standardize column names and merge data
      names(IDTP)[names(IDTP) == "Spectrum ID"] <- "Scan"
      dat <- merge(MS2TP, IDTP, by = "Scan", all = T)
      names(dat)[names(dat)=="Protein accession"] <- "Protein.Descriptions"
      dat$Mass <- as.numeric(dat$Mass)
      dat$Identification <- ifelse(!is.na(as.character(dat$Protein.Descriptions)), "IDed", "Not IDed")
      return(dat)
    }
  })
  
  #--------------------------------------------------------------------------
  # SECTION 7: DATA PROCESSING AND PREPARATION FOR PLOTTING
  #--------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------
  # 7.1: Final data processing and formatting for visualization
  #--------------------------------------------------------------------------
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
      lfiles <- ThresholdCleaning(lfiles, input$IntensityThresh)
      
      if (filetype$BioPharma == 0 & filetype$ProMex == 0) {
        # Only RoWinPro/Bruker/TopPIC files: rename columns and combine
        l <- list()
        for (i in seq_along(lfiles)) {
          l[[i]] <- RenameBioPharma(lfiles[[i]])
        }
        names(l) <- names(lfiles)
        lfiles <- l
        return(RBindList(lfiles))
        
      } else if (filetype$RoWinPro == 0) {
        # Only BioPharma or ProMex files: standardize and combine
        lfiles <- lapply(lfiles, function(x) {
          RenameBioPharma(x)
        })
        return(RBindList(lfiles))
        
      } else {
        # Mixed file types: standardize but keep separate
        lfiles <- lapply(lfiles, function(x) {
          RenameBioPharma(x) # ProMex files will have an empty RT column
        })
        return(lfiles)
      }
    }
  }
  
  # ====================
  # PLOT RANGE MANAGEMENT
  # ====================
  # Reactive values to store the current plot ranges for zooming functionality
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  # Observer for the "De-Zoom" button - restores previous zoom level
  observeEvent(input$DeZoom, {
    ranges$x <- oldranges$x
    ranges$y <- oldranges$y
  })
  
  # Observer for the "Total De-Zoom" button - resets to full data range
  observeEvent(input$TotalDeZoom, {
    # Scenario 1: Only MS data loaded (no MS/MS identification data)
    if (is.null(filedataMS2())) {
      if (class(filedata()) != "list") { # Single data table
        # For ProMex and BioPharma formats, use column 5 for RT, column 2 for mass
        if (filetype$ProMex > 0 | filetype$BioPharma > 0) {
          ranges$x <- c(0, range(filedata()[,5])[2])  # RT range starting from 0
          ranges$y <- range(filedata()[,2])           # Mass range
        } else { # For other formats, use column 1 for RT
          ranges$x <- c(0, range(filedata()[,1])[2])  # RT range starting from 0
          ranges$y <- range(filedata()[,2])           # Mass range
        }
      } else { # Multiple data tables (mixed file types)
        # Extract retention time values from all possible sources
        x1 <- sapply(filedata(), function(z) {
          z$RT
        })
        x2 <- sapply(filedata(), function(z) {
          z$PeakStart
        })
        x3 <- sapply(filedata(), function(z) {
          z$PeakStop
        })
        # Combine all time values and remove NAs
        x <- c(unlist(x1), unlist(x2), unlist(x3))
        x <- x[!is.na(x)]
        # Extract mass values
        y <- sapply(filedata(), function(z) {
          z$Mass
        })
        ranges$x <- c(0, range(x)[2])  # Combined RT range
        ranges$y <- range(y)           # Combined mass range
      }
    } else if (is.null(filedata())) { # Scenario 2: Only MS/MS data loaded (no MS trace)
      if (input$PDPFModeCheck == "PD") { # Proteome Discoverer format
        ranges$x <- c(0, range(filedataMS2()$MS2file$RT.in.min)[2])
        ranges$y <- range(filedataMS2()$MS2file$Precursor.MHplus.in.Da)
      } else { # TopPIC format
        ranges$x <- c(0, range(filedataMS2TP()$MS2file$RT)[2])
        ranges$y <- range(filedataMS2TP()$MS2file$Mass)
      }
    } else { # Scenario 3: Both MS and MS/MS data loaded (overlay mode)
      # Calculate combined ranges from both data sources
      if (filetype$ProMex > 0 | filetype$BioPharma > 0) {
        # Combine RT and mass ranges from both MS and MS/MS data
        ranges$x <- c(0, range(c(filedataMS2()$MS2file$RT.in.min, filedata()[,5]))[2])
        ranges$y <- range(c(filedataMS2()$MS2file$Precursor.MHplus.in.Da, filedata()[,2]))
      } else if (input$PDPFModeCheck == "PD") { # Proteome Discoverer overlay
        ranges$x <- c(0, range(c(filedataMS2()$MS2file$RT.in.min, filedata()[,1]))[2])
        ranges$y <- range(c(filedataMS2()$MS2file$Precursor.MHplus.in.Da, filedata()[,2]))
      } else if (input$PDPFModeCheck == "TP") { # TopPIC overlay
        ranges$x <- c(0, range(c(filedataMS2TP()$MS2file$RT, filedata()[,1]))[2])
        ranges$y <- range(c(filedataMS2()$MS2file$Mass, filedata()[,2]))
      }
    }
  })
  
  # Reactive values to store previous zoom ranges for "De-Zoom" functionality
  oldranges <- reactiveValues(x = NULL, y = NULL)
  
  # Observer for interactive plot selection events (brush/lasso selection for zooming)
  observeEvent(event_data("plotly_selected"), {
    # Store current ranges before updating (for De-Zoom functionality)
    oldranges$x <- ranges$x
    oldranges$y <- ranges$y
    
    # Only process selection if DataPoints mode is enabled
    if (input$DataPoints) {
      newdata <- event_data("plotly_selected")
      
      # Process valid selection data
      if (!is.null(newdata) & class(newdata)=="data.frame") {
        
        # Case 1: Single table with ProMex/BioPharma format
        if (class(filedata()) != "list" & (filetype$ProMex > 0 | filetype$BioPharma > 0)) {
          # For ProMex/BioPharma: use peak start/stop boundaries for X range
          ranges$x <- c(min(filedata()[filedata()[,5] >= min(newdata$x),4]), 
                        max(filedata()[filedata()[,4] <= max(newdata$x),5]))
          ranges$y <- range(newdata$y)  
          
          # Case 2: Multiple file types (mixed ProMex/BioPharma with others)
        } else if (class(filedata()) == "list") {
          tab <- filedata()
          
          # Extract ProMex/BioPharma data
          if (sum(ftype()=="ProMex" | ftype()=="BioPharma") > 1) {
            tab2 <- RBindList(tab[ftype()=="ProMex" | ftype()=="BioPharma"])
          } else {
            tab2 <- tab[ftype()=="ProMex" | ftype()=="BioPharma"][[1]]
          }
          
          # Extract other format data
          if (sum(!(ftype()=="ProMex" | ftype()=="BioPharma")) > 1) {
            tab <- RBindList(tab[!(ftype()=="ProMex" | ftype()=="BioPharma")])
          } else {
            tab <- tab[!(ftype()=="ProMex" | ftype()=="BioPharma")][[1]]
          }
          
          # Calculate combined X range from both data types
          minx <- min(tab2[tab2[,5] >= min(newdata$x),4])
          minx <- min(minx, min(tab$RT))
          maxx <- max(tab2[tab2[,4] <= max(newdata$x),5])
          maxx <- max(maxx, max(tab$RT))
          ranges$x <- c(minx, maxx)
          ranges$y <- range(newdata$y)  
          
          # Case 3: Standard single table format
        } else {
          ranges$x <- range(newdata$x)
          ranges$y <- range(newdata$y)      
        }
      } else {
        newdata <- NULL
      }
    }
  })
  
  # ====================
  # DYNAMIC RANGE SELECTION
  # ====================
  
  # Function to determine appropriate plot ranges based on current zoom state and data
  defineranges <- function() {
    # Use current zoom ranges if they exist
    if (!is.null(ranges$x) & !is.null(ranges$y)) {
      rangesx <- ranges$x
      rangesy <- ranges$y
      
      # Scenario 1: Only MS data loaded (no MS/MS identification data)
    } else if (is.null(filedataMS2())) {
      if (class(filedata()) != "list") { # Single data table
        rangesx <- range(filedata()[,1])  # RT range
        rangesy <- range(filedata()[,2])  # Mass range
      } else { # Multiple data tables (mixed file types)
        # Extract retention time values from all possible sources
        x1 <- sapply(filedata(), function(z) {
          as.numeric(z$RT)
        })
        x2 <- sapply(filedata(), function(z) {
          as.numeric(z$PeakStart)
        })
        x3 <- sapply(filedata(), function(z) {
          as.numeric(z$PeakStop)
        })
        # Combine all time values and remove NAs
        x <- c(unlist(x1), unlist(x2), unlist(x3))
        x <- x[!is.na(x)]
        # Extract mass values
        y <- sapply(filedata(), function(z) {
          as.numeric(z$Mass)
        })
        rangesx <- range(x)  # Combined RT range
        rangesy <- range(y)  # Combined mass range
      }
      
      # Scenario 2: Only MS/MS data loaded (no MS trace) OR MS trace disabled
    } else if (is.null(filedata()) | input$MSTrace == FALSE) {
      rangesx <- range(filedataMS2()$MS2file$RT.in.min)
      rangesy <- range(filedataMS2()$MS2file$Precursor.MHplus.in.Da)
      
      # Scenario 3: Both MS and MS/MS data loaded (overlay mode)
    } else {
      if (filetype$ProMex > 0 | filetype$BioPharma > 0) {
        # Combine RT and mass ranges from both MS and MS/MS data
        rangesx <- range(c(filedataMS2()$MS2file$RT.in.min, filedata()[,5]))
        rangesy <- range(c(filedataMS2()$MS2file$Precursor.MHplus.in.Da, filedata()[,2]))
      }  else if (input$PDPFModeCheck == "PD") { # Proteome Discoverer overlay
        ranges$x <- range(c(filedataMS2()$MS2file$RT.in.min, filedata()$RT))
        ranges$y <- range(c(filedataMS2()$MS2file$Precursor.MHplus.in.Da, filedata()$Mass))
      } else if (input$PDPFModeCheck == "TP") { # TopPIC overlay
        ranges$x <- range(c(filedataMS2TP()$RT, filedata()[,1]))[2]
        ranges$y <- range(c(filedataMS2TP()$Mass, filedata()[,2]))
      }
    }
    return(list(rangesx, rangesy))
  }
  
  # ====================
  # SIGNAL SUMMARIZATION
  # ====================
  ## Calculate sum of signal in the plotting window
  summarized_signal <- reactiveVal(NULL)
  
  observeEvent(input$CalcSum, {
    
    # Get current plot ranges
    rangesx <- defineranges()[[1]]
    rangesy <- defineranges()[[2]]
    
    # Make table
    int_tab <- filedata()
    indexes <- int_tab$Mass <= rangesy[2] & int_tab$Mass >= rangesy[1] &
      int_tab$RT <= rangesx[2] & int_tab$RT >= rangesx[1]
    input_summarization <- data.frame(
      Sum_intensities = sum(int_tab$intensity[indexes], na.rm = T),
      Prop_sum_intensities = sum(int_tab$intensity[indexes], na.rm = T) / sum(int_tab$intensity, na.rm = T),
      Point_count = length(int_tab$intensity[indexes][!is.na(int_tab$intensity[indexes])]),
      Mass_min = rangesy[1],
      Mass_max = rangesy[2],
      RT_min = rangesx[1],
      RT_max = rangesx[2],
      Percentage_signal_threshold = input$IntensityThresh
    )
    summarized_signal(t(input_summarization))
    
  })
  
  output$SumSignal <- renderTable({ # return table with sum of signal
    if (!is.null(summarized_signal())) {
      summarized_signal()
    } else {
      return(invisible())
    }
  }, rownames = T, colnames = F)
  
  observeEvent(input$CalcSumReset, { # empty table when press reset
    summarized_signal(NULL)
  })
  
  ## To mask button when bars are plotted instead of points:
  output$show_calcsum <- reactive({
    if (linput() == 1) {
      if (!is.null(ftype())) {
        !(ftype() %in% c("BioPharma", "ProMex"))
      } else {
        TRUE
      }
    } else {
      FALSE
    }
  })
  outputOptions(output, "show_calcsum", suspendWhenHidden = FALSE)
  
  # ====================
  # PLOT GENERATION
  # ====================
  
  # Main plotting function - generates the 2D mass spectrometry visualization
  plotInput1 <- function(){
    # Validate user inputs for point size and intensity threshold
    validate(
      need(input$pch <= 10 & input$pch >= 0.1, "Please define a size of points between 0.1 and 10.")
    )
    validate(
      need(input$IntensityThresh <= 100 & input$IntensityThresh >= 0.1, "Please define a threshold between 0.1 and 100. This value corresponds to the proportion of the points in the data set that you want to upload (highest intensities).")
    )
    
    # Return NULL if no data is loaded
    if (is.null(filedata()) & is.null(filedataMS2()) & is.null(filedataMS2PF()) & is.null(filedataMS2TP())) {
      return(NULL)
    } else {
      # Generate plot if data is available
      if (!is.null(linput())) {
        # Get current plot ranges
        rangesx <- defineranges()[[1]]
        rangesy <- defineranges()[[2]]
        
        # Plot generation for RoWinPro/Bruker/TopPIC formats only
        if (filetype$BioPharma == 0 & filetype$ProMex == 0) {
          gtab <- filedata()
          
          # Multi-file comparison mode
          if (linput() >= 2) {
            g <- ggplot() + 
              geom_point(data = gtab, aes(x = RT, y = Mass, col = File, text = paste(RT, "min\n", Mass, "Da\nSignal:", intensity, "\n", File)), alpha = 0.7, size = input$pch) +
              geom_text(parse = TRUE) +
              coord_cartesian(xlim = rangesx, ylim = rangesy, expand = TRUE) +
              theme_bw() + 
              scale_colour_brewer(palette = colval()) + 
              ylab("Molecular Weight (Da)") + 
              xlab("Retention time (min)")
            
            # Single file mode with intensity coloring
          } else {
            g <- ggplot() + 
              geom_point(data = gtab, aes(x = RT, y = Mass, col = log10(intensity), text = paste(RT, "min\n", Mass, "Da\nSignal:", intensity)), alpha = 0.7, size = input$pch) +
              coord_cartesian(xlim = rangesx, ylim = rangesy, expand = TRUE) +
              theme_bw() + 
              scale_colour_distiller(palette = colval()) + 
              ylab("Molecular Weight (Da)") + 
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
              ylab("Molecular Weight (Da)") + 
              xlab("Retention time (min)")
          } else {
            g <- ggplot() + 
              geom_segment(data = gtab, aes(x = PeakStart, y = Mass, xend = PeakStop, yend = Mass, col = log10(intensity), text = paste("Start:", PeakStart, "min\n", "Stop:", PeakStop, "min\n", Mass, "Da\nSignal:", intensity)), alpha = 0.7, size = input$pch) +
              geom_point(data = gtab, aes(x = RT, y = Mass, col = log10(intensity), text = paste("Start:", PeakStart, "min\n", "Stop:", PeakStop, "min\n", Mass, "Da\nSignal:", intensity)), alpha = 0) +
              coord_cartesian(xlim = rangesx, ylim = rangesy, expand = TRUE) + 
              theme_bw() + 
              scale_colour_distiller(palette = colval()) + 
              ylab("Molecular Weight (Da)") + 
              xlab("Retention time (min)")
          }
        } else { # several types of input format
          gtabRWP <- RBindList(filedata()[ftype()=="RoWinPro" | ftype()=="Bruker" | ftype()=="TopPic"])
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
            geom_point(data = gtabRWP, aes(y = Mass, x = RT, col = File, text = paste(RT, "min\n", Mass, "Da\nSignal:", intensity, "\n", File)), size = input$pch) + 
            ylab("Molecular Weight (Da)") + 
            xlab("Retention time (min)")
        }
      }
      if (is.null(filedataMS2()) & is.null(filedataMS2PF()) & is.null(filedataMS2TP())) {
        return(g)
      } else { # MS2 files loaded
        # When in MS2 mode: overlay of the MS2 values:
        if (input$MSModeCheck == "MS2" & (!is.null(filedataMS2()) | !is.null(filedataMS2PF()) | !is.null(filedataMS2TP()))) {
          if (input$PDPFModeCheck == "PD" | testfileinput() == 3) {
            PSM <- filedataMS2()$PSM
            MS2 <- filedataMS2()$MS2
            PSM$ID <- paste0(PSM$Spectrum.File, "|", PSM$RT, "|", PSM$Precursor.MHplus.in.Da)
            MS2$ID <- paste0(MS2$Spectrum.File, "|", MS2$RT, "|", MS2$Precursor.MHplus.in.Da)
            print(MS2[MS2$ID %in% PSM$ID,])
            
            # Retrieve protein IDs in the MS2 table:
            MS2$Master.Protein.Descriptions <- PSM$Master.Protein.Descriptions[match(MS2$ID, PSM$ID)]
            print(head(MS2))
            print(MS2$Master.Protein.Descriptions[!is.na(MS2$Master.Protein.Descriptions)])
            
            # Make table for plot:
            gtabMS2 <- MS2[,c("RT.in.min", "Precursor.MHplus.in.Da", "Master.Protein.Descriptions")]
            gtabMS2$Identification <- ifelse(!is.na(gtabMS2$Master.Protein.Descriptions), "IDed", "Not IDed")
            # print(head(gtabMS2))
            # print(gtabMS2$Master.Protein.Descriptions)
            
            # Action button:
            if (input$HideMSMS == TRUE) {
              gtabMS2 <- gtabMS2[gtabMS2$Identification == "IDed",]
            }
            
            # names(gtabMS2)[3] <- "intensity"
            gtabMS2 <- gtabMS2[order(gtabMS2$Identification, decreasing = T),]
            
            if (!is.null(input$SelectProt) & input$PDPFModeCheck == "PD") {
              vec <- unique(gtabMS2$Master.Protein.Descriptions[gtabMS2$Master.Protein.Descriptions %in% input$SelectProt])
              vec <- vec[!is.na(vec)]
              getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
            }
            names(gtabMS2)[names(gtabMS2)=="Master.Protein.Descriptions"] <- "Protein.Descriptions"
            
            if (is.null(filedata0()) | input$MSTrace == FALSE) { # No MS trace
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
  
  
  # Plotly output if DataPoints == T
  output$plot1 <- renderPlotly({
    validate(
      need(!is.null(plotInput1()), '')
    )
    if (input$DataPoints == TRUE) {
      g <- plotInput1() 
      g <- g + 
        theme(legend.title = element_blank()) 
      p <- ggplotly(g, tooltip = "text", height = 800) %>%
        layout(dragmode = "select", 
               xaxis=list(fixedrange=TRUE),
               yaxis=list(fixedrange=TRUE),
               title = "") %>%
        config(displayModeBar = F) %>%
        layout(margin = list(l = 110, b = 40, r = 3, t = 10, pad = -2))
      # Remove IDed and Not IDed from the legend:
      if (input$MSModeCheck == "MS2") {
        #   vec <- sapply(p$x$data, function(x) {
        #     grepl("ID", x$name)
        #   })
        # p <- style(p, showlegend = FALSE, traces = which(vec))
        # Remove all legends:
        p <- style(p, showlegend = FALSE, traces = seq_len(length(p$x$data)))
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
      if (!is.null(linput())) {
        if ((linput() > 1 & input$MSModeCheck == "MS")) {
          plotInput1() + 
            theme(legend.direction ="vertical", legend.position="bottom") +
            guides(fill=guide_legend(ncol=2))
        } else if (nProtSelection() > 0 & input$MSModeCheck == "MS2") {
          plotInput1() + 
            theme(legend.direction ="vertical", legend.position="bottom") +
            guides(fill=guide_legend(ncol=2))
        } else {
          plotInput1()+ 
            theme(legend.direction ="vertical", legend.position="right") 
        }
      } else if (nProtSelection() > 0 & input$MSModeCheck == "MS2") {
        plotInput1() + 
          theme(legend.direction ="vertical", legend.position="bottom") +
          guides(fill=guide_legend(ncol=2))
      } else {
        plotInput1()+ 
          theme(legend.direction ="vertical", legend.position="right") 
      }
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
        summarized_signal(NULL)
      }
    }
  })
  
  # For export:
  # ====================
  # PLOT DOWNLOAD HANDLERS
  # ====================
  
  # PDF download handler - exports the current plot as a PDF file
  output$Download <- downloadHandler(
    # Dynamic filename generation based on loaded files and current date
    filename = function(){
      if (is.null(InputFilesMS2())) {
        paste0("VisioProt-MS_", substring(InputFileMS()$name, first = 1, last = (nchar(InputFileMS()$name)-4)), "_", Sys.Date(), ".pdf")
      } else {
        paste0("VisioProt-MS_", substring(InputFilesMS2()[[1]]$name, first = 1, last = (nchar(InputFilesMS2()[[1]]$name)-4)), "_", Sys.Date(), ".pdf")
      }
    },
    content = function(file) {
      # Dynamic height calculation based on legend requirements
      h <- 8  # Base height
      if (nProtSelection() > 0) {  # Protein selection mode
        if (!is.null(filedata0()) & input$MSTrace) {
          h <- h + 1.6  # Extra height for MS trace overlay
        }
        if (nProtSelection() > 5) {
          h <- h+(nProtSelection()*0.152)  # Additional height per protein
        }
      } else if (linput() > 1) {  # Multiple file comparison mode
        h <- h+(linput()*0.152)  # Additional height per file
      }
      
      # Adjust legend position based on plot type
      if ((linput() > 1 & input$MSModeCheck == "MS")|(nProtSelection() > 0 & input$MSModeCheck == "MS2")) {
        g <- plotInput1() + 
          theme(legend.direction ="vertical", legend.position="bottom")
      } else {
        g <- plotInput1()+ 
          theme(legend.direction ="vertical", legend.position="right") 
      }
      ggsave(file, plot = g, device = "pdf", width = 10, height = h)
    })
  
  # PNG download handler - exports the current plot as a PNG file
  output$Download1 <- downloadHandler(
    # Dynamic filename generation for PNG export
    filename = function(){
      if (is.null(InputFilesMS2())) {
        paste0("VisioProt-MS_", substring(InputFileMS()$name, first = 1, last = (nchar(InputFileMS()$name)-4)), "_", Sys.Date(), ".png")
      } else {
        paste0("VisioProt-MS_", substring(InputFilesMS2()[[1]]$name, first = 1, last = (nchar(InputFilesMS2()[[1]]$name)-4)), "_", Sys.Date(), ".png")
      }
    },
    content = function(file) {
      # Dynamic height calculation for PNG (in pixels)
      h <- 800  # Base height in pixels
      if (nProtSelection() > 0) {  # Protein selection mode
        if (!is.null(filedata0()) & input$MSTrace) {
          h <- h + 160  # Extra height for MS trace overlay
        }
        if (nProtSelection() > 5) {
          h <- h+(nProtSelection()*15.2)  # Additional height per protein
        }
      } else if (linput() > 1) {  # Multiple file comparison mode
        h <- h+(linput()*15.2)  # Additional height per file
      }
      
      # PNG device configuration with resolution settings
      device <- function(..., width, height) {
        grDevices::png(..., width = 1000, height = h, res = 120)
      }
      
      # Generate plot with appropriate legend positioning
      if ((linput() > 1 & input$MSModeCheck == "MS")|(nProtSelection() > 0 & input$MSModeCheck == "MS2")) {
        g <- plotInput1() + 
          theme(legend.direction ="vertical", legend.position="bottom")
      } else {
        g <- plotInput1()+ 
          theme(legend.direction ="vertical", legend.position="right") 
      }
      ggsave(file, plot = g, device = device)
    })
  
  # SVG download handler - exports the current plot as a vector SVG file
  output$Download2 <- downloadHandler(
    filename = function(){
      if (is.null(InputFilesMS2())) {
        paste0("VisioProt-MS_", substring(InputFileMS()$name, first = 1, last = (nchar(InputFileMS()$name)-4)), "_", Sys.Date(), ".svg")
      } else {
        paste0("VisioProt-MS_", substring(InputFilesMS2()[[1]]$name, first = 1, last = (nchar(InputFilesMS2()[[1]]$name)-4)), "_", Sys.Date(), ".svg")
      }
    },
    content = function(file) {
      h <- 8
      if (nProtSelection() > 0) {
        if (!is.null(filedata0()) & input$MSTrace) {
          h <- h + 1.6
        }
        if (nProtSelection() > 5) {
          h <- h+(nProtSelection()*0.152)
        }
      } else if (is.null(linput())) {
        if (linput() > 1) {
          h <- h+(linput()*0.152)
        }
      }
      device <- function(..., width, height) {
        grDevices::svg(..., width = 10, height = h)
      }
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
      ggsave(file, plot = g, device = device)
    })
  
  # tsv download handler - exports the current MS trace as a table
  output$Download3 <- downloadHandler( 
    filename = function(){
      validate(need(is.null(InputFilesMS2()),
                    "Only the MS trace will be saved in the exported table."
      ))
      paste0("VisioProt-MS_", substring(InputFileMS()$name, first = 1, last = (nchar(InputFileMS()$name)-4)), "_", Sys.Date(), ".tsv")
    },
    content = function(file) {
      gtab <- filedata()
      write.table(gtab, file = file, sep = "\t", row.names = F)
    }
  )
  
} 




# End of server function