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


# App:
############################################################################
## UI:
############################################################################

ui <- fluidPage(
  titlePanel("VisioProt-MS"),
  # Control side bar:
  sidebarLayout( 
    sidebarPanel(
      # File selection:
      fileInput("file", "Select input file(s):",  
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv"),
                multiple = T
      ),
      actionButton("TestFile1", "Single test file"),
      actionButton("TestFile2", "Multiple test files"),
      checkboxInput("DataPoints", "Show data labels (slower)", FALSE), # To switch between ggplot and plotly.
      # Selection of the colour scales. This depends on the number of input files:
      uiOutput("colourUI"),
      # Parameters for the plot:
      numericInput("pch", label = "Point size:", value = 1, min = 0.01, step = 0.1),
      numericInput("IntensityThresh", label = "Plotting threshold\n(Percentage points with the highest intensity):", value = 0.2, min = 0, max = 1, step = 0.1),
      # Information regarding how to zoom (depends on the plotting type):
      htmlOutput("ZoomParam"),
      
      # For debugging:
      textOutput("info"),
      
      br(),
      actionButton("DeZoom", "Unzoom", style='padding:8px; font-size:150%'),
      br(),
      br(),
      # Buttons for download:
      downloadButton("Download", "Download .pdf"),
      downloadButton("Download1", "Download .png")
    ),
    # Main panel for plotting (output different in function of the checkbox DataPoints):
    mainPanel(
      uiOutput("plotUI")
    )
  ),
  # Footer
  tabsetPanel(
    tabPanel(
      HTML('<footer><font size="0.8">copyright 2017 - CNRS - All rights reserved - VisioProt-MS V1.1</font></footer>')
    )
  )
)

## Server:
############################################################################

server <- function(input, output, clientData, session) {
  
  # Prevents error message when slow to render:
  
  
  # Text output to describe how to zoom (dependent on the checkbox DataPoints):
  output$ZoomParam <- renderUI({
    if (input$DataPoints) {
      HTML("<h4>Zoom in: select the ranges of interest.<br/>Zoom out: Click on the unzoom button.</h4>")
    } else {
      HTML("<h4>Zoom in: select the ranges of interest and double click.<br/>Zoom out: double click.</h4>")
    }
  })
  
  # Number of input file(s) from the same type:
  linput <- reactiveVal(1)
  
  # Define the number of input files to update UI in function:
  output$colourUI <- renderUI({
    if (linput() == 1) {
      selectInput("colourscale", "Colour scale:", # for continuous scales
                  c("Spectral" = "Spectral",
                    "Red/yellow/blue" = "RdYlBu", 
                    "Red/yellow/green" = "RdYlGn",
                    "yellow to red" = "YlOrRd"
                  ))
    } else {
      selectInput("colourscale", "Colour scale:", # for discete scales 
                  c("Set1" = "Set1",
                    "Dark2" = "Dark2", 
                    "Paired" = "Paired",
                    "Accent" = "Accent",
                    "Set2" = "Set2",
                    "Set3" = "Set3"
                  ))
    }
  })
  
  # Test files input:
  testfileinput <- reactiveVal(0) # 0: no test file; 1: single file; 2: multiple file
  
  filetype <- reactiveValues(RoWinPro = 0, BioPharma = 0) # Number of files of each type. Bruker files fall into the "RoWinPro" category once recognised and opened properly.
  
  observeEvent(input$TestFile1, {
    linput(1)
    testfileinput(1)
    filetype$RoWinPro <- 1
    filetype$BioPharma <- 0
  })
  
  observeEvent(input$TestFile2, {
    linput(4)
    testfileinput(2)
    filetype$RoWinPro <- 4
    filetype$BioPharma <- 0
  })
  
  observeEvent(input$file, { # Return to 0 when uploading new file
    testfileinput(0)
    if (!is.null(input$file)) {
      l <- list()
      l2 <- list()
      for(i in 1:nrow(input$file)){
        l[[i]] <- substr(readLines(input$file[i, 'datapath'])[1], 0, 17) == "Monoisotopic Mass"  | (substr(readLines(input$file[i, 'datapath'])[1], 0, 12) == "Protein Name") # TRUE if Biopharma
        l2[[i]] <- substr(readLines(input$file[i, 'datapath'])[2], 0, 13) == "Compound Name" # TRUE if Bruker
      }
      l <- unlist(l)
      l2 <- unlist(l2)
      filetype$RoWinPro <- length(l[l==F | l2==T]) # Bruker files too
      filetype$BioPharma <- length(l[l==T & l2==F])
      linput(max(as.numeric(table(l))))
    } else {
      filetype$RoWinPro <- 0
      filetype$BioPharma <- 0
    }
  })
  ftype <- reactive({
    if (is.null(input$file) & testfileinput() == 0) {
      return(NULL)
    } else {
      l <- list()
      for(i in 1:nrow(input$file)){
        val <- (substr(readLines(input$file[i, 'datapath'])[1], 0, 17) == "Monoisotopic Mass") | (substr(readLines(input$file[i, 'datapath'])[1], 0, 12) == "Protein Name") # T for BioPharma, F for RoWinPro
        val2 <- substr(readLines(input$file[i, 'datapath'])[2], 0, 13) == "Compound Name" # TRUE if Bruker
        val <- ifelse(val,  "BioPharma", "RoWinPro")
        val[val2] <- "Bruker"
        # Check it is the correct input format:
        if (val == "BioPharma") { # Find "Apex RT" in BioPharma files
          if (!grepl("Apex RT", readLines(input$file[i, 'datapath'])[1])) {
            val <- "DoNotApply"
            validate(
              need(val!="DoNotApply", "Incorrect input format")
            )
          }
        }
        if (val == "RoWinPro") { # Four columns in RoWinPro files
          if (length(as.numeric(gregexpr("\t", readLines(input$file[i, 'datapath'])[1])[[1]]))!=3) {
            val <- "DoNotApply"
            validate(
              need(val!="DoNotApply", "Incorrect input format")
            )
          }
        }
        if (val == "Bruker") { # First line has no "," and second line contains the column header " RT / min"
          if (!(!grepl(",", readLines(input$file[i, 'datapath'])[1]) & grepl(" RT", readLines(input$file[i, 'datapath'])[2]))) {
            val <- "DoNotApply"
            validate(
              need(val!="DoNotApply", "Incorrect input format")
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
    
    # Warning if trying to plot several types of data AND several files:
    validate(
      need(!((filetype$RoWinPro >= 2 & filetype$BioPharma >= 1) | (filetype$RoWinPro >= 1 & filetype$BioPharma >= 2)), "Can only input one file per type of format for comparison")
    )
    
    if (testfileinput() == 0) { # no input test file
      if (is.null(input$file)) {
        # User has not uploaded a file yet
        return(NULL)
      } else {
        lfiles <- list()
        for(i in 1:nrow(input$file)){
          if (ftype()[i] == "BioPharma") { # If the file is from Thermo BioPharma
            if (substr(readLines(input$file[i, 'datapath'])[1], 0, 17) == "Monoisotopic Mass") { # No IDs
              lfiles[[i]] <- read.table(input$file[i, 'datapath'], sep = "\t", header = T)
              lfiles[[i]] <- lfiles[[i]][,c("Apex.RT", "Monoisotopic.Mass", "Sum.Intensity", "Start.Time..min.", "Stop.Time..min.")] # Map the columns as in RoWinPro format, but with apex RT, start and stop instead of all the points of the peak.
            }
            if (substr(readLines(input$file[i, 'datapath'])[1], 0, 12) == "Protein Name") { # IDs
              lfiles[[i]] <- read.table(input$file[i, 'datapath'], sep = "\t", header = T)
              lfiles[[i]] <- lfiles[[i]][,c("Apex.RT", "Monoisotopic.Mass", "Sum.Intensity", "Start.Time..min.", "Stop.Time..min.")] # Map the columns as in RoWinPro format, but with apex RT, start and stop instead of all the points of the peak.
            }
          } else if (ftype()[i] == "RoWinPro")  { # RoWinPro output
            lfiles[[i]] <- read.table(input$file[i, 'datapath'], sep = "\t", header = F)
            lfiles[[i]] <- cbind(lfiles[[i]][,1:3], " Temp1" = rep(NA, nrow(lfiles[[i]])), "Temp2" = rep(NA, nrow(lfiles[[i]]))) # add one more column to allow row binding later on 
          } else if (ftype()[i] == "Bruker")  { # Bruker output
            lfiles[[i]] <- read.table(input$file[i, 'datapath'], sep = ",", header = F, skip = 2)
            lfiles[[i]] <- cbind(lfiles[[i]][,2:4], " Temp1" = rep(NA, nrow(lfiles[[i]])), "Temp2" = rep(NA, nrow(lfiles[[i]]))) # add one more column to allow row binding later on 
          }
        }
        names(lfiles) <- input$file$name
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
    } else { # Multiple file tests
      infile <- list.files("files/Multiple/", pattern = ".csv", full.names = T)
      lfiles <- list()
      for(i in 1:length(infile)){
        lfiles[[i]] <- read.table(infile[i], sep = "\t", header = F)
        lfiles[[i]] <- cbind(lfiles[[i]][,1:3], " Temp1" = rep(NA, nrow(lfiles[[i]])), "Temp2" = rep(NA, nrow(lfiles[[i]]))) # add one more column to allow row binding later on 
      }
      names(lfiles) <- c("test data 1", "test data 2", "test data 3", "test data 4")
      return(lfiles)   
      testfileinput(0)
    }
  })
  
  # Create table for plotting:
  filedata <- function() {
    if (is.null(filedata0())) {
      return(NULL)
    } else {
      lfiles <- filedata0()
      lfiles <- ThresholdCleaning(lfiles, input$IntensityThresh)
      if (filetype$BioPharma == 0) { # Only RoWinPro files
        l <- list()
        for (i in seq_along(lfiles)) {
          l[[i]] <- RenameBioPharma(lfiles[[i]])
        }
        names(l) <- names(lfiles)
        lfiles <- l
        return(RBindList(lfiles))
      } else if (filetype$RoWinPro == 0) { # Only BioPharma
        lfiles <- lapply(lfiles, function(x) {
          RenameBioPharma(x)
        })
        return(RBindList(lfiles))
      } else { # More than one type of files
        lfiles <- lapply(lfiles, function(x) {
          RenameBioPharma(x)
        })
        return(lfiles)
      }
    }
  }
  
  # For zoomable plot:
  ranges <- reactiveValues(x = NULL, y = NULL)
  observeEvent(input$DeZoom, {
    ranges$x <- NULL
    ranges$y <- NULL
  })
  observeEvent(event_data("plotly_selected"), {
    if (input$DataPoints) {
      newdata <- event_data("plotly_selected")
      if (!is.null(newdata) & class(newdata)=="data.frame") {
        ranges$x <- range(newdata$x)
        ranges$y <- range(newdata$y)       
      } else {
        newdata <- NULL
      }
    }
  })
  
  
  # Plot:
  # For export:
  ## ggplot for the option DataPoints == T:
  defineranges <- function() {
    if (!is.null(ranges$x) & !is.null(ranges$y)) {
      rangesx <- ranges$x
      rangesy <- ranges$y
      return(list(rangesx, rangesy))
    } else {
      if (filetype$BioPharma == 0 | filetype$RoWinPro == 0) { # one table
        rangesx <- range(filedata()[,1])
        rangesy <- range(filedata()[,2])
        return(list(rangesx, rangesy))
      } else  { # two tables because two types of files
        x <- c(filedata()[[1]][,1], filedata()[[2]][,1])
        y <- c(filedata()[[1]][,2], filedata()[[2]][,2])
        rangesx <- range(x)
        rangesy <- range(y)
        return(list(rangesx, rangesy))
      }
    }
  }
  
  plotInput1 <- function(){
    if (input$DataPoints == F | is.null(filedata())) {
      return(NULL)
    } else {
      rangesx <- defineranges()[[1]]
      rangesy <- defineranges()[[2]]
      if (filetype$BioPharma == 0) { # Only RoWinPro
        gtab <- filedata()
        if (linput() >= 2) { # if comparing several plots
          g <- ggplot(gtab, aes(x = RT, y = Mass, col = File, text = paste0("Intensity: ", intensity))) + 
            geom_point(alpha = 0.7, size = input$pch) +
            coord_cartesian(xlim = rangesx, ylim = rangesy, expand = TRUE) +
            theme_bw() + 
            scale_colour_brewer(palette = input$colourscale) + 
            ylab("Protein mass (Da)") + 
            xlab("Retention time (min)")
        } else { # only one plot, no overlay
          g <- ggplot(gtab, aes(x = RT, y = Mass, col = log10(intensity), text = paste0("Intensity: ", intensity))) + 
            geom_point(alpha = 0.7, size = input$pch) +
            coord_cartesian(xlim = rangesx, ylim = rangesy, expand = TRUE) +
            theme_bw() + 
            scale_colour_distiller(palette = input$colourscale) + 
            ylab("Protein mass (Da)") + 
            xlab("Retention time (min)")
        }
      } else if (filetype$RoWinPro == 0) { # Only type BioPharma
        gtab <- filedata()
        # Define the ranges for margins in the plot:
        rangesyB <- c(min(gtab$PeakStart, na.rm = T) - 0.05*min(gtab$PeakStart, na.rm = T), max(gtab$PeakStop, na.rm = T) + 0.05*max(gtab$PeakStop, na.rm = T))
        
        if (linput() >= 2) { # if comparing several plots
          g <- ggplot(gtab, aes(y = RT, x = Mass, col = File, ymin = PeakStart, ymax = PeakStop, text = paste0("Intensity: ", intensity))) + 
            geom_pointrange(alpha = 0.7, size = input$pch) +
            coord_flip(xlim = rangesy, ylim = rangesyB) +
            theme_bw() + 
            scale_colour_brewer(palette = input$colourscale) + 
            xlab("Protein mass (Da)") + 
            ylab("Retention time (min)")
        } else {
          g <- ggplot(gtab, aes(y = RT, x = Mass, ymin = PeakStart, ymax = PeakStop, col = log10(intensity), text = paste0("Intensity: ", intensity))) + 
            geom_pointrange(alpha = 0.7, size = input$pch) +
            coord_flip(xlim = rangesy, ylim = rangesyB) +
            theme_bw() + 
            scale_colour_distiller(palette = input$colourscale) + 
            xlab("Protein mass (Da)") + 
            ylab("Retention time (min)")
        }
      } else { # several types of input format
        gtabRWP <- filedata()[ftype()=="RoWinPro" | ftype()=="Bruker"][[1]]
        gtabBP <- filedata()[ftype()=="BioPharma"][[1]]
        
        # Define the ranges for margins in the plot:
        rangesyB <- c(min(gtabBP$PeakStart, na.rm = T) - 0.05*min(gtabBP$PeakStart, na.rm = T), max(gtabBP$PeakStop, na.rm = T) + 0.05*max(gtabBP$PeakStop, na.rm = T))
        rangesyB <- c(min(rangesyB[1], rangesx[1]), max(rangesyB[2], rangesx[2]))
        
        g <- ggplot() + 
          geom_pointrange(data = gtabBP, aes(y = RT, x = Mass, col = log10(intensity), ymin = PeakStart, ymax = PeakStop), size = input$pch, alpha = 0.7) + 
          coord_flip(xlim = rangesy, ylim = rangesyB) +
          theme_bw() + 
          scale_colour_distiller(palette = input$colourscale) + 
          geom_point(data = gtabRWP, aes(y = RT, x = Mass, col = log10(intensity))) + 
          xlab("Protein mass (Da)") + 
          ylab("Retention time (min)")
      }
      return(g)
    }
  }
  
  ## plotly for the option DataPoints == F:
  plotInput2 <- function(){
    if (is.null(filedata()) | input$DataPoints == T) {
      return(NULL)
    } else {
      rangesx <- defineranges()[[1]]
      rangesy <- defineranges()[[2]]
      if (filetype$BioPharma == 0) { # Only one type of plot: RoWinPro
        gtab <- filedata()
        if (linput() >= 2) { # For plotting multiple plots.
          g <- ggplot(gtab, aes(x = RT, y = Mass, col = File)) + 
            geom_point(alpha = 0.7, size = input$pch) + 
            coord_cartesian(xlim = rangesx, ylim = rangesy, expand = TRUE) + 
            theme_bw() + 
            scale_colour_brewer(palette = input$colourscale) + 
            ylab("Protein mass (Da)") + 
            xlab("Retention time (min)")
        } else { # For simple plot.
          g <- ggplot(gtab, aes(x = RT, y = Mass, col = log10(intensity))) + 
            geom_point(alpha = 0.7, size = input$pch) + 
            coord_cartesian(xlim = rangesx, ylim = rangesy, expand = TRUE) + 
            theme_bw() + 
            scale_colour_distiller(palette = input$colourscale) + 
            ylab("Protein mass (Da)") + 
            xlab("Retention time (min)")
        }
      } else if (filetype$RoWinPro == 0) { # Only one type of plot: BioPharma
        gtab <- filedata()
        # Define the ranges for margins in the plot:
        rangesyB <- c(min(gtab$PeakStart, na.rm = T) - 0.05*min(gtab$PeakStart, na.rm = T), max(gtab$PeakStop, na.rm = T) + 0.05*max(gtab$PeakStop, na.rm = T))
        
        if (linput() >= 2) { # For plotting multiple plots.
          g <- ggplot(gtab, aes(y = RT, x = Mass, col = File)) + 
            geom_pointrange(aes(ymin = PeakStart, ymax = PeakStop), alpha = 0.7, size = input$pch) + 
            coord_flip(xlim = rangesy, ylim = rangesyB, expand = TRUE) +
            theme_bw() + 
            scale_colour_brewer(palette = input$colourscale) + 
            xlab("Protein mass (Da)") + 
            ylab("Retention time (min)")
        } else { # For simple plot.
          g <- ggplot(gtab, aes(y = RT, x = Mass, col = log10(intensity))) + 
            geom_pointrange(aes(ymin = PeakStart, ymax = PeakStop), alpha = 0.7, size = input$pch) + 
            coord_flip(xlim = rangesy, ylim = rangesyB, expand = TRUE) +
            theme_bw() + 
            scale_colour_distiller(palette = input$colourscale) + 
            xlab("Protein mass (Da)") + 
            ylab("Retention time (min)")
        }
      } else { # two types of plot
        gtabRWP <- filedata()[ftype()=="RoWinPro" | ftype()=="Bruker"][[1]]
        gtabBP <- filedata()[ftype()=="BioPharma"][[1]]
        
        # Define the ranges for margins in the plot:
        rangesyB <- c(min(gtabBP$PeakStart, na.rm = T) - 0.05*min(gtabBP$PeakStart, na.rm = T), max(gtabBP$PeakStop, na.rm = T) + 0.05*max(gtabBP$PeakStop, na.rm = T))
        rangesyB <- c(min(rangesyB[1], rangesx[1]), max(rangesyB[2], rangesx[2]))
        
        g <- ggplot() + 
          geom_pointrange(data = gtabBP, aes(y = RT, x = Mass, col = log10(intensity), ymin = PeakStart, ymax = PeakStop), size = input$pch, alpha = 0.7) + 
          coord_flip(xlim = rangesy, ylim = rangesyB, expand = TRUE) + 
          theme_bw() + 
          scale_colour_distiller(palette = input$colourscale) + 
          geom_point(data = gtabRWP, aes(y = RT, x = Mass, col = log10(intensity))) + 
          xlab("Protein mass (Da)") + 
          ylab("Retention time (min)")
      }
      return(g)
    }
  }
  
  
  # Plotly output if DataPoints == T
  output$plot1 <- renderPlotly({
    validate(
      need(!is.null(plotInput1()), '')
    )
    if (!is.null(plotInput1)) {p <- ggplotly(plotInput1()) %>%
      layout(height = 800, dragmode = "select") %>%
      config(displayModeBar = F) %>%
      layout(xaxis=list(fixedrange=TRUE)) %>%
      layout(yaxis=list(fixedrange=TRUE))
    } else {
      plotly_empty()
    }
  })
  
  # Plotly output if DataPoints == F
  output$plot2 <- renderPlot({
    validate(
      need(!is.null(plotInput2()), '')
    )
    plotInput2()  
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
  output$Download <- downloadHandler(
    filename = "VisioProt-MS_Output.pdf",
    content = function(file) {
      device <- function(..., width, height) {
        grDevices::pdf(..., width = 10, height = 8)
      }
      if (input$DataPoints) {
        ggsave(file, plot = plotInput1(), device = device)
      } else {
        ggsave(file, plot = plotInput2(), device = device)
      }
    })
  output$Download1 <- downloadHandler(
    filename = "VisioProt-MS_Output.png",
    content = function(file) {
      device <- function(..., width, height) {
        grDevices::png(..., width = 1000, height = 800, res = 120)
      }
      if (input$DataPoints) {
        ggsave(file, plot = plotInput1(), device = device)
      } else {
        ggsave(file, plot = plotInput2(), device = device)
      }
    })
  
  # For debugging:
  output$info <- renderText({
   print(unlist(defineranges()))
  })
}

############################################################################

shinyApp(ui = ui, server = server)

############################################################################

# rsconnect::deployApp("T:/RRelatedWork/VisioProt-MS")
