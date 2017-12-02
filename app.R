############################################################################
# Copyright CNRS 2017
# Contributor : Marie Locard-Paulet (21/11/2017) [marie.locard@ipbs.fr]
# This software is a computer program whose purpose is to visualize and inspect deconvoluted MS 2D maps.
# This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software. You can  use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info". 
# As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive licensors have only limited liability. In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or reproducing the software by the user in light of its specific status of free software,that may mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data to be ensured and, more generally, to use and operate it in the  same conditions as regards security. 
# The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
############################################################################




############################################################################
# Packages:

library(shiny)
# devtools::install_github('hadley/ggplot2')
library(ggplot2)
library(plotly)
library(dplyr)

############################################################################
# Functions:

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

############################################################################
# App:

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
      #textOutput("info"),
      
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
      HTML('<footer><font size="0.8">copyright 2017 - CNRS - All rights reserved - VisioProt-MS V1.0</font></footer>')
    )
  )
)

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
  
  # Number of input file(s) :
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
  observeEvent(input$TestFile1, {
          linput(1)
          testfileinput(1)
  })
  observeEvent(input$TestFile2, {
          linput(2)
          testfileinput(2)
  })
  observeEvent(input$file, { # Return to 0 when uploading new file
          testfileinput(0)
  })
  
  # Input the data table:
  filedata0 <- reactive({
          #This function is repsonsible for loading in the selected file
          if (testfileinput() == 0) { # no input file
                  infile <- input$file
                  if (is.null(input$file)) {
                          # User has not uploaded a file yet
                          return(NULL)
                  } else {
                          lfiles <- list()
                          for(i in 1:nrow(input$file)){
                                  lfiles[[i]] <- read.table(input$file[i, 'datapath'], sep = "\t", header = F)
                          }
                          names(lfiles) <- input$file$name
                          linput(length(lfiles))
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
          } else {
                  infile <- list.files("files/Multiple/", pattern = ".csv", full.names = T)
                  lfiles <- list()
                  for(i in 1:length(infile)){
                          lfiles[[i]] <- read.table(infile[i], sep = "\t", header = F)
                  }
                  names(lfiles) <- c("test data 1", "test data 2", "test data 3", "test data 4")
                  return(lfiles)   
                  testfileinput(0)
          }
  })
  
  
  
  
  
  # Create table for plotting:
  filedata <- function() {
    if (!is.null(filedata0())) {
      lfiles <- filedata0()
      lfiles <- lapply(seq_along(lfiles), function(x) {
        temp <- cbind(lfiles[[x]], "File" = rep(names(lfiles)[x], nrow(lfiles[[x]])))
        temp <- temp[order(temp[,3], decreasing = T),]
        thresh <- floor(input$IntensityThresh * nrow(temp))
        temp <- temp[c(1:thresh),]
        temp
      })
      if (length(lfiles) == 1) {
        gtab <- lfiles[[1]]
      } else {
        gtab <- RBindList(lfiles)
      }
      gtab <- gtab[!is.na(gtab[,3]),]
      names(gtab)[3] <- "intensity"
      names(gtab)[2] <- "Mass"
      names(gtab)[1] <- "RT"
      gtab
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
  defineranges <- function(){
    if (!is.null(filedata())) {
      gtab <- filedata()
      if (!is.null(ranges$x) & !is.null(ranges$y)) {
        rangesx <- ranges$x
        rangesy <- ranges$y
      } else {
        rangesx <- range(gtab$RT)
        rangesy <- range(gtab$Mass)
      }
      list(rangesx, rangesy)
    }
  }
  
  plotInput1 <- function(){
    if(input$DataPoints) {
      gtab <- filedata()
      rangesx <- defineranges()[[1]]
      rangesy <- defineranges()[[2]]
      if (linput() == 2) { # if comparing several plots
        ggplot(gtab, aes(x = RT, y = Mass, col = File, text = paste0("Intensity: ", intensity))) + 
          geom_point(alpha = 0.7, size = input$pch) +
          coord_cartesian(xlim = rangesx, ylim = rangesy, expand = TRUE) +
          theme_bw() + 
          scale_colour_brewer(palette = input$colourscale) + 
          ylab("Protein mass (Da)") + 
          xlab("Retention time (min)")
      } else {
        ggplot(gtab, aes(x = RT, y = Mass, col = log10(intensity), text = paste0("Intensity: ", intensity))) + 
          geom_point(alpha = 0.7, size = input$pch) +
          coord_cartesian(xlim = rangesx, ylim = rangesy, expand = TRUE) +
          theme_bw() + 
          scale_colour_distiller(palette = input$colourscale) + 
          ylab("Protein mass (Da)") + 
          xlab("Retention time (min)")
      }
    }
  }
  
  ## plotly for the option DataPoints == F:
  plotInput2 <- function(){
    if (input$DataPoints == F) {
      rangesx <- defineranges()[[1]]
      rangesy <- defineranges()[[2]]
      if (!is.null(filedata())) {
        gtab <- filedata()
        if (linput() == 2) { # For plotting multiple plots.
          ggplot(gtab, aes(x = RT, y = Mass, col = File)) + 
            geom_point(alpha = 0.8, size = input$pch) + 
            coord_cartesian(xlim = rangesx, ylim = rangesy, expand = TRUE) + 
            theme_bw() + 
            scale_colour_brewer(palette = input$colourscale) + 
            ylab("Protein mass (Da)") + 
            xlab("Retention time (min)")
        } else { # For simple plot.
          ggplot(gtab, aes(x = RT, y = Mass, col = log10(intensity))) + 
            geom_point(alpha = 0.8, size = input$pch) + 
            coord_cartesian(xlim = rangesx, ylim = rangesy, expand = TRUE) + 
            theme_bw() + 
            scale_colour_distiller(palette = input$colourscale) + 
            ylab("Protein mass (Da)") + 
            xlab("Retention time (min)")
        }
      } else {
        plotly_empty()
      }
    }
  }
  
  # Plotly output if DataPoints == T
  output$plot1 <- renderPlotly({
    validate(
      need(!is.null(plotInput1()), 'Plotting...')
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
  #output$info <- renderText({
  #  unlist(defineranges())
  #})
}

############################################################################

shinyApp(ui = ui, server = server)

############################################################################

# rsconnect::deployApp('//GANDALF/mloca/RRelatedWork/VisioProt-MS_V0.0')
