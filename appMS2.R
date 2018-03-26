library(shiny)
library(ggplot2)
#library(Cairo)   # For nicer ggplot2 output when deployed on Linux
# setwd("P:/RRelatedWork/Rowinpro")

ui <- fluidPage(
  titlePanel("Plot from RoWinPro files"),
  sidebarLayout(
    sidebarPanel(
      column(5,
             fileInput("MS2file", "Choose MS2 File:", 
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
             ),
             actionButton("Button0", "Hide/show MSMS withouth ID"),
             numericInput("pch", label = "Point size:", value = 3, min = 0),
             numericInput("IntensityThresh", label = "Minimum intensity threshold:", value = 0),
             h4("To zoom in: select the ranges and double click"),
             h4("To zoom out: double click."),
             h4(textOutput("info")),
             downloadButton("Download", "Download .pdf")
      ),
      column(7,
             uiOutput("CheckboxProt")
      )
    ),
    mainPanel(
      plotOutput("plot",
                 #hover = "plot_hover",
                 click = "plot_click",
                 dblclick = "plot_dblclick",
                 brush = brushOpts(
                   id = "plot_brush",
                   resetOnNew = TRUE)
      )
    )
  )
)

server <- function(input, output, clientData, session) {
  #This function is repsonsible for loading in the selected file
  filedataMS2 <- reactive({
    infile <- input$MS2file
    if (is.null(infile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    read.table(infile$datapath, sep = "\t", header = T)
  })
  filedataPSM <- reactive({
    infile <- input$PSMfile
    if (is.null(infile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    read.table(infile$datapath, sep = "\t", header = T)
  })
  # For zoomable plot:
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  # List of proteins to select:
  output$CheckboxProt <- renderUI({
    proteins <- sort(unique(filedataPSM()$Master.Protein.Descriptions))
    checkboxGroupInput("Proteins", "Choose Protein", proteins)
  })
  
  output$plot <- renderPlot({
    if (!is.null(filedataMS2()) & !is.null(filedataPSM())) {
      # Create mapping ID:
      PSM <- filedataPSM()
      MS2 <- filedataMS2()
      PSM$ID <- paste0(PSM$Spectrum.File, "|", PSM$First.Scan)
      MS2$ID <- paste0(MS2$Spectrum.File, "|", MS2$First.Scan)
      # Retrieve protein IDs in the MS2 table:
      MS2$Master.Protein.Descriptions <- PSM$Master.Protein.Descriptions[match(MS2$ID, PSM$ID)]
      
      # Plot:
      gtab <- MS2[,c("RT.in.min", "Precursor.MHplus.in.Da", "Precursor.Intensity", "Master.Protein.Descriptions")]
      gtab$Identification <- ifelse(!is.na(gtab$Master.Protein.Descriptions), "IDed", "NoID")
      
      # Action button:
      if (input$Button0 %% 2 == 1) {
        gtab <- gtab[gtab$Identification == "IDed",]
      }
      
      names(gtab)[3] <- "intensity"
      gtab <- gtab[gtab[,3]>=input$IntensityThresh,]
      gtab <- gtab[order(gtab$Identification, decreasing = T),]

      ggplot(data = gtab[!(gtab$Master.Protein.Descriptions %in% input$Proteins),], aes(x = RT.in.min, y = Precursor.MHplus.in.Da, shape = Identification)) + 
        geom_point(alpha = 0.8, size = input$pch, col = "grey30") + 
        coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)  + 
        theme_bw() + 
        scale_shape_manual(values = c(16, 1)) + 
        ylab("Protein mass (Da)") + 
        xlab("Retention time (min)") + 
        geom_point(data = gtab[gtab$Master.Protein.Descriptions %in% input$Proteins,], aes(x = RT.in.min, y = Precursor.MHplus.in.Da, col = Master.Protein.Descriptions), size = input$pch, alpha = 0.8) 
      }
  }, height = 800)
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$plot_dblclick, {
    brush <- input$plot_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
  output$info <- renderText({
    if(is.null(input$plot_click)) {
      "Click on the plot to get the mass and retention time."
    } else {
      paste0(
        "RT(min)= ", round(input$plot_click$x, 1), "; Mass(Da)=", round(input$plot_click$y, 1), "\n"
      )
    }
  })
  
  # For export:
  plotInput <- function(){
    if (!is.null(filedataMS2()) & !is.null(filedataPSM())) {
      # Create mapping ID:
      PSM <- filedataPSM()
      MS2 <- filedataMS2()
      PSM$ID <- paste0(PSM$Spectrum.File, "|", PSM$First.Scan)
      MS2$ID <- paste0(MS2$Spectrum.File, "|", MS2$First.Scan)
      # Retrieve protein IDs in the MS2 table:
      MS2$Master.Protein.Descriptions <- PSM$Master.Protein.Descriptions[match(MS2$ID, PSM$ID)]
      
      # Plot:
      gtab <- MS2[,c("RT.in.min", "Precursor.MHplus.in.Da", "Precursor.Intensity", "Master.Protein.Descriptions")]
      gtab$Identification <- ifelse(!is.na(gtab$Master.Protein.Descriptions), "IDed", "NoID")
      
      # Action button:
      if (input$Button0 %% 2 == 1) {
        gtab <- gtab[gtab$Identification == "IDed",]
      }
      
      names(gtab)[3] <- "intensity"
      gtab <- gtab[gtab[,3]>=input$IntensityThresh,]
      gtab <- gtab[order(gtab$Identification, decreasing = T),]
      
      ggplot(data = gtab[!(gtab$Master.Protein.Descriptions %in% input$Proteins),], aes(x = RT.in.min, y = Precursor.MHplus.in.Da, shape = Identification)) + 
        geom_point(alpha = 0.8, size = input$pch, col = "grey30") + 
        coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)  + 
        theme_bw() + 
        scale_shape_manual(values = c(16, 1)) + 
        ylab("Protein mass (Da)") + 
        xlab("Retention time (min)") + 
        geom_point(data = gtab[gtab$Master.Protein.Descriptions %in% input$Proteins,], aes(x = RT.in.min, y = Precursor.MHplus.in.Da, col = Master.Protein.Descriptions), size = input$pch, alpha = 0.8) + 
        theme(legend.text=element_text(size=7))
    }
  }
  
  output$Download <- downloadHandler(
    filename = "ShinyAppMLPOutput.pdf",
    content = function(file) {
      device <- function(..., width, height) {
        grDevices::pdf(..., width = 11.69, height = 8.27)
      }
      ggsave(file, plot = plotInput(), device = device)
    })
  
}

shinyApp(ui = ui, server = server)

# rsconnect::deployApp('P:/RRelatedWork/RoWinProTestPD')

# sessionInfo()
