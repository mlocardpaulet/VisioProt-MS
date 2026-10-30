############################################################################
# VisioProt-MS: Interactive 2D maps from intact protein mass spectrometry
# Copyright CNRS 2017
# Contributor : Marie Locard-Paulet (21/11/2017) [marie.locard@ipbs.fr]
# 
# DESCRIPTION:
# This software is a computer program whose purpose is to visualize and inspect 
# deconvoluted MS 2D maps from intact protein mass spectrometry experiments.
# It supports various deconvolution software outputs including:
# - RoWinPro and Bruker DataAnalysis
# - BioPharma Finder (Thermo Fisher Scientific)
# - ProMex (Pacific Northwest National Laboratory)
# - TopPIC (University of California, San Diego)
#
# The application can operate in two modes:
# 1. MS mode: Visualization of deconvoluted MS1 data only
# 2. MS/MS mode: Overlay of top-down MS/MS identification results on MS1 data
#
# LICENSE:
# This software is governed by the CeCILL license under French law and abiding 
# by the rules of distribution of free software. You can use, modify and/or 
# redistribute the software under the terms of the CeCILL license as circulated 
# by CEA, CNRS and INRIA at the following URL "http://www.cecill.info". 
# As a counterpart to the access to the source code and rights to copy, modify 
# and redistribute granted by the license, users are provided only with a limited 
# warranty and the software's author, the holder of the economic rights, and the 
# successive licensors have only limited liability. In this respect, the user's 
# attention is drawn to the risks associated with loading, using, modifying and/or 
# developing or reproducing the software by the user in light of its specific 
# status of free software, that may mean that it is complicated to manipulate, 
# and that also therefore means that it is reserved for developers and experienced 
# professionals having in-depth computer knowledge. Users are therefore encouraged 
# to load and test the software's suitability as regards their requirements in 
# conditions enabling the security of their systems and/or data to be ensured and, 
# more generally, to use and operate it in the same conditions as regards security. 
# The fact that you are presently reading this means that you have had knowledge 
# of the CeCILL license and that you accept its terms.
############################################################################

############################################################################
# REQUIRED PACKAGES AND GLOBAL CONFIGURATION
############################################################################

# Core Shiny framework for building interactive web applications
library(shiny)

# Advanced plotting libraries
# devtools::install_github('hadley/ggplot2')  # Use development version if needed
library(ggplot2)    # Grammar of graphics for static plots
library(plotly)     # Interactive web-based data visualization
library(dplyr)      # Data manipulation and transformation

# Visualization and UI enhancement packages
library(RColorBrewer)  # Color palettes for data visualization
library(shinyBS)       # Bootstrap components for Shiny (tooltips, modals, etc.)
library(data.table)    # High-performance data manipulation and reading

# Import file
library(rhdf5)  # For Unidec input compatibility

# Global application settings
options(shiny.maxRequestSize=90*1024^2) # Set maximum upload size to 90MB for large MS files

############################################################################
# UI
############################################################################
source("ui.R")

############################################################################
# SERVER
############################################################################
source("server.R")

# ====================
# APPLICATION LAUNCH
# ====================
# Launch the Shiny application with the defined UI and server components
shinyApp(ui = ui, server = server)

# ====================
# DEPLOYMENT CONFIGURATION
# ====================
# Deployment command for shinyapps.io or RStudio Connect (commented out)
# rsconnect::deployApp("T:/RRelatedWork/VisioProt-MS")

