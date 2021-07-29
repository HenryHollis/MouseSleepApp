library(shiny)
library(tidyverse)
library(readxl)
library(pracma)

shinyUI(fluidPage(
    titlePanel("Upload xlsx files"),
    sidebarLayout(
        sidebarPanel(
            fileInput("file","Upload the file", multiple = TRUE,  accept = c(".xlsx")), # fileinput() function is used to get the file upload contorl option
            #helpText("Default max. file size is 5MB"),
            #helpText("Select the read.table parameters below"),
            #checkboxInput(inputId = 'header', label = 'Header', value = TRUE),
            #checkboxInput(inputId = "stringAsFactors", "stringAsFactors", FALSE),
            #radioButtons(inputId = 'sep', label = 'Separator', choices = c(Comma=',',Semicolon=';',Tab='\t', Space=''), selected = ','),
            uiOutput("selectfile"),
            tags$hr(),
            sliderInput("num", "Cell to plot:",min = 2, max = 400, step=1,value=c(2)),
            tags$hr(),
            sliderInput("thresholdFactor", "Scale Threshold:",min = 0, max = 3, step=.05, value=c(1)),
            tags$hr(),
            sliderInput("lengthCutoff", "Minimum time units of spike:",min = 0, max = 100, step=1,value=c(10)),
            tags$hr(),
            downloadButton("download", "Download File(s)")
            
            
        ),
        mainPanel(
            uiOutput("tb")
            
        )
        
    )
))
