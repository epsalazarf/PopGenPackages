# AUTO PCA PLOTTER: CSV (Shiny)
# Shiny: app.R
# Author: Pavel Salazar-Fernandez
# Developed at: LANGEBIO - Mexico
# Last Edit: August 11 2016

# Requirements:
# - CSV file with the PCA values from fine STRUCTURE

# Pipeline:
# 1. Reads an csv file from a chosen directory.
# 2. Identifies names and populations from the file.
# 3. Generates a PCA plot with color-coded populations.

# Features:
# - Plot types: Can select between points or tags for the plot.
# - Select Population: Shows only a particular population.
# - Emphasize Population: Highlights points or tags for a chosen population.
# - Show/Hide Legend: Displays all values from the selected category.
# - Interactive Zoom: select an area and double click to zoom in, double click again to zoom out.
# - Area Info: brush an area to see the selected points info and coordinates.

#<START>
# Load required libraries 
library(shiny)

#<INPUT>
# Choose data file
data <- read.csv(choose.files(default = "", caption = "Select CSV file", multi=F))
#</INPUT>

#<PREPARATIONS>

# Load PCA files
evec <- as.numeric(data[1,-1])
evpct <- round(evec/sum(evec) * 100,2)
data <- data[-1,1:21]
rownames(data) <- NULL

# Read samples and populations
inds <- as.character(data[,1])
nsamples <- length(inds)
pops <- unique(substr(inds,0,3))
data$POP <- substr(inds,0,3)
pops <- unique(data$POP)

ncomps <- table(sapply(data, class))["numeric"]
PC1Col <- (match("numeric",sapply(data, class)))
PCnames<- colnames(data[,PC1Col:(PC1Col+ncomps-1)])

PCnames
#</PREPARATIONS>

#<UI>
ui <- fluidPage(    
  # Page Title
  img(src="MorLabLogo.jpg", style = "float:right"),
  titlePanel("Auto PCA Plotter: CSV Data"),
  helpText("Author: Pavel Salazar-Fernandez | Developed at: LANGEBIO (MX)"),
  hr(),    
  # Sidebar
  sidebarLayout(      
    
    # Input Panels
    sidebarPanel(width=3,
                 h4("Settings"),
                 textInput("plottitle", label="Title", value = ""),
                 numericInput("PCa", label = "First Component (X)", value = 1,
                              min = 1, max = ncomps, step = 1),
                 numericInput("PCb", label = "Second Component (Y)", value = 2,
                              min = 1, max = ncomps, step = 1),
                 selectizeInput("pops", label = "Populations displayed:", 
                                choices= pops, selected = NULL,
                                options= list(maxItems = length(pops)-1,
                                              placeholder = 'Select population(s)',
                                              onInitialize = I('function() { this.setValue(""); }'))
                 ),
                 selectizeInput("pope", label = "Populations emphasis:", 
                                choices= pops, selected = NULL,
                                options= list(maxItems = 1,
                                              placeholder = 'Choose population',
                                              onInitialize = I('function() { this.setValue(""); }'))
                 ),
                 hr(),
                 radioButtons("type", label = "Type:",
                              choices = list("Points" = 1, "Text" = 2), 
                              selected = 1),
                 checkboxInput("leg", label = "Legend", value = FALSE),
                 hr(),
                 h5("Points Info"),
                 verbatimTextOutput("brshinfo")
    ),
    
    # Plotting Area
    mainPanel(width=9,
              plotOutput("PCAPlot", width = "960px", height= "960px", dblclick = "dclk", brush = brushOpts(id= "brsh", resetOnNew = TRUE))  
    )
  )
)
#</UI>

#<SERVER>
server <- function(input, output) {
  
  #<REACTIVES>
  # Inputs
  sub.pops <- reactive(as.vector(data[data$POP %in% input$pops,"Label"]))
  pope.idn <- reactive(as.vector(data[data$POP == input$pope,"Label"]))
  #grouping <- reactive(if(input$flds != "")as.character(unique(popinfo[,input$flds])) else(""))
  PCaCol <- reactive(PC1Col + input$PCa - 1)
  PCbCol <- reactive(PC1Col + input$PCb - 1)
  pct.varPCa <- reactive(evpct[input$PCa])
  pct.varPCb <- reactive(evpct[input$PCb])
  ranges <- reactiveValues(x = NULL, y = NULL)
  #</REACTIVES>
  
  #<OBSERVERS>
  observeEvent(input$dclk, {
    brush <- input$brsh
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  #</OBSERVERS>
  
  #<PROCESSING>
  hiplot <- reactive({
    sbP.idx <- which(data$Label %in% sub.pops())
    pope.idx <- which(data$Label %in% pope.idn())
    groups.color <- substr(c(rainbow(length(pops)),"#000000"),0,7)
    names(groups.color) <- c(pops,"(Emphasis)")
    data$Colors <- groups.color[data$POP]

    if(length(pope.idx) > 0){groups.color[length(groups.color)] <- input$pope}
    plot(data[,c(PCaCol(),PCbCol())], col= data$Colors, pch= 20, cex.main = 1.5,
         xlab= paste("PC",input$PCa," (",pct.varPCa(),"%)",sep=""),
         ylab= paste("PC",input$PCb," (",pct.varPCb(),"%)",sep=""), 
         main= ifelse(input$plottitle=="","PCA Plot",input$plottitle),
         type= ifelse(is.null(input$pops) && input$type == 1,"p","n"), cex=3,
         xlim= ranges$x,
         ylim= ranges$y
    )
    abline(v= 0, h= 0, lty= 2, col= "grey")
    if(input$leg) legend("topright", names(groups.color), ncol= 1, col= groups.color, pch= 20,
                         pt.cex= 3, bty= "n")
    if (input$type==1) {
      points(data[sbP.idx,c(PCaCol(),PCbCol())], cex= 2, pch= 19, col= data$Colors[sbP.idx])
      points(data[pope.idx,c(PCaCol(),PCbCol())], cex= 2.5, pch= 23, bg="#000000", col= "#FFFFFF")
    }
    else {
      if (is.null(input$pops)) text(data[,c(PCaCol(),PCbCol())], labels= data$Label, col= data$Colors, cex= 0.9)
      else text(data[sbP.idx,c(PCaCol(),PCbCol())], labels=data[sbP.idx,"Label"], cex= 1, col= data$Colors[sbP.idx])
      if (length(pope.idx) > 0) text(data[pope.idx,c(PCaCol(),PCbCol())], labels= data[pope.idx,"Label"], font= 2, cex=1, col= "#000000")
    }
  })
  #</PROCESSING>
  
  #<OUTPUT>
  output$PCAPlot <- renderPlot(hiplot())
  output$brshinfo <- renderPrint({data[rownames(brushedPoints(data[,],input$brsh, xvar= PCaCol(), yvar= PCbCol())),c("Label","POP",PCnames[c(PCaCol()-1,PCbCol()-1)])]})
  #</OUTPUT>
}
#</SERVER>

#<APP>
shinyApp(ui = ui, server = server)
#</APP>

#<END>