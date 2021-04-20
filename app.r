### Volcano Plot
library(ggplot2)
library(extrafont)
require(grid)
library(shiny)
library(shinyBS)
library(shinyjs)
library(colourpicker)
library(metricsgraphics)
library(RColorBrewer)
library(scales)
library(ggrepel)
#library(DT)
library(plotly)
options(shiny.usecairo=TRUE)

ui = fluidPage(
  titlePanel("Volcano Plot" ),
  sidebarLayout(
    sidebarPanel(      # uploading file
         fileInput("csvFile", "Choose List File", accept = c("text/csv",
                    "text/comma-separated-values,text/plain", ".csv")
        ),
    "Note: The input file should be text file (comma, tab separated), You can download all datasets for the different comparisons", a(href="https://github.com/Vilchezlab/shiny-Volcano-Plot","here"),".",
     
      br(),
     
     br(),
      # colname of your geneName, foldchange and pvalue
      fluidRow(
#        column(4,
#               textInput("geneName", "column name of your gene ID", value = "geneName") ),
#        column(4,
 #              textInput("logFC", "column name of your logFC", value = "logFC")),
#        column(4,
 #              textInput("pvalue", "column name of your pvalue", value = "pvalue"))
      ), br(),
      
      # select Protein 
      textInput("SP", label = "Gene Name from your List"),
      br(),
      
      # select fold change and pvalue 
      fluidRow(
        column(6,
               numericInput("logFC1", "log2 Fold Change threshold 1",
                            value = 1, min = 0, max = 2, step = 0.1) ),
        column(6,
               numericInput("logFC2", "log2 Fold Change threshold 2",
                            value = 2, min = 2, max =10, step = 0.2) )
         ),
      
      fluidRow(
        column(4,
               numericInput("pvalue1", "p value threshold 1",
                            value = 0.05, min = 0.01, max = 0.1, step = 0.01) ),
        column(4,
               numericInput("pvalue2", "p vlaue threshold 2",
                            value = 1e-2, min= 1e-5, max = 1e-3, step = 1e-4) ),
        )  ),
    
    # Show a plot 
    mainPanel(
      tabsetPanel(
        tabPanel("volcano output", plotOutput("volcanoImage") 
        ),
        tabPanel("data table", dataTableOutput("inputdata"))
      )  )
  ))


server <- function(input, output, session) {
  
  inputdf <- reactive({
    inFile <- input$csvFile
    
    if (is.null(inFile))
      return(NULL)
    
    df <- read.csv(inFile$datapath)
    pos <- which(colnames(df) %in% c(input$geneName, input$log2FC, input$pvalue))
    colnames(df)[pos] <- c("geneName","log2FC","pvalue")
    df
  })
  
  
  output$inputdata <- renderDataTable({
  # Inpute file
    inputdf()
      })
    output$volcanoImage <- renderPlot({
    if (! is.null(inputdf())){
      Vplot(inputdf(),
      
                  selectprotein = input$SP,
                  log2FC1 = input$logFC1,
                  log2FC2 = input$logFC2,
                  pval1 = input$pvalue1,
                  pval2 = input$pvalue2,
            ) }
    })
 
  Vplot <- function(df, 
                           selectprotein = "", 
                           circlesize=2.0,
                           log2FC1 = 1.5, log2FC2 = 2.0,
                           pval1 = 0.05, pval2 = 1e-2
  ){
    require(grid)
    mycol <- c("darkgreen","chocolate4","blueviolet","#223D6C","#D20A13",
               "#088247","#58CDD9","#7A142C","#5D90BA","#431A3D",
              "#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
   
    x <- subset(df, ! is.na(df$pvalue))
    
   # set color for different logFC
    n1 <- length(x[, 1])
    cols <- rep("grey", n1)
    names(cols)<- rownames(x)
    
    cols[x$pvalue < pval1 & x$log2FC > log2FC1]<- "red"
    cols[x$pvalue < pval2 & x$log2FC > log2FC2]<- "blue"
    cols[x$pvalue < pval1 & x$log2FC < -log2FC1]<- "blue"
    cols[x$pvalue < pval2 & x$log2FC < -log2FC2]<- "red"
    color_transparent <- adjustcolor(cols, alpha.f = 0.5)
    x$color_transparent <- color_transparent
    
    # size  for different logFC
    n1 <- length(x[, 1])
    size <- rep(1, n1)
    size[x$pvalue < pval1 & x$log2FC > log2FC1]<- 2
    size[x$pvalue < pval2 & x$log2FC > log2FC2]<- 4
    size[x$pvalue < pval1 & x$log2FC < -log2FC1]<- 2
    size[x$pvalue < pval2 & x$log2FC < -log2FC2]<- 4
    
    # ggplot
    p1 <- ggplot(data=x, aes(log2FC, pvalue)) +
      geom_point(alpha = 0.6, size = size, col = cols)+ #colour = x$color_transparent) +
      labs(x="Log2FoldChange", y="-LogP-value") + scale_x_continuous(
        breaks = c(-10, -5, -log2FC1, 0, log2FC1, 5, 10), 
        labels = c(-10, -5, -log2FC1, 0, log2FC1, 5, 10),  )  +
      geom_vline(xintercept = c(-log2FC1, log2FC1), color="grey", linetype="longdash", lwd = 0.5) + 
      geom_hline(yintercept = -log10(pval1), color="grey", linetype="longdash", lwd = 0.5) +
      theme_bw()+ theme(panel.grid=element_blank())
    
    # add threshold line
    p1 <- p1 + 
      geom_vline(xintercept = c(-log2FC2, log2FC2), color="green", linetype="longdash", lwd = 0.5) +
      geom_hline(yintercept = -log10(pval2), color="green", linetype="longdash", lwd = 0.5)
    
    #show Proteins Name
    if (nchar(selectprotein) < 1){
      return(p1) 
    }else{

      protein <- strsplit(selectprotein, ",")[[1]]
      selectprotein <- x[x$geneName == protein,]
      p2 <- p1 + 
        # point the search Protein
        geom_point(data = selectprotein, 
                            alpha = 2, size = circlesize, shape = 2, color = "blue") +
        scale_color_manual(values = cols) + 
        geom_text_repel(data = selectprotein, aes(label=geneName), show.legend = FALSE, 
                                 size = 5, box.padding = unit(0.35, "lines"), 
                                 point.padding = unit(0.3, "lines")) #+
         
      return(p2)
      }
    }
} 
shinyApp(ui, server)
#library(rsconnect)
#deployApp()