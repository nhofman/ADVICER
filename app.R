#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
#

library(shiny)
#library(shinyBS)
#library(shinyWidgets)
library(DT)
library(VennDiagram)
library(UpSetR)
library(upsetjs)
library(ggplot2)
library(plotly)
library(pheatmap)
library(gtools)
library(openxlsx)
library(vcfR)
library(tidytext)
library(dplyr)
library(heatmaply)
library(shinycssloaders)

plotHeatmap <- function(x, row_subset = NA, distMethod = "euclidean", clusterMethod = "complete", clrn = 1, clcn = 1, setWidth = F,
                        rowClust = T, colClust = T, fontsize_r = 0.8, fontsize_c = 10, annCol = NA, annRow = NA, border_col = "grey60", plot.fig = T,
                        legend.limit = 1, filter_col = NA, annotation_colors = NA, break_step = 0.1, display_numbers = F, hover = NULL, file = NA){
  if(is.na(row_subset[1])){
    xx <- x
  }else{
    xx <- x[row_subset,, drop = F]
  }
  #xx[xx==0] <- 0.000001
  #xx <- na.omit(xx)
  if (nrow(xx)<2) {
    rowClust = F
  }
  if(!is.na(filter_col)){
    xx <- xx[,which(colMaxs(as.matrix(abs(xx)))>filter_col)]
  }
  xx.unlist <- as.numeric(unlist(xx))
  color <- vector()
  if(-1 %in% sign(xx.unlist)){
    color_down <- colorRampPalette(c("blue", "white"))(length(seq(quantile(xx.unlist, na.rm = TRUE, probs = (1-legend.limit)), 0, break_step))) #blue(n=length(breakSeq)-1) #, low="blue", mid = "white", high="red")
    color <- c(color,color_down)
  }
  if(1 %in% sign(xx.unlist)){
    color_up <- colorRampPalette(c("white", "red"))(length(seq(0, quantile(xx.unlist, na.rm = TRUE, probs = legend.limit), break_step)))
    color <- c(color,color_up)
  }
  #if((nrow(xx)>1 | ncol(xx)>1)){
  if(setWidth){
    #p <- p %>% layout(width=ncol(xx)*50, height = nrow(xx)*10)
    if(!is.na(file)){
      p <- heatmaply(xx, Colv = F, Rowv = rowClust, colors = color, custom_hovertext = hover, plot_method = "plotly", show_dendrogram = F, width = ncol(xx)*70, height = nrow(xx)*30, fontsize_col = fontsize_c, fontsize_row = fontsize_r, file = file) 
    }else{
      p <- heatmaply(xx, Colv = F, Rowv = rowClust, colors = color, custom_hovertext = hover, plot_method = "plotly", show_dendrogram = F, width = ncol(xx)*70, height = nrow(xx)*30, fontsize_col = fontsize_c, fontsize_row = fontsize_r) 
    }
  }else{
    if(!is.na(file)){
      p <- heatmaply(xx, Colv = F, Rowv = rowClust, colors = color, custom_hovertext = hover, plot_method = "plotly", show_dendrogram = F, fontsize_col = fontsize_c, fontsize_row = fontsize_r, file = file) 
    }else{
      p <- heatmaply(xx, Colv = F, Rowv = rowClust, colors = color, custom_hovertext = hover, plot_method = "plotly", show_dendrogram = F, fontsize_col = fontsize_c, fontsize_row = fontsize_r) 
    }
    p <- p %>% config(displayModeBar = TRUE)
    #}
    #p <- p %>% config(toImageButtonOptions=list(format="png", width="None", height="None"))
    #pheatmap(xx, cluster_cols=colClust, cluster_rows=rowClust, clustering_distance_rows = distMethod, clustering_distance_cols = distMethod,
    #         clustering_method = clusterMethod, annotation_col=annCol, annotation_row = annRow, 
    #         breaks = breakSeq, color = color, annotation_colors = annotation_colors,
    #         fontsize_row = fontsize_row, fontsize_col = fontsize_col, cutree_rows = clrn, 
    #         border_color = border_col, display_numbers = display_numbers)
  }
  return(p)
}


createLink <- function(val) {
  sprintf('<a href="https://www.ncbi.nlm.nih.gov/gene?term=(%s%%5BGene%%20Name%%5D)%%20AND%%20human%%5BOrganism%%5D" target="_blank">Info</a>', val)
}

plotExpression <- function(expr.df, padj.df, gene){
  lfc.gene <- data.frame(t(expr.df[gene,]))
  colnames(lfc.gene) <- c("symbol")
  lfc.gene$padj <- t(padj.df[gene,])
  lfc.gene$sig <- ifelse(lfc.gene$padj<0.0001,"***",ifelse(lfc.gene$padj<0.001,"**",ifelse(lfc.gene$padj<0.01,"*",NA)))
  lfc.gene$group <- sub("_.*","",rownames(lfc.gene))
  lfc.gene$time <- sub(".*_(.*_.*)","\\1",rownames(lfc.gene)) #factor(sub(".*_(.*)","\\1",rownames(lfc.gene)), levels = mixedsort(unique(sub(".*_(.*)","\\1",rownames(lfc.gene)))))
  lfc.gene$time <- ifelse(lfc.gene$time=="Mock_BPL", "Virus BPL:Mock BPL", ifelse(grepl("BPL",lfc.gene$time), "Virus 24h:Virus BPL", sub("Mock_","",lfc.gene$time)))
  lfc.min <- min(lfc.gene$symbol)
  lfc.max <- max(lfc.gene$symbol) + 0.25
  p1 <- ggplot(lfc.gene[grepl("h$",lfc.gene$time),], aes(x=factor(time, levels = mixedsort(unique(time))), y=symbol, group=group, color=group)) + geom_line() + geom_point() +
    labs(title = paste("Expression profile of", gene), y = "Log2FoldChange", x = "Time", color = "Virus", caption = "*: padj < 0.01    **: padj < 0.001    ***: padj < 0.0001") + 
    geom_text(aes(label=sig), nudge_y = 0.1, show.legend = FALSE) + ylim(c(lfc.min, lfc.max)) +
    theme(plot.caption = element_text(hjust = 0, size = 15, family = "Helvetica"), legend.text = element_text(size = 15, family = "Helvetica"), 
          legend.title = element_text(size = 15, face = "bold", family = "Helvetica"), axis.text = element_text(size = 15),
          axis.title = element_text(size = 15, face = "bold", family = "Helvetica"), plot.title = element_text(size = 20, family = "Helvetica", face = "bold")) +
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  p2 <-  ggplot(lfc.gene[grepl("BPL",lfc.gene$time),], aes(x=time, y=symbol, group=group, fill=group)) + geom_col(position = "dodge") +
    labs(title = paste("Expression profile of", gene), y = "Log2FoldChange", x = "", fill = "Virus", caption = "*: padj < 0.01    **: padj < 0.001    ***: padj < 0.0001") + 
    geom_text(aes(label=sig), position=position_dodge(width=0.9), vjust=-0.25, show.legend = F) + ylim(c(lfc.min, lfc.max)) +
    theme(plot.caption = element_text(hjust = 0, size = 15, family = "Helvetica"), legend.text = element_text(size = 15, family = "Helvetica"), 
          legend.title = element_text(size = 15, face = "bold", family = "Helvetica"), axis.text = element_text(size = 15),
          axis.title = element_text(size = 15, face = "bold", family = "Helvetica"), plot.title = element_text(size = 20, family = "Helvetica", face = "bold")) +
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  #p <- c(p1, p2)
  #gridExtra::grid.arrange(p1, p2, ncol = 1)
  #margin = theme(plot.margin = unit(c(2,2,2,2), "cm"))
  p <- gridExtra::arrangeGrob(p1, p2, ncol = 1, heights = c(10, 10))
  dev.off()
  return(p)
}

readRenviron(".Renviron")
if(Sys.getenv("DATADIR") != ""){
  shinyOptions(filedir = Sys.getenv("DATADIR"))
}
Sys.setenv("PATH" = paste(Sys.getenv("PATH"), "/root/miniconda3/bin", sep = .Platform$path.sep))
#addResourcePath("imgResources", "Documents/Virus_project/virus-shiny-app/")

#setwd("/tmp/")

datasetInput <- NULL
# Define UI for application that draws a histogram
ui <- fluidPage(
  tags$head(
    tags$style(
      HTML(".shiny-notification {
              height: 100px;
              width: 800px;
              position:fixed;
              top: calc(50% - 50px);;
              left: calc(50% - 400px);;
              font-size: 28px;
            }
            .tooltip > .tooltip-inner {
                width: 400px;
                color: white;
                background-color: grey;
            }
            .tooltip.right > .tooltip-arrow {
              border-right: 5px solid grey;
            }
            #help{
                  border-radius: 15px;
                  background-color: rgb(245, 245, 245, 0);
                  border-color: rgb(245, 245, 245, 0);
                  
                  }
            "
      ),
      "#shiny-modal img { max-width: 100%; }",
      "#uploadText{color: black;
                                 font-size: 20px;
                                 font-style: bold;
                                 }"
    )
  ),
  
  # Application title
  titlePanel("Differential Gene Expression Analysis"),
  
  # Sidebar with a slider input for number of bins 
  navbarPage(
    "", 
    tabPanel("Download Data",
             #sidebarPanel(
             #    fileInput("file", label = h5(br("File input")), multiple = T),
             #),
             sidebarPanel(
               #textOutput("uploadText"),
               #br(),
               div(style = "text-align:right", actionButton("help0", label = div(strong("Help"), icon("question")))),
               htmlOutput("download")
               #tableOutput("table"),
               #h2("Uploaded files"),
             ),
             mainPanel(
               downloadButton("filesDown", "Download selected files")
               #textOutput("debug", container = pre)
             )
    ),
    tabPanel("Volcano Plot / MA-Plot",
             sidebarPanel(
               div(style = "text-align:right", actionButton("help1", label = div(strong("Help"), icon("question")))),
               htmlOutput("fileSelect"),
               radioButtons("plotTypeSingle", label = h5(strong("Displayed plot")), choices = c("Volcano plot", "MA-plot"), selected = "Volcano plot"),
               sliderInput("LFCsingle",
                           h5(strong("LFC cutoff:")),
                           min = 0,
                           max = 10,
                           value = 1, 
                           step = 0.5),
               br()
               #downloadButton("downSingle", "Download")
             ),
             mainPanel(
               #downloadButton("downPlot", "Download plot"),
               #actionButton("export", "Export plot"),
               conditionalPanel(
                 condition = 'input.plotTypeSingle == "Volcano plot"',
                 plotlyOutput("volcanoPlot", height = "500px"),
                 br(),
                 dataTableOutput("clickVP")
               ),
               conditionalPanel(
                 condition = 'input.plotTypeSingle == "MA-plot"',
                 plotlyOutput("MAPlot", height = "500px"),
                 br(),
                 dataTableOutput("clickMA")
               )
             )),
    tabPanel("Compare Time Points",
             sidebarPanel(
               div(style = "text-align:right", actionButton("help2", label = div(strong("Help"), icon("question")))),
               radioButtons("plotType", label = h5(strong("Displayed plot")), choiceNames = c("Venn diagram", "UpSet plot"), choiceValues = c("Venn", "UpSet"), selected = "Venn"),
               selectInput("select_v", label = h5(strong("Select virus")), 
                           choices = list("H1N1", "H5N1", "RVFV", "SFSV", "RSV", "NiV", "EBOV", "MARV", "LASV"), 
                           selected = NULL),
               htmlOutput("selectTime"),
               sliderInput("LFC",
                           h5(strong("LFC cutoff:")),
                           min = 0,
                           max = 10,
                           value = 1, 
                           step = 0.5),
               actionButton("addHeat", "Show heatmap"),
               #downloadButton("downHeatTime", "Download heatmap"),
               br(),
               br(),
               h5(strong("Download table")),
               downloadButton("down_vXLSX", "Download as xlsx", class = "butt1"),
               tags$head(tags$style(".butt1{
                                    border: 2px solid black;
                                    font-size: 16px;}")),
               br(),
               downloadButton("down_vCSV", "Download as csv")
             ),
             mainPanel(
               conditionalPanel(
                 condition = 'input.plotType == "UpSet"',
                 upsetjsOutput("upset"), #, click = "upsetClick"),
                 br(),
                 DT::dataTableOutput("clickedElements"),
                 br(),
                 br()
                 #plotOutput("heatmapTimeUp", height = "750px")
               ),
               conditionalPanel(
                 condition = 'input.plotType == "Venn"',
                 upsetjsOutput("upsetVenn"), #, click = "upsetClick"),
                 br(),
                 DT::dataTableOutput("clickedElementsVenn"),
                 br(),
                 br()
                 #plotOutput("heatmapTimeVenn", height = "750px")
               ),
               #conditionalPanel("output.show", plotlyOutput("heatmapTime", height = "750px")),
               plotlyOutput("heatmapTime", height = "750px") %>% withSpinner(type = 5, color.background = "white", color = "grey", size = 1.5)
             )
    ),
    tabPanel("Gene Expression",
             sidebarPanel(
               div(style = "text-align:right", actionButton("help3", label = div(strong("Help"), icon("question")))),
               checkboxGroupInput("selectVirus", label = h5(strong("Select viruses")), 
                                  choices = list("H1N1", "H5N1", "RVFV", "SFSV", "RSV", "NiV", "EBOV", "MARV", "LASV"), 
                                  selected = NULL),
               textInput("gene", label = h5(strong("Gene symbol")), value = NULL, placeholder = "e.g. TNF or cxcl2"),
               h5(strong("Download plot")),
               downloadButton("downGeneXpng","Download as png"),
               downloadButton("downGeneXsvg","Download as svg")
             ), 
             mainPanel(
               #textOutput("debug", container = pre),
               #verbatimTextOutput("geneX"),
               plotOutput("geneX", height = 1000),
               br()
               #tableOutput("sig")
             )),
    tabPanel("Virus Comparison",
             sidebarPanel(
               div(style = "text-align:right", actionButton("help4", label = div(strong("Help"), icon("question")))),
               checkboxGroupInput("select", label = h5(strong("Select viruses")), 
                                  choices = list("H1N1", "H5N1", "RVFV", "SFSV", "RSV", "NiV", "EBOV", "MARV", "LASV"), 
                                  selected = NULL), #c("H1N1", "H5N1", "MERS", "CoV229E", "RVFV", "SFSV", "RSV", "NIV", "EBOV")),
               htmlOutput("selectCond"),
               em("(BPL samples were taken after 24h)"),
               #checkboxGroupInput("selectCond", label = h5(strong(" Select conditions to include")), 
               #                   choices = list("")),
               #choices = c("infected"="Mock_.*h", "inactivated control"="24h_vs_.*BPL", "uninfected control"="Mock_BPL")),
               br(),
               br(),
               sliderInput("LFCall",
                           h5(strong("LFC cutoff:")),
                           min = 0,
                           max = 10,
                           value = 1, 
                           step = 0.5),
               sliderInput("limit",
                           h5(strong("Number of intersections:")),
                           min = 0,
                           max = 100,
                           value = 20, 
                           step = 5),
               #br(),
               br(),
               h5(strong("Download table (includes LFC and padj)")),
               downloadButton("downallXLSX", "Download as xlsx", class = "butt"),
               tags$head(tags$style(".butt{
                                    border: 2px solid black;
                                    font-size: 16px;}")),
               br(),
               downloadButton("downallCSV", "Download as csv")
             ),
             mainPanel(
               upsetjsOutput("upsetAll"), #, click = "upsetClick"),
               br(),
               DT::dataTableOutput("clickedElementsAll")
             )
    ),
    tabPanel("Heatmap Virus Comparison",
             sidebarPanel(
               div(style = "text-align:right", actionButton("help5", label = div(strong("Help"), icon("question")))),
               h5("Please select an intersection in the tab 'Virus Comparison' to show the corresponding heatmap.")
             ),
             mainPanel(
               #downloadButton("downHeatAll", "Download heatmap"),
               plotlyOutput("heatmapAll", height = 1500) %>% withSpinner(type = 5, color.background = "white", color = "grey", size = 2)
             )),
    tabPanel("SNP Analysis",
             sidebarPanel(
               div(style = "text-align:right", actionButton("help6", label = div(strong("Help"), icon("question")))),
               selectInput("virusSNP", label = h5(strong("Select virus")), 
                           choices = list("H1N1", "H5N1", "RVFV", "SFSV", "RSV", "NiV", "EBOV", "MARV", "LASV"), 
                           selected = NULL)
             ),
             mainPanel(
               #plotOutput("heatSNP", height = "600px", click = "SNPclick"),
               #DT::dataTableOutput("SNPdata")
               plotlyOutput("heatSNP", height = "1000px")
               #verbatimTextOutput("SNPdata")
             ))
  )
)

# Define directory containing data
filedir <- getShinyOption("filedir")
files.list <- list.files(filedir, pattern = ".*.csv", full.names = T)
vcf.files <- list.files(filedir, pattern = ".*[1|2].vcf", full.names = T)

# Define server logic 
server = function(input, output, session) {
  
  # Define reactive values
  data <- reactiveValues(d=NULL, df=NULL, dM=NULL, dfM=NULL, heat=FALSE, heat_data = NULL, heat_hover = NULL, dt=NULL, snp=NULL, all.df=NULL, virus.df=NULL, plotX=NULL)
  
  # Read data
  withProgress(message = 'PROCESSING DATA...', detail = "This may take a while...", value = 0,{
    datasetInput <- lapply(files.list, function(x){
      f <- read.csv(x, header = TRUE, stringsAsFactors = F, row.names = 1)
      incProgress(1/length(files.list))
      return(f)
    })
    #Sys.sleep(5)
    names(datasetInput) <- sub("deseq2_results_(.*).csv","\\1",basename(files.list))
    names(datasetInput) <- sub("Vs","vs", names(datasetInput))
    datasetInput <- datasetInput[mixedorder(sub("vs_Mock_","",names(datasetInput)))]
    vcf.list <- lapply(vcf.files, read.vcfR)
    names(vcf.list) <- sub(".vcf", "", basename(vcf.files))
  })
  
  # round numeric values
  #datasetInput <- lapply(datasetInput, function(x){
  #  x[,c(5:ncol(x))] <- signif(x[,c(5:ncol(x))], 4)
  #  return(x)
  #})
  
  # Help text
  observeEvent(input$help0,{
    showModal(modalDialog(
      title = tags$div("Help for ", tags$b("Download Data")),
      size = "l",
      "Results of the differential gene expression analysis with DESeq2 in table format. For each gene the table contains:",
      tags$div(
        tags$ul(
          tags$li("UNIPROT: Uniprot identifier"),
          tags$li("GENENAME: gene name"),
          tags$li("PATH: KEGG pathways associated with a gene"),
          tags$li("baseMean: base mean of normalized read counts"),
          tags$li("log2FoldChange: log2FoldChange (LFC) for conditions treated vs untreated"),
          tags$li("lfcSE: standard error"),
          tags$li("pvalue: p-value for LFC (describes the probability to get this LFC by chance)"),
          tags$li("padj: Benjamini-Hochberg adjusted p-value (p-value adjustment for multiple hypothesis testing)"),
          tags$li("log10(padj): logarithmic padj used in Volcano plot"),
          tags$li("normalized_[sample-name]: normalized read counts for every replicate")
        )),
      "These data are used as a basis for all analyses and visualizations in the following tabs.",
      easyClose = TRUE
    ))
  })
  
  output$uploadText <- renderText({
    if(!is.null(datasetInput)){
      print("Uploaded data:")
    }
  })
  
  output$table <- renderTable(
    data.frame(names(datasetInput))
    , colnames = F, rownames = F)
  
  output$download <- renderUI({
    checkboxGroupInput("filesToDown", "Choose files to download", as.list(c("all",names(datasetInput))))
  })
  
  output$filesDown <- downloadHandler(
    filename = function(){
      "data.zip"
    },
    content = function(f){
      files <- NULL
      if("all" %in% input$filesToDown){
        files <- files.list
      }else{
        for (i in input$filesToDown){
          #write each sheet to a csv file, save the name
          fileName <- files.list[grep(sub("vs","Vs",i),files.list)] #paste(i,".csv",sep = "")
          files <- c(files,fileName)
        }
      }
      zip(f, files)
    },
    contentType = "csv"
  )
  
  ## Volcano Plot / MA-Plot
  
  # Help text
  observeEvent(input$help1,{
    showModal(modalDialog(
      title = tags$div("Help for ", tags$b("Volcano Plot/MA-Plot")),
      size = "l",
      tags$div("Volcano plots and MA-plots for uploaded files. Each file contains the results of the differential gene expression analysis with DEseq2.
      The sidebar shows the available options to adapt the plots:"),
      img(src="volcano_sidebar.png"),
      tags$div("The Volcano plot shows the log2FoldChange (LFC) and the negative log10 of the adjusted p-value (padj) for the chosen comparison. Significantly differentially expressed genes are colored blue (down-regulated) or red (up-regulated). A gene is considered differentially expressed if: LFC > ", tags$i("LFC cutoff") ,"and padj < 0.05 and a normalized gene count of at least one sample ≥ 10."),
      HTML("<br><br>"),
      img(src="Volcano_plot.png"),
      HTML("<br><br>"),
      tags$div("The modebar in the upper right corner of the plot shows:"),
      HTML("<br><br>"),
      img(src="volcano_modebar.png"),
      HTML("<br><br>"),
      tags$div("Further information for one or multiple points, like gene symbol, LFC and padj, can be shown by mouse over."),
      img(src="volcano_plot_select.png"),
      tags$div("The selected genes are shown in tabular form along with their LFCs and padjs. A link to the NCBI database is provided for each gene. The link opens in a new tab outside the application. The table can be searched and all columns can be sorted in ascending and descending order."),
      img(src="Volcano_table.png"),
      easyClose = TRUE
    ))
  })
  
  output$fileSelect <- renderUI({
    selectInput("fileToPlot", label = "Choose file", choices = names(datasetInput))
  })
  
  # download results        
  output$downSingle <- downloadHandler(
    filename = function(){
      return(paste0(input$fileToPlot, "_selected", ".csv"))},
    content = function(f){
      write.csv(data.frame("Symbol"=unlist(input$upsetAll_click$elems)), f)
    },
    contentType = "csv")
  
  # set plotly_click to NULL on input change
  observeEvent(input$fileToPlot,{
    data$d <- NULL
    data$dV <- NULL
    data$dM <- NULL
    data$dfM <- NULL
  })
  
  # save clicked elements to vector
  observeEvent(event_data("plotly_selected", source = "V"), {
    req(datasetInput)
    req(input$fileToPlot)
    data$d <- event_data("plotly_selected", source = "V")
    a <- data$dV
    #a <- c(a, data$d$customdata)
    data$dV <- a
  })
  
  observeEvent(event_data("plotly_selected", source = "M"), {
    req(datasetInput)
    req(input$fileToPlot)
    data$dM <- event_data("plotly_selected", source = "M")
    #a <- data$dfM
    #a <- c(a, data$dM$customdata)
    #data$dfM <- a
  })
  
  # display table of clicked elements and associated LFC and padj in Volcano plot
  output$clickVP <- renderDataTable({
    d <- data$d
    d.V <- data$dV
    if(is.null(d)){
      return(NULL)
    }else{
      df <- datasetInput[[input$fileToPlot]]
      df <- df[df$SYMBOL %in% d$customdata,c("SYMBOL","log2FoldChange","padj")]
      df[,c(2:ncol(df))] <- signif(df[,c(2:ncol(df))], 4)
      df$NCBI <- createLink(df$SYMBOL)
      return(df)
    }
  }, escape = F)
  
  # plot Volcano plot of selected sample
  output$volcanoPlot <- renderPlotly({
    p <- plotVolcano()
    p %>% config(toImageButtonOptions=list(format="png", width=800, height=600, scale=2))
  })
  
  plotVolcano <- reactive({
    if(is.null(datasetInput) | is.null(input$fileToPlot))
      return(NULL)
    file.data <- datasetInput[[input$fileToPlot]]
    file.data$col <- ifelse(file.data$log2FoldChange > input$LFCsingle & file.data$padj < 0.05 & apply(file.data[,grep("normalized", colnames(file.data))],1,max) >= 10, "up", 
                            ifelse(file.data$log2FoldChange < -(input$LFCsingle) & file.data$padj < 0.05 & apply(file.data[,grep("normalized", colnames(file.data))],1,max) >= 10, "down", "not significant"))
    p <- plot_ly(file.data, type = "scatter", x = ~log2FoldChange, y = ~-log10(padj), color = ~col, colors = c("up"="red","down"="steelblue3","not significant"="black"), 
                 mode = "markers", marker = list(size = 5), customdata = ~SYMBOL,
                 text = ~paste("Gene: ", SYMBOL, '<br>LFC: ', signif(log2FoldChange, 4), '<br>padj: ', signif(padj, 4)),
                 source = "V")
    p <- p %>% layout(legend = list(orientation = "h", y = -0.2), dragmode = "select") %>% config(displayModeBar = TRUE)
    return(p)
  })
  
  # display table of clicked elements and associated LFC and padj in MA-plot
  output$clickMA <- renderDataTable({
    d <- data$dM
    d.MA <- data$dfM
    if(is.null(d)){
      return(NULL)
    }else{
      df <- datasetInput[[input$fileToPlot]]
      df <- df[df$SYMBOL %in% d$customdata,c("SYMBOL","log2FoldChange","padj")]
      df[,c(2:ncol(df))] <- signif(df[,c(2:ncol(df))], 4)
      df$NCBI <- createLink(df$SYMBOL)
      return(df)
    }
  }, escape = F)
  
  # plot MA-plot of selected sample
  output$MAPlot <- renderPlotly({
    p <- plotMA()
    p %>% config(toImageButtonOptions=list(format="png", width=800, height=600, scale=2))
  })
  
  plotMA <- reactive({
    if(is.null(datasetInput) | is.null(input$fileToPlot)){
      return(NULL)
    }else{
      file.data <- datasetInput[[input$fileToPlot]]
      file.data$col <- ifelse(file.data$log2FoldChange > input$LFCsingle & file.data$padj < 0.05 & apply(file.data[,grep("normalized", colnames(file.data))],1,max) >= 10, "up", 
                              ifelse(file.data$log2FoldChange < -(input$LFCsingle) & file.data$padj < 0.05 & apply(file.data[,grep("normalized", colnames(file.data))],1,max) >= 10, "down", "not significant"))
      file.data$baseMeanSample <- rowMeans(file.data[,grep("normalized", colnames(file.data))])
      p <- plot_ly(file.data, type = "scatter", x = ~log2(baseMean), y = ~log2FoldChange, color = ~col, colors =c("up"="red","down"="steelblue3","not significant"="black"),
                   mode = "markers", marker = list(size = 5), customdata = ~SYMBOL,
                   text = ~paste("Gene: ", SYMBOL, '<br>LFC: ', signif(log2FoldChange, 4), '<br>padj: ', signif(padj, 4)),
                   source = "M")
      p <- p %>% layout(legend = list(orientation = "h", y = -0.2), dragmode = "select") %>% config(displayModeBar = TRUE)
      return(p)
    }
  })
  
  ## Compare time points of defined virus
  
  # Help text 
  observeEvent(input$help2,{
    showModal(modalDialog(
      title = tags$div("Help for ", tags$b("Compare Time Points")),
      size = "l",
      tags$div("This tab is meant to compare the differential gene expression between different time points of infection. Each gene list contains all genes differentially expressed with |LFC| > ", tags$i("LFC cutoff")," and padj < 0.05 and a normalized gene count of at least one sample ≥ 10. Intersections can either be shown as Venn diagram or UpSet plot ", tags$a("(info)", href="https://jku-vds-lab.at/tools/upset/", target="_blank"), ". The sidebar shows the available options to adapt the plots:"),
      img(src="time_sidebar.png"),
      "The Venn diagram shows the intersections of genes from the selected gene lists. For reasons of clarity, a maximum number of 5 lists can be displayed.",
      img(src="time_venn.png"),
      tags$div("Clicking an intersection will display underlying genes in a sortable table along with log2FoldChange (LFC), adjusted p-value (padj) and a link to the NCBI database. The table can be downloaded as xlsx or csv file."),
      HTML("<br><br>"),
      img(src="time_venn_select.png"),
      HTML("<br><br>"),
      tags$div("The UpSet plot shows the intersections of genes in a matrix-like layout. This approach allows to compare a larger number of sets. Each bar shows an intersection of genes with the corresponding number. In the respective column below you can see the samples involved, marked with black dots. For example in the plot below 4467 genes are differentially expressed (see bar) after 6h, 12h and 24h (see black dots in column below)."),
      img(src="time_upset.png"),
      tags$div("Genes in an intersection can be listed in a sortable table along with LFC, padj and a link to the NCBI database by clicking on the appropriate bar. The table can be downloaded as xlsx or csv file."),
      HTML("<br><br>"),
      img(src="time_upset_select.png"),
      HTML("<br><br>"),
      tags$div("The selected genes can also be shown in a heatmap by clicking on the button in the sidebar:"),
      img(src="time_heatmap.png"),
      HTML("<br><br>"),
      tags$div("The modebar in the upper right corner of the heatmap shows the following functions, that allow the user to interact with the plot:"),
      HTML("<br><br>"),
      img(src="heatmap_modebar.png"),
      HTML("<br><br>"),
      tags$div("The user is able to zoom into the heatmap by dragging a box over the area of interest. To get information about the underlying data mouseover a specific cell."),
      HTML("<br><br>"),
      img(src="time_heatmap_zoom.png"),
      img(src="time_heatmap_zoom2.png"),
      easyClose = TRUE
    ))
  })
  
  # download genes of selected plot area as csv or xlsx
  output$down_vCSV <- downloadHandler(
    filename = function(){
      if(input$plotType == "UpSet"){
        return(paste0(input$select_v, "_", sub("&","_",input$upset_click$name), ".csv"))
      }else{
        return(paste0(input$select_v, "_", sub("&","_",input$upsetVenn_click$name), ".csv"))
      }
    },
    content = function(f){
      #if(input$plotType == "UpSet"){
      #write.csv(data.frame("Symbol"=unlist(input$upset_click$elems)), f, row.names = F)
      write.csv(virus.df(), f, row.names = F)
      #}
    },
    contentType = "csv")
  
  output$down_vXLSX <- downloadHandler(
    filename = function(){
      if(input$plotType == "UpSet"){
        return(paste0(input$select_v, "_", sub("&","_",input$upset_click$name), ".xlsx"))
      }else{
        return(paste0(input$select_v, "_", sub("&","_",input$upsetVenn_click$name), ".xlsx"))
      }
    },
    content = function(f){
      write.xlsx(virus.df(), f)
    })
  
  # download heatmap
  output$downHeatTime <- downloadHandler(
    filename = function(){
      return(paste0(input$select_v, ".pdf"))},
    content = function(f){
      data.df <- data$heat_data
      hover <- data$heat_hover
      p <- plotHeatmap(data.df, colClust = F, border_col = NA, fontsize_r = 10, hover = hover, file = f)
      #orca(f, p, format = "pdf", more_args = c("--disable-gpu", "--enable-webgl"))
    }
  )
  
  # select times to plot
  output$selectTime <- renderUI({
    checkboxGroupInput("time", label = "Choose time points (for Venn: max. 5)",
                       choiceNames = unique(sapply(names(datasetInput)[grep(input$select_v, names(datasetInput))], 
                                                   function(x){
                                                     #s <- sub("^[^_]*_","",x)
                                                     #s <- ifelse(grepl("Mock",x), x, sub("vs_.*_","vs_Virus_",x)) #ifelse(grepl("Mock",x), sub("_Mock_.*","_Mock",s), s)
                                                     s <- sub("_vs_", " : ", x)
                                                     #s <- sub("48h$", "48h (HCV only)", s)
                                                     return(s)
                                                   }, USE.NAMES = F)),
                       choiceValues = sub(".*vs","vs",names(datasetInput)[grep(input$select_v, names(datasetInput))]))
  })
  
  
  # set heatmap to NULL at selection of new virus
  observe({
    input$select_v
    input$plotType
    data$heat <- FALSE   
    data$dt <- NULL
  })
  
  # clear data table and heatmap if plot type changes
  observeEvent({
    input$upset_click$elems
  },
  {
    data$heat <- FALSE
    data$dt <- 1
    #data_heat_time()
  })
  
  observeEvent({
    input$upsetVenn_click$elems
  },
  {
    data$heat <- FALSE
    data$dt <- 1
    #data_heat_time()
  })
  
  # add heatmap
  observeEvent(input$addHeat, {
    data$heat <- TRUE
    #data_heat_time()
  })
  
  output$show <- reactive({
    print(data$heat)
    return(data$heat)
  })
  
  output$heatmapTime <- renderPlotly({
    if(!data$heat){
      data$heat_data <- NULL
      data$heat_hover <- NULL
      #return(NULL)
    }else{
      data.df <- virus.df()
      if(is.null(data.df)){
        data$heat_data <- NULL
        data$heat_hover <- NULL
        #return(NULL)
      }
      data.rows <- data.df$SYMBOL
      data.df <- data.df[,-c(1,ncol(data.df)), drop=F]
      if(nrow(data.df)>1){
        data.df <- apply(data.df,2,as.numeric)
      }else{
        data.df <- as.data.frame(t(apply(data.df,2,as.numeric)))
      }
      rownames(data.df) <- data.rows
      data.df <- data.df[,c(TRUE, FALSE), drop=F]
      colnames(data.df) <- sub(".LFC","",colnames(data.df))
      hover <- Reduce(function(x,y)merge(x,y,by="SYMBOL",all=T), sapply(colnames(data.df), function(s){
        s.gr <- sub("_.*","",strsplit(s,"_vs_")[[1]])
        dataset <- datasetInput[[s]]
        dataset[,c(5:ncol(dataset))] <- signif(dataset[,c(5:ncol(dataset))], 4) 
        dataset <- dataset[dataset$SYMBOL %in% rownames(data.df), ]
        dataset$c1 <- apply(dataset[ , grep(paste0("normalized.*", s.gr[1]), colnames(dataset)), drop = F] , 1 , paste , collapse = "; " )
        dataset$c2 <- apply(dataset[ , grep(paste0("normalized.*", s.gr[2]), colnames(dataset))] , 1 , paste , collapse = "; " )
        df.tmp <- data.frame(dataset[,"SYMBOL"], sprintf("Gene: %s<br />Sample: %s<br />LFC: %s<br />padj: %s<br />%s: %s<br />%s: %s", 
                                                         dataset[,"SYMBOL"], s, dataset[,"log2FoldChange"], dataset[,"padj"],
                                                         s.gr[1], dataset[,"c1"], s.gr[2], dataset[,"c2"]))
        colnames(df.tmp) <- c("SYMBOL", s)
        return(df.tmp)
      }, simplify = F)
      )
      rownames(hover) <- hover$SYMBOL
      hover <- subset(hover, select=-c(SYMBOL))
      hover <- hover[order(rownames(hover)),]
      data.df <- data.df[order(rownames(data.df)),]
      data$heat_data <- data.df
      data$heat_hover <- hover
      plotHeatmap(data.df, colClust = F, rowClust = T, border_col = NA, fontsize_r = 10, hover = hover)
    }
    
    
  })
  
  virus.df <- reactive({
    virus_id <- NULL
    if(input$plotType=="UpSet"){
      virus_id <- unlist(input$upset_click$elems)
    }
    if(input$plotType=="Venn"){
      virus_id <- unlist(input$upsetVenn_click$elems)
    }
    if(is.null(input$select_v) | is.null(input$time) | is.null(virus_id)){
      return(NULL)
    }else{
      data.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL", all = T),lapply(names(datasetInput[grep(paste0(input$select_v, ".*_", input$time, collapse = "|"), names(datasetInput))]), function(x){
        y <- datasetInput[[x]]
        y <- y[y$SYMBOL %in% virus_id, c("SYMBOL","log2FoldChange","padj"), drop=F]
        colnames(y) <- c("SYMBOL", paste0(x,".LFC"), paste0(x,".padj"))
        return(y)
      })
      )
      data.df <- data.df[,c(1,mixedorder(colnames(data.df)[-1])+1),drop=F]
      data.df$Link <- createLink(data.df$SYMBOL)
      return(data.df)
    }
  })
  
  # display table with elements in clicked plot area of UpSet plot
  output$clickedElements <- renderDataTable({
    dt <- data$dt
    if(is.null(dt)){
      data.frame("SYMBOL"=NULL, "LinkToNCBI"=NULL)
    }else{
      df <- virus.df()
      df[,c(2:(ncol(df)-1))] <- signif(df[,c(2:(ncol(df)-1))], 4)
      return(df)
    }
  }, escape = F)
  
  # plot UpSet plot of selected virus und times
  output$upset <- renderUpsetjs({
    if(is.null(input$select_v) | is.null(input$time)){
      return(NULL)
    }else{
      id.list <- lapply(datasetInput[grep(paste0(input$select_v, ".*_", input$time, collapse = "|"), names(datasetInput))], function(x){
        return(x[which(abs(x$log2FoldChange) > input$LFC & x$padj < 0.05 & apply(x[,grep("normalized", colnames(x))],1,max) >= 10),"SYMBOL"])
      })
      id.list <- id.list[sapply(id.list, length) > 0]
      names(id.list) <- sub(".*_vs_", "", names(id.list))
      names(id.list) <- sub("Mock_", "", names(id.list))
      if(length(id.list)>0){
          upsetjs(sizingPolicy = upsetjsSizingPolicy(padding = 10)) %>% upsetjs::fromList(id.list) %>% generateDistinctIntersections() %>%
          chartFontSizes(font.family = "Helvetica", chart.label = "18px", set.label = "16px", bar.label = "16px", axis.tick = "16px", export.label = "13px") %>% 
          chartLayout(padding = 40, bar.padding = 0.2) %>% interactiveChart() %>% chartProps(exportButtons=list(share=FALSE, vega=FALSE, dump=FALSE))
      }
    }
  })
  
  # display table with elements in clicked plot area of Venn plot
  output$clickedElementsVenn <- renderDataTable({
    dt <- data$dt
    if(is.null(dt)){
      data.frame("SYMBOL"=NULL, "LinkToNCBI"=NULL)
    }else{
      df <- virus.df()
      df[,c(2:(ncol(df)-1))] <- signif(df[,c(2:(ncol(df)-1))], 4)
      return(df)
    }
  }, escape = F)
  
  # plot Venn diagram of selected virus und times
  output$upsetVenn <- renderUpsetjs({
    if(is.null(input$select_v) | is.null(input$time)){
      return(NULL)
    }else{
      id.list <- lapply(datasetInput[grep(paste0(input$select_v, ".*_", input$time, collapse = "|"), names(datasetInput))], function(x){
        return(x[which(abs(x$log2FoldChange) > input$LFC & x$padj < 0.05 & apply(x[,grep("normalized", colnames(x))],1,max) >= 10),"SYMBOL"])
      })
      id.list <- id.list[sapply(id.list, length) > 0]
      names(id.list) <- sub(".*_vs_", "", names(id.list))
      names(id.list) <- sub("Mock_", "", names(id.list))
      if(length(id.list)>0){
        upsetjsVennDiagram() %>% upsetjs::fromList(id.list) %>% chartProps(exportButtons=list(share=FALSE, vega=FALSE, dump=FALSE)) %>%
          chartFontSizes(font.family = "Helvetica", chart.label = "18px", set.label = "16px", bar.label = "16px", axis.tick = "16px", export.label = "13px", value.label = "14px") %>% interactiveChart()
      }
    }
  })
  
  ## Gene expression plot
  
  # Help text
  observeEvent(input$help3,{
    showModal(modalDialog(
      title = tags$div("Help for ", tags$b("Gene Expression")),
      size = "l",
      tags$div("This tab is meant to explore the expression profile of genes over time.",
               "The sidebar shows the available options to adapt the plots:"),
      img(src="expression_sidebar.png"),
      HTML("<br><br>"),
      tags$div("The plot shows the expression of the selected gene for the chosen viruses as log2FoldChange (LFC) over time. The stars represent the significance of the LFC."),
      HTML("<br><br>"),
      img(src="expression_plot_upper.png"),
      HTML("<br><br>"),
      tags$div("The plot shows the expression of the selected gene as LFC. The left hand side shows the 24 h time point of the selected viruses compared to the inactivated virus (Virus BPL), whereas the right hand side shows the inactivated virus compared to the uninfected control (Mock BPL)."),
      HTML("<br><br>"),
      img(src="expression_plot_lower.png"),
      easyClose = TRUE
    ))
  })
  
  output$debug <- renderPrint({
    print("Installed packages:")
    print(installed.packages()[,c(1,3)])
  })
  
  # download gene expression plot
  output$downGeneXpng <- downloadHandler(
    filename = function(){
      return(paste0(paste(input$selectVirus, collapse = "_"), "_", toupper(input$gene), ".png"))},
    content = function(f){
      ggsave(f, data$plotX, "png", width = 8, height = 10, dpi = 400)
    }
  )
  
  output$downGeneXsvg <- downloadHandler(
    filename = function(){
      return(paste0(paste(input$selectVirus, collapse = "_"), "_", toupper(input$gene), ".svg"))},
    content = function(f){
      svg(f, width = 8, height = 10, family = "Helvetica")
      plot(data$plotX)
      dev.off()
      #ggsave(f, data$plotX, "svg", width = 8, height = 10, dpi = 400)
    }
  )
  
  # plot gene expression of selected viruses over time
  output$geneX <- renderPlot({
    p <- data$plotX
    if(is.null(p)){
      return(NULL)
    }else{
      plot(p)
    }
  })
  
  plot_geneX <- reactive({
    list(input$selectVirus, input$gene)
  })
  
  observeEvent(plot_geneX(), {
    if(is.null(input$selectVirus) | is.null(input$gene)){
      return(NULL)
    }else{
      lfc.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL",all=T),lapply(names(datasetInput[grep(paste(input$selectVirus, collapse = "|"), names(datasetInput))]), function(x){
        x.df <- data.frame(datasetInput[[x]]$SYMBOL,datasetInput[[x]]$log2FoldChange)
        colnames(x.df) <- c("SYMBOL",x)
        return(x.df)}))
      rownames(lfc.df) <- lfc.df$SYMBOL
      lfc.df <- lfc.df[,-1]
      padj.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL",all=T),lapply(names(datasetInput[grep(paste(input$selectVirus, collapse = "|"), names(datasetInput))]), function(x){
        x.df <- data.frame(datasetInput[[x]]$SYMBOL,datasetInput[[x]]$padj)
        colnames(x.df) <- c("SYMBOL",x)
        return(x.df)}))
      rownames(padj.df) <- padj.df$SYMBOL
      padj.df <- padj.df[,-1]
      data$plotX <- plotExpression(expr.df = lfc.df, padj.df = padj.df, gene = toupper(input$gene))
    }
  })
  
  ## Virus comparison
  
  # select time points to include in comparison
  output$selectCond <- renderUI({
    checkboxGroupInput("cond", label = "Choose time points", 
                       choiceNames = unique(sapply(unique(sub("[^_]*_","Virus_",names(datasetInput))), 
                                                   function(x){
                                                     s <- ifelse(grepl("Mock",x), x, sub("vs_.*_","vs_Virus_",x))
                                                     s <- sub("_vs_", " : ", s)
                                                     s <- sub("48h$", "48h (HCV only)", s)
                                                     return(s)
                                                   }, USE.NAMES = F)),
                       choiceValues = unique(sapply(unique(sub(".*vs","vs",names(datasetInput))), function(x) if(grepl("Mock",x)){return(x)}else{return(sub("_.*_","_Virus_",x))}, USE.NAMES = F)))
  })
  
  # Help text
  observeEvent(input$help4,{
    showModal(modalDialog(
      title = tags$div("Help for ", tags$b("Virus Comparison")),
      size = "l",
      tags$div("This tab is meant to compare the genes that are differentially expressed after the infection with different viruses.",
               "The sidebar shows the available options to adapt the plots:"),
      img(src="viruses_sidebar.png"),
      HTML("<br><br>"),
      tags$div("The UpSet plot shows the intersections of genes in a matrix-like layout ", tags$a("(info)", href = "https://jku-vds-lab.at/tools/upset/", target="_blank"), ".",
               "This approach allows to compare a large number of sets. Each bar shows an intersection of genes with the corresponding number. In the respective column below you can see the samples involved, marked with black dots.",
               "For each virus all genes differentially expressed in at least one of the selected conditions are pooled. A gene is considered differentially expressed if: log2FoldChange (LFC) > ",tags$i("LFC cutoff")," and adjusted p-value (padj) < 0.05 and a normalized gene count of at least one sample ≥ 10."),
      img(src="viruses_upset.png"),
      tags$div("Genes in an intersection can be listed in a sortable table along with a link to the NCBI database by clicking on the appropriate bar. The table can be downloaded in xlsx or csv format and then also contains the LFC and padj for all selected conditions.",
               "The selected bar in the example (colored orange) shows, that there are 218 genes that are diff. expressed in at least one of the included time points (see bar) in MARV, EBOV and NIV (see black dots in column below)."),
      HTML("<br><br>"),
      img(src="viruses_upset_select.png"),
      HTML("<br><br>"),
      tags$div("To display a selected intersection as heatmap click on the respective intersection, wait until the table is updated and switch to the tab 'Heatmap Virus Comparison'."),
      easyClose = TRUE
    ))
  })
  
  all.df <- reactive({
    gene_id <- unlist(input$upsetAll_click$elems)
    cond <- ifelse(grepl("Virus",input$cond), sub("Virus", "(?!Mock).*",input$cond), input$cond)
    data.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL", all = T),lapply(names(datasetInput[grep(paste(as.vector(outer(input$select, cond, paste, sep=".*")), collapse = "|"), names(datasetInput), perl = T)]), function(x){
      y <- datasetInput[[x]]
      y <- y[y$SYMBOL %in% gene_id, c("SYMBOL","log2FoldChange", "padj")]
      colnames(y) <- c("SYMBOL", paste0(x, ".LFC"), paste0(x,".padj"))
      return(y)
    })
    )
    data.df <- data.df[,c(1,mixedorder(colnames(data.df)[-1])+1),drop=F]
    data.df$Link <- createLink(data.df$SYMBOL)
    return(data.df)
  })
  
  # download genes of selected plot area as csv or xlsx
  output$downallCSV <- downloadHandler(
    filename = function(){
      return(paste0(gsub("&","_",unlist(input$upsetAll_click$name)), ".csv"))},
    content = function(f){
      write.csv(all.df(), f, row.names = F)
    },
    contentType = "csv")
  
  output$downallXLSX <- downloadHandler(
    filename = function(){
      return(paste0(gsub("&","_",unlist(input$upsetAll_click$name)), ".xlsx"))},
    content = function(f){
      write.xlsx(all.df(), f)
    })
  
  # display table of genes in selected area
  output$clickedElementsAll <- renderDataTable({
    df <- data.frame("Symbol"=unlist(input$upsetAll_click$elems), "Link" = createLink(unlist(input$upsetAll_click$elems)))
    if(length(unlist(input$upsetAll_click$name)) > 0){
      colnames(df)[1] <- gsub("&"," & ",unlist(input$upsetAll_click$name))
    }
    return(df)
  }, escape = F)
  
  # plot UpSet plot with selected viruses and conditions
  output$upsetAll <- renderUpsetjs({
    if(is.null(input$select) | is.null(input$cond)){
      return(NULL)
    }else{
      cond <- ifelse(grepl("Virus",input$cond), sub("Virus", "(?!Mock).*",input$cond), input$cond)
      id.list <- sapply(input$select,function(virus){
        unique(unlist(lapply(datasetInput[grep(paste(virus, cond, sep = ".*", collapse = "|"), names(datasetInput), perl = T)], function(x){
          return(x[which(abs(x$log2FoldChange) > input$LFCall & x$padj < 0.05 & apply(x[,grep("normalized", colnames(x))],1,max) >= 10),"SYMBOL"])
        })))
      }, USE.NAMES = T, simplify = F)
      id.list <- id.list[sapply(id.list, length) > 0]
      upsetjs() %>% fromList(id.list) %>% generateDistinctIntersections(limit = input$limit) %>%  chartLayout(padding = 40, bar.padding = 0.2) %>%
      chartFontSizes(font.family = "Helvetica", chart.label = "18px", set.label = "16px", bar.label = "16px", axis.tick = "16px", export.label = "13px") %>% 
      interactiveChart()  %>% chartProps(exportButtons=list(share=FALSE, vega=FALSE, dump=FALSE))
    }
  })
  
  # download heatmap
  output$downHeatAll <- downloadHandler(
    filename = function(){
      return(paste0(paste(input$select, collapse = "_"), ".png"))},
    content = function(f){
      png(filename = f)
      print(plot_heat_all())
      dev.off()
    }
  )
  
  observeEvent(input$help5,{
    showModal(modalDialog(
      title = tags$div("Help for ", tags$b("Heatmap Virus Comparison")),
      size = "l",
      tags$div("This tab shows the corresponding heatmap after an intersection is selected in the tab 'Virus Comparison'."),
      img(src="viruses_heatmap.png"),
      tags$div("The modebar in the upper right corner of the plot shows the following functions, that allow the user to interact with the plot: "),
      HTML("<br><br>"),
      img(src="heatmap_modebar.png"),
      HTML("<br><br>"),
      tags$div("The user is able to zoom into the heatmap by dragging a box over the area of interest. To get information about the underlying data mouseover a specific cell."),
      HTML("<br><br>"),
      img(src="viruses_heatmap_zoom.png"),
      img(src="viruses_heatmap_zoom2.png"),
      easyClose = TRUE
    ))
  })
  
  # plot heatmap of genes in selected area
  output$heatmapAll <- renderPlotly({
    plot_heat_all()
  })
  
  plot_heat_all <- reactive({
    if(is.null(input$select) | is.null(input$cond) | is.null(input$upsetAll_click$elems)){
      return(NULL)
    }else{
      data <- all.df()
      data.rows <- data$SYMBOL
      data <- data[,-c(1,ncol(data))]
      data <- data[,c(TRUE, FALSE)]
      colnames(data) <- sub(".LFC","",colnames(data))
      data <- apply(data,2,as.numeric)
      rownames(data) <- data.rows
      hover <- Reduce(function(x,y)merge(x,y,by="SYMBOL",all=T), sapply(colnames(data), function(s){
        s.gr <- sub("_.*","",strsplit(s,"_vs_")[[1]])
        dataset <- datasetInput[[s]]
        dataset[,c(5:ncol(dataset))] <- signif(dataset[,c(5:ncol(dataset))], 4) 
        dataset <- dataset[dataset$SYMBOL %in% rownames(data), ]
        dataset$c1 <- apply(dataset[ , grep(paste0("normalized.*", s.gr[1]), colnames(dataset)), drop = F] , 1 , paste , collapse = "; " )
        dataset$c2 <- apply(dataset[ , grep(paste0("normalized.*", s.gr[2]), colnames(dataset))] , 1 , paste , collapse = "; " )
        df.tmp <- data.frame(dataset[,"SYMBOL"], sprintf("Gene: %s<br />Sample: %s<br />LFC: %s<br />padj: %s<br />%s: %s<br />%s: %s", 
                                                         dataset[,"SYMBOL"], s, dataset[,"log2FoldChange"], dataset[,"padj"],
                                                         s.gr[1], dataset[,"c1"], s.gr[2], dataset[,"c2"]))
        colnames(df.tmp) <- c("SYMBOL", s)
        return(df.tmp)
      }, simplify = F)
      )
      rownames(hover) <- hover$SYMBOL
      hover <- subset(hover, select=-c(SYMBOL))
      hover <- hover[order(rownames(hover)),]
      data <- data[order(rownames(data)),]
      plotHeatmap(data, colClust = F, rowClust = T, border_col = NA, hover = hover, setWidth = T, fontsize_r = 12, fontsize_c = 12)
    }
  })
  
  ## SNP analysis
  
  # Help text
  observeEvent(input$help6,{
    showModal(modalDialog(
      title = tags$div("Help for ", tags$b("SNP Analysis")),
      size = "l",
      tags$div("This tab is meant to explore the results of a variant analysis."),
      HTML("<br>"),
      img(src="snp_sidebar.png"),
      HTML("<br><br>"),
      tags$div("The modebar in the upper right corner of the plot shows the following functions, that allow the user to interact with the plot: "),
      HTML("<br>"),
      img(src="snp_modebar.png"),
      HTML("<br><br>"),
      tags$div("The plot shows the frequency of single nucleotide polymorphisms (SNP) and insertions and deletions (INDEL) of a chosen virus over time. The user can get further information about a SNP/INDEL by hovering over a cell. The tooltip includes information about the SNP position on the genome, the nucleotide in the reference genome and the alternative nucleotide(s) as well as the frequency of a SNP and the read depth at the SNP position."),
      img(src="snp_mouseover.png"),
      easyClose = TRUE
    ))
  })
  
  # add vertical line to plot
  vline <- function(x = 0, color = "black") {
    list(
      type = "line", 
      y0 = 0, 
      y1 = 1, 
      yref = "paper",
      x0 = x, 
      x1 = x, 
      line = list(color = color)
    )
  }
  
  # plot SNPs for selected  virus as heatmap
  output$heatSNP <- renderPlotly({
    v <- input$virusSNP
    v.df <- Reduce(rbind, sapply(names(vcf.list[grep(v, names(vcf.list))]), function(n){
      x <- vcf.list[[n]]
      x.df <- as.data.frame(x@fix, stringsAsFactors = F)[, c("CHROM", "POS", "REF", "ALT")]
      if(nrow(x.df)>0){
        x.df$AF <- extract.info(x, "AF", as.numeric = T)
        x.df$DP <- extract.info(x, "DP", as.numeric = F)
        x.df$Sample <- n
        colnames(x.df) <- c("Segment", "POS", "REF", "ALT", "AF", "DP", "Sample")
        x.df$snp <- paste(x.df$POS, x.df$ALT, sep = "_")
        return(x.df)
      }
    }, simplify = F))
    
    x_break <- which(sub("_[1|2]$","",mixedsort(unique(v.df$Sample)))[-1] != sub("_[1|2]$","",mixedsort(unique(v.df$Sample)))[-length(sub("_[1|2]$","",mixedsort(unique(v.df$Sample))))])
    v.sum <- v.df %>% group_by(Segment) %>% summarise(n_distinct(snp))
    height <- ifelse(sum(v.sum[,2])>30, sum(v.sum[,2])*30, sum(v.sum[,2])*50)
    palette <- colorRampPalette(c("white","yellow", "orange", "red", "darkred"))
    
    # create plotly heatmap for each segment
    p <- lapply(mixedsort(unique(v.df$Segment)), function(x){
      v.df.mod <- v.df[v.df$Segment==x,]
      v.df.mod$snp <- factor(v.df.mod$snp, levels = mixedsort(unique(v.df.mod$snp)))
      v.df.mod$Sample <- factor(v.df.mod$Sample, levels = mixedsort(unique(v.df$Sample)))
      
      v.df.split <- split(v.df.mod[,c("snp","AF")], v.df.mod$Sample)
      v.df.mat <- Reduce(function(x,y)merge(x,y,by="snp",all=T,no.dups=F),v.df.split)
      colnames(v.df.mat) <- c("snp",names(v.df.split))
      rownames(v.df.mat) <- v.df.mat$snp
      v.df.mat <- as.matrix(v.df.mat[,-1])
      v.df.mat[is.na(v.df.mat)] <- 0
      
      
      #create matrix with hoverinfo per sample and snp
      hover <- matrix(ncol = length(levels(v.df.mod$Sample)), nrow = length(unique(v.df.mod$snp)))
      colnames(hover) <- levels(v.df.mod$Sample)
      rownames(hover) <- mixedsort(unique(v.df.mod$snp))
      for(i in unique(v.df.mod$snp)){
        for(s in unique(v.df.mod$Sample)){
          if(nrow(v.df.mod[v.df.mod$Sample==s & v.df.mod$snp==i, ]) > 0){
            hover[i,s] <- sprintf("Sample: %s<br />Position: %s<br />Reference: %s<br />Alternative: %s<br />Frequency: %s<br />Depth: %s", 
                                  s, v.df.mod[v.df.mod$Sample==s & v.df.mod$snp==i, "POS"], v.df.mod[v.df.mod$Sample==s & v.df.mod$snp==i, "REF"], v.df.mod[v.df.mod$Sample==s & v.df.mod$snp==i, "ALT"], 
                                  v.df.mod[v.df.mod$Sample==s & v.df.mod$snp==i, "AF"], v.df.mod[v.df.mod$Sample==s & v.df.mod$snp==i, "DP"])
          }else{
            hover[i,s] <- NA
          }
        }
      }
      
      plot_ly(x=colnames(v.df.mat), y=rownames(v.df.mat), z = v.df.mat, type = "heatmap", colors = "Greys", showlegend = F, zmin = 0, zmax = 1, zauto = F,
              text = hover, hoverinfo = "text") %>%  #hovertemplate = "x : %{x}\ny : %{y}\nDepth : %{customdata}<extra></extra>") %>% 
        layout(shapes=lapply(x_break-0.5, vline)) %>% config(displayModeBar = TRUE) %>% 
        add_annotations(
          text = x,
          x = 0.5,
          y = 1,
          yref = "paper",
          xref = "paper",
          xanchor = "left",
          yanchor = "bottom",
          showarrow = FALSE,
          font = list(size = 15)
        )
    })
    subplot(p, shareX = T, nrows = length(unique(v.df$Segment)))
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
