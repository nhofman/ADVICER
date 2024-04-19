# ADVICER - Analysis Dashboard for Virus-Induced Cell Response based on RNA-seq data
# 
# Shiny application for the comparative visualization of virus-induced differential gene expression
#
# 

library(shiny)
library(DT)
library(VennDiagram)
library(UpSetR)
library(upsetjs)
library(ggplot2)
library(plotly)
library(gtools)
library(openxlsx)
library(vcfR)
library(tidytext)
library(dplyr)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(GetoptLong)
library(S4Vectors)

# Create interactive heatmap with InteractiveComplexHeatmap
plotHeatmap <- function(x, row_subset = NA, distMethod = "euclidean", clusterMethod = "complete", clrn = 1, clcn = 1, setWidth = F, 
                        rowClust = T, colClust = T, fontsize_r = 0.8, fontsize_c = 10, annCol = NA, annRow = NA, border_col = "grey60", plot.fig = T, 
                        legend.cut = 1, filter_col = NA, annotation_colors = NA, break_step = 0.1, display_numbers = F, hover = NULL, file = NA, 
                        legend.limit.up = NA, legend.limit.down = NA){
  if(is.na(row_subset[1])){
    xx <- x
  }else{
    xx <- x[row_subset, , drop = F]
  }
  #xx[xx==0] <- 0.000001
  #xx <- na.omit(xx)
  if (nrow(xx)<2) {
    rowClust <- F
  }
  if(!is.na(filter_col)){
    xx <- xx[, which(colMaxs(as.matrix(abs(xx)))>filter_col)]
  }
  xx.unlist <- as.numeric(unlist(xx))
  
  # Calculate legend limits, if not defined
  if(is.na(legend.limit.up)){
    legend.limit.up <- quantile(unlist(xx), na.rm = TRUE, probs = legend.cut)
  }
  if(is.na(legend.limit.down)){
    legend.limit.down <- quantile(unlist(xx), na.rm = TRUE, probs = (1-legend.cut))  
  }
  
  color <- vector()
  legend.limit <- max(legend.limit.up, abs(legend.limit.down))
  legend.limit.round <- ceiling(legend.limit)
  
  # Color function
  col_fun <- circlize::colorRamp2(c(-legend.limit, 0, legend.limit), c("blue", "white", "red"), space = "LAB")
  
  p <- Heatmap(as.matrix(xx), cluster_columns = colClust, cluster_rows = rowClust, col = col_fun, show_row_dend = F, 
               heatmap_legend_param = list(title = "", at = c(-legend.limit.round, -legend.limit.round/2, 0, legend.limit.round/2, legend.limit.round)), 
               row_names_gp = gpar(family = "Helvetica"), column_names_gp = gpar(family = "Helvetica"))
  return(p)
}

# Create Link to NCBI
createLink <- function(val) {
  sprintf('<a href="https://www.ncbi.nlm.nih.gov/gene?term=(%s%%5BGene%%20Name%%5D)%%20AND%%20human%%5BOrganism%%5D" target="_blank">Info</a>', val)
}

# plot expression course of specific gene
plotExpression <- function(expr.df, padj.df, gene){
  lfc.gene <- data.frame(t(expr.df[gene, ]))
  lfc.gene <- merge(lfc.gene, t(padj.df[gene, ]), by = "row.names")
  rownames(lfc.gene) <- lfc.gene[, 1]
  lfc.gene <- lfc.gene[, -1]
  colnames(lfc.gene) <- c("lfc", "padj")
  #lfc.gene$padj[is.na(lfc.gene$padj)] <- 10
  #lfc.gene <- lfc.gene[complete.cases(lfc.gene), ]
  lfc.gene$sig <- ifelse(lfc.gene$padj<0.0001, "***", ifelse(lfc.gene$padj<0.001, "**", ifelse(lfc.gene$padj<0.05, "*", "")))
  lfc.gene$group <- sub("_.*", "", rownames(lfc.gene))
  lfc.gene$time <- sub(".*_(.*_.*)", "\\1", rownames(lfc.gene)) 
  lfc.gene$time <- ifelse(lfc.gene$time=="Mock_BPL", "Virus BPL:Mock BPL", ifelse(grepl("BPL", lfc.gene$time), "Virus 24h:Virus BPL", sub("Mock_", "", lfc.gene$time)))
  lfc.min <- ifelse(min(lfc.gene$lfc, na.rm = T)>0, 0, min(lfc.gene$lfc, na.rm = T)-0.25)
  lfc.max <- max(lfc.gene$lfc, na.rm = T) + 0.25
  p1 <- ggplot(lfc.gene[grepl("h$", lfc.gene$time), ], aes(x=factor(time, levels = mixedsort(unique(time))), y=lfc, group=group, color=group)) + 
    geom_line(na.rm = T) + geom_point(na.rm = T) +
    labs(title = paste("Expression profile of", gene), y = "Log2FoldChange", x = "Time", color = "Virus", caption = "*: padj < 0.05    **: padj < 0.001    ***: padj < 0.0001") + 
    geom_text(aes(label=sig), nudge_y = 0.1, show.legend = FALSE, na.rm = T) + ylim(c(lfc.min, lfc.max)) + scale_color_manual(values = color) +
    theme(plot.caption = element_text(hjust = 0, size = 15, family = "Helvetica"), legend.text = element_text(size = 15, family = "Helvetica"), 
          legend.title = element_text(size = 15, face = "bold", family = "Helvetica"), axis.text = element_text(size = 15), 
          axis.title = element_text(size = 15, face = "bold", family = "Helvetica"), plot.title = element_text(size = 20, family = "Helvetica", face = "bold")) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
  p2 <-  ggplot(lfc.gene[grepl("BPL", lfc.gene$time), ], aes(x=time, y=lfc, group=group, fill=group)) + 
    geom_col(position = "dodge", na.rm = T) +
    labs(title = paste("Expression profile of", gene), y = "Log2FoldChange", x = "", fill = "Virus", caption = "*: padj < 0.05    **: padj < 0.001    ***: padj < 0.0001") + 
    geom_text(aes(label=sig), position=position_dodge(width=0.9), vjust=-0.25, show.legend = F, na.rm = T) + ylim(c(lfc.min, lfc.max)) + scale_fill_manual(values = color) +
    theme(plot.caption = element_text(hjust = 0, size = 15, family = "Helvetica"), legend.text = element_text(size = 15, family = "Helvetica"), 
          legend.title = element_text(size = 15, face = "bold", family = "Helvetica"), axis.text = element_text(size = 15), 
          axis.title = element_text(size = 15, face = "bold", family = "Helvetica"), plot.title = element_text(size = 20, family = "Helvetica", face = "bold")) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
  p <- gridExtra::arrangeGrob(p1, p2, ncol = 1, heights = c(10, 10))
  #dev.off()
  return(p)
}

# read Renviron file and extract name of data folder
readRenviron(".Renviron")
if(Sys.getenv("DATADIR") != ""){
  shinyOptions(filedir = Sys.getenv("DATADIR"))
}
Sys.setenv("PATH" = paste(Sys.getenv("PATH"), "/root/miniconda3/bin", sep = .Platform$path.sep))

datasetInput <- NULL
# Define UI for application
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
            #geneNotFound{
              font-size: 20px;
              font-style: bold;
           }
           .butt{
              border: 2px solid black;
              font-size: 16px;
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
  titlePanel(div(HTML("<b>ADVICER</b> - <b>A</b>nalysis <b>D</b>ashboard for <b>V</b>irus-<b>I</b>nduced <b>CE</b>ll <b>R</b>esponse based on RNA-Seq data")), windowTitle = "ADVICER"), 
  # Define tabs
  navbarPage(
    "", 
    tabPanel("Download Data", 
             sidebarPanel(
               div(style = "text-align:right", actionButton("help0", label = div(strong("Help"), icon("question")))), 
               htmlOutput("download")
             ), 
             mainPanel(
               downloadButton("filesDownXlsx", HTML("Download selected file(s) in Excel format <b>(xlsx)</b>"), class = "butt"),
               downloadButton("filesDown", HTML("Download selected file(s) in comma separated values format <b>(csv)</b>"), class = "butt")
             )
    ), 
    tabPanel("Volcano Plot / MA Plot", 
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
               br(),
               downloadButton("downPlotSelect_xlsx", "Download table as xlsx", class = "butt"),
               downloadButton("downPlotSelect_csv", "Download table as csv")
             ), 
             mainPanel(
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
               br(), 
               actionButton("addHeat2", "Generate heatmap", class = "butt"), 
               br(), 
               br(), 
               h5(strong("Download table")), 
               downloadButton("down_vXLSX", "Download as xlsx", class = "butt"), 
               #br(), 
               downloadButton("down_vCSV", "Download as csv")
             ), 
             mainPanel(
               conditionalPanel(
                 condition = 'input.plotType == "UpSet"', 
                 upsetjsOutput("upset"), #, click = "upsetClick"), 
                 br(), 
                 br(), 
                 DT::dataTableOutput("clickedElements"), 
                 br(), 
                 br()
               ), 
               conditionalPanel(
                 condition = 'input.plotType == "Venn"', 
                 upsetjsOutput("upsetVenn"), #, click = "upsetClick"), 
                 br(), 
                 textOutput("headerVenn"), 
                 br(), 
                 DT::dataTableOutput("clickedElementsVenn"), 
                 br(), 
                 br()
               ), 
               htmlOutput("heatmapTimeComplex")
             )
    ), 
    tabPanel("Gene Expression", 
             sidebarPanel(
               div(style = "text-align:right", actionButton("help3", label = div(strong("Help"), icon("question")))), 
               checkboxGroupInput("selectVirus", label = h5(strong("Select viruses")), 
                                  choices = list("H1N1", "H5N1", "RVFV", "SFSV", "RSV", "NiV", "EBOV", "MARV", "LASV"), 
                                  selected = NULL), 
               textInput("gene", label = h5(strong("Gene symbol")), value = "", placeholder = "e.g. TNF or cxcl2"), 
               h5(strong("Download plot")), 
               downloadButton("downGeneXpng", "Download as png"), 
               downloadButton("downGeneXsvg", "Download as svg")
             ), 
             mainPanel(
               plotOutput("geneX", height = 1000)
             )), 
    tabPanel("Virus Comparison", 
             sidebarPanel(
               div(style = "text-align:right", actionButton("help4", label = div(strong("Help"), icon("question")))), 
               checkboxGroupInput("select", label = h5(strong("Select viruses")), 
                                  choices = list("H1N1", "H5N1", "RVFV", "SFSV", "RSV", "NiV", "EBOV", "MARV", "LASV"), 
                                  selected = NULL), #c("H1N1", "H5N1", "MERS", "CoV229E", "RVFV", "SFSV", "RSV", "NIV", "EBOV")), 
               htmlOutput("selectCond"), 
               em("(BPL samples were taken after 24h)"), 
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
               br(), 
               actionButton("addHeatAll", "Generate heatmap", class = "butt"), 
               br(), 
               br(), 
               h5(strong("Download table (includes LFC and padj)")), 
               downloadButton("downallXLSX", "Download as xlsx", class = "butt"), 
               #br(), 
               downloadButton("downallCSV", "Download as csv")
             ), 
             mainPanel(
               upsetjsOutput("upsetAll"), #, click = "upsetClick"), 
               br(), 
               DT::dataTableOutput("clickedElementsAll"), 
               br(), 
               htmlOutput("heatmapAll")
             )
    ), 
    tabPanel("SNP Analysis", 
             sidebarPanel(
               div(style = "text-align:right", actionButton("help6", label = div(strong("Help"), icon("question")))), 
               selectInput("virusSNP", label = h5(strong("Select virus")), 
                           choices = list("H1N1", "H5N1", "RVFV", "SFSV", "RSV", "NiV", "EBOV", "MARV", "LASV"), 
                           selected = NULL)
             ), 
             mainPanel(
               plotlyOutput("heatSNP", height = "1000px")
             ))
  )
)

# Define directory containing data
filedir <- getShinyOption("filedir")
files.list <- list.files(filedir, pattern = ".*.csv", full.names = T)
files.list.xlsx <- list.files(filedir, pattern = ".*.xlsx", full.names = T)
vcf.files <- list.files(filedir, pattern = ".*[1|2].vcf", full.names = T)

# Define font family
#family = "Helvetica"

# read colors for viruses
color.df <- read.table("virus_colors.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
color <- color.df[, 2]
names(color) <- color.df[, 1]

# Define server logic 
server = function(input, output, session) {
  
  session$onSessionEnded(function(){
    stopApp()
  })
  
  # Define reactive values
  data <- reactiveValues(select_points=NULL, heat=FALSE, heat_data = NULL, heat_hover = NULL, dt=NULL, snp=NULL, all.df=NULL, all.dt = 1, 
                         virus.df=NULL, plotX=NULL, heatPlty=FALSE, heatCplx=FALSE, heatmap_id=0, heatAll=FALSE, heatmap_id_all=0)
  
  # Read data from files.list
  withProgress(message = 'PROCESSING DATA...', detail = "This may take a while...", value = 0, {
    datasetInput <- lapply(files.list, function(x){
      f <- read.csv(x, header = TRUE, stringsAsFactors = F, row.names = 1)
      incProgress(1/length(files.list))
      return(f)
    })
    #Sys.sleep(5)
    names(datasetInput) <- sub("deseq2_results_(.*).csv", "\\1", basename(files.list))
    datasetInput <- datasetInput[mixedorder(sub("vs_Mock_", "", names(datasetInput)))]
    vcf.list <- lapply(vcf.files, read.vcfR)
    names(vcf.list) <- sub(".vcf", "", basename(vcf.files))
  })
  
  # Help text
  observeEvent(input$help0, {
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
  
  # List all uploaded files
  output$uploadText <- renderText({
    if(!is.null(datasetInput)){
      print("Uploaded data:")
    }
  })
  
  # List all files in list
  output$table <- renderTable(
    data.frame(names(datasetInput))
    , colnames = F, rownames = F)
  
  output$download <- renderUI({
    checkboxGroupInput("filesToDown", "Choose files to download", as.list(c("all", names(datasetInput))))
  })
  
  # Download the selected files
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
          fileName <- files.list[grep(i, files.list)] 
          files <- c(files, fileName)
        }
      }
      zip(f, files)
    }, 
    contentType = "csv"
  )
  
  output$filesDownXlsx <- downloadHandler(
    filename = function(){
      "data_xlsx.zip"
    }, 
    content = function(f){
      files <- NULL
      if("all" %in% input$filesToDown){
        files <- files.list.xlsx
      }else{
        for (i in input$filesToDown){
          fileName <- files.list.xlsx[grep(i, files.list.xlsx)] 
          files <- c(files, fileName)
        }
      }
      zip(f, files)
    }, 
    #contentType = "csv"
  )
  
  ## Volcano Plot / MA Plot
  
  # Help text
  observeEvent(input$help1, {
    showModal(modalDialog(
      title = tags$div("Help for ", tags$b("Volcano Plot / MA Plot")), 
      size = "l", 
      tags$div("Volcano plots and MA-plots for uploaded files. Each file contains the results of the differential gene expression analysis with DEseq2.
      The sidebar shows the available options to adapt the plots:"), 
      img(src="volcano_sidebar.jpg"), 
      tags$div("The Volcano plot shows the log2FoldChange (LFC) and the negative log10 of the adjusted p-value (padj) for the chosen comparison. Significantly differentially expressed genes are colored blue (down-regulated) or red (up-regulated). A gene is considered differentially expressed if: LFC > ", tags$i("LFC cutoff") , "and padj < 0.05 and a normalized gene count of at least one sample ≥ 10."), 
      HTML("<br><br>"), 
      img(src="volcano_plot.png"), 
      HTML("<br><br>"), 
      tags$div("The modebar in the upper right corner of the plot shows:"), 
      HTML("<br><br>"), 
      img(src="volcano_modebar.jpg"), 
      HTML("<br><br>"), 
      tags$div("Further information for one or multiple points, like gene symbol, LFC and padj, can be shown by mouse over."), 
      tags$div("The selected genes are shown in tabular form along with their LFCs and padjs. A link to the NCBI database is provided for each gene. The link opens in a new tab outside the application. The table can be searched and all columns can be sorted in ascending and descending order. It can be downloaded in xlsx or csv file format."), 
      img(src="volcano_plot_select.jpg"), 
      easyClose = TRUE
    ))
  })
  
  # render file selection
  output$fileSelect <- renderUI({
    selectInput("fileToPlot", label = "Choose file", choices = names(datasetInput))
  })
  
  # download results        
  output$downPlotSelect_csv <- downloadHandler(
    filename = function(){
      return(paste0(input$fileToPlot, "_selected", ".csv"))}, 
    content = function(f){
      df <- select.df()
      df$LinkToNCBI <- sub(".+\\\"(.+?)\\\".+", "\\1", df$LinkToNCBI)
      write.csv(df, f, row.names = F)
    }, 
    contentType = "csv")
  
  output$downPlotSelect_xlsx <- downloadHandler(
    filename = function(){
      return(paste0(input$fileToPlot, "_selected", ".xlsx"))}, 
    content = function(f){
      df <- select.df()
      df$LinkToNCBI <- sub(".+\\\"(.+?)\\\".+", "\\1", df$LinkToNCBI)
      write.xlsx(df, f)
    }, 
    contentType = NULL)
  
  # set plotly_click to NULL or 1 on input change
  observeEvent(input$fileToPlot, {
    data$select_points <- NULL
  })
  
  observeEvent(input$plotTypeSingle, {
    data$select_points <- NULL
  })
  
  observeEvent(event_data("plotly_selected"), {
    data$select_points <- 1
  })
  
  # create table of selected elements
  select.df <- reactive({
    selected <- event_data("plotly_selected")
    df <- datasetInput[[input$fileToPlot]]
    df <- df[df$SYMBOL %in% selected$customdata, c("SYMBOL", "log2FoldChange", "padj")]
    df[, c(2:ncol(df))] <- signif(df[, c(2:ncol(df))], 4)
    df$LinkToNCBI <- createLink(df$SYMBOL)
    return(df)
  })
  
  # display table of clicked elements and associated LFC and padj in Volcano plot
  output$clickVP <- renderDataTable({
    req(!is.null(data$select_points))
    select.df()
  }, escape = F, rownames = F)
  
  # plot Volcano plot of selected sample
  output$volcanoPlot <- renderPlotly({
    p <- plotVolcano()
    p %>% config(toImageButtonOptions=list(format="svg", width=800, height=600, scale=2))
  })
  
  plotVolcano <- reactive({
    if(is.null(datasetInput) | is.null(input$fileToPlot))
      return(NULL)
    file.data <- datasetInput[[input$fileToPlot]]
    file.data$col <- ifelse(file.data$log2FoldChange > input$LFCsingle & file.data$padj < 0.05 & apply(file.data[, grep("normalized", colnames(file.data))], 1, max) >= 10, "up", 
                            ifelse(file.data$log2FoldChange < -(input$LFCsingle) & file.data$padj < 0.05 & apply(file.data[, grep("normalized", colnames(file.data))], 1, max) >= 10, "down", "not significant"))
    p <- plot_ly(file.data, type = "scatter", x = ~log2FoldChange, y = ~-log10(padj), color = ~col, colors = c("up"="red", "down"="steelblue3", "not significant"="black"), 
                 mode = "markers", marker = list(size = 5), customdata = ~SYMBOL, 
                 text = ~paste("Gene: ", SYMBOL, '<br>LFC: ', signif(log2FoldChange, 4), '<br>padj: ', signif(padj, 4)))#, 
    #source = "V")
    p <- p %>% layout(legend = list(orientation = "h", y = -0.2), xaxis = list(exponentformat = "none"), dragmode = "select") %>% config(displayModeBar = TRUE)  
    #p <- p %>% event_register("plotly_selected")
    return(p)
  })
  
  # display table of clicked elements and associated LFC and padj in MA-plot
  output$clickMA <- renderDataTable({
    req(!is.null(data$select_points))
    select.df()
  }, escape = F, rownames = F)
  
  # plot MA-plot of selected sample
  output$MAPlot <- renderPlotly({
    p <- plotMA()
    p %>% config(toImageButtonOptions=list(format="svg", width=800, height=600, scale=2))
  })
  
  plotMA <- reactive({
    if(is.null(datasetInput) | is.null(input$fileToPlot)){
      return(NULL)
    }else{
      file.data <- datasetInput[[input$fileToPlot]]
      file.data$col <- ifelse(file.data$log2FoldChange > input$LFCsingle & file.data$padj < 0.05 & apply(file.data[, grep("normalized", colnames(file.data))], 1, max) >= 10, "up", 
                              ifelse(file.data$log2FoldChange < -(input$LFCsingle) & file.data$padj < 0.05 & apply(file.data[, grep("normalized", colnames(file.data))], 1, max) >= 10, "down", "not significant"))
      file.data$baseMeanSample <- rowMeans(file.data[, grep("normalized", colnames(file.data))])
      p <- plot_ly(file.data, type = "scatter", x = ~log2(baseMean), y = ~log2FoldChange, color = ~col, colors =c("up"="red", "down"="steelblue3", "not significant"="black"), 
                   mode = "markers", marker = list(size = 5), customdata = ~SYMBOL, 
                   text = ~paste("Gene: ", SYMBOL, '<br>LFC: ', signif(log2FoldChange, 4), '<br>padj: ', signif(padj, 4)))#, 
      #source = "M")
      p <- p %>% layout(legend = list(orientation = "h", y = -0.2), yaxis = list(exponentformat = "none"), dragmode = "select") %>% config(displayModeBar = TRUE)
      return(p)
    }
  })
  
  ## Compare time points of defined virus
  
  # Help text 
  observeEvent(input$help2, {
    showModal(modalDialog(
      title = tags$div("Help for ", tags$b("Compare Time Points")), 
      size = "l", 
      tags$div("This tab is meant to compare the differential gene expression between different time points of infection. Each gene list contains all genes differentially expressed with |LFC| > ", tags$i("LFC cutoff"), " and padj < 0.05 and a normalized gene count of at least one sample ≥ 10. Intersections can either be shown as Venn diagram or UpSet plot ", tags$a("(info)", href="https://jku-vds-lab.at/tools/upset/", target="_blank"), ". The sidebar shows the available options to adapt the plots:"), 
      HTML("<br>"), 
      img(src="time_sidebar.jpg"), 
      HTML("<br><br>"), 
      "The Venn diagram shows the intersections of genes from the selected gene lists. For reasons of clarity, a maximum number of 5 lists can be displayed.", 
      img(src="time_venn.jpg"), 
      tags$div("Clicking an intersection will display underlying genes in a sortable table along with log2FoldChange (LFC), adjusted p-value (padj) and a link to the NCBI database. The table can be downloaded as xlsx or csv file."), 
      #HTML("<br><br>"), 
      img(src="time_venn_select.jpg"), 
      HTML("<br><br>"), 
      tags$div("The UpSet plot shows the intersections of genes in a matrix-like layout. This approach allows to compare a larger number of sets. Each bar shows an intersection of genes with the corresponding number. In the respective column below you can see the samples involved, marked with black dots. For example in the plot below 263 genes are differentially expressed (see bar) after 6h and 12h (see black dots in column below)."), 
      tags$div("Genes in an intersection can be listed in a sortable table along with LFC, padj and a link to the NCBI database by clicking on the appropriate bar. The table can be downloaded as xlsx or csv file."), 
      #img(src="time_upset.jpg"), 
      img(src="time_upset_select.jpg"), 
      HTML("<br><br>"), 
      #HTML("<br><br>"), 
      tags$div("The selected genes can also be shown in a heatmap by clicking on the button in the sidebar. To get information about the underlying data click on a specific cell."), 
      img(src="time_heatmap_select.jpg"), 
      HTML("<br><br>"), 
      tags$div("The modebar beneath the heatmap shows the following functions, that allow the user to interact with and download the plot:"), 
      #HTML("<br><br>"), 
      img(src="heatmap_modebar.jpg"), 
      #HTML("<br><br>"), 
      tags$div("The user is able to zoom into the heatmap by dragging a box over the area of interest. The new heatmap is drawn in the box 'Selected sub-heatmap' next to the original heatmap."), 
      HTML("<br>"), 
      img(src="time_heatmap_zoom.jpg"), 
      HTML("<br><br>"), 
      tags$div("The user can also interact with this sub-heatmap using the buttons beneath the box. There the user can configure or download the plot as well as export it as table. To zoom further into the current sub-heatmap, the user can make it interactive by clicking the corresponding button in 'Configure sub-heatmap' and zoom into the heatmap as previously described."), 
      HTML("<br>"), 
      img(src="time_heatmap_zoom2.jpg"), 
      easyClose = TRUE
    ))
  })
  
  # download genes of selected plot area as csv or xlsx
  output$down_vCSV <- downloadHandler(
    filename = function(){
      if(input$plotType == "UpSet"){
        return(paste0(input$select_v, "_", sub("&", "_", input$upset_click$name), ".csv"))
      }else{
        return(paste0(input$select_v, "_", sub("&", "_", input$upsetVenn_click$name), ".csv"))
      }
    }, 
    content = function(f){
      df <- virus.df()
      df$LinkToNCBI <- sub(".+\\\"(.+?)\\\".+", "\\1", df$LinkToNCBI)
      write.csv(df, f, row.names = F)
    }, 
    contentType = "csv")
  
  output$down_vXLSX <- downloadHandler(
    filename = function(){
      if(input$plotType == "UpSet"){
        return(paste0(input$select_v, "_", sub("&", "_", input$upset_click$name), ".xlsx"))
      }else{
        return(paste0(input$select_v, "_", sub("&", "_", input$upsetVenn_click$name), ".xlsx"))
      }
    }, 
    content = function(f){
      df <- virus.df()
      df$LinkToNCBI <- sub(".+\\\"(.+?)\\\".+", "\\1", df$LinkToNCBI)
      write.xlsx(df, f)
    })
  
  # select times to plot
  output$selectTime <- renderUI({
    checkboxGroupInput("time", label = "Choose time points (for Venn: max. 5)", 
                       choiceNames = unique(sapply(names(datasetInput)[grep(input$select_v, names(datasetInput))], 
                                                   function(x){
                                                     s <- sub("_vs_", " : ", x)
                                                     return(s)
                                                   }, USE.NAMES = F)), 
                       choiceValues = sub(".*vs", "vs", names(datasetInput)[grep(input$select_v, names(datasetInput))]))
  })
  
  # set heatmap to NULL at selection of new virus or plot type or time point
  observe({
    input$select_v
    input$plotType
    input$time
    data$heatCplx <- FALSE
    data$dt <- NULL
  })
  
  # clear data table and heatmaps if intersection changes
  observe({
    input$upsetVenn_click$elems
    input$upset_click$elems
    #data$heatCplx <- FALSE
    data$dt <- 1
  })
  
  # set heatmap to TRUE if button is clicked 
  observe({
    if(input$addHeat2){
      data$heatCplx <- T
    }
  })
  
  # if heatCplx TRUE: plot heatmap with InteractiveCompelxHeatmap
  # if heatCplx FALSE: remove heatmap
  observe({
    if(data$heatCplx){
      data$heatmap_id <- sample(1:1000, 1)
      df <- data_heatmap()
      if(!is.null(df)){
        ht <- plotHeatmap(df, colClust = F, rowClust = T)
        ht <- draw(ht)
        dev.off()
        InteractiveComplexHeatmapWidget(input, output, session, ht, output_id = "heatmapTimeComplex", layout = "1-(2|3)", compact = FALSE, 
                                        heatmap_id = paste0("ht", data$heatmap_id), output_ui_float = F, output_ui = htmlOutput("info"), 
                                        click_action = click_action, close_button = T, width1 = 600, height1 = 750, title1 = "Heatmap of genes shown in table", 
                                        height2 = 500)
      }
    }else{
      removeUI(paste0("#ht", data$heatmap_id, "_heatmap_widget *"), multiple = T)
    }
  }) #%>% bindEvent(input$addHeat2, ignoreInit = T)
  
  click_action <- function(df, output) {
    output$info <- renderUI({
      if(is.null(df)) { 
        DataFrame()
      } else {
        data.df <- data_heatmap()
        gene <- rownames(data.df[df$row_index, , drop=F])
        sample <- colnames(data.df[, df$column_index, drop=F])
        sample.parts <- strsplit(sample, "_vs_")[[1]]
        dataset <- datasetInput[[sample]]
        datainfo <- dataset[dataset$SYMBOL==gene, , drop=F]
        datainfo[, c(5:ncol(datainfo))] <- signif(datainfo[, c(5:ncol(datainfo))], 4) 
        datainfo$GENENAME <- ifelse(is.na(datainfo$GENENAME), "", datainfo$GENENAME)
        virus.norm <- apply(datainfo[, grep(sample.parts[1], colnames(datainfo)), drop=F], 1, paste, collapse="; ")
        mock.norm <- apply(datainfo[, grep(sub("_", ".*", sample.parts[2]), colnames(datainfo)), drop=F], 1, paste, collapse="; ")
        HTML(qq("<p style='background-color:white;color:black;border:1px solid black;padding:5px;font-size:14px;'>Sample: @{sample}<br />Gene: @{gene}<br />
                Genename: @{datainfo$GENENAME}<br />LFC: @{datainfo$log2FoldChange}<br />padj: @{datainfo$padj}<br />
                Normalized counts @{sample.parts[1]}: @{virus.norm}<br />Normalized counts @{sample.parts[2]}: @{mock.norm} </p>"))
      }
    })
  }
  
  # create dataframe for heatmap
  data_heatmap <- reactive({ 
    data.df <- virus.df()
    if(is.null(data.df)){
      return(NULL)
    }else{
      data.rows <- data.df$SYMBOL
      data.df <- data.df[, -c(1, ncol(data.df)), drop=F]
      if(nrow(data.df)>1){
        data.df <- apply(data.df, 2, as.numeric)
      }else{
        data.df <- as.data.frame(t(apply(data.df, 2, as.numeric)))
      }
      rownames(data.df) <- data.rows
      data.df <- data.df[, c(TRUE, FALSE), drop=F]
      colnames(data.df) <- sub(".LFC", "", colnames(data.df))
      return(data.df)
    }
  }) 
  
  # create dataframe with LFC and padj for selection
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
      data.df <- Reduce(function(x, y)merge(x, y, by="SYMBOL", all = T), lapply(names(datasetInput[grep(paste0(input$select_v, ".*_", input$time, collapse = "|"), names(datasetInput))]), function(x){
        y <- datasetInput[[x]]
        y <- y[y$SYMBOL %in% virus_id, c("SYMBOL", "log2FoldChange", "padj"), drop=F]
        colnames(y) <- c("SYMBOL", paste0(x, ".LFC"), paste0(x, ".padj"))
        return(y)
      })
      )
      data.df <- data.df[, c(1, mixedorder(colnames(data.df)[-1])+1), drop=F]
      data.df$LinkToNCBI <- createLink(data.df$SYMBOL)
      return(data.df)
    }
  })
  
  # display table with elements in clicked plot area of UpSet plot
  output$clickedElements <- renderDataTable({
    dt <- data$dt
    df <- virus.df()
    if(is.null(df) | is.null(dt)){
      return(data.frame("SYMBOL"=character(), "LinkToNCBI"=character()))
    }else{
      df[, c(2:(ncol(df)-1))] <- signif(df[, c(2:(ncol(df)-1))], 4)
      return(df)
    }
  }, escape = F, rownames = F)
  
  # plot UpSet plot of selected virus and times
  output$upset <- renderUpsetjs({
    if(is.null(input$select_v) | is.null(input$time)){
      return(NULL)
    }else{
      id.list <- lapply(datasetInput[grep(paste0(input$select_v, ".*_", input$time, collapse = "|"), names(datasetInput))], function(x){
        return(x[which(abs(x$log2FoldChange) > input$LFC & x$padj < 0.05 & apply(x[, grep("normalized", colnames(x))], 1, max) >= 10), "SYMBOL"])
      })
      names(id.list) <- sub(".*_vs_", "", names(id.list))
      names(id.list) <- sub("Mock_", "", names(id.list))
      upsetjs(sizingPolicy = upsetjsSizingPolicy(padding = 10)) %>% upsetjs::fromList(id.list) %>% generateDistinctIntersections(empty = T) %>%
        chartFontSizes(font.family = "Helvetica", chart.label = "16px", set.label = "14px", bar.label = "14px", axis.tick = "14px", export.label = "13px") %>% 
        chartLayout(padding = 40, bar.padding = 0.2) %>% chartTheme(selection.color = "orange") %>% 
        chartProps(exportButtons=list(share=FALSE, vega=FALSE, dump=FALSE)) %>% interactiveChart("click") 
    }
  })
  
  # display table with elements in clicked plot area of Venn plot
  output$clickedElementsVenn <- renderDataTable({
    dt <- data$dt
    df <- virus.df()
    if(is.null(df) | is.null(dt)){
      return(data.frame("SYMBOL"=character(), "LinkToNCBI"=character()))
    }else{
      df[, c(2:(ncol(df)-1))] <- signif(df[, c(2:(ncol(df)-1))], 4)
      return(df)
    }
  }, escape = F, rownames = F)
  
  # plot Venn diagram of selected virus und times
  output$upsetVenn <- renderUpsetjs({
    if(is.null(input$select_v) | is.null(input$time)){
      return(NULL)
    }else{
      id.list <- lapply(datasetInput[grep(paste0(input$select_v, ".*_", input$time, collapse = "|"), names(datasetInput))], function(x){
        return(x[which(abs(x$log2FoldChange) > input$LFC & x$padj < 0.05 & apply(x[, grep("normalized", colnames(x))], 1, max) >= 10), "SYMBOL"])
      })
      names(id.list) <- sub(".*_vs_", "", names(id.list))
      names(id.list) <- sub("Mock_", "", names(id.list))
      upsetjsVennDiagram() %>% fromList(id.list) %>% 
        chartFontSizes(font.family = "Helvetica", chart.label = "16px", set.label = "18px", bar.label = "14px", axis.tick = "14px", export.label = "13px") %>% 
        chartTheme(selection.color = "orange") %>% interactiveChart("click") %>% chartProps(exportButtons=list(share=FALSE, vega=FALSE, dump=FALSE))
    }
  })
  
  ## Gene expression plot
  
  # Help text
  observeEvent(input$help3, {
    showModal(modalDialog(
      title = tags$div("Help for ", tags$b("Gene Expression")), 
      size = "l", 
      tags$div("This tab is meant to explore the expression profile of genes over time.", 
               "The sidebar shows the available options to adapt the plots:"), 
      img(src="expression_sidebar.jpg"), 
      HTML("<br><br>"), 
      tags$div("The plot shows the expression of the selected gene for the chosen viruses as log2FoldChange (LFC) over time. The asterisks represent the significance of the LFC."), 
      HTML("<br><br>"), 
      img(src="expression_plot_upper.png"), 
      HTML("<br><br>"), 
      tags$div("The plot shows the expression of the selected gene as LFC. The left hand side shows the 24 h time point of the selected viruses compared to the inactivated virus (Virus BPL), whereas the right hand side shows the inactivated virus compared to the uninfected control (Mock BPL)."), 
      HTML("<br><br>"), 
      img(src="expression_plot_lower.png"), 
      easyClose = TRUE
    ))
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
      ggsave(f, data$plotX, "svg", width = 8, height = 10, dpi = 400)
    }
  )
  
  # plot gene expression of selected viruses over time
  output$geneX <- renderPlot({
    if(toupper(input$gene) %in% rownames(data_lfc())){
      data$plotX <- plotExpression(expr.df = data_lfc(), padj.df = data_padj(), gene = toupper(input$gene))
      plot(data$plotX)
    }
  })
  
  # Show notification, if gene not found  
  observe({
    if(input$gene != "" & !toupper(input$gene) %in% rownames(data_lfc())){
      id <- showNotification("Gene not found!", duration = NULL, type = "warning", id = "abc")
    }
    if(toupper(input$gene) %in% rownames(data_lfc())){
      removeNotification("abc")
    }
  })
  
  # create dataframe of LFCs for selected virus(es)
  data_lfc <- reactive({
    if(is.null(input$selectVirus)){
      return(NULL)
    }else{
      lfc.df <- Reduce(function(x, y)merge(x, y, by="SYMBOL", all=T), lapply(names(datasetInput[grep(paste(input$selectVirus, collapse = "|"), names(datasetInput))]), function(x){
        x.df <- data.frame(datasetInput[[x]]$SYMBOL, datasetInput[[x]]$log2FoldChange)
        colnames(x.df) <- c("SYMBOL", x)
        return(x.df)}))
      rownames(lfc.df) <- lfc.df$SYMBOL
      lfc.df <- lfc.df[, -1]
      return(lfc.df)
    }
  })
  
  # create dataframe of padjs for selected virus(es)
  data_padj <- reactive({
    if(is.null(input$selectVirus)){
      return(NULL)
    }else{
      padj.df <- Reduce(function(x, y)merge(x, y, by="SYMBOL", all=T), lapply(names(datasetInput[grep(paste(input$selectVirus, collapse = "|"), names(datasetInput))]), function(x){
        x.df <- data.frame(datasetInput[[x]]$SYMBOL, datasetInput[[x]]$padj)
        colnames(x.df) <- c("SYMBOL", x)
        return(x.df)}))
      rownames(padj.df) <- padj.df$SYMBOL
      padj.df <- padj.df[, -1]
      return(padj.df)
    }
  })
  
  ## Virus comparison
  
  # select time points to include in comparison
  output$selectCond <- renderUI({
    checkboxGroupInput("cond", label = "Choose time points", 
                       choiceNames = unique(sapply(unique(sub("[^_]*_", "Virus_", names(datasetInput))), 
                                                   function(x){
                                                     s <- ifelse(grepl("Mock", x), x, sub("vs_.*_", "vs_Virus_", x))
                                                     s <- sub("_vs_", " : ", s)
                                                     s <- sub("48h$", "48h (HCV only)", s)
                                                     return(s)
                                                   }, USE.NAMES = F)), 
                       choiceValues = unique(sapply(unique(sub(".*vs", "vs", names(datasetInput))), function(x) if(grepl("Mock", x)){return(x)}else{return(sub("_.*_", "_Virus_", x))}, USE.NAMES = F)))
  })
  
  # Help text
  observeEvent(input$help4, {
    showModal(modalDialog(
      title = tags$div("Help for ", tags$b("Virus Comparison")), 
      size = "l", 
      tags$div("This tab is meant to compare the genes that are differentially expressed after the infection with different viruses.", 
               "The sidebar shows the available options to adapt the plots:"), 
      HTML("<br>"), 
      img(src="viruses_sidebar.jpg"), 
      HTML("<br><br>"), 
      tags$div("The UpSet plot shows the intersections of genes in a matrix-like layout ", tags$a("(info)", href = "https://jku-vds-lab.at/tools/upset/", target="_blank"), ".", 
               "This approach allows to compare a large number of sets. Each bar shows an intersection of genes with the corresponding number. In the respective column below you can see the samples involved, marked with black dots.", 
               "For each virus all genes differentially expressed in at least one of the selected conditions are pooled. A gene is considered differentially expressed if: log2FoldChange (LFC) > ", tags$i("LFC cutoff"), " and adjusted p-value (padj) < 0.05 and a normalized gene count of at least one sample ≥ 10."), 
      img(src="viruses_upset.jpg"), 
      tags$div("Genes in an intersection can be listed in a sortable table along with a link to the NCBI database by clicking on the appropriate bar. The table can be downloaded in xlsx or csv format and then also contains the LFC and padj for all selected conditions.", 
               "The selected bar in the example (colored orange) shows, that there are 318 genes that are diff. expressed in at least one of the included time points (see sidebar) in NiV and RSV (see black dots in column below)."), 
      #HTML("<br><br>"), 
      img(src="viruses_upset_select.jpg"), 
      #HTML("<br><br>"), 
      tags$div("To display a selected intersection as heatmap click on the button in the sidebar. The user can get further information about the underlying data by clicking on a specific cell."), 
      img(src="viruses_heatmap_select.jpg"), 
      HTML("<br><br>"), 
      tags$div("The modebar beneath the heatmap allows the user to interact with and download the plot:"), 
      img(src="heatmap_modebar.jpg"), 
      tags$div("The user is able to zoom into the heatmap by dragging a box over the area of interest. The new heatmap is drawn in the box 'Selected sub-heatmap' next to the original heatmap."), 
      HTML("<br><br>"), 
      img(src="viruses_heatmap_zoom.jpg"), 
      tags$div("The user can also interact with this sub-heatmap using the buttons beneath the box. There the user can configure or download the plot as well as export it as table. To zoom further into the current sub-heatmap, the user can make it interactive by clicking the corresponding button in 'Configure sub-heatmap' and zoom into the heatmap as previously described."), 
      img(src="viruses_heatmap_zoom2.jpg"), 
      easyClose = TRUE
    ))
  })
  
  # download genes of selected plot area as csv or xlsx
  output$downallCSV <- downloadHandler(
    filename = function(){
      return(paste0(gsub("&", "_", unlist(input$upsetAll_click$name)), ".csv"))}, 
    content = function(f){
      df <- all.df()
      df$LinkToNCBI <- sub(".+\\\"(.+?)\\\".+", "\\1", df$LinkToNCBI)
      write.csv(df, f, row.names = F)
    }, 
    contentType = "csv")
  
  output$downallXLSX <- downloadHandler(
    filename = function(){
      return(paste0(gsub("&", "_", unlist(input$upsetAll_click$name)), ".xlsx"))}, 
    content = function(f){
      df <- all.df()
      df$LinkToNCBI <- sub(".+\\\"(.+?)\\\".+", "\\1", df$LinkToNCBI)
      write.xlsx(df, f)
    })
  
  # create dataframe of LFC and padj for the selected area
  all.df <- reactive({
    gene_id <- unlist(input$upsetAll_click$elems)
    cond <- ifelse(grepl("Virus", input$cond), sub("Virus", "(?!Mock).*", input$cond), input$cond)
    data.df <- Reduce(function(x, y)merge(x, y, by="SYMBOL", all = T), lapply(names(datasetInput[grep(paste(as.vector(outer(input$select, cond, paste, sep=".*")), collapse = "|"), names(datasetInput), perl = T)]), function(x){
      y <- datasetInput[[x]]
      y <- y[y$SYMBOL %in% gene_id, c("SYMBOL", "log2FoldChange", "padj")]
      colnames(y) <- c("SYMBOL", paste0(x, ".LFC"), paste0(x, ".padj"))
      return(y)
    })
    )
    data.df <- data.df[, c(1, mixedorder(colnames(data.df)[-1])+1), drop=F]
    data.df$LinkToNCBI <- createLink(data.df$SYMBOL)
    return(data.df)
  })
  
  # display table of genes in selected area
  output$clickedElementsAll <- renderDataTable({
    genes <- unlist(input$upsetAll_click$elems)
    if(is.null(data$all.dt) | is.null(genes) | is.null(input$select) | is.null(input$cond)){
      df <- data.frame("Symbol" = character(), "LinkToNCBI" = character())
    }else{
      df <- data.frame("Symbol" = genes, "LinkToNCBI" = createLink(unlist(input$upsetAll_click$elems)))
    }
    return(df)
  }, escape = F, rownames = F)
  
  # plot UpSet plot with selected viruses and conditions
  output$upsetAll <- renderUpsetjs({
    if(is.null(input$select) | is.null(input$cond)){
      return(NULL)
    }else{
      cond <- ifelse(grepl("Virus", input$cond), sub("Virus", "(?!Mock).*", input$cond), input$cond)
      id.list <- sapply(input$select, function(virus){
        #print(names(datasetInput[grep(paste(virus, cond, sep = ".*", collapse = "|"), names(datasetInput))]))
        unique(unlist(lapply(datasetInput[grep(paste(virus, cond, sep = ".*", collapse = "|"), names(datasetInput), perl = T)], function(x){
          return(x[which(abs(x$log2FoldChange) > input$LFCall & x$padj < 0.05 & apply(x[, grep("normalized", colnames(x))], 1, max) >= 10), "SYMBOL"])
        })))
      }, USE.NAMES = T, simplify = F)
      upsetjs() %>% fromList(id.list) %>% generateDistinctIntersections(limit = input$limit, empty = T) %>%  chartLayout(padding = 40, bar.padding = 0.2) %>%
        chartFontSizes(font.family = "Helvetica", chart.label = "16px", set.label = "14px", bar.label = "14px", axis.tick = "14px", export.label = "13px") %>% 
        interactiveChart("click") %>% chartProps(exportButtons=list(share=FALSE, vega=FALSE, dump=FALSE))
    }
  })
  
  # set heatmap to TRUE if button is clicked 
  observe({
    if(input$addHeatAll){
      data$heatAll <- T
    }
  })
  
  # remove heatmap if input changes
  observe({
    input$select
    input$cond
    data$heatAll <- F
    data$all.dt <- NULL
  })
  
  # render table if an intersection is clicked
  observe({
    input$upsetAll_click$elems
    data$all.dt <- 1
  })
  
  # create dataframe for heatmap of genes in table
  data_heat_all <- reactive({
    if(is.null(input$select) | is.null(input$cond) | is.null(input$upsetAll_click$elems)){
      return(NULL)
    }else{
      data <- all.df()
      data.rows <- data$SYMBOL
      data <- data[, -c(1, ncol(data)), drop = F]
      data <- data[, c(TRUE, FALSE), drop = F]
      colnames(data) <- sub(".LFC", "", colnames(data))
      if(nrow(data)>1){
        data <- apply(data, 2, as.numeric)
      }else{
        data <- as.data.frame(t(apply(data, 2, as.numeric)))
      }
      rownames(data) <- data.rows
      data <- data[order(rownames(data)), ]
      return(data)
    }
  })
  
  # click action function to create hover data
  click_action_all <- function(df, output) {
    output$infoAll <- renderUI({
      if(is.null(df)) { # have not clicked or brushed into the heatmap body
        DataFrame()
      } else {
        data.df <- data_heat_all()
        gene <- rownames(data.df[df$row_index, , drop=F])
        sample <- colnames(data.df[, df$column_index, drop=F])
        sample.parts <- strsplit(sample, "_vs_")[[1]]
        dataset <- datasetInput[[sample]]
        datainfo <- dataset[dataset$SYMBOL==gene, , drop=F]
        datainfo[, c(5:ncol(datainfo))] <- signif(datainfo[, c(5:ncol(datainfo))], 4) 
        datainfo$GENENAME <- ifelse(is.na(datainfo$GENENAME), "", datainfo$GENENAME)
        virus.norm <- apply(datainfo[, grep(sample.parts[1], colnames(datainfo)), drop=F], 1, paste, collapse="; ")
        mock.norm <- apply(datainfo[, grep(sub("_", ".*", sample.parts[2]), colnames(datainfo)), drop=F], 1, paste, collapse="; ")
        HTML(qq("<p style='background-color:white;color:black;border:1px solid black;padding:5px;font-size:14px;'>Sample: @{sample}<br />Gene: @{gene}<br />
                Genename: @{datainfo$GENENAME}<br />LFC: @{datainfo$log2FoldChange}<br />padj: @{datainfo$padj}<br />
                Normalized counts @{sample.parts[1]}: @{virus.norm}<br />Normalized counts @{sample.parts[2]}: @{mock.norm} </p>"))
      }
    })
  }
  
  # plot heatmap
  observe({
    if(data$heatAll){
      data$heatmap_id_all <- sample(1:1000, 1)
      df <- data_heat_all()
      if(!is.null(df)){
        ht <- plotHeatmap(df, colClust = F, rowClust = T)
        ht <- draw(ht)
        dev.off()
        InteractiveComplexHeatmapWidget(input, output, session, ht, output_id = "heatmapAll", action = "click", compact = FALSE, 
                                        heatmap_id = paste0("ht", data$heatmap_id_all), output_ui_float = F, output_ui = htmlOutput("infoAll"), 
                                        click_action = click_action_all, close_button = T, height2 = 500, layout = "1-(2|3)", 
                                        width1 = 600, height1 = 750, title1 = "Heatmap of genes shown in table")
      }
    }else{
      removeUI(paste0("#ht", data$heatmap_id_all, "_heatmap_widget *"), multiple = T)
    }
  }) 
  
  ## SNP analysis
  
  # Help text
  observeEvent(input$help6, {
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
      img(src="snp_mouseover.jpg"), 
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
        x.df$snp <- paste(x.df$POS, x.df$REF, x.df$ALT, sep = "_")
        return(x.df)
      }
    }, simplify = F))
    
    x_break <- which(sub("_[1|2]$", "", mixedsort(unique(v.df$Sample)))[-1] != sub("_[1|2]$", "", mixedsort(unique(v.df$Sample)))[-length(sub("_[1|2]$", "", mixedsort(unique(v.df$Sample))))])
    
    # create plotly heatmap for each segment
    p <- lapply(mixedsort(unique(v.df$Segment)), function(x){
      v.df.mod <- v.df[v.df$Segment==x, ]
      v.df.mod$snp <- factor(v.df.mod$snp, levels = mixedsort(unique(v.df.mod$snp)))
      v.df.mod$Sample <- factor(v.df.mod$Sample, levels = mixedsort(unique(v.df$Sample)))
      
      v.df.split <- split(v.df.mod[, c("snp", "AF")], v.df.mod$Sample)
      v.df.mat <- Reduce(function(x, y)merge(x, y, by="snp", all=T, no.dups=F), lapply(names(v.df.split), function(x){
        df.tmp <- v.df.split[[x]]
        colnames(df.tmp) <- c("snp", x)
        #df.tmp$id <- make.unique(as.character(df.tmp$snp))
        return(df.tmp)
      }))
      rownames(v.df.mat) <- v.df.mat$snp
      v.df.mat <- as.matrix(v.df.mat[, -1])
      v.df.mat[is.na(v.df.mat)] <- 0
      
      
      #create matrix with hoverinfo per sample and snp
      hover <- matrix(ncol = length(levels(v.df.mod$Sample)), nrow = length(unique(v.df.mod$snp)))
      colnames(hover) <- levels(v.df.mod$Sample)
      rownames(hover) <- mixedsort(unique(v.df.mod$snp))
      for(i in unique(v.df.mod$snp)){
        for(s in unique(v.df.mod$Sample)){
          if(nrow(v.df.mod[v.df.mod$Sample==s & v.df.mod$snp==i, ]) > 0){
            hover[i, s] <- sprintf("Sample: %s<br />Position: %s<br />Reference: %s<br />Alternative: %s<br />Frequency: %s<br />Depth: %s", 
                                   s, v.df.mod[v.df.mod$Sample==s & v.df.mod$snp==i, "POS"], v.df.mod[v.df.mod$Sample==s & v.df.mod$snp==i, "REF"], v.df.mod[v.df.mod$Sample==s & v.df.mod$snp==i, "ALT"], 
                                   v.df.mod[v.df.mod$Sample==s & v.df.mod$snp==i, "AF"], v.df.mod[v.df.mod$Sample==s & v.df.mod$snp==i, "DP"])
          }else{
            hover[i, s] <- NA
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
