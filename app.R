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
  #print(seq(quantile(xx.unlist, na.rm = TRUE, probs = (1-legend.limit)), quantile(xx.unlist, na.rm = TRUE, probs = legend.limit), break_step))
  #breakSeq <- seq(quantile(xx.unlist, na.rm = TRUE, probs = (1-legend.limit)), quantile(xx.unlist, na.rm = TRUE, probs = legend.limit), break_step)
  #if(length(breakSeq)==1){
  #    if(sign(breakSeq)==1){
  #        breakSeq <- seq(0,breakSeq,break_step)
  #    }else{
  #        breakSeq <- seq(breakSeq,0,break_step)
  #    }
  #}
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
  if((nrow(xx)>1 | ncol(xx)>1)){
    print(paste("File: ", file))
    if(setWidth){
      #p <- p %>% layout(width=ncol(xx)*50, height = nrow(xx)*10)
      if(!is.na(file)){
        p <- heatmaply(xx, Colv = F, colors = color, custom_hovertext = hover, plot_method = "plotly", width = ncol(xx)*70, height = nrow(xx)*30, fontsize_col = fontsize_c, fontsize_row = fontsize_r, file = file) 
      }else{
        p <- heatmaply(xx, Colv = F, colors = color, custom_hovertext = hover, plot_method = "plotly", width = ncol(xx)*70, height = nrow(xx)*30, fontsize_col = fontsize_c, fontsize_row = fontsize_r) 
      }
    }else{
      if(!is.na(file)){
        p <- heatmaply(xx, Colv = F, colors = color, custom_hovertext = hover, plot_method = "plotly", fontsize_col = fontsize_c, fontsize_row = fontsize_r, file = file) 
      }else{
        p <- heatmaply(xx, Colv = F, colors = color, custom_hovertext = hover, plot_method = "plotly", fontsize_col = fontsize_c, fontsize_row = fontsize_r) 
      }
    }
    p <- p %>% config(toImageButtonOptions=list(format="png", width=800, height=600, scale=2))
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
  #lfc.gene$time <- factor(sub(".*_(.*)\\..*","\\1",rownames(lfc.gene)), levels = c("BPL","3h","6h","12h","24h","48h"))
  #lfc.gene$time <- factor(sub(".*_(.*)","\\1",rownames(lfc.gene)), levels = c("BPL","3h","6h","12h","24h","48h"))
  lfc.gene$time <- sub(".*_(.*_.*)","\\1",rownames(lfc.gene)) #factor(sub(".*_(.*)","\\1",rownames(lfc.gene)), levels = mixedsort(unique(sub(".*_(.*)","\\1",rownames(lfc.gene)))))
  lfc.gene$time <- ifelse(lfc.gene$time=="Mock_BPL", "Virus BPL:Mock BPL", ifelse(grepl("BPL",lfc.gene$time), "Virus 24h:Virus BPL", sub("Mock_","",lfc.gene$time)))
  lfc.min <- min(lfc.gene$symbol)
  lfc.max <- max(lfc.gene$symbol) + 0.25
  p1 <- ggplot(lfc.gene[grepl("h$",lfc.gene$time),], aes(x=factor(time, levels = mixedsort(unique(time))), y=symbol, group=group, color=group)) + geom_line() + geom_point() +
    labs(title = paste("Expression profile of", gene), y = "Log2FoldChange", x = "Time", color = "Virus", caption = "*: padj < 0.01    **: padj < 0.001    ***: padj < 0.0001") + 
    geom_text(aes(label=sig), nudge_y = 0.1, show.legend = FALSE) + ylim(c(lfc.min, lfc.max)) +
    theme(plot.caption = element_text(hjust = 0, size = 15, family = "sans"), legend.text = element_text(size = 15, family = "sans"), 
          legend.title = element_text(size = 15, face = "bold", family = "sans"), axis.text = element_text(size = 15),
          axis.title = element_text(size = 15, face = "bold", family = "sans"), plot.title = element_text(size = 20, family = "sans", face = "bold")) +
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  p2 <-  ggplot(lfc.gene[grepl("BPL",lfc.gene$time),], aes(x=time, y=symbol, group=group, fill=group)) + geom_col(position = "dodge") +
    labs(title = paste("Expression profile of", gene), y = "Log2FoldChange", x = "Time", fill = "Virus", caption = "*: padj < 0.01    **: padj < 0.001    ***: padj < 0.0001") + 
    geom_text(aes(label=sig), position=position_dodge(width=0.9), vjust=-0.25, show.legend = F) + ylim(c(lfc.min, lfc.max)) +
    theme(plot.caption = element_text(hjust = 0, size = 15, family = "sans"), legend.text = element_text(size = 15, family = "sans"), 
          legend.title = element_text(size = 15, face = "bold", family = "sans"), axis.text = element_text(size = 15),
          axis.title = element_text(size = 15, face = "bold", family = "sans"), plot.title = element_text(size = 20, family = "sans", face = "bold")) +
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  #p <- c(p1, p2)
  #gridExtra::grid.arrange(p1, p2, ncol = 1)
  #margin = theme(plot.margin = unit(c(2,2,2,2), "cm"))
  gridExtra::grid.arrange(p1, p2, ncol = 1, heights = c(10, 10))
}

readRenviron(".Renviron")
if(Sys.getenv("DATADIR") != ""){
  shinyOptions(filedir = Sys.getenv("DATADIR"))
}

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
           "
      ),
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
    "DGE analysis", 
    tabPanel("Download data",
             #sidebarPanel(
             #    fileInput("file", label = h5(br("File input")), multiple = T),
             #),
             sidebarPanel(
               #textOutput("uploadText"),
               #br(),
               htmlOutput("download"),
               #tableOutput("table"),
               #h2("Uploaded files"),
             ),
             mainPanel(
               downloadButton("filesDown", "Download selected files"),
               textOutput("debug", container = pre)
             )
    ),
    tabPanel("Volcano Plot / MAPlot",
             sidebarPanel(
               htmlOutput("fileSelect"),
               radioButtons("plotTypeSingle", label = h5(strong("Displayed Plot")), choices = c("Volcano", "MAPlot"), selected = "Volcano"),
               sliderInput("LFCsingle",
                           h5(strong("LFC cutoff:")),
                           min = 0,
                           max = 10,
                           value = 1, 
                           step = 0.5),
               br(),
               #downloadButton("downSingle", "Download")
             ),
             mainPanel(
               #downloadButton("downPlot", "Download plot"),
               actionButton("export", "Export plot"),
               conditionalPanel(
                 condition = 'input.plotTypeSingle == "Volcano"',
                 #d.V <- vector(),
                 plotlyOutput("volcanoPlot", height = "500px"),
                 br(),
                 dataTableOutput("clickVP")
               ),
               conditionalPanel(
                 condition = 'input.plotTypeSingle == "MAPlot"',
                 plotlyOutput("MAPlot", height = "500px"),
                 br(),
                 dataTableOutput("clickMA")
               )
             )),
    tabPanel("Compare time points",
             sidebarPanel(
               radioButtons("plotType", label = h5(strong("Displayed Plot")), choices = c("Venn", "UpSet"), selected = "Venn"),
               selectInput("select_v", label = h5(strong("Select virus")), 
                           choices = list("H1N1", "H5N1", "MERS", "CoV229E", "RVFV", "SFSV", "RSV", "NIV", "EBOV", "MARV", "HCV", "LASV"), 
                           selected = NULL),
               htmlOutput("selectTime"),
               sliderInput("LFC",
                           h5(strong("LFC cutoff:")),
                           min = 0,
                           max = 10,
                           value = 1, 
                           step = 0.5),
               actionButton("addHeat", "Show heatmap"),
               downloadButton("downHeatTime", "Download heatmap"),
               br(),
               br(),
               h5(strong("Download table")),
               downloadButton("down_vCSV", "Download as csv"),
               br(),
               downloadButton("down_vXLSX", "Download as xlsx")
             ),
             mainPanel(
               conditionalPanel(
                 condition = 'input.plotType == "UpSet"',
                 upsetjsOutput("upset"), #, click = "upsetClick"),
                 br(),
                 DT::dataTableOutput("clickedElements"),
                 br(),
                 br(),
                 #plotOutput("heatmapTimeUp", height = "750px")
               ),
               conditionalPanel(
                 condition = 'input.plotType == "Venn"',
                 upsetjsOutput("upsetVenn"), #, click = "upsetClick"),
                 br(),
                 DT::dataTableOutput("clickedElementsVenn"),
                 br(),
                 br(),
                 #plotOutput("heatmapTimeVenn", height = "750px")
               ),
               plotlyOutput("heatmapTime", height = "750px"),
               
             )
    ),
    tabPanel("Gene expression",
             sidebarPanel(
               checkboxGroupInput("selectVirus", label = h5(strong("Select viruses")), 
                                  choices = list("H1N1", "H5N1", "MERS", "CoV229E", "RVFV", "SFSV", "RSV", "NIV", "EBOV", "MARV", "HCV", "LASV"), 
                                  selected = NULL),
               textInput("gene", label = h3("Gene input"), value = NULL),
               downloadButton("downGeneX","Download plot"),
             ), 
             mainPanel(
               plotOutput("geneX", height = 1000),
               br(),
               #tableOutput("sig")
             )),
    tabPanel("Virus Comparison",
             sidebarPanel(
               checkboxGroupInput("select", label = h5(strong("Select viruses")), 
                                  choices = list("H1N1", "H5N1", "MERS", "CoV229E", "RVFV", "SFSV", "RSV", "NIV", "EBOV", "MARV", "HCV", "LASV"), 
                                  selected = NULL), #c("H1N1", "H5N1", "MERS", "CoV229E", "RVFV", "SFSV", "RSV", "NIV", "EBOV")),
               htmlOutput("selectCond"),
               em("(BPL samples were taken after 24h)"),
               #checkboxGroupInput("selectCond", label = h5(strong(" Select conditions to include")), 
               #                   choices = list("")),
               #choices = c("infected"="Mock_.*h", "inactivated control"="24h_Vs_.*BPL", "uninfected control"="Mock_BPL")),
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
               h5(strong("Download table")),
               downloadButton("downallCSV", "Download as csv"),
               br(),
               downloadButton("downallXLSX", "Download as xlsx")
             ),
             mainPanel(
               upsetjsOutput("upsetAll"), #, click = "upsetClick"),
               br(),
               DT::dataTableOutput("clickedElementsAll"),
             )
    ),
    tabPanel("Heatmap",
             mainPanel(
               #downloadButton("downHeatAll", "Download heatmap"),
               plotlyOutput("heatmapAll", height = 1500)
               
             )),
    tabPanel("SNP analysis",
             sidebarPanel({
               selectInput("virusSNP", label = h5(strong("Select virus")), 
                           choices = list("H1N1", "H5N1", "MERS", "CoV229E", "RVFV", "SFSV", "RSV", "NIV", "EBOV", "MARV", "HCV", "LASV"), 
                           selected = NULL)
             }),
             mainPanel(
               #plotOutput("heatSNP", height = "600px", click = "SNPclick"),
               #DT::dataTableOutput("SNPdata")
               plotlyOutput("heatSNP", height = "1000px"),
               #verbatimTextOutput("SNPdata")
             ))
  )
)

# Define directory containing data
#filedir <- "/data"
#filedir <- "/home/nina/Documents/Virus_project/analyses/host/deseq2_new/deseq2_comparisons_shrunken/data/"
filedir <- getShinyOption("filedir")
#filedir <- "/vol/sfb1021/SFB1021_Virus/dge_analyses_new/analyses/host/deseq2/deseq2_comparisons_shrunken/" #"/home/nina/Documents/Virus_project/analyses/host/deseq2_stranded/csv/"
files.list <- list.files(filedir, pattern = "*.csv", full.names = T)
#filedir.vcf <- getShinyOption("filedirvcf")
#filedir.vcf <- "/vol/sfb1021/SFB1021_Virus/dge_analyses_new/variant_analyses/variant_calling/" 
#filedir <- "/home/nina/Documents/Virus_project/variant_calling_new/" #"/data" 
vcf.files <- list.files(filedir, pattern = ".*[1|2].vcf", full.names = T)

# Define server logic 
server = function(input, output, session) {
  
  # Define reactive values
  data <- reactiveValues(d=NULL, df=NULL, dM=NULL, dfM=NULL, heat=NULL, heat_data = NULL, heat_hover = NULL, dt=NULL, snp=NULL, all.df=NULL, virus.df=NULL)
  
  # Read data
  withProgress(message = 'PROCESSING DATA...', detail = "This may take a while...", value = 0,{
    datasetInput <- lapply(files.list, function(x){
      f <- read.csv(x, header = TRUE, stringsAsFactors = F, row.names = 1)
      incProgress(1/length(files.list))
      return(f)
    })
    #Sys.sleep(5)
    names(datasetInput) <- sub("deseq2_results_(.*).csv","\\1",basename(files.list))
    
    vcf.list <- lapply(vcf.files, read.vcfR)
    names(vcf.list) <- sub(".vcf", "", basename(vcf.files))
  })
  
  datasetInput <- lapply(datasetInput, function(x){
    x[,c(5:ncol(x))] <- signif(x[,c(5:ncol(x))], 4)
    return(x)
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
      if(input$filesToDown=="all"){
        files <- files.list
      }else{
        for (i in input$filesToDown){
          #write each sheet to a csv file, save the name
          fileName <- files.list[grep(i,files.list)] #paste(i,".csv",sep = "")
          #write.table(datasetInput[i],fileName,sep = ',', row.names = F, col.names = T)
          files <- c(files,fileName)
        }
      }
      zip(f, files)
    },
    contentType = "csv"
  )
  
  output$debug <- renderPrint({
    orca_sys <- Sys.which("orca") 
    #orca_help <- base::system("orca -h", intern = T) 
    orca_help <- processx::run("orca", "-h")
    cat(paste("Orca\n"))
    cat(paste("Sys.which:", orca_sys))
    cat(paste("\nOrca help:", orca_help$stdout))
  })
  
  ## Volcano Plot / MAPlot
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
  observeEvent(event_data("plotly_click", source = "V"), {
    req(datasetInput)
    req(input$fileToPlot)
    data$d <- event_data("plotly_click", source = "V")
    a <- data$dV
    a <- c(a, data$d$customdata)
    data$dV <- a
  })
  
  observeEvent(event_data("plotly_click", source = "M"), {
    req(datasetInput)
    req(input$fileToPlot)
    data$dM <- event_data("plotly_click", source = "M")
    a <- data$dfM
    a <- c(a, data$dM$customdata)
    data$dfM <- a
  })
  
  # display table of clicked elements and associated LFC and padj in Volcano plot
  output$clickVP <- renderDataTable({
    d <- data$d
    d.V <- data$dV
    if(is.null(d)){
      return(NULL)
    }else{
      df <- datasetInput[[input$fileToPlot]]
      df <- df[df$SYMBOL %in% d.V,c("SYMBOL","log2FoldChange","padj")]
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
    p <- ggplot(file.data, aes(x = log2FoldChange, y = -log10(padj), color = col, text = SYMBOL)) +
      geom_point(aes(size = 0.4)) +
      xlab("log fold change") +
      ylab("-log10(P-value)") + scale_color_manual(values = c("yellow", "purple", "black")) + theme(legend.position = "bottom", legend.direction = "horizontal")
    #p %>% ggplotly(tooltip = c("SYMBOL", "log2FoldChange", "padj")) %>% layout(legend = list(orientation = "h"))#, x = 0.4, y = -0.2))
    p <- plot_ly(file.data, type = "scatter", x = ~log2FoldChange, y = ~-log10(padj), color = ~col, colors = c("up"="red","down"="blue","not significant"="black"), 
                 mode = "markers", marker = list(size = 5), customdata = ~SYMBOL,
                 text = ~paste("Gene: ", SYMBOL, '<br>LFC: ', log2FoldChange, '<br>padj: ', padj),
                 source = "V")
    p <- p %>% layout(legend = list(orientation = "h", y = -0.2))
    return(p)
  })
  
  # display table of clicked elements and associated LFC and padj in MA plot
  output$clickMA <- renderDataTable({
    d <- data$dM
    d.MA <- data$dfM
    if(is.null(d)){
      return(NULL)
    }else{
      df <- datasetInput[[input$fileToPlot]]
      df <- df[df$SYMBOL %in% d.MA,c("SYMBOL","log2FoldChange","padj")]
      df$NCBI <- createLink(df$SYMBOL)
      return(df)
    }
  }, escape = F)
  
  # plot MA plot of selected sample
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
      p <- ggplot(file.data, aes(x = log2FoldChange, y = -log(padj), color = col, text = SYMBOL)) +
        geom_point(aes(size = 0.4)) +
        xlab("log fold change") + ylab("-log10(P-value)") + scale_color_manual(values = c("yellow", "purple", "black")) + 
        theme(legend.position = "bottom", legend.direction = "horizontal", axis.title = element_text(size = 16, face = "bold"), axis.text = element_text(size = 16, face = "bold"))
      #p %>% ggplotly(tooltip = c("SYMBOL", "log2FoldChange", "padj")) %>% layout(legend = list(orientation = "h"))#, x = 0.4, y = -0.2))
      p <- plot_ly(file.data, type = "scatter", x = ~log2(baseMean), y = ~log2FoldChange, color = ~col, colors =c("up"="red","down"="blue","not significant"="black"),
                   mode = "markers", marker = list(size = 5), customdata = ~SYMBOL,
                   text = ~paste("Gene: ", SYMBOL, '<br>LFC: ', log2FoldChange, '<br>padj: ', padj),
                   source = "M")
      p <- p %>% layout(legend = list(orientation = "h", y = -0.2))
      return(p)
    }
  })
  
  output$downPlot <- downloadHandler(
    filename = function(){
      return(paste0(input$plotTypeSingle,"_",input$fileToPlot,".pdf"))
    },
    content = function(f){
      print(f)
      if(input$plotTypeSingle=="MAPlot"){
        orca(plotMA(), f,"pdf")
      }else{
        orca(plotVolcano(), file = f, format = "pdf")
      } 
    }, contentType = "application/pdf"
  )
  
  observeEvent(input$export, {
    f <- paste0(input$plotTypeSingle,"_",input$fileToPlot,".pdf")
    if(input$plotTypeSingle=="MAPlot"){
      orca(plotMA(), f,"pdf", more_args = c("--disable-gpu", "--enable-webgl"))
    }else{
      orca(plotVolcano(), file = f, format = "pdf", more_args = c("--disable-gpu", "--enable-webgl"))
    } 
  })
  
  ## Compare time points of defined virus
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
      #if(!is.null(data.df)){
      plotHeatmap(data.df, colClust = F, border_col = NA, fontsize_r = 10, hover = hover, file = f)
      #}else{
      #  return(NULL)
      #}
      #orca(f, p, format = "pdf", more_args = c("--disable-gpu", "--enable-webgl"))
    }
  )
  
  # select times to plot
  output$selectTime <- renderUI({
    checkboxGroupInput("time", label = "Choose times to plot (for Venn: max. 5)",
                       choiceNames = unique(sapply(names(datasetInput)[grep(input$select_v, names(datasetInput))], 
                                                   function(x){
                                                     s <- sub("^[^_]*_","",x)
                                                     s <- ifelse(grepl("Mock",x), sub("_Mock_.*","_Mock",s), s)
                                                     #  s <- sub("_Vs_", " : ", s)
                                                     #s <- sub("48h$", "48h (HCV only)", s)
                                                     return(s)
                                                   }, USE.NAMES = F)),
                       choiceValues = sub(".*Vs","Vs",names(datasetInput)[grep(input$select_v, names(datasetInput))]))
  })
  
  # set heatmap to NULL at selecion of new virus
  observeEvent(input$select_v, {
    data$heat <- NULL   
    data$dt <- NULL
  })
  
  # add heatnmap
  observeEvent(input$addHeat, {
    data$heat <- 1
    data_heat_time()
  })
  
  output$heatmapTime <- renderPlotly({
    data.df <- data$heat_data
    hover <- data$heat_hover
    if(!is.null(data.df)){
      p <- plotHeatmap(data.df, colClust = F, border_col = NA, fontsize_r = 10, hover = hover)
      print(p)
    }else{
      return(NULL)
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
      #print(paste0(input$select_v, ".*_", input$time, collapse = "|"))
      data.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL", all = T),lapply(names(datasetInput[grep(paste0(input$select_v, ".*_", input$time, collapse = "|"), names(datasetInput))]), function(x){
        y <- datasetInput[[x]]
        #y$log2FoldChange <- ifelse(y$padj<0.05 & apply(y[,grep("normalized", colnames(y))],1,max) >= 10, y$log2FoldChange, "-")
        y <- y[y$SYMBOL %in% virus_id, c("SYMBOL","log2FoldChange","padj")]
        colnames(y) <- c("SYMBOL", paste0(x,".LFC"), paste0(x,".padj"))
        return(y)
      })
      )
      data.df <- data.df[,c(1,mixedorder(colnames(data.df)[-1])+1),drop=F]
      data.df$Link <- createLink(data.df$SYMBOL)
      #data$virus.df <- data.df
      return(data.df)
    }
  })
  
  # plot heatmap
  plot_heat_time <- reactive({
    data.df <- data$heat_data
    hover <- data$heat_hover
    if(!is.null(data.df)){
      p <- plotHeatmap(data.df, colClust = F, border_col = NA, fontsize_r = 10, hover = hover)
      return(p)
    }else{
      return(NULL)
    }
  })
  
  data_heat_time <- reactive({
    h <- data$heat
    if(is.null(h)){
      return(NULL)
    }else{
      data.df <- virus.df()
      data.rows <- data.df$SYMBOL
      data.df <- data.df[,-c(1,ncol(data.df))]
      data.df <- apply(data.df,2,as.numeric)
      rownames(data.df) <- data.rows
      data.df <- data.df[,c(TRUE, FALSE)]
      colnames(data.df) <- sub(".LFC","",colnames(data.df))
      hover <- matrix(ncol = ncol(data.df), nrow = nrow(data.df))
      colnames(hover) <- colnames(data.df)
      rownames(hover) <- rownames(data.df)
      for(i in rownames(hover)){
        for(s in colnames(hover)){
          dataset <- datasetInput[[s]]
          s.gr <- sub("_.*","",strsplit(s,"_Vs_")[[1]])
          count.1 <- paste(dataset[dataset$SYMBOL==i, grep(paste0("normalized.*", s.gr[1]), colnames(dataset))], collapse = "; ")
          count.2 <- paste(dataset[dataset$SYMBOL==i, grep(paste0("normalized.*", s.gr[2]), colnames(dataset))], collapse = "; ")
          if(i %in% dataset$SYMBOL){
            hover[i,s] <- sprintf("Gene: %s<br />Sample: %s<br />LFC: %s<br />padj: %s<br />%s: %s<br />%s: %s", 
                                  i, s, dataset[dataset$SYMBOL==i,"log2FoldChange"], dataset[dataset$SYMBOL==i,"padj"],
                                  s.gr[1], count.1, s.gr[2], count.2)
            #paste(colnames(df), df, sep = ":", collapse = "\n"))
          }else{
            hover[i,s] <- sprintf("Gene: %s<br />Sample: %s<br />LFC: %s<br />padj: %s<br />", 
                                  i, s, NA, NA)
          }
        }
      }
      data$heat_data <- data.df
      data$heat_hover <- hover
      #return(list(data=data.df, hover=hover.df))
    }
  })
  
  # clear data table if plot type changes
  observeEvent(input$plotType,{
    data$dt <- NULL
    data$heat <- NULL
  })
  
  observeEvent(input$upset_click$elems,{
    data$dt <- 1
    #return(unlist(input$upset_click$elems))
  })
  
  observeEvent(input$upsetVenn_click$elems,{
    data$dt <- 1
    #return(unlist(input$upsetVenn_click$elems))
  })
  
  
  # display table with elements in clicked plot area of UpSet plot
  output$clickedElements <- renderDataTable({
    #if(input$plotType=="UpSet"){
    dt <- data$dt
    if(is.null(dt)){
      data.frame("SYMBOL"=NULL, "LinkToNCBI"=NULL)
    }else{
      virus.df()
      #data.frame("Symbol"=unlist(input$upset_click$elems), "LinkToNCBI" = createLink(unlist(input$upset_click$elems)))
    }
  }, escape = F)
  
  # plot UpSet plot of selected virus und times
  output$upset <- renderUpsetjs({
    if(is.null(input$select_v) | is.null(input$time)){
      return(NULL)
    }else{
      id.list <- lapply(datasetInput[grep(paste0(input$select_v, ".*_", input$time, collapse = "|"), names(datasetInput))], function(x){
        return(x[which(abs(x$log2FoldChange) > input$LFC & x$padj < 0.05, apply(x[,grep("normalized", colnames(x))],1,max) >= 10),"SYMBOL"])
      })
      id.list <- id.list[sapply(id.list, length) > 0]
      names(id.list) <- sub(".*_Vs_", "", names(id.list))
      names(id.list) <- sub("Mock_", "", names(id.list))
      upsetjs() %>% fromList(id.list) %>% generateDistinctIntersections %>% chartFontSizes(font.family = "sans", set.label = "14px", bar.label = "14px", axis.tick = "12px") %>% interactiveChart()
    }
  })
  
  
  # display table with elements in clicked plot area of Venn plot
  output$clickedElementsVenn <- renderDataTable({
    dt <- data$dt
    if(is.null(dt)){
      return(NULL)
    }else{
      virus.df()
    }
  }, escape = F)
  
  # plot Venn diagram of selected virus und times
  output$upsetVenn <- renderUpsetjs({
    if(is.null(input$select_v) | is.null(input$time)){
      return(NULL)
    }else{
      id.list <- lapply(datasetInput[grep(paste0(input$select_v, ".*_", input$time, collapse = "|"), names(datasetInput))], function(x){
        return(x[which(abs(x$log2FoldChange) > input$LFC & x$padj < 0.05, apply(x[,grep("normalized", colnames(x))],1,max) >= 10),"SYMBOL"])
      })
      id.list <- id.list[sapply(id.list, length) > 0]
      names(id.list) <- sub(".*_Vs_", "", names(id.list))
      names(id.list) <- sub("Mock_", "", names(id.list))
      upsetjsVennDiagram() %>% fromList(id.list) %>% chartFontSizes(font.family = "sans", set.label = "14px", value.label = "14px", axis.tick = "12px") %>% interactiveChart()
    }
  })
  
  ## Gene expression plot
  # download gene expression plot
  output$downGeneX <- downloadHandler(
    filename = function(){
      return(paste0(paste(input$selectVirus, collapse = "_"), "_", toupper(input$gene), ".png"))},
    content = function(f){
      ggsave(f, plot_geneX(), "png", width = 8, height = 10, dpi = 400)
    }
  )
  
  # plot gene expression of selected viruses over time
  output$geneX <- renderPlot({
    print(plot_geneX())
  })
  
  plot_geneX <- reactive({
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
      return(plotExpression(expr.df = lfc.df, padj.df = padj.df, gene = toupper(input$gene)))
    }
  })
  
  ## Virus comparison
  
  # select time points to include in comparison
  output$selectCond <- renderUI({
    checkboxGroupInput("cond", label = "Choose times to include", 
                       choiceNames = unique(sapply(unique(sub("[^_]*_","Virus_",names(datasetInput))), 
                                                   function(x){
                                                     s <- ifelse(grepl("Mock",x), x, sub("Vs_.*_","Vs_Virus_",x))
                                                     s <- sub("_Vs_", " : ", s)
                                                     s <- sub("48h$", "48h (HCV only)", s)
                                                     return(s)
                                                   }, USE.NAMES = F)),
                       choiceValues = unique(sapply(unique(sub(".*Vs","Vs",names(datasetInput))), function(x) if(grepl("Mock",x)){return(x)}else{return(sub("_.*_","_Virus_",x))}, USE.NAMES = F)))
  })
  
  all.df <- reactive({
    #data$all.df
    gene_id <- unlist(input$upsetAll_click$elems)
    cond <- ifelse(grepl("Virus",input$cond), sub("Virus", "(?!Mock).*",input$cond), input$cond)
    data.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL", all = T),lapply(names(datasetInput[grep(paste(as.vector(outer(input$select, cond, paste, sep=".*")), collapse = "|"), names(datasetInput), perl = T)]), function(x){
      y <- datasetInput[[x]]
      #y$log2FoldChange <- ifelse(y$padj<0.05 & apply(y[,grep("normalized", colnames(y))],1,max) >= 10, y$log2FoldChange, "-")
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
      return(paste0(paste(input$select, collapse = "_"), ".csv"))},
    content = function(f){
      #write.csv(data.frame("Symbol"=unlist(input$upsetAll_click$elems)), f)
      write.csv(all.df(), f, row.names = F)
    },
    contentType = "csv")
  
  output$downallXLSX <- downloadHandler(
    filename = function(){
      return(paste0(paste(input$select, collapse = "_"), ".xlsx"))},
    content = function(f){
      write.xlsx(all.df(), f)
    })
  
  # display table of genes in selected area
  output$clickedElementsAll <- renderDataTable({
    data.frame("Symbol"=unlist(input$upsetAll_click$elems), "Link" = createLink(unlist(input$upsetAll_click$elems)))
  }, escape = F)
  
  # plot UpSet plot with selected viruses and conditions
  output$upsetAll <- renderUpsetjs({
    if(is.null(input$select) | is.null(input$cond)){
      return(NULL)
    }else{
      cond <- ifelse(grepl("Virus",input$cond), sub("Virus", "(?!Mock).*",input$cond), input$cond)
      id.list <- sapply(input$select,function(virus){
        #print(names(datasetInput[grep(paste(virus, cond, sep = ".*", collapse = "|"), names(datasetInput))]))
        unique(unlist(lapply(datasetInput[grep(paste(virus, cond, sep = ".*", collapse = "|"), names(datasetInput), perl = T)], function(x){
          return(x[which(abs(x$log2FoldChange) > input$LFCall & x$padj < 0.05 & apply(x[,grep("normalized", colnames(x))],1,max) >= 10),"SYMBOL"])
        })))
      }, USE.NAMES = T, simplify = F)
      id.list <- id.list[sapply(id.list, length) > 0]
      upsetjs() %>% fromList(id.list) %>% generateDistinctIntersections(limit = input$limit) %>% chartFontSizes(font.family = "sans", set.label = "14px", bar.label = "14px", axis.tick = "12px") %>% interactiveChart()
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
  
  # plot heatmap of genes in selected area
  output$heatmapAll <- renderPlotly({
    plot_heat_all()
  })
  
  plot_heat_all <- reactive({
    if(is.null(input$select) | is.null(input$cond) | is.null(input$upsetAll_click$elems)){
      return(NULL)
    }else{
      # cond <- ifelse(grepl("Virus",input$cond), sub("Virus", "(?!Mock).*",input$cond), input$cond)
      # #   paste(as.vector(outer(virus.levels, cond, paste, sep=".*")), collapse = "|")
      # 
      # data.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL", all = T),lapply(names(datasetInput[grep(paste(as.vector(outer(input$select, cond, paste, sep=".*")), collapse = "|"), names(datasetInput), perl = T)]), function(x){
      #     y <- datasetInput[[x]]
      #     y <- y[unlist(input$upsetAll_click$elems), c("SYMBOL","log2FoldChange")]
      #     colnames(y) <- c("SYMBOL", x)
      #     return(y)
      # })
      # )
      # #print(data.df$SYMBOL[duplicated(data.df$SYMBOL)])
      # data.df <- data.df[rowSums(is.na(data.df))!=ncol(data.df),,drop=F]
      # rownames(data.df) <- data.df$SYMBOL
      # data.df <- data.df[,-1,drop=F]
      # data.df <- data.df[,mixedorder(colnames(data.df)),drop=F]
      data <- all.df()
      data.rows <- data$SYMBOL
      data <- data[,-c(1,ncol(data))]
      data <- data[,c(TRUE, FALSE)]
      colnames(data) <- sub(".LFC","",colnames(data))
      #for(i in 1:ncol(data)){
      #    data[,i] <- sub("-", NA, data[,i])
      #}
      #data[data=="-"] <- NA
      data <- apply(data,2,as.numeric)
      rownames(data) <- data.rows
      hover <- matrix(ncol = ncol(data), nrow = nrow(data))
      colnames(hover) <- colnames(data)
      rownames(hover) <- rownames(data)
      for(i in rownames(hover)){
        for(s in colnames(hover)){
          dataset <- datasetInput[[s]]
          s.gr <- sub("_.*","",strsplit(s,"_Vs_")[[1]])
          count.1 <- paste(dataset[dataset$SYMBOL==i, grep(paste0("normalized.*", s.gr[1]), colnames(dataset))], collapse = "; ")
          count.2 <- paste(dataset[dataset$SYMBOL==i, grep(paste0("normalized.*", s.gr[2]), colnames(dataset))], collapse = "; ")
          if(i %in% dataset$SYMBOL){
            hover[i,s] <- sprintf("Gene: %s<br />Sample: %s<br />LFC: %s<br />padj: %s<br />%s: %s<br />%s: %s", 
                                  i, s, dataset[dataset$SYMBOL==i,"log2FoldChange"], dataset[dataset$SYMBOL==i,"padj"],
                                  s.gr[1], count.1, s.gr[2], count.2)
            #paste(colnames(df), df, sep = ":", collapse = "\n"))
          }else{
            hover[i,s] <- sprintf("Gene: %s<br />Sample: %s<br />LFC: %s<br />padj: %s<br />", 
                                  i, s, NA, NA)
          }
        }
      }
      plotHeatmap(data, colClust = F, rowClust = T, border_col = NA, hover = hover, setWidth = T, fontsize_r = 12, fontsize_c = 12)
    }
  })
  
  ## SNP analysis
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
        layout(shapes=lapply(x_break-0.5, vline)) %>% 
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
