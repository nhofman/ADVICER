#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
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

plotHeatmap <- function(x, row_subset = NA, distMethod = "euclidean", clusterMethod = "complete", clrn = 1, clcn = 1,
                        rowClust = T, colClust = T, fontsize_row = 0.8, fontsize_col = 10, annCol = NA, annRow = NA, border_col = "grey60", plot.fig = T,
                        legend.limit = 1, filter_col = NA, annotation_colors = NA, break_step = 0.1, display_numbers = F){
    if(is.na(row_subset[1])){
        xx <- x
    }else{
        xx <- x[row_subset,, drop = F]
    }
    #xx[xx==0] <- 0.000001
    xx <- na.omit(xx)
    if (!is.na(clrn) & !is.na(clcn)) {
        # set the custom distance and clustering functions
        hclustfunc <- function(x) hclust(x, method = clusterMethod)
        distfunc <- function(x) dist(x, method = distMethod)
        # perform clustering on rows and columns
        cl.row <- hclustfunc(distfunc(xx))
        cl.col <- hclustfunc(distfunc(t(xx)))
        
        # extract cluster assignments; i.e. k=8 (rows) k=5 (columns)
        gr.row <- cutree(cl.row, clrn)
        gr.col <- cutree(cl.col, clcn)
    }else{
        # set default values
        gr.row <- rep(1, nrow(xx))
        gr.col <- rep(1, ncol(xx))
    }
    if(!is.na(filter_col)){
        xx <- xx[,which(colMaxs(as.matrix(abs(xx)))>filter_col)]
    }
    
    breakSeq <- seq(quantile(unlist(xx), na.rm = TRUE, probs = (1-legend.limit)), quantile(unlist(xx), na.rm = TRUE, probs = legend.limit), break_step)
    if(length(breakSeq)==1){
        if(sign(breakSeq)==1){
            breakSeq <- seq(0,breakSeq,break_step)
        }else{
            breakSeq <- seq(breakSeq,0,break_step)
        }
    }
    color <- vector()
    if(-1%in%sign(unlist(xx))){
        color_down <- colorRampPalette(c("green", "white"))(length(seq(quantile(unlist(xx), na.rm = TRUE, probs = (1-legend.limit)), 0, break_step))) #blue(n=length(breakSeq)-1) #, low="blue", mid = "white", high="red")
        color <- c(color,color_down)
    }
    if(1%in%sign(unlist(xx))){
        color_up <- colorRampPalette(c("white", "red"))(length(seq(0, quantile(unlist(xx), na.rm = TRUE, probs = legend.limit), break_step)))
        color <- c(color,color_up)
    }
    if((nrow(xx)>1 | ncol(xx)>1)){
        pheatmap(xx, cluster_cols=colClust, cluster_rows=rowClust, clustering_distance_rows = distMethod, clustering_distance_cols = distMethod,
                 clustering_method = clusterMethod, annotation_col=annCol, annotation_row = annRow, 
                 breaks = breakSeq, color = color, annotation_colors = annotation_colors,
                 fontsize_row = fontsize_row, fontsize_col = fontsize_col, cutree_rows = clrn, 
                 border_color = border_col, display_numbers = display_numbers)
    }
}


createLink <- function(val) {
    sprintf('<a href="https://www.ncbi.nlm.nih.gov/gene?term=(%s%%5BGene%%20Name%%5D)%%20AND%%20human%%5BOrganism%%5D">Info</a>', val)
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
    lfc.gene$time <- ifelse(lfc.gene$time=="Mock_BPL", "uninfected", ifelse(grepl("BPL",lfc.gene$time), "inactivated", sub("Mock_","",lfc.gene$time)))
    lfc.min <- min(lfc.gene$symbol)
    lfc.max <- max(lfc.gene$symbol) + 0.25
    p1 <- ggplot(lfc.gene[grepl("h",lfc.gene$time),], aes(x=factor(time, levels = mixedsort(unique(time))), y=symbol, group=group, color=group)) + geom_line() + geom_point() +
        labs(title = paste("Expression profile of", gene), y = "Log2FoldChange", x = "Time", color = "Virus", caption = "*: padj < 0.01    **: padj < 0.001    ***: padj < 0.0001") + 
        geom_text(aes(label=sig), nudge_y = 0.1, show.legend = FALSE) + ylim(c(lfc.min, lfc.max)) +
        theme(plot.caption = element_text(hjust = 0, size = 15, family = "sans"), legend.text = element_text(size = 15, family = "sans"), 
              legend.title = element_text(size = 15, face = "bold", family = "sans"), axis.text = element_text(size = 15),
              axis.title = element_text(size = 15, face = "bold", family = "sans"), plot.title = element_text(size = 20, family = "sans", face = "bold")) +
        theme(plot.margin = unit(c(1,1,1,1), "cm"))
    p2 <-  ggplot(lfc.gene[!grepl("h",lfc.gene$time),], aes(x=time, y=symbol, group=group, fill=group)) + geom_col(position = "dodge") +
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

datasetInput <- NULL
#d.V <- vector()
#d.MA <- vector()

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
        tabPanel("Upload",
                 #sidebarPanel(
                 #    fileInput("file", label = h5(br("File input")), multiple = T),
                 #),
                 mainPanel(
                     textOutput("uploadText"),
                     br(),
                     tableOutput("table"),
                     #h2("Uploaded files"),
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
                     plotOutput("heatmapTime", height = "750px"),
                     
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
                     checkboxGroupInput("selectCond", label = h5(strong(" Select conditions to include")), 
                                        choices = c("infected"="Mock_.*h", "inactivated control"="24h_Vs_.*BPL", "uninfected control"="Mock_BPL")),
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
                     downloadButton("downHeatAll", "Download heatmap"),
                     plotOutput("heatmapAll", height = 800)
                     
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
filedir <- "/home/nina/Documents/Virus_project/analyses/host/deseq2_new/deseq2_comparisons_shrunken/csv/"
files.list <- list.files(filedir, full.names = T)
filedir.vcf <- "/home/nina/Documents/Virus_project/variant_calling_new/"
vcf.files <- list.files(filedir.vcf, pattern = ".*[1|2].vcf", full.names = T)

# Define server logic 
server = function(input, output, session) {
    
    # Define reactive values
    data <- reactiveValues(d=NULL, df=NULL, dM=NULL, dfM=NULL, heat=NULL, dt=NULL, snp=NULL)
    
    # Read data
    withProgress(message = 'PROCESSING DATA...', detail = "This may take a while...", value = 0,{
        datasetInput <- lapply(files.list, function(x){
            print(x)
            f <- read.csv(x, header = TRUE, stringsAsFactors = F, row.names = 1)
            incProgress(1/length(files.list))
            return(f)
        })
        #Sys.sleep(5)
            names(datasetInput) <- sub("deseq2_results_(.*).csv","\\1",basename(files.list))
            
            vcf.list <- lapply(vcf.files, read.vcfR)
            names(vcf.list) <- sub(".vcf", "", basename(vcf.files))
    })
    
    output$uploadText <- renderText({
        if(!is.null(datasetInput)){
            print("Uploaded data:")
        }
    })
    
    output$table <- renderTable(
        data.frame(names(datasetInput))
        , colnames = F, rownames = F)
    
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
            df <- df[d.V,c("SYMBOL","log2FoldChange","padj")]
            df$NCBI <- createLink(df$SYMBOL)
            return(df)
        }
    }, escape = F)
    
    # plot Volcano plot of selected sample
    output$volcanoPlot <- renderPlotly({
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
        p <- plot_ly(file.data, type = "scatter", x = ~log2FoldChange, y = ~-log10(padj), color = ~col, colors = c("up"="red","down"="green","not significant"="black"), 
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
            df <- df[d.MA,c("SYMBOL","log2FoldChange","padj")]
            df$NCBI <- createLink(df$SYMBOL)
            return(df)
        }
    }, escape = F)
    
    # plot MA plot of selected sample
    output$MAPlot <- renderPlotly({
        if(is.null(datasetInput) | is.null(input$fileToPlot))
            return(NULL)
        file.data <- datasetInput[[input$fileToPlot]]
        file.data$col <- ifelse(file.data$log2FoldChange > input$LFCsingle & file.data$padj < 0.05 & apply(file.data[,grep("normalized", colnames(file.data))],1,max) >= 10, "up", 
                                ifelse(file.data$log2FoldChange < -(input$LFCsingle) & file.data$padj < 0.05 & apply(file.data[,grep("normalized", colnames(file.data))],1,max) >= 10, "down", "not significant"))
        p <- ggplot(file.data, aes(x = log2FoldChange, y = -log(padj), color = col, text = SYMBOL)) +
            geom_point(aes(size = 0.4)) +
            xlab("log fold change") +
            ylab("-log10(P-value)") + scale_color_manual(values = c("yellow", "purple", "black")) + theme(legend.position = "bottom", legend.direction = "horizontal")
        #p %>% ggplotly(tooltip = c("SYMBOL", "log2FoldChange", "padj")) %>% layout(legend = list(orientation = "h"))#, x = 0.4, y = -0.2))
        p <- plot_ly(file.data, type = "scatter", x = ~log2(baseMean), y = ~log2FoldChange, color = ~col, colors =c("up"="red","down"="green","not significant"="black"),
                     mode = "markers", marker = list(size = 5), customdata = ~SYMBOL,
                     text = ~paste("Gene: ", SYMBOL, '<br>LFC: ', log2FoldChange, '<br>padj: ', padj),
                     source = "M")
        p %>% layout(legend = list(orientation = "h", y = -0.2))
    })
    
    ## Compare time points of defined virus
    # download genes of selected plot area as csv or xlsx
    output$down_vCSV <- downloadHandler(
        filename = function(){
            if(input$plotType == "UpSet"){
                return(paste0(paste(input$upset_click$name, collapse = "_"), ".csv"))
            }else{
                return(paste0(paste(input$upsetVenn_click$name, collapse = "_"), ".csv"))
            }
        },
        content = function(f){
            if(input$plotType == "UpSet"){
                write.csv(data.frame("Symbol"=unlist(input$upset_click$elems)), f, row.names = F)
            }else{
                write.csv(data.frame("Symbol"=unlist(input$upsetVenn_click$elems)), f, row.names = F)
            }
        },
        contentType = "csv")
    
    output$down_vXLSX <- downloadHandler(
        filename = function(){
            if(input$plotType == "UpSet"){
                return(paste0(paste(input$upset_click$name, collapse = "_"), ".xlsx"))
            }else{
                return(paste0(paste(input$upsetVenn_click$name, collapse = "_"), ".xlsx"))
            }
        },
        content = function(f){
            if(input$plotType == "UpSet"){
                write.xlsx(data.frame("Symbol"=unlist(input$upset_click$elems)), f, row.names = F)
            }else{
                write.xlsx(data.frame("Symbol"=unlist(input$upsetVenn_click$elems)), f, row.names = F)
            }
        })
    
    # select times to plot
    output$selectTime <- renderUI({
        selectInput("time", label = "Choose times to plot", choices = sub(".*Vs","Vs",names(datasetInput)[grep(input$select_v, names(datasetInput))]), multiple = TRUE)
    })
    
    # set heatmap to NULL at selecion of new virus
    observeEvent(input$select_v, {
        data$heat <- NULL   
        data$dt <- NULL
    })
    
    # add heatnmap
    observeEvent(input$addHeat, {
        data$heat <- 1
    })
    
    output$heatmapTime <- renderPlot({
        print(plot_heat_time())
    })
    
    # download heatmap
    output$downHeatTime <- downloadHandler(
        filename = function(){
            return(paste0(input$select_v, ".png"))},
        content = function(f){
            png(filename = f)
            print(plot_heat_time())
            dev.off()
        }
    )
    
    # plot heatmap
    plot_heat_time <- reactive({
        h <- data$heat
        click <- NULL
        if(input$plotType=="UpSet"){
            click <- input$upset_click$elems
        }
        if(input$plotType=="Venn"){
            click <- input$upsetVenn_click$elems
        }
        if(is.null(input$select_v) | is.null(input$time) | is.null(click) | is.null(h)){
            return(NULL)
        }else{
            data.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL", all = T),lapply(names(datasetInput[grep(paste0(input$select_v, ".*_", input$time, collapse = "|"), names(datasetInput))]), function(x){
                y <- datasetInput[[x]]
                y <- y[unlist(click), c("SYMBOL","log2FoldChange")]
                colnames(y) <- c("SYMBOL", x)
                return(y)
            })
            )
            data.df <- data.df[rowSums(is.na(data.df))!=ncol(data.df),]
            rownames(data.df) <- data.df$SYMBOL
            data.df <- data.df[,-1]
            data.df <- data.df[,mixedorder(colnames(data.df))]
            plotHeatmap(data.df, colClust = F, border_col = NA, fontsize_row = 10)
        }
    })
    
    # clear data table if plot type changes
    observeEvent(input$plotType,{
        data$dt <- NULL
    })
    
    observeEvent(input$upset_click$elems,{
        data$dt <- 1
    })
    
    observeEvent(input$upsetVenn_click$elems,{
        data$dt <- 1
    })
    
    # display table with elements in clicked plot area of UpSet plot
    output$clickedElements <- renderDataTable({
        #if(input$plotType=="UpSet"){
        dt <- data$dt
        if(is.null(dt)){
            data.frame("SYMBOL"=NULL, "LinkToNCBI"=NULL)
        }else{
            data.frame("Symbol"=unlist(input$upset_click$elems), "LinkToNCBI" = createLink(unlist(input$upset_click$elems)))
        }
        #if(input$plotType=="Venn"){
        #    data.frame("Symbol"=unlist(input$upsetVenn_click$elems), "LinkToNCBI" = createLink(unlist(input$upsetVenn_click$elems)))
        #}
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
            #names(id.list) <- sub("Vs","Vs<br>",names(id.list))
            upsetjs() %>% fromList(id.list) %>% generateDistinctIntersections %>% chartFontSizes(font.family = "sans", set.label = "8px") %>% interactiveChart()
        }
    })
    
    # display table with elements in clicked plot area of Venn plot
    output$clickedElementsVenn <- renderDataTable({
        dt <- data$dt
        if(is.null(dt)){
            return(NULL)
        }else{
            data.frame("Symbol"=unlist(input$upsetVenn_click$elems), "LinkToNCBI" = createLink(unlist(input$upsetVenn_click$elems)))
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
            #names(id.list) <- sub("Vs","Vs<br>",names(id.list))
            upsetjsVennDiagram() %>% fromList(id.list) %>% interactiveChart()
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
    # download genes of selected plot area as csv or xlsx
    output$downallCSV <- downloadHandler(
        filename = function(){
            return(paste0(paste(input$select, collapse = "_"), ".csv"))},
        content = function(f){
            write.csv(data.frame("Symbol"=unlist(input$upsetAll_click$elems)), f)
        },
        contentType = "csv")
    
    output$downallXLSX <- downloadHandler(
        filename = function(){
            return(paste0(paste(input$select, collapse = "_"), ".xlsx"))},
        content = function(f){
            write.xlsx(data.frame("Symbol"=unlist(input$upsetAll_click$elems)), f)
        })
    
    # display table of genes in selected area
    output$clickedElementsAll <- renderDataTable({
        data.frame("Symbol"=unlist(input$upsetAll_click$elems), "Link" = createLink(unlist(input$upsetAll_click$elems)))
    }, escape = F)
    
    # plot UpSet plot with selected viruses and conditions
    output$upsetAll <- renderUpsetjs({
        if(is.null(input$select) | is.null(input$selectCond)){
            return(NULL)
        }else{
            id.list <- sapply(input$select,function(virus){
                unique(unlist(lapply(datasetInput[grep(paste(virus, input$selectCond, sep = ".*"), names(datasetInput))], function(x){
                    return(x[which(abs(x$log2FoldChange) > input$LFCall & x$padj < 0.05 & apply(x[,grep("normalized", colnames(x))],1,max) >= 10),"SYMBOL"])
                })))
            }, USE.NAMES = T, simplify = F)
            id.list <- id.list[sapply(id.list, length) > 0]
            upsetjs() %>% fromList(id.list) %>% generateDistinctIntersections(limit = input$limit) %>% interactiveChart()
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
    output$heatmapAll <- renderPlot({
        print(plot_heat_all())
    })
    
    plot_heat_all <- reactive({
        if(is.null(input$select) | is.null(input$selectCond) | is.null(input$upsetAll_click$elems)){
            return(NULL)
        }else{
            data.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL", all = T),lapply(names(datasetInput[grep(paste(input$select, input$selectCond, collapse = "|", sep = ".*"), names(datasetInput))]), function(x){
                y <- datasetInput[[x]]
                y <- y[unlist(input$upsetAll_click$elems), c("SYMBOL","log2FoldChange")]
                colnames(y) <- c("SYMBOL", x)
                return(y)
            })
            )
            #print(data.df$SYMBOL[duplicated(data.df$SYMBOL)])
            data.df <- data.df[rowSums(is.na(data.df))!=ncol(data.df),]
            rownames(data.df) <- data.df$SYMBOL
            data.df <- data.df[,-1]
            data.df <- data.df[,mixedorder(colnames(data.df))]
            plotHeatmap(data.df, colClust = F, border_col = NA, fontsize_row = 12)
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
            
            plot_ly(x=colnames(v.df.mat), y=rownames(v.df.mat), z=v.df.mat, type = "heatmap", colors = palette(50), showlegend = F, zmin = 0, zmax = 1, zauto = F,
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
