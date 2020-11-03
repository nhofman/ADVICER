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
        color_down <- colorRampPalette(c("blue", "white"))(length(seq(quantile(unlist(xx), na.rm = TRUE, probs = (1-legend.limit)), 0, break_step))) #blue(n=length(breakSeq)-1) #, low="blue", mid = "white", high="red")
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
    lfc.gene$key <- ifelse(lfc.gene$padj<0.0001,"padj<0.0001",ifelse(lfc.gene$padj<0.001,"padj<0.001",ifelse(lfc.gene$padj<0.01,"padj<0.01",NA)))
    lfc.gene$group <- sub("_.*","",rownames(lfc.gene))
    #lfc.gene$time <- factor(sub(".*_(.*)\\..*","\\1",rownames(lfc.gene)), levels = c("BPL","3h","6h","12h","24h","48h"))
    #lfc.gene$time <- factor(sub(".*_(.*)","\\1",rownames(lfc.gene)), levels = c("BPL","3h","6h","12h","24h","48h"))
    lfc.gene$time <- factor(sub(".*_(.*)","\\1",rownames(lfc.gene)), levels = mixedsort(unique(sub(".*_(.*)","\\1",rownames(lfc.gene)))))
    ggplot(lfc.gene, aes(x=time, y=symbol, group=group, color=group)) + geom_line() + geom_point() +
        labs(title = paste("Expression profile of", gene), y = "Log2FoldChange", x = "Time", color = "Virus", caption = "*: padj < 0.01    **: padj < 0.001    ***: padj < 0.0001") + 
        geom_text(aes(label=sig), nudge_y = 0.1, show.legend = NA) + 
        theme(plot.caption = element_text(hjust = 0, size = 15, family = "sans"), legend.text = element_text(size = 15, family = "sans"), 
              legend.title = element_text(size = 15, face = "bold", family = "sans"), axis.text = element_text(size = 15),
              axis.title = element_text(size = 15, face = "bold", family = "sans"), plot.title = element_text(size = 20, family = "sans", face = "bold")) 
}

datasetInput <- NULL
#d.V <- vector()
#d.MA <- vector()

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("Differential Gene Expression Analysis"),
    
    # Sidebar with a slider input for number of bins 
    navbarPage(
        "DGE analysis", 
                    tabPanel("Upload",
                             sidebarPanel(
                                 fileInput("file", label = h5(br("File input")), multiple = T),
                             ),
                             mainPanel(
                                 h3("Uploaded files"),
                                 tableOutput("table")
                             )
                    ),
                    tabPanel("Volcano Plot / MAPlot",
                             sidebarPanel(
                                 selectInput("fileToPlot", label = "Choose file", choices = names(datasetInput)),
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
                                             selected = "CoV229E"),
                                 sliderInput("LFC",
                                             h5(strong("LFC cutoff:")),
                                             min = 0,
                                             max = 10,
                                             value = 1, 
                                             step = 0.5),
                                 actionButton("addHeat", "Show heatmap"),
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
                                     plotOutput("heatmapTimeUp", height = "750px")
                                 ),
                                 conditionalPanel(
                                     condition = 'input.plotType == "Venn"',
                                     upsetjsOutput("upsetVenn"), #, click = "upsetClick"),
                                     br(),
                                     DT::dataTableOutput("clickedElementsVenn"),
                                     br(),
                                     br(),
                                     plotOutput("heatmapTimeVenn", height = "750px")
                                 ),
                                 
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
                                 plotOutput("geneX"),
                                 br(),
                                 #tableOutput("sig")
                             )),
                    tabPanel("Virus Comparison",
                             sidebarPanel(
                                 checkboxGroupInput("select", label = h5(strong("Select viruses")), 
                                                    choices = list("H1N1", "H5N1", "MERS", "CoV229E", "RVFV", "SFSV", "RSV", "NIV", "EBOV", "MARV", "HCV", "LASV"), 
                                                    selected = NULL), #c("H1N1", "H5N1", "MERS", "CoV229E", "RVFV", "SFSV", "RSV", "NIV", "EBOV")),
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
                                 
                             ))
        )
    )

# Define server logic required to draw a histogram
server = function(input, output, session) {
    datasetInput <- reactive({
        infile <- input$file
        if(is.null(infile))
            return(NULL)
        files.list <- list()
        for(i in 1:length(input$file[,1])){
            files.list[[i]] <- read.csv(input$file[[i,'datapath']], header = TRUE, stringsAsFactors = F, row.names = 1)
            print(input$file[[i,'name']])
            names(files.list)[i] <- sub("deseq2_results_(.*).csv","\\1",input$file[[i,'name']])
        }
        return(files.list)
    })
    
    output$table = renderTable(
        data.frame(names(datasetInput()))
        , colnames = F, rownames = F)
    
    
    output$fileToPlot <- renderUI({
        selectInput("file", label = "Choose file", choices = names(datasetInput()))
    })
    
    observe({
        updateSelectInput(
            session, "fileToPlot",
            choices = names(datasetInput())
        )
    })
    
    output$downSingle <- downloadHandler(
        filename = function(){
            return(paste0(input$fileToPlot, "_selected", ".csv"))},
        content = function(f){
            write.csv(data.frame("Symbol"=unlist(input$upsetAll_click$elems)), f)
        },
        contentType = "csv")
    
    data <- reactiveValues(d=NULL, df=NULL, dM=NULL, dfM=NULL, heat=NULL)
    
    observeEvent(input$fileToPlot,{
        data$d <- NULL
        data$dV <- NULL
        data$dM <- NULL
        data$dfM <- NULL
    })
    
    observeEvent(event_data("plotly_click", source = "V"), {
        req(datasetInput())
        req(input$fileToPlot)
        data$d <- event_data("plotly_click", source = "V")
        a <- data$dV
        a <- c(a, data$d$customdata)
        data$dV <- a
    })
    
    observeEvent(event_data("plotly_click", source = "M"), {
        req(datasetInput())
        req(input$fileToPlot)
        data$dM <- event_data("plotly_click", source = "M")
        a <- data$dfM
        a <- c(a, data$dM$customdata)
        data$dfM <- a
    })
    
    output$clickVP <- renderDataTable({
        d <- data$d
        d.V <- data$dV
        if(is.null(d)){
            return(NULL)
        }else{
            df <- datasetInput()[[input$fileToPlot]]
            df <- df[d.V,c("SYMBOL","log2FoldChange","padj")]
            df$NCBI <- createLink(df$SYMBOL)
            return(df)
        }
    }, escape = F)
    
    output$volcanoPlot <- renderPlotly({
        if(is.null(datasetInput()) | is.null(input$fileToPlot))
            return(NULL)
        file.data <- datasetInput()[[input$fileToPlot]]
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
    
    output$clickMA <- renderDataTable({
        d <- data$dM
        d.MA <- data$dfM
        if(is.null(d)){
            return(NULL)
        }else{
            df <- datasetInput()[[input$fileToPlot]]
            df <- df[d.MA,c("SYMBOL","log2FoldChange","padj")]
            df$NCBI <- createLink(df$SYMBOL)
            return(df)
        }
    }, escape = F)
    
    output$MAPlot <- renderPlotly({
        if(is.null(datasetInput()) | is.null(input$fileToPlot))
            return(NULL)
        file.data <- datasetInput()[[input$fileToPlot]]
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
    
    observeEvent(input$select_v, {
        data$heat <- NULL   
    })
    
    observeEvent(input$addHeat, {
        data$heat <- 1
    })
    
    
    output$heatmapTimeUp <- renderPlot({
        h <- data$heat
        if(is.null(input$select_v) | is.null(input$upset_click$elems) | is.null(h)){
            return(NULL)
        }else{
            data.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL", all = T),lapply(names(datasetInput()[grep(paste(input$select_v, collapse = "|"), names(datasetInput()))]), function(x){
                y <- datasetInput()[[x]]
                y <- y[unlist(input$upset_click$elems), c("SYMBOL","log2FoldChange")]
                colnames(y) <- c("SYMBOL", x)
                return(y)
            })
            )
            data.df <- data.df[rowSums(is.na(data.df))!=ncol(data.df),]
            rownames(data.df) <- data.df$SYMBOL
            data.df <- data.df[,-1]
            plotHeatmap(data.df, colClust = F, border_col = NA, fontsize_row = 12)
        }
    })
    
    output$heatmapTimeVenn <- renderPlot({
        h <- data$heat
        if(is.null(input$select_v) | is.null(input$upsetVenn_click$elems) | is.null(h)){
            return(NULL)
        }else{
            data.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL", all = T),lapply(names(datasetInput()[grep(paste(input$select_v, collapse = "|"), names(datasetInput()))]), function(x){
                y <- datasetInput()[[x]]
                y <- y[unlist(input$upsetVenn_click$elems), c("SYMBOL","log2FoldChange")]
                colnames(y) <- c("SYMBOL", x)
                return(y)
            })
            )
            data.df <- data.df[rowSums(is.na(data.df))!=ncol(data.df),]
            rownames(data.df) <- data.df$SYMBOL
            data.df <- data.df[,-1]
            plotHeatmap(data.df, colClust = F, border_col = NA, fontsize_row = 12)
        }
    }) 
    
    output$clickedElements <- renderDataTable({
        data.frame("Symbol"=unlist(input$upset_click$elems), "LinkToNCBI" = createLink(unlist(input$upset_click$elems)))
    }, escape = F)
    
    output$upset <- renderUpsetjs({
        id.list <- lapply(datasetInput()[grep(input$select_v, names(datasetInput()))], function(x){
            return(x[which(abs(x$log2FoldChange) > input$LFC & x$padj < 0.05, apply(x[,grep("normalized", colnames(x))],1,max) >= 10),"SYMBOL"])
        })
        id.list <- id.list[sapply(id.list, length) > 0]
        #names(id.list) <- sub("Vs","Vs<br>",names(id.list))
        upsetjs() %>% fromList(id.list) %>% generateDistinctIntersections %>% chartFontSizes(font.family = "sans", set.label = "8px") %>% interactiveChart()
    })
    
    output$clickedElementsVenn <- renderDataTable({
        data.frame("Symbol"=unlist(input$upsetVenn_click$elems), "LinkToNCBI" = createLink(unlist(input$upsetVenn_click$elems)))
    }, escape = F)
    
    output$upsetVenn <- renderUpsetjs({
        #if(!input$plotType=="Venn"){
        id.list <- lapply(datasetInput()[grep(input$select_v, names(datasetInput()))], function(x){
            return(x[which(abs(x$log2FoldChange) > input$LFC & x$padj < 0.05, apply(x[,grep("normalized", colnames(x))],1,max) >= 10),"SYMBOL"])
        })
        id.list <- id.list[sapply(id.list, length) > 0]
        #names(id.list) <- sub("Vs","Vs<br>",names(id.list))
        upsetjsVennDiagram() %>% fromList(id.list) %>% interactiveChart()
        #}
    })
    
    output$downGeneX <- downloadHandler(
        filename = function(){
            return(paste0(paste(input$selectVirus, collapse = "_"), "_", toupper(input$gene), ".png"))},
        content = function(f){
            ggsave(f, plot_geneX(), "png")
        }
    )
    
    output$geneX <- renderPlot({
        print(plot_geneX())
    })
    
    plot_geneX <- reactive({
        if(is.null(input$selectVirus) | is.null(input$gene)){
            return(NULL)
        }else{
            lfc.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL",all=T),lapply(names(datasetInput()[grep(paste(input$selectVirus, collapse = "|"), names(datasetInput()))]), function(x){
                x.df <- data.frame(datasetInput()[[x]]$SYMBOL,datasetInput()[[x]]$log2FoldChange)
                colnames(x.df) <- c("SYMBOL",x)
                return(x.df)}))
            rownames(lfc.df) <- lfc.df$SYMBOL
            lfc.df <- lfc.df[,-1]
            padj.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL",all=T),lapply(names(datasetInput()[grep(paste(input$selectVirus, collapse = "|"), names(datasetInput()))]), function(x){
                x.df <- data.frame(datasetInput()[[x]]$SYMBOL,datasetInput()[[x]]$padj)
                colnames(x.df) <- c("SYMBOL",x)
                return(x.df)}))
            rownames(padj.df) <- padj.df$SYMBOL
            padj.df <- padj.df[,-1]
            return(plotExpression(expr.df = lfc.df, padj.df = padj.df, gene = toupper(input$gene)))
        }
    })
    
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
    
    output$clickedElementsAll <- renderDataTable({
        data.frame("Symbol"=unlist(input$upsetAll_click$elems), "Link" = createLink(unlist(input$upsetAll_click$elems)))
    }, escape = F)
    
    output$upsetAll <- renderUpsetjs({
        if(is.null(input$select)){
            return(NULL)
        }else{
            id.list <- sapply(input$select,function(virus){
                unique(unlist(lapply(datasetInput()[grep(virus, names(datasetInput()))], function(x){
                    return(x[which(abs(x$log2FoldChange) > input$LFCall & x$padj < 0.05 & apply(x[,grep("normalized", colnames(x))],1,max) >= 10),"SYMBOL"])
                })))
            }, USE.NAMES = T, simplify = F)
            id.list <- id.list[sapply(id.list, length) > 0]
            upsetjs() %>% fromList(id.list) %>% generateDistinctIntersections(limit = input$limit) %>% interactiveChart()
        }
    })
    
    output$downHeatAll <- downloadHandler(
        filename = function(){
            return(paste0(paste(input$select, collapse = "_"), ".png"))},
        content = function(f){
            png(filename = f)
            print(plot_heat_all())
            dev.off()
        }
    )
    
    output$heatmapAll <- renderPlot({
        print(plot_heat_all())
    })
    
    plot_heat_all <- reactive({
        if(is.null(input$select) | is.null(input$upsetAll_click$elems)){
            return(NULL)
        }else{
            data.df <- Reduce(function(x,y)merge(x,y,by="SYMBOL", all = T),lapply(names(datasetInput()[grep(paste(input$select, collapse = "|"), names(datasetInput()))]), function(x){
                y <- datasetInput()[[x]]
                y <- y[unlist(input$upsetAll_click$elems), c("SYMBOL","log2FoldChange")]
                colnames(y) <- c("SYMBOL", x)
                return(y)
            })
            )
            #print(data.df$SYMBOL[duplicated(data.df$SYMBOL)])
            data.df <- data.df[rowSums(is.na(data.df))!=ncol(data.df),]
            rownames(data.df) <- data.df$SYMBOL
            data.df <- data.df[,-1]
            plotHeatmap(data.df, colClust = F, border_col = NA, fontsize_row = 12)
        }
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
