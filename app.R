library(shiny)
library(MALDIHDX)
library(DT)
library(shinyjs)
library(ggplot2)


options(shiny.maxRequestSize = 30*1024^2)



peptide.identifications <- import.identifications()

TP <- unique(pep10[,2])

defSNR = 5

#L <- lapply(peptide.features,function(x) DFtoSpec(x))
DefaultAllCents <- lapply(isolate(peptide.features$data)  ,function(x) mainCentNewMod2(x))
#DefMEMTable <- MEMHDXall2(DefaultAllCents)

# Define UI
ui <- fluidPage(theme=shinytheme("cerulean"),
                
                
                
                useShinyjs(),
                
                navbarPage(
                
                            strong("MALDIHDX"),
                           
                          tabPanel(icon=icon("home"),"About"),
                           
                           
                           tabPanel(icon=icon("clipboard"),"Tutorial",
                                    
                                    fluidRow(column(12,titlePanel("Tutorial"),align="center")),
                                    
                                    hr(),
                                    
                                    fluidRow(column(12,img(src = "dodecane.gif", height = 140, width = 400),align="center")),hr(),
                                    
                                    fluidRow(column=12,htmlOutput("video"),align="center")
                                    
                           ),
                           
                           
                           
                           
                           
                           tabPanel(icon=icon("cogs"),"Analysis",
                                    
                                    # App title
                                    titlePanel("Centroid calculation and validation of HDX-MS Experiments"),
                                    
                                    # Main panel for displaying outputs ----
                                    mainPanel(
                                        
                                        # Output: Tabset
                                        tabsetPanel(type = "tabs",
                                                    
                                                    tabPanel("Data Import",
                                                             
                                                             sidebarLayout(
                                                                 
                                                                 sidebarPanel(
                                                                     
                                                                     helpText("Please upload a Zip file containing
                                          Scaffold file and Mass Spectra"),
                                                                     
                                                                     
                                                                     fileInput("FileInput", "Choose folder",
                                                                               multiple = TRUE, buttonLabel = "Browse...",placeholder = "No file selected"
                                                                     ),
                                                                     
                                                                     actionButton("FileImport","Import all Spectra")
                                                                 ),
                                                                 
                                                                 mainPanel()
                                                             )  
                                                    ),  
                                                    
                                                    tabPanel("Identifications",
                                                             
                                                             dataTableOutput("table")),
                                                    
                                                    tabPanel("Centroid Plots",
                                                             
                                                             fluidRow(
                                                                 
                                                                 column(8,
                                                                        
                                                                        plotOutput("plot"),
                                                                        hr(),
                                                                        
                                                                        fluidRow(column(6,
                                                                                        
                                                                                        selectInput("var", label = "Select peptide",
                                                                                                    choices = peptide.identifications[,3],
                                                                                                    selected = peptide.identifications[1,3]),
                                                                                        
                                                                                        selectInput("varRep", label = "Select Replicate",
                                                                                                    choices = c("1","2","3"),
                                                                                                    selected = "1")),
                                                                                 
                                                                                 column(6,
                                                                                        
                                                                                        selectInput("varTime", label = "Select timepoint",
                                                                                                    choices = TP,
                                                                                                    selected = TP[1]),
                                                                                        
                                                                                        selectInput("varBound", label = "Select State",
                                                                                                    choices = c("Bound", "Unbound"),
                                                                                                    selected = "Unbound"))
                                                                                 
                                                                                 
                                                                        ) 
                                                                 ),
                                                                 
                                                                 column(4,br(),
                                                                        
                                                                        h4("Centroid Parameters"), hr(),
                                                                        
                                                                        sliderInput('SNR','Signal to Noise',min=1,max=10,
                                                                                    value = 5, step =0.1),
                                                                        actionButton("resSNR","Reset to default"),br(),br(),
                                                                        
                                                                        sliderInput('BPI','% Base Peak Intensity',min=0.1,max=95,
                                                                                    value = 50, step =1),
                                                                        actionButton("resBPI","Reset to default"),br(),br(),hr(),br(),
                                                                        
                                                                        fluidRow(column(12,
                                                                                        
                                                                                        tableOutput("tableCent"),actionButton("expCent","Export centroid to output table"),span(textOutput("Done"), style="color:blue")  
                                                                                        
                                                                        ))
                                                                 )
                                                             )
                                                    ),
                                                    
                                                    
                                                    tabPanel("Uptake Plots",
                                                             
                                                             sidebarLayout(
                                                                 
                                                                 sidebarPanel(
                                                                     
                                                                     helpText("Select peptide to view Deuterium uptake plots"),
                                                                     
                                                                     selectInput("varPep", label = "Select peptide",
                                                                                 choices = peptide.identifications[,3],
                                                                                 selected = peptide.identifications[1,3])
                                                                     
                                                                 ),
                                                                 
                                                                 mainPanel(plotOutput("plotUptake"))
                                                             ) 
                                                             
                                                             
                                                    ),
                                                    
                                                    
                                                    tabPanel("Output table",
                                                             
                                                             dataTableOutput("MEMHDXTable")
                                                    ),
                                                    tabPanel("Test Table",
                                                             dataTableOutput("Testtable")
                                                    )
                                                    
                                                    
                                                    
                                        )
                                    ) 
                                    
                           ),tabPanel(icon=icon("book"),"Publications"),tabPanel(icon=icon("star"),"News"))
)



# Define server logic
server <- function(input, output) {
    
    
    peptide.features <- reactiveValues()
    AllCentReact <- reactiveValues()
    
    
    
    observeEvent(input$FileImport,
                 
                 
                 {
                     withProgress(message="Importing and analysing data...",value=0,{
                         
                         
                         setwd("data/HDX220318")
                         peptide.features$data <- importNew()
                         incProgress(1/2,"Calculating centroids based on default parameters..")
                         AllCentReact$data <- lapply(peptide.features$data  ,function(x) mainCentNewMod2(x))
                         setwd("..")
                         setwd("..")
                         
                         incProgress(1/2,"Complete")
                         
                     })
                     
                 }
    )
    
    L <- reactive({
        
        lapply(peptide.features$data,function(x) DFtoSpec(x))
        
    })
    
    
    #DefaultAllCents <- lapply(peptide.features,function(x) mainCentNewMod2(x))
    #DefMEMTable <- MEMHDXall2(DefaultAllCents)
    
    
    observeEvent(input$FileInput,
                 {
                     
                     
                     setwd("data")
                     file.remove(list.files())
                     unlink(list.files(),recursive=T)
                     infile <- input$FileInput
                     if(is.null(infile))
                         return(NULL)
                     unzip(infile$datapath)
                     setwd("..")
                 }   
    )
    
    
    
    
    
    curRow <- reactive({
        
        state <- switch(input$varBound,"Unbound" = "UNBOUND","Bound" = "BOUND")
        index <- which(c(DefMEMTable$Exposure == input$varTime & DefMEMTable$Replicate==input$varRep&DefMEMTable$State == state & DefMEMTable$Sequence==input$var))
        return(index)
    })
    
    heading <- reactive({
        paste(input$var,"Timepoint",input$varTime,"\n Replicate",
              input$varRep,"-",input$varBound)
    })
    
    
    
    
    
    
    observeEvent(c(input$var,input$varRep,input$resSNR,input$varBound,input$varTime), {
        
        reset("SNR")
        
    })
    
    observeEvent(c(input$var,input$varRep,input$resBPI,input$varBound,input$varTime), {
        
        reset("BPI")
        output$Done <- renderText("Export manual centroid") 
        
    })
    
    
    MEMTable <- reactiveValues(data=DefMEMTable) 
    
    MEMTable2 <- reactive({
        MEMTable
    })
    
    #AllCentReact <- reactiveValues(data = DefaultAllCents)
    
    #AllCentReact <- reactive({
    
    #   lapply(peptide.features$data  ,function(x) mainCentNewMod2(x))
    
    #})
    
    
    observeEvent(input$expCent,{
        temp <- MEMTable$data
        temp[curRow(),9] <- Centroid()
        MEMTable$data <- temp
        text <- paste("Row",as.character(curRow()),"Centroid was exported")
        output$Done <- renderText(text)
        
        
        PepNumber <- match(input$var,peptide.identifications[,3])
        stateUp <- switch(input$varBound,"Unbound" = "A","Bound" = "B")
        stateRepUp <- paste(stateUp,input$varRep,sep="")
        
        subNo <- switch(stateRepUp,"A1"=1,"A2"=2,"A3"=3,"B1"=4,"B2"=5,"B3"=6)
        
        
        tempCent <- AllCentReact$data[[PepNumber]][[subNo]]
        tempCent[tempCent$'time (min)'==input$varTime,2] <- Centroid()
        AllCentReact$data[[PepNumber]][[subNo]] <- tempCent
        
    } )
    
    
    
    
    
    
    
    CurSpec <- reactive({
        
        PepNo <- match(input$var,peptide.identifications[,3])
        
        PepSpectra <- L()[[PepNo]]
        
        state <- switch(input$varBound,"Unbound" = "A","Bound" = "B")
        
        stateRep <- paste(state,input$varRep,sep="")
        
        cond <- lapply(PepSpectra,function(x)
            metaData(x)$StateRep == stateRep & metaData(x)$time == input$varTime)
        
        CurSpec <- PepSpectra[unlist(cond)]
        
        return(CurSpec)
        
    })
    
    CurPepUptake <- reactive({
        
        match(input$varPep,peptide.identifications[,3])
        
    })
    
    peaks <- reactive({
        peakPick(CurSpec()[[1]],SNR=input$SNR)
    })
    
    Interp <- reactive({
        if(length(peaks()) > 1){
            linearInterp(peaks())
        }
    })
    
    Width <- reactive({
        if(length(peaks()) > 1){
            widthFinder(peaks(),CurSpec()[[1]],(input$BPI)/100)
        }
    })
    
    Centroid <- reactive({
        centroidCalc(Width(),CurSpec()[[1]])
    })
    
    CentTable <- reactive({
        data.frame(cbind(Centroid = Centroid(),Width = (Width()[2]-Width()[1])))
    })
    
    
    
    
    
    
    
    
    
    
    output$table <- renderDataTable({
        peptide.identifications},rownames=T)
    
    output$tableCent <- renderTable({
        CentTable()
    },bordered = T)
    
    output$MEMHDXTable <- renderDataTable(server=FALSE,{
        datatable(MEMTable$data, extensions = 'Buttons'
                  , options = list( 
                      dom = "Blfrtip"
                      , buttons = 
                          list( list(
                              extend = "collection"
                              , buttons = c("csv")
                              , text = "Download"
                          ))),rownames = FALSE)%>%
            formatStyle(columns=9,backgroundColor = styleEqual(levels=NA,values = 'red'))},rownames=T
    )
    
    output$plot <- renderPlot({
        
        plot(CurSpec()[[1]], main = heading())
        points(peaks(),col="red",pch=4)
        plotLin(Interp(),peaks())
        if(length(Width()>1)){
            plotWidth(Width())
        }
        centroidPlot(Centroid(),CurSpec()[[1]])
        
    },height = 400, width = 400)
    
    
    output$plotUptake <- renderPlot({
        
        
        PlotUptakeCompare(CurPepUptake(),all.cents = AllCentReact$data[[CurPepUptake()]],times=TP)
        
    })
    
    output$Testtable <- renderDataTable({
        datatable(AllCentReact$data[[1]][[4]])
    })
    
    output$video <- renderUI({
        tags$iframe(width="560", height="315", src="https://www.youtube.com/embed/DOClPhUJWcY",frameborder="0", allow="autoplay; encrypted-media")
    })
    
}

# Create Shiny app ----
shinyApp(ui, server)