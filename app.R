library(shiny)
library(MALDIHDX)
library(DT)
library(shinyjs)

options(shiny.maxRequestSize = 30*1024^2)


peptide.identifications <- import.identifications()

TP <- unique(pep10[,2])

defSNR = 5

DefaultAllCents <- lapply(peptide.features,function(x) mainCentNewMod2(x))
DefMEMTable <- MEMHDXall2(DefaultAllCents)

# Define UI
ui <- fluidPage(
    
    useShinyjs(),
    
    # App title
    titlePanel("MALDIHDX"),
    
    # Main panel for displaying outputs ----
    mainPanel(
        
        # Output: Tabset
        tabsetPanel(type = "tabs",
                    
                    tabPanel("Data Import",
                             
                           sidebarLayout(
                               
                             sidebarPanel(
                                 
                                 helpText("Please upload a Zip file containing
                                          Scaffold file and Mass Spectra"),
                                 
                                 
                                 fileInput("file1", "Choose folder",
                                           multiple = TRUE, buttonLabel = "Browse...",placeholder = "No file selected"
                                       )
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
                                        
                                        sliderInput('BPI','% Base Peak Intensity',min=0.1,max=100,
                                                    value = 50, step =1),
                                        actionButton("resBPI","Reset to default"),br(),br(),hr(),br(),
                                        
                                        fluidRow(column(12,
                                                        
                                                        tableOutput("tableCent"),actionButton("expCent","Export centroid to output table")  
                                                        
                                                        ))
                                 )
                             )
                    ),
                    
                    
                    tabPanel("Uptake Plots",
                             
                             dataTableOutput("table2")),
                    
                    
                    tabPanel("Output table",
                             
                             dataTableOutput("MEMHDXTable"))
                    
        )
    ) 
)



# Define server logic
server <- function(input, output) {
    
    
    curRow <- reactive({
        
        state <- switch(input$varBound,"Unbound" = "UNBOUND","Bound" = "BOUND")
        index <- which(c(DefMEMTable$Exposure == input$varTime & DefMEMTable$Replicate==input$varRep&DefMEMTable$State == state & DefMEMTable$Sequence==input$var))
        return(index)
    })
    
    
   

    
    
    heading <- reactive({
        paste(input$var,"Timepoint",input$varTime,"\n Replicate",
              input$varRep,"-",input$varBound)
    })
    
  
    observeEvent(input$resSNR, {
        
        reset("SNR")
            
    })
    
    observeEvent(input$resBPI, {
        
        reset("BPI")
        
    })

    
    observeEvent(input$expCent, {
        
        
        
        
    })
    
    
    
    
    CurSpec <- reactive({
        
        PepNo <- match(input$var,peptide.identifications[,3])
        
        PepSpectra <- L[[PepNo]]
        
        state <- switch(input$varBound,"Unbound" = "A","Bound" = "B")
        
        stateRep <- paste(state,input$varRep,sep="")
        
        cond <- lapply(PepSpectra,function(x)
            metaData(x)$StateRep == stateRep & metaData(x)$time == input$varTime )
        
        CurSpec <- PepSpectra[unlist(cond)]
        
        return(CurSpec)
        
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
    
    output$MEMHDXTable <- renderDataTable({
        datatable(DefMEMTable)%>%
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
    
}

# Create Shiny app ----
shinyApp(ui, server)