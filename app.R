library(shiny)
library(stringr)
library(MALDIHDX)
library(DT)
library(shinyjs)
library(ggplot2)
library(shinythemes)
library(MALDIquantForeign)
library(reshape2)

options(shiny.maxRequestSize = *1024^2)

defSNR = 5

# Define UI
ui <- fluidPage(theme=shinytheme("cerulean"),
                
                useShinyjs(),
                
                navbarPage(
                    
                    strong("MALDIHDX"),
                    
                    tabPanel(icon=icon("home"),"About",
                             
                             strong(h3("MALDIHDX: Semi-automated centroid analysis for HDX-MS data",align="center")),
                             em(h5("developed by Sam MacIntyre and Thomas Nebl (CSIRO)",align="center")),
                             div(em("contact us at:", a("sam.macintyre@csiro.au")),align="center"),
                             br(),
                             
                             p("This tool allows users to perform a semi-automated workflow to analyse, validate and visualize large HDX-MS datasets.
                               The input file is a zipped folder containing the ",strong("Raw MS Data"),"and an",strong("Identifications File"),
                               "(which has specified required data fields) Output includes centroid plots for each sample (with manual editing based on available parameters),
                               comparative deuterium uptake plots for each peptide and a downloadable",tags$a(href="http://memhdx.c3bi.pasteur.fr/", "MEMHDX")," compatible table.",align="center",style="font-family: verdana",style="font-size: 40%"),
                             
                             img(src="csiro-logo.jpg",height=60,width=60,style="display: block; margin-left: auto; margin-right: auto;"), br(),
                             
                             fluidRow(column(5,strong(h4("A. Centroid Plot")),img(src="Capture.PNG",height=350,width=400),
                                             
                                             fluidRow(column(12,strong(h4("B. Uptake Plot")),img(src="Capture2.PNG",height=350,width=400)))
                                             
                                             ),column(7,strong(h4("C. MEMHDX compatible Output Table")),img(src="Capture3.PNG",height=500,width=800)))
                             
                             ),
                    
                    
                    tabPanel(icon=icon("clipboard"),"Tutorial",
                             
                             fluidRow(column(12,h3("How to use MALDIHDX"),
                                             align="center")),
                             
                             hr(),
                             
                             fluidRow(column(12,img(src = 
                                                        "TutorialGraphicDraft.png",height=260,width=750)
                                             ,align="center")),hr(),
                             
                             fluidRow(column=12,htmlOutput("video"),
                                      align="center")
                             
                    ),
                    
                    
                    
                    
                    
                    tabPanel(icon=icon("cogs"),"Analysis",
                             
                             # App title
                             h3("Centroid calculation and 
                                        validation of HDX-MS Experiments"),
                             
                             # Main panel for displaying outputs ----
                             mainPanel(
                                 
                                 # Output: Tabset
                                 tabsetPanel(id="tabs",
                                             
                                             tabPanel("Data Import",br(),
                                                      
                                                      
                                                      
                                                      
                                                      sidebarPanel(
                                                          
                                                          
                                                          fileInput
                                                          ("FileInput", "Choose folder",
                                                              
                                                              multiple = TRUE, buttonLabel = "Browse...",placeholder = "No file selected"
                                                          ),
                                                          
                                                          
                                                          actionButton("FileImport","Import all Spectra"),
                                                          
                                                          downloadLink("TESTDATA",label="  Download test data set here")
                                                      ),
                                                      
                                                      
                                                      
                                                      mainPanel(
                                                          
                                                          helpText(p(strong("
                                                              MALDIHDX")," requires a",strong("zip file")," containing raw mass spectrometry data with the following format:",br(),br(),
                                                                  div(em("STATE_TIMEPOINT_REPLICATE e.g. A_300_1, A_300_2, B_3600_3"),style="color:blue"),br(),br(),
                                                              "Note that ",span("STATE",style="color:blue")," can only have value A or B currently where A = Untreated/Unbound and B = Treated/Bound",br(),br(),
                                                              span("TIMEPOINT",style="color:blue")," should be in minutes to conform to the MEMHDX input",br(),br(),
                                                              "Secondly, a",strong(" .csv (comma separated values)")," file is required with the following fields:"),br(),

                                                                  strong("Sequence:")," Peptide sequence",br(),br(),
                                                              strong("Observed:")," Observed monoisotopic mass of peptide",br(),br(),
                                                              strong("Charge:")," Peptide charge (MALDIHDX currently only handles singly charged peptides)",br(),br(),
                                                              strong("Start:")," Peptide Start position on the protein",br(),br(),
                                                              strong("Stop:")," Peptide End position on the protein
                                                              "))
                                                      
                                                  )
                                                      
                                                      
                                             
                                             
                                             ,  
                                             
                                             tabPanel("Identifications"
                                                      ,
                                                       
        
                                                          DT::dataTableOutput(
                                                              "table")),
                                             
                                             tabPanel("Centroid Plots",br(),
                                                      
                                                      fluidRow(
                                                          
                                                          column(8,
                                                                 
                                                                 
                                                                 plotOutput("plot"),
                                                                 hr(),
                                                                 
                                                                 
                                                                 fluidRow(column(6,
                                                                                 
                                                                                 
                                                                                 
                                                                                 selectInput("var", label = "Select peptide",
                                                                                             
                                                                                             ""),
                                                                                 
                                                                                 
                                                                                 
                                                                                 selectInput("varRep", label = "Select Replicate",
                                                                                             
                                                                                             "")),
                                                                          
                                                                          
                                                                          
                                                                          column(6,
                                                                                 
                                                                                 
                                                                                 
                                                                                 selectInput("varTime", label = "Select timepoint",
                                                                                             
                                                                                             ""),
                                                                                 
                                                                                 
                                                                                 
                                                                                 selectInput("varBound", label = "Select State",
                                                                                             
                                                                                             choices = c("Bound", "Unbound"),
                                                                                             
                                                                                             selected = "Unbound"))
                                                                          
                                                                          
                                                                          
                                                                          
                                                                 )
                                                          ),
                                                          
                                                          column(4,
                                                                 
                                                                 h4(
                                                                     "Centroid Parameters"), hr(),
                                                                 
                                                                 
                                                                 sliderInput('SNR','Signal to Noise',min=1,max=10,
                                                                             
                                                                             value = 5, step =0.1),
                                                                 
                                                                 actionButton("resSNR","Reset to default"),br(),br(),
                                                                 
                                                                 
                                                                 sliderInput('BPI','% Base Peak Intensity',min=0.1,max=95,
                                                                             
                                                                             value = 50, step =1),
                                                                 
                                                                 actionButton("resBPI","Reset to default"),br(),br(),hr(),br(),
                                                                 
                                                                 
                                                                 fluidRow(column(12,
                                                                                 
                                                                                 
                                                                                 
                                                                                 tableOutput("tableCent"),actionButton("expCent","Export centroid to 
                                                                                                                       output table"),span(textOutput("Done"), style="color:blue")  
                                                                                 
                                                                                 
                                                                                 ))
                                                          )
                                                      )
                                                          ),
                                             
                                             
                                             tabPanel("Uptake Plots",br(),
                                                      
                                                      sidebarLayout(
                                                          
                                                          sidebarPanel(
                                                              
                                                              helpText(
                                                                  "Select peptide to view Deuterium uptake plots"),
                                                              
                                                              
                                                              selectInput("varPep", label = "Select peptide",
                                                                          
                                                                          "")
                                                              
                                                          ),
                                                          
                                                          mainPanel(
                                                              plotOutput(
                                                                  "plotUptake"))
                                                      ) 
                                                      
                                                      
                                             ),
                                             
                                             
                                             tabPanel("Output table",br(),
                                                      
                                                      DT::dataTableOutput(
                                                          "MEMHDXTable")
                                             )
                                             
                                             
                                             
                                             )
                                 ) 
                             
                             ),tabPanel(icon=icon("book"),"Publications",h3("MALDIHDX publications",
                                                                            align="center"),hr(),
                                        h4("Papers"),br(),
                                        p("MEMHDX paper from The Institut Pasteur: "),
                                        strong("MEMHDX: An interactive tool to expedite the statistical validation and visualization of large HDX-MS datasets"),
                                        em("Hourdel V, Volant S, O'Brien DP, Chenal A, Chamot-Rooke J, Dillies MA, Brier S, 2016 Jul 13"),tags$a(href="https://watermark.silverchair.com/btw420.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAcUwggHBBgkqhkiG9w0BBwagggGyMIIBrgIBADCCAacGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQM-xmiVS26AhuCBnpBAgEQgIIBeH7hjPdBfnz4wWufJ6NhyBYWf1cXxUS_OwyAXrwFnaW6TcZXBlcsUmwoKgQe8QQLkCiFse2HvU-wEZ_upPkpb1rwYG_RKkJt9XA-MWziaanlW7m2mINqgw8IafYuYTh_mIqhdkjPH8IHctg4CQnV90aUqtppjWLGoWbvBpFbWCphg7CMBlaq7uh-ZNPFGpe9v-IM5qK2HushIXnP5srzazPcwydeFZ-wvwDi-pz4yXELu-fjRTc4lIg-haOgyIHfu1hToX65p9gX8QuKuL3g70BxLulsYGcyc2DX_52KdZSLOBmKCTm3sCl_NtOwQhc0ivqMfHp9YrIa_rQYjqJfs9AEXo71G5AqToBBnFDGewqd07yDvMbpO5cpzEcgwBMBZ__LyEDZDEPqKi_-3eMpDbzmmeEruBnK87i_AbpglRUK3i-MWZGyzhQbVZz_uQjH9rwev7VxUDZTURkaVa_miKHle40zHQBzbYut-y9LwbIpCmJDt_1Uh-I", "Full Text"),
                                        br(),br(),h4("Posters and Presentations"),"- IMSC, 2018:",br(),strong("Development and validation of semi-automatic MALDI-HDX sample preparation and data analysis tools"),br(),em("S. Macintyre, T. Nebl")),
                    tabPanel(icon=icon("star"),"News",h3("Latest news in MALDIHDX",hr(),align="center")
                             
                             
                             
                             ))
                )



# Define server logic
server <- function(input, output,session) {
    
    FLAG <- reactiveValues(data=FALSE)
    
    
    # Pre-hides computation based tabs to prevent user seeing null related
    # errors
    hideTab(inputId = "tabs", target = "Centroid Plots")
    hideTab(inputId = "tabs", target = "Uptake Plots")
    hideTab(inputId = "tabs", target = "Identifications")
    hideTab(inputId = "tabs", target = "Output table")
    
    # Shows computation-dependent tabs after computation has been flagged as 
    # complete 
    observe({
        
        if(FLAG$data==TRUE){
            showTab(inputId = "tabs", target = "Centroid Plots")
            showTab(inputId = "tabs", target = "Uptake Plots")
            showTab(inputId = "tabs", target = "Identifications")
            showTab(inputId = "tabs", target = "Output table")
        }
        
        else{
            return()
        }
    }
    )
    
    
    
    
    
    
    # NULL initialises reactiveValues objects for important variables
    peptide.features <- reactiveValues()
    AllCentReact <- reactiveValues()
    MEMTable <- reactiveValues() 
    peptide.identifications <- reactiveValues()
    TP <- reactiveValues()
    
    
    
    # Handles import of spectra and identifications.
    # Also initialises peptide.features, the reactive centroid table and MEMHDX 
    # output table
    observeEvent(input$FileImport,
                 
                 
                 {
                     
                     
                     withProgress(message="Importing and analysing data...",
                                  value=0,{
                                    
                                      
                                      peptide.identifications$data <- 
                                          import.identifications(path=paste("data/",list.files("data"),sep=""))
                                      
                                      
                                      peptide.features$data <- importNew(ids=
                                                                             peptide.identifications$data,path=paste("data/",list.files("data"),sep=""))
                                      
                                      incProgress(1/2,"Calculating centroids based on 
                                                  default parameters..")
                                      AllCentReact$data <- lapply(peptide.features$data  ,
                                                                  function(x) mainCentNewMod2(x))
                                      
                                      MEMTable$data <- MEMHDXall2(AllCentReact$data,Idents=
                                                                      peptide.identifications$data)
                                      
                                      TP$data <- as.numeric(unique(peptide.features$data[[1
                                                                                          ]][,2]))
                                      
                                      
                                   
                                      
                                      incProgress(1/2,"Complete")
                                      FLAG$data <- TRUE
                                  })
                     
                 }
                     )
    
    # Reactive variable which takes peptide.feature's list of matrices format 
    # and converts it to a list of spectra
    L <- reactive({
        
        lapply(peptide.features$data,function(x) DFtoSpec(x))
        
    })
    
    # Handles .zip file upload and unzipping in "data" folder
    observeEvent(input$FileInput,
                 {
                     
                     
                     file.remove(paste("data/",list.files("data"),sep=""))
                     unlink(paste("data/",list.files("data"),sep=""),recursive = T)
                     infile <- input$FileInput
                     if(is.null(infile))
                         return(NULL)
                     unzip(infile$datapath,exdir = "data")
                     
                 }   
    )
    
    # Reactive index to track the corresponding row of MEMHDX table 
    curRow <- reactive({
        
        state <- switch(input$varBound,"Unbound" = "UNBOUND","Bound" = "BOUND")
        index <- which(c(MEMTable$data$Exposure == input$varTime & MEMTable$
                             data$Replicate==input$varRep&MEMTable$data$State == state & MEMTable$
                             data$Sequence==input$var))
        return(index)
    })
    
    # Reactive heading to be used in output plots
    heading <- reactive({
        paste(input$var,"Timepoint",input$varTime,"\n Replicate",
              input$varRep,"-",input$varBound)
    })
    
    # Resets "SNR" slider if the current spectrum is modified
    observeEvent(c(input$var,input$varRep,input$resSNR,input$varBound,input$
                       varTime), {
                           
                           reset("SNR")
                           
                       })
    
    # Resets "BPI" slider if the current spectrum is modified
    observeEvent(c(input$var,input$varRep,input$resBPI,input$varBound,input$
                       varTime), {
                           
                           reset("BPI")
                           output$Done <- renderText("Export manual centroid") 
                           
                       })
    
    # Modifies MEMHDX output table with an updated centroid if "Export Centroid"
    # is activated
    observeEvent(input$expCent,
                 
                 {
                     temp <- MEMTable$data
                     temp[curRow(),9] <- Centroid()
                     MEMTable$data <- temp
                     text <- paste("Row",as.character(curRow()),"Centroid was exported")
                     output$Done <- renderText(text)
                     
                     
                     PepNumber <- match(input$var,peptide.identifications$data[,4])
                     stateUp <- switch(input$varBound,"Unbound" = "A","Bound" = "B")
                     stateRepUp <- paste(stateUp,input$varRep,sep="")
                     
                     subNo <- switch(stateRepUp,"A1"=1,"A2"=2,"A3"=3,"B1"=4,"B2"=5,"B3"=6)
                     
                     
                     tempCent <- AllCentReact$data[[PepNumber]][[subNo]]
                     tempCent[tempCent$'time (min)'==input$varTime,2] <- Centroid()
                     AllCentReact$data[[PepNumber]][[subNo]] <- tempCent
                 } 
    )
    
    # Selects the correct spectrum based on input parameters (Bound/Unbound,
    # Timepoint, Peptide etc.)
    CurSpec <- reactive({
        
        PepNo <- match(input$var,peptide.identifications$data[,4])
        
        PepSpectra <- L()[[PepNo]]
        
        state <- switch(input$varBound,"Unbound" = "A","Bound" = "B")
        
        stateRep <- paste(state,input$varRep,sep="")
        
        cond <- lapply(PepSpectra,function(x)
            metaData(x)$StateRep == stateRep & metaData(x)$time == input$
                varTime)
        
        CurSpec <- PepSpectra[unlist(cond)]
        
        return(CurSpec)
        
    })
    
    # Reactive index to denote current Uptake plot 
    CurPepUptake <- reactive({
        
        match(input$varPep,peptide.identifications$data[,4])
        
    })
    
    
    # Reactive peaks list
    peaks <- reactive({
        peakPick(CurSpec()[[1]],SNR=input$SNR)
    })
    
    
    # Reactive Linear Interpolation between peaks
    Interp <- reactive({
        if(length(peaks()) > 1){
            linearInterp(peaks())
        }
    })
    
    
    # Reactive width value
    Width <- reactive({
        if(length(peaks()) > 1){
            widthFinder(peaks(),CurSpec()[[1]],(input$BPI)/100)
        }
    })
    
    # Reactive centroid value
    Centroid <- reactive({
        centroidCalc(Width(),CurSpec()[[1]])
    })
    
    
    # Reactive centroid table displaying width and centroid
    CentTable <- reactive({
        data.frame(cbind(Centroid = Centroid(),Width = (Width()[2]-Width()[1
                                                                           ])))
    })
    
    
    
    
    # Allows timepoint input list to be dependent on user input data
    observe({
        
        updateSelectInput(
            session,
            "varTime",
            choices=TP$data)
        
    })
    
    # Allows petide input list to be dependent on user input data
    observe({
        updateSelectInput(
            session,
            "var",
            choices=peptide.identifications$data[,4])
        
    })
    
    
    # Allows petide input list to be dependent on user input data
    observe({
        updateSelectInput(
            session,
            "varPep",
            choices=peptide.identifications$data[,4])
        
    })
    
    observe({
        updateSelectInput(
            session,
            "varRep",
            choices=unique(MEMTable$data[,8]))
        
    })
    
    
    
    # Output: Peptide identifications table
    output$table <- DT::renderDataTable({
        peptide.identifications$data},rownames=T)
    
    
    # Output: Table displaying current centroid value and width
    output$tableCent <- renderTable({
        CentTable()
    },bordered = T)
    
    
    
    # Output: MEMHDX formatted table
    output$MEMHDXTable <- DT::renderDataTable(server=FALSE,{
        datatable(MEMTable$data, extensions = 'Buttons'
                  , options = list( 
                      dom = "Blfrtip"
                      , buttons = 
                          list( list(
                              extend = "collection"
                              , buttons = c("csv")
                              , text = "Download"
                          ))),rownames = FALSE)%>%
            formatStyle(columns=9,backgroundColor = styleEqual(levels=NA,
                                                               values = 'red'))},rownames=T
    )
    
    
    # Output: Main centroid plot
    output$plot <- renderPlot({
        
        plot(CurSpec()[[1]], main = heading())
        points(peaks(),col="red",pch=4)
        plotLin(Interp(),peaks())
        if(length(Width()>1)){
            plotWidth(Width())
        }
        centroidPlot(Centroid(),CurSpec()[[1]])
        
    },height = 400, width = 400)
    
    # Output: Uptake plots
    output$plotUptake <- renderPlot({
        
        PlotUptakeCompare(CurPepUptake(),all.cents = AllCentReact$data[[
            CurPepUptake()]],times=TP,pep.ids = peptide.identifications$data)
        
    })
    
    # Tutorial video: embedded from Youtube
    output$video <- renderUI({
        tags$iframe(width="560", height="315", src="https://www.youtube.com/embed/3g9-RZWAm9w", frameborder="0", allow="autoplay; encrypted-media")
    })
    
    output$TESTDATA <- downloadHandler(
        
        filename=function(){
            paste("MALDIHDXtest-data","zip",sep=".")
        },
        
        content = function(file){
            file.copy("TESTDATA/test-data5.zip",file)
        },
        contentType = "application/zip"
    )
    
}

# Create Shiny app ----
shinyApp(ui, server)
