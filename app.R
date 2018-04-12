library(shiny)
library(MALDIHDX)
library(DT)


peptide.identifications <- import.identifications()
TP <- unique(pep10[,2])




# Define UI
ui <- fluidPage(
    
    # App title
    titlePanel("MALDIHDX"),
    
    # Main panel for displaying outputs ----
    mainPanel(
        
        # Output: Tabset
        tabsetPanel(type = "tabs",
                    
                    tabPanel("Identifications",
                             
                             dataTableOutput("table")),
                    
                    tabPanel("Centroid Plots",
                             
                             fluidRow(
                                 
                                 column(8,
                                        
                                        plotOutput("plot"),
                                        hr(),
                                        
                                        fluidRow(column(6,
                                                        
                                                        selectInput("var", label = "Select your peptide",
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
                                                    value = 5, step =0.1),br(),
                                        
                                        sliderInput('BPI','% Base Peak Intensity',min=0.1,max=100,
                                                    value = 50, step =1)
                                 )
                             )
                    ),
                    
                    
                    tabPanel("Uptake Plots",
                             
                             dataTableOutput("table2"))
                    
        )
    ) 
)



# Define server logic
server <- function(input, output) {
    
    
    heading <- reactive({
        paste(input$var,"Timepoint",input$varTime,"\n Replicate",
              input$varRep,"-",input$varBound)
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
    
    
    
    output$table <- renderDataTable({
        peptide.identifications},rownames=T)
    
    
    output$plot <- renderPlot({
        
        plot(CurSpec()[[1]], main = heading())
        
    },height = 400, width = 400)
    
}

# Create Shiny app ----
shinyApp(ui, server)