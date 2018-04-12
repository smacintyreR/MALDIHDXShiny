library(shiny)
library(MALDIHDX)
library(DT)


peptide.identifications <- import.identifications()
TP <- unique(pep10[,2])

Pep10Spectra <- DFtoSpec(pep10)


# Define UI for random distribution app ----
ui <- fluidPage(

  # App title ----
  titlePanel("MALDIHDX"),

  # Sidebar layout with input and output definitions ----


    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("Identifications",dataTableOutput("table")),
                  tabPanel("Centroid Plots",
                        
                        
                               sidebarPanel(
                                   helpText("Investigate centroid plots for bound/unbound peptides"),
                                   
                                   selectInput("var", 
                                               label = "Select your peptide",
                                               choices = c("Peptide 1", "Peptide 2",
                                                           "Peptide 3", "Peptide 4"),
                                               selected = "Peptide 1"),
                                   
                                   selectInput("varRep", 
                                               label = "Select Replicate",
                                               choices = c("1","2","3"),
                                               selected = "1"),
                                   
                                   selectInput("varTime", 
                                               label = "Select your timepoint",
                                               choices = TP,
                                               selected = TP[1]),
                                   
                                   selectInput("varBound", 
                                               label = "Select State",
                                               choices = c("Bound", "Unbound"),
                                               selected = "Unbound")
                               ),
                           
                           mainPanel(
                               
                               textOutput("text"),
                               
                               plotOutput("plot")
                           )
                               
                               
                             
                           
                           
                           
                           
                           
                           
                           
                           ),
                  tabPanel("Uptake Plots",dataTableOutput("table2"))
      
                  )
    )
  
)
  


# Define server logic for random distribution app ----
server <- function(input, output) {
    


    
    CurSpec <- reactive({
        
   
        
        state <- switch(input$varBound,"Unbound" = "A","Bound" = "B")
        
        stateRep <- paste(state,input$varRep,sep="")
        
        cond <- lapply(Pep10Spectra,function(x) metaData(x)$StateRep == stateRep & metaData(x)$time == input$varTime )
        
        CurSpec <- Pep10Spectra[unlist(cond)]
        
        return(CurSpec)
        
    })
    
    
    
    output$table <- renderDataTable({
       peptide.identifications
    },rownames=T)
    
    output$plot <- renderPlot({
    
     
        plot(CurSpec()[[1]])
 
   
    },height = 600, width = 600)
    
    
    output$table2 <- renderDataTable({
       
    },rownames=F)
    
    output$text <- renderText({ 
        paste(input$var,"Timepoint",input$varTime,"Replicate",input$varRep,"-",input$varBound)
    })

}

# Create Shiny app ----
shinyApp(ui, server)
