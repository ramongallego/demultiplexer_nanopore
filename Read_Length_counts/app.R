#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(see)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Length Distribution and demultiplexing success"),

    # Fluid Row layout
    sidebarLayout(
        sidebarPanel(
               fileInput(inputId = "length.file",
                         label= "Select plate length file",
                        multiple = F),
               fileInput(inputId = "well.file",
                         label= "Select Well length file",
                         multiple = F),
               
        ),

        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(
          tabPanel("First Step demultiplex",
                   plotOutput("Histogram"),
                   
                   plotOutput("Barplot"))           ,
          tabPanel("Second Step",
                   plotOutput("Well.Histogram"),
                   
                   plotOutput("Well.Barplot")),
          tabPanel("Well Heatmap",
                   plotOutput("Well.Heatmap")
           
          
        )))
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  options(shiny.maxRequestSize=300*1024^2) 
  
    output$Histogram <- renderPlot({
        # Generate Lengths histograms
      path <- input$length.file
      basedir <- dirname(path$datapath)
      df <- read_table(path$datapath,
                                col_names = c("seq","length", "plate")) %>%
                       mutate(plate = str_remove(plate, "_round1.fastq"),
                              plate = basename(plate))


         df %>%
           ggplot ( aes(y = length, x = plate, fill = plate)) +
           geom_violindot(fill_dots = "black") +
           theme_modern()+
           scale_fill_material_d() +
           # facet_wrap(~version, nrow = 2) +
           guides(fill = "none")
      })

      output$Barplot <- renderPlot({
         path <- input$length.file
           df <- read_table(path$datapath,
                       col_names = c("seq","length", "plate")) %>%
            mutate(plate = str_remove(plate, "_round1.fastq"),
                   plate = basename(plate))
        df %>%
          ggplot( aes(x = plate, fill= length>1000)) +
          theme_modern(legend.position=c(.9,.75))+
          geom_bar(position = "stack")
      })

       output$Well.Histogram <- renderPlot({
         path <- input$well.file
          df <- read_table(file = path$datapath,
                          col_names = c("seq","length", "plate")) %>%
           mutate(plate = basename(plate)) %>%
           separate(plate, into = c("plate", "well"), sep = "_Well_")


         df %>%
           ggplot ( aes(y = length, x = plate, fill = plate)) +
           geom_violindot(fill_dots = "black") +
           theme_modern()+
           scale_fill_material_d() +
           # facet_wrap(~version, nrow = 2) +
           guides(fill = "none")

       })

       output$Well.Barplot <- renderPlot({
         path <- input$well.file
         df <- read_table(file = path$datapath,
                          col_names = c("seq","length", "plate")) %>%
           mutate(plate = basename(plate)) %>%
           separate(plate, into = c("plate", "well"), sep = "_Well_")

         df %>%
           ggplot( aes(x = plate, fill= length>1000)) +
           theme_modern(legend.position=c(.9,.75))+
           geom_bar(position = "stack")
       })
       output$Well.Heatmap <- renderPlot({
         path <- input$well.file
         df <- read_table(file = path$datapath,
                          col_names = c("seq","length", "plate")) %>%
           mutate(plate = basename(plate)) %>%
           separate(plate, into = c("plate", "well"), sep = "_Well_") %>% 
           separate(well, into = c("Row", "Column"), sep = 1, convert = T) %>% 
           mutate(Column = str_extract(Column, "^[:digit:]+"),
                  Column = fct_reorder(Column, as.numeric(Column))) %>% 
           group_by(plate, Row, Column) %>% 
           tally(name = "nReads")
         
         df %>%
           ggplot( aes(x = Column, y = Row, fill= nReads)) +
           theme_modern()+
           geom_raster()+
           facet_wrap(~plate)+
           scale_fill_viridis_c() +
           scale_y_discrete(limits = rev)
       })
    }


# Run the application 
shinyApp(ui = ui, server = server)
