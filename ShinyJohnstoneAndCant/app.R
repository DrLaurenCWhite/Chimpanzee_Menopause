library(shiny)
library(data.table)
library(ggplot2)
source("./GetRel.R")

ui <- fluidPage(
  fluidRow(
    tags$h3(tags$strong("Relatedness Asymmetries Over Age"), style='padding-left: 15px'),
    tags$p(" As presented in Johnstone and Cant (2010)",
           tags$a(href="https://royalsocietypublishing.org/doi/10.1098/rspb.2010.0988", "https://doi.org/10.1098/rspb.2010.0988"),
           style='padding-left: 15px'),
  ),
  fluidRow(
    column(4,
           wellPanel(sliderInput(inputId="DF", label="Female dispersal rate (df)", value=0.85, min=0, max=1),
                     sliderInput(inputId="DM", label="Male dispersal rate (dm)", value=0.15, min=0, max=1),
                     sliderInput(inputId="NF", label="Number of adult female in group (nf)", value=3, min=1, max=100),
                     sliderInput(inputId="NM", label="Number of adult male in group (nm)", value=3, min=1, max=100),
                     sliderInput(inputId="u", label="u", value=0.1, min=0, max=1),
                     sliderInput(inputId="m", label="Proportion of local breeding (m)", value=1, min=0, max=1),
                     #selectInput(inputId="m", label="Proportion of local breeding (m)", selected=1, choices=c(0,1)),
                     actionButton("reset1", "Ape Case", width="100%"),
                     actionButton("reset2", "Whale Case", width="100%"),
                     actionButton("reset3", "Ngogo Case - Community", width="100%"),
                     actionButton("reset4", "Ngogo Case - Subgroup", width="100%"),
           )),
    column(8,
           plotOutput("lines"),
           tags$br(),
           tags$br(),
           plotOutput("CB"))
  )
)

?plotOutput

server <- function(input, output, session) {
  
  observeEvent(input$reset1,{
    updateSliderInput(session, 'DF',value = 0.85)
    updateSliderInput(session, 'DM',value = 0.15)
    updateSliderInput(session, 'NF',value = 3)
    updateSliderInput(session, 'NM',value = 3)
    updateSelectInput(session, 'm', selected = 1)
    updateSliderInput(session, 'u',value = 0.1)
  })
  
  observeEvent(input$reset2,{
    updateSliderInput(session, 'DF',value = 0.15)
    updateSliderInput(session, 'DM',value = 0.15)
    updateSliderInput(session, 'NF',value = 3)
    updateSliderInput(session, 'NM',value = 3)
    updateSelectInput(session, 'm', selected = 0)
    updateSliderInput(session, 'u',value = 0.1)
  })
  
  observeEvent(input$reset3,{
    updateSliderInput(session, 'DF',value = 0.5)
    updateSliderInput(session, 'DM',value = 0)
    updateSliderInput(session, 'NF',value = 66)
    updateSliderInput(session, 'NM',value = 44)
    updateSelectInput(session, 'm', selected = 1)
    updateSliderInput(session, 'u',value = 0.076)
  })
  
  observeEvent(input$reset4,{
    updateSliderInput(session, 'DF',value = 0.9)
    updateSliderInput(session, 'DM',value = 0.1)
    updateSliderInput(session, 'NF',value = 18)
    updateSliderInput(session, 'NM',value = 12)
    updateSelectInput(session, 'm', selected = 1)
    updateSliderInput(session, 'u',value = 0.076)
  })
  
  Rel <- reactive( { 
    fox=GetRel(input$DF, input$DM, input$NF, input$NM, as.numeric(input$m), input$u)
    fox} )
  
  Cost <- reactive( {
    df <- as.data.table(Rel())
    trot=GetCB(input$DF, input$DM, input$NF, input$NM, as.numeric(input$m), df)
    trot
  })
  
  output$lines <- renderPlot( {
    df <- as.data.frame(Rel())
    
    cols = c("RFF"="#CC0C00FF", "RMM"="#5C88DAFF", "RFM"="#FFCD00FF", "RMF"="#84BD00FF")
    p=ggplot(data=df, aes(x=Ns))
    p<- p + geom_line(aes(x=Ns, y=rff, col="RFF"), lwd=2) + geom_line(aes(x=Ns, y=rmf, col="RMF"), lwd=2)+
      geom_line(aes(x=Ns, y=rmm, col="RMM"), lwd=2)+
      geom_line(aes(x=Ns, y=rfm, col="RFM") , lwd=2) +
      scale_colour_manual(name="Key:", values=cols, labels=c(
        "Average relatedness of a female of age x to all other adult females",
        "Average relatedness of a male of age x to all adult females",
        "Average relatedness of a female of age x to all adult males",
        "Average relatedness of a male of age x to all other adult males"
        
      )) +
      guides(col=guide_legend(nrow=4,byrow=TRUE)) +
      xlab("age (relative to mean generation time)") + ylab("Relatedness") +
      theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
            panel.background = element_blank(),
            axis.line = element_line(colour="black"),
            panel.grid.major = element_line(colour="#f0f0f0"),
            axis.text = element_text(size=12),
            axis.title = element_text(size=15),
            legend.text = element_text(size=15),
            legend.title = element_blank(),
            legend.position = "bottom",
            legend.justification = "left")
    print(p)
  } )
  
  output$CB <- renderPlot( {
    df <- as.data.table(Cost())
    
    cols = c("Help"="pink", "Harm"="blue")
    p=ggplot(data=df, aes(x=Ns, y=res))
    if (nrow(df[res<0])>0 & nrow(df[res>0])>0){
      p <- p + geom_line(lwd=2) + xlab('age (relative to mean generation time)') + ylab("c/b") +
        geom_hline(yintercept=0, lwd=1) + 
        geom_area(data=df[res<0], aes(fill="Harm"), alpha=0.5) +
        geom_area(data=df[res>0], aes(fill="Help"), alpha=0.5) +
        scale_fill_manual(values=cols) +
        theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
              panel.background = element_blank(),
              axis.line = element_line(colour="black"),
              panel.grid.major = element_line(colour="#f0f0f0"),
              axis.text = element_text(size=12),
              axis.title = element_text(size=15),
              legend.text = element_text(size=15),
              legend.title = element_blank(),
              legend.position="bottom")
      print(p)
    }
    if (nrow(df[res<0])==0){
      p <- p + geom_line(lwd=2) + xlab('age (relative to mean generation time)') + ylab("c/b") +
        geom_hline(yintercept=0, lwd=1) + 
        geom_area(aes(fill="Help"), alpha=0.5) +
        scale_fill_manual(values=cols) +
        theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
              panel.background = element_blank(),
              axis.line = element_line(colour="black"),
              panel.grid.major = element_line(colour="#f0f0f0"),
              axis.text = element_text(size=12),
              axis.title = element_text(size=15),
              legend.text = element_text(size=15),
              legend.title = element_blank(),
              legend.position="bottom")
      print(p)
    }
    if (nrow(df[res>0])==0){
      p <- p + geom_line(lwd=2) + xlab('age (relative to mean generation time)') + ylab("c/b") +
        geom_hline(yintercept=0, lwd=1) + 
        geom_area(aes(fill="Harm"), alpha=0.5) +
        scale_fill_manual(values=cols) +
        theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
              panel.background = element_blank(),
              axis.line = element_line(colour="black"),
              panel.grid.major = element_line(colour="#f0f0f0"),
              axis.text = element_text(size=12),
              axis.title = element_text(size=15),
              legend.text = element_text(size=15),
              legend.title = element_blank(),
              legend.position="bottom")
      print(p)
    }
  } )
  
}

shinyApp(ui=ui, server=server)

