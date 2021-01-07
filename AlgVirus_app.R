######################################################################################################
rm(list = ls()) # clean the R history
######################################################################################################
#libraries that I need
######################################################################################################

library(tidyverse)
library(deSolve)       # solving the differential equation
library(ggplot2)       # ploting with ggplot2
library(shiny)

######################################################################################################
#the initial values of N=nutrients, C=Chlorella, P=parasite or virus
# C has five values. We assume that we have evolution from the type 1 to type 2, from type 2 to type 3
# from type 3 to type 4 and from type 4 to type 5. the type is a general reistant host which can not be infected
# by any virus. More info about the nature of the model come to the attached leaflet
######################################################################################################

yini= c(N=80, C=c(1,0,0,0,0), P=c(0.1,0,0,0))

######################################################################################################
#  PARAMETERS
######################################################################################################

Ni = 80          # nutrients
V = 0.33         # chemostat volume

###############################################
# algae-chlorella parameters
###############################################

Kc = 4.3         # Half saturation constant
Cm = 0.437       # Critical Chlorella concentration
xc = 0.05        # Algae conversion efficiency
Wc = 20.0        # N content in 10^9 Chlorella cells
e = 1.0          # Assimilation efficiency
epsilon=10^-3    # mutation rate
# algae evolutionary matrix

M_C =  matrix(c(1-epsilon, epsilon/2, 0,0,0,epsilon, (1-epsilon), (epsilon/2), 0, 0, 0, epsilon/2, (1-epsilon), epsilon/2, 0, 0, 0 , epsilon/2, (1-epsilon), epsilon,0,0,0, epsilon/2, 1-epsilon),
              nrow=5,
              ncol=5,
              byrow=TRUE) 


###############################################
# virus-parasite parameters
###############################################

phi = 9*10^-2    # virus adsorption rate
beta = 50        # virus burst size
# virus evolutioanry matrix
M_P =  matrix(c(1-epsilon, epsilon/2,0,0,epsilon, (1-epsilon),(epsilon/2), 0, 0, epsilon/2, (1-epsilon), epsilon, 0,0,  epsilon/2, (1-epsilon)),
              nrow=4,
              ncol=4,
              byrow=TRUE) 
# infection matrix
A= matrix(c(1,1,1,1,0,1,1,1,0,0,1,1,0,0,0,1,0,0,0,0),
          nrow=5,
          ncol=4,
          byrow=TRUE) 

######################################################################################################
# Function
######################################################################################################

pp <- function(t, yini, parameters) {
    with(as.list(c(yini, parameters)), {
        C=yini[paste("C",1:5, sep = "")]
        P=yini[paste("P",1:4, sep = "")]
        
        dN = D*(V*Ni - N) - sum(C*(bc*Wc*(N/V))/(e*(Kc + (N/V))))  #change in nutrients over time
        
        dC = M_C %*%(xc*C*bc*Wc*(N/V)/(e*(Kc + (N/V)))) - D*C  - (phi * (A %*% P))*C -d*C #changes in the alga types over time
        
        dP = (M_P * beta) %*% (phi*(t(A)%*%(C)*P)) - (phi*(t(A)%*%(C))*P) - D*P   #changes in the viral types over time
        
        return(list(c(dN,dC,dP)))  
    })
}

######################################################################################################
######################################################################################################
#     \  :  /
#  `. __/ \__ .'
#  _ _\     /_ _      # Construct the shiny app - ui= input$model, output$graph1, output$graph2 
#     /_   _\           - server = reactive(input$model), output$graph1, output$graph2, session
#   .'  \ /  `.
#     /  :  \    
######################################################################################################
######################################################################################################

ui <- fluidPage(                                                            # start the user interface
    titlePanel("Abiotic stress and Alage-Virus eco-evolutionary dynamics"),   # title- add the title of the user interface
    sidebarLayout(                                                            # I want to have a sidebar layout
        sidebarPanel(                                                           # with a side bar
            h3("Model options"),                                           # a third level header which is called model options
            h4("Initial values"),                                        # choose the initial values of the model
            textInput("C1", "Ancestral alga", value=1),                 # after your titles add your inputs 
            textInput("P1", "Ancestral virus", value=1),                # after your titles add your inputs 
            textInput("time","Number of steps to simulate", value=500), # thats your other input
            
            h4("Parameters"),
            sliderInput("D", label = "Dilution rate, D",                 # this time I do not use text but a sliderInput
                        min = 0.1, max = 0.8, value = 0.3, step = 0.1),  # one sliderInput
            sliderInput("d", label = "Host death rate, d",
                        min = 0, max = 0.6, value = 0, step = 0.05),     # second sliderInput
            sliderInput("bc", label= "Lowest cost of resitance",
                        min=0.5, max=0.75, value=0.7, step=0.02)),     # and like this you finish the sidebarLayout and panel
        mainPanel(
            plotOutput("dynamics"),    # first plot as an output
            plotOutput("tradeoff"))                          # second plot as an output
    ))                      


server <- function(input, output, session){
    
    
    dataInput <- reactive({              #input the data first
        yini= c(N=80, C=c(as.numeric(input$C1),0,0,0,0), P=c(as.numeric(input$P1),0,0,0)) # this input is connected with the user interface and affects the outcom
        bc <- sort(seq(as.numeric(input$bc), as.numeric(input$bc)+0.08,0.02), decreasing = TRUE)
        print(bc)
        
        parameters <- list(Ni, V, D=as.numeric(input$D), Kc, Cm, xc, bc=as.numeric(bc), Wc, e, epsilon, M_C, phi, beta, M_P, A,d=as.numeric(input$d))
        print(parameters)
        time <- seq(0, as.numeric(input$time), 0.1)
        
        ode(y = yini, times = time, func = pp, parms = parameters, method="lsoda", atol=10^-15,rtol=10^-15) 
        
    })
    
    ###################
    #plot the outputs
    ###################
    
    output$dynamics <- renderPlot({
        data <- as.data.frame(dataInput())
        algae=rowSums(data[,c("C1","C2","C3","C4","C5")])
        virus=rowSums(data[,c("P1","P2","P3","P4")])
        data=cbind(data,algae,virus)
        ggplot(data, aes(x=time)) +
            geom_line(aes(y=log(algae), color="algae"), size=1.5) +
            geom_line(aes(y=log(virus), color="virus"), size=1.5) +
            scale_color_manual(values=c("algae"= "#36BA45", "virus"= "#48BBFC"), name="species")+
            labs(x="time", y="log (population size)", title="Alga - Virus population dynamics") +
            theme_bw()+
            theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold", family="Avenir"), 
                  legend.title.align=0.5,
                  #panel.grid.major  = element_line(colour = "#F0F0F2", size=0.5),
                  axis.ticks.x = element_line(colour = "#333333"),
                  axis.ticks.y = element_line(colour = "#333333"),
                  axis.ticks.length =  unit(0.26, "cm"), 
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank(),
                  strip.background =element_rect(fill="white", colour="#333333"))
    })
    
    output$tradeoff <- renderPlot({
        data <- as.data.frame(dataInput())
        bc <- sort(seq(as.numeric(input$bc), as.numeric(input$bc)+0.08,0.02), decreasing = TRUE)
        trade.off <- as.data.frame(bc)
        type=c("C1","C2","C3","C4","C5")
        type=as.data.frame(type)
        data <- cbindX(data,trade.off,type)
        data <- data %>% 
            drop_na()
        ggplot(data, aes(type,bc,group=1)) +
            geom_point() +
            geom_line() +
            labs(x="alga types", y="growth rate", title="Trade off of the alga types") +
            theme_bw()+
            theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold", family="Avenir"), 
                  legend.title.align=0.5,
                  #panel.grid.major  = element_line(colour = "#F0F0F2", size=0.5),
                  axis.ticks.x = element_line(colour = "#333333"),
                  axis.ticks.y = element_line(colour = "#333333"),
                  axis.ticks.length =  unit(0.26, "cm"), 
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank(),
                  strip.background =element_rect(fill="white", colour="#333333"))
    })
}
