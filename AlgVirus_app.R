# After modifying code and running app locally, you need to set up the account
# Visit https://www.shinyapps.io and get token then run following in console:
# rsconnect::setAccountInfo(name='lennonlab', token='F8514597B861B8FD199A346586E3EF58', secret='owstvIzgllZN3hxqTIXzXKOEwB7Mr0ayiPcOus6R')
# Then run following in console to deploy app
# rsconnect::deployApp('/Users/lennonj/GitHub/community-modules/bacteria-phage-eco/')
# install.packages("rsconnect")
# 

#######################################################################################################
# Alg-Virus model in R
######################################################################################################
rm(list = ls()) # clean the R history
######################################################################################################
#libraries that I need
######################################################################################################

library(tidyverse)
library(deSolve)       # solving the differential equation
library(ggplot2)       # ploting with ggplot2
library(vegan)         # best package for calculating diversity
library(reshape2)      # data manipulation
library(gt)            # allows you to make gt tables 
library(glue)          # allows you to glue strings togehter
library(patchwork)     # allow you to combine plots
library("gdata")       # allows you to bind rows with diffent numbwer of rows
library("wesanderson") # best color palettes around
library(changepoint)   # crucial for idenifying the transient time
library(shiny)


######################################################################################################
#  the time period of the differential equation
######################################################################################################

time <- seq(0,1000, by = 0.01)
time

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

D = 0.3          # dilution rate
Ni = 80          # nutrients
V = 0.33         # chemostat volume
d = 0.1          # death rate caused by the environment

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

bc=c(0.70,0.68,0.66,0.64,0.62)

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
       sliderInput("x", label= "Lowest cost of resitance",
                   min=0.5, max=0.75, value=0.7, step=0.02)),     # and like this you finish the sidebarLayout and panel
    mainPanel(
    plotOutput("dynamics"),    # first plot as an output
    plotOutput("trade.off"))                          # second plot as an output
  ))                      


server <- function(input, output, session){
  
  
  dataInput <- reactive({              #input the data first
  yini= c(N=80, C=c(as.numeric(input$C1),0,0,0,0), P=c(as.numeric(input$P1),0,0,0)) # this input is connected with the user interface and affects the outcome
  x= as.numeric(input$x)
  theta=replicate(1,{seq(x,x+0.08,0.02)})
  theta = apply(theta,2,sort,decreasing=T) # this is cool because your sort number in descending order
  theta =c(t(theta)) # here I make it correctly: I create a sequence of 5 numbers starting from the highest to lowest
  bc =theta  
  print(bc)
  parameters <- list(Ni, V, D=as.numeric(input$D), Kc, Cm, xc, bc, Wc, e, epsilon, M_C, phi, beta, M_P, A,d=as.numeric(input$d))
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
     geom_line(aes(y=log(algae), color="algae")) +
     geom_line(aes(y=log(virus), color="virus")) +
     scale_color_manual(values=c("algae"= "#36BA45", "virus"= "#48BBFC"), name="species")+
     labs(x="time", y="log(population size)", title="Alga - Virus population dynamics")
 })

 output$trade.off <- renderPlot({
   data <- as.data.frame(dataInput())
   trade.off <- as.data.frame(bc)
   type=c("C1","C2","C3","C4","C5")
   type=as.data.frame(type)
   data <- cbindX(data,trade.off,type)
   data <- data %>% 
     drop_na()
    ggplot(data, aes(type,bc,group=1)) +
    geom_point() +
    geom_line() +
    labs(x="alga types", y="growth rate", title="Trade off of the alga types")
     })
}

shinyApp(ui, server)


# out=as.data.frame(out)        # convert the output to a data frame
# out$algae= rowSums(out[3:7])  # algae will be the sum of all alga types
# out$virus= rowSums(out[8:11]) # virus will be the sum of all viral typ
# 
 trade.off <- as.data.frame(bc)
 type=c("C1","C2","C3","C4","C5")
 type=as.data.frame(type)
 data <- cbindX(data,trade.off,type)
