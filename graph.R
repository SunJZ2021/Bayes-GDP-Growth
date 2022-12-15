graph <- function(country1,country2,stat=TRUE,g.prior=FALSE){
  
  suppressMessages(library(ggplot2))
  suppressMessages(library(patchwork))
  suppressMessages(library(gridExtra))
  
  if (!((country1 %in% codes$Code) & (country2 %in% codes$Code))) 
    stop("No data or wrong codes")
  
  Tag <- array()
  Tag[1] <- codes$Country[codes$Code==country1]
  Tag[2] <- codes$Country[codes$Code==country2]
  
  Data <- array()
  if (stat==TRUE){
    Data[1] <- paste("MCout_st",country1,sep=".")
    Data[2] <- paste("MCout_st",country2,sep=".")
  }else{
    Data[1] <- paste("MCout_ns",country1,sep=".")
    Data[2] <- paste("MCout_ns",country2,sep=".")
  }
  
  if (g.prior==TRUE){
    for (j in seq(1,2)){
      Data[j] <- paste(Data[j],".g",sep="")
    }
  }
  
  for (k in seq(1,2)){
    
    df <- get(Data[k])
    
    for (i in c("beta","alpha","lambda")){
      
      title <- switch(i,
                      "beta" = "Slope",
                      "alpha" = "Intercept",
                      "lambda" = "Precision")
      
      figure.name <- paste(i,as.character(k),sep=".")
      
      assign(figure.name, ggplot(
        data.frame(vs=df[[i]],
                   t=seq(0,length(df[[i]])-1,1)),
        aes(x=t,y=vs)) +
          geom_line(color="deepskyblue4") +
          theme(legend.position="none",
                panel.grid.minor = 
                  element_line(
                    colour="lightgrey",
                    linewidth=0.35),
                panel.background = element_blank(),
                axis.line = element_line(
                  linewidth = 0.7,
                  colour = "black",
                  linetype=1)) +
          labs(x="Draw", y="") +
          ggtitle(paste(Tag[k],title,sep=": ")) +
          scale_x_continuous(expand = c(0,0)))
    }
  }
  (beta.1|(alpha.1/lambda.1))/(beta.2|(alpha.2/lambda.2))
}