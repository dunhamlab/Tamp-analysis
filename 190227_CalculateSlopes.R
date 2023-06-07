#!/usr/bin/env Rscript
# Modified version of Anna's CalculateSlope script to define function calculateSlope
library(splines)
# A function to return the slopes estimate for a particular gene for all the biological replicates
# X is your full data set
## gene is an argument for the gene we are interested in
# If plot==TRUE, it will plot all the scatterplots for each biological replicate
# the lines are the fitted values from the model
# If filter==TRUE, it will return NA for a particular gene that has less than 3 biological replicates
calculateSlope <- function(X,gene,dir,exp,plot=FALSE,filter=TRUE){
  
  slope <- NULL
  tempdata <- X[which(X$gene==gene),]
  uniquemer <- unique(tempdata$mer)
  
  if(filter==TRUE && length(uniquemer)<=3) return(NA)
  
  # par sets number of col and row blocks in the graph
  if(plot==TRUE){
    pdf(file=paste(dir,'/',gene,'_',exp,'.pdf',sep=""), width=11, height=8.5)
    par(mfrow=c(4,4))
  }
  
  i <- 1
  for(mer in uniquemer){
    tempdata2 <- tempdata[which(tempdata$mer==mer),]
    
    # Check how many replicates there are, need to fit separate intercept
    if(length(unique(tempdata2$replicate))==1){
      
      m1 <- lm(growth~generation,data=tempdata2)
      m2 <- lm(growth~pmin(generation,quantile(generation,0.5))+pmax(generation,quantile(generation,0.5)),data=tempdata2)
      # or lm(y~x+I((x-knot)*(x>knot)))		 
      
      # Test if piecewise linear model is preferable 
      checkpvalue <- anova(m1,m2)[[6]][2]	 	  	
      # Record the slope estimate 
      #	  	 slope <- c(slope,m1$coefficients[2])
      slope <- c(slope,ifelse(checkpvalue>=0.05,m1$coefficients[2],m2$coefficients[2]))
      
      if(plot==TRUE){
        if(i %% 16 ==1 && i!=1 ){
          par(mfrow=c(4,4))
          
        }	
        # plot the scatterplot			
        plot(tempdata2$generation,tempdata2$growth,xlab="generation",ylab="frequency",main=mer)
        if(checkpvalue>=0.05){
          # plot the fitted value
          lines(tempdata2$generation,predict(m1),col="black")	
        }
        
        else if(checkpvalue<0.05){
          # plot the fitted value
          lines(tempdata2$generation,predict(m2),col="black")	
        }
        
      }	
      
    }
    
    # If there are more than one technical replicates	  
    else if(length(unique(tempdata2$replicate))!=1){
      
      m1 <- lm(growth~generation+replicate,data=tempdata2)	  		  
      m2 <- lm(growth~pmin(generation,quantile(generation,0.5))+pmax(generation,quantile(generation,0.5))+replicate,data=tempdata2)
      # Test if piecewise linear model is preferable 
      checkpvalue <- anova(m1,m2)[[6]][2]	 	  	
      # Record the slope estimate 
      slope <- c(slope,ifelse(checkpvalue>=0.05,m1$coefficients[2],m2$coefficients[2]))
      
      
      if(plot==TRUE){
        if(i %% 16 ==1 && i!=1){
          par(mfrow=c(4,4))
          
        }
        replicate<-as.numeric(tempdata2$replicate)
        tempreplicate <- rep(unique(tempdata2$replicate),each=((round(max(tempdata2$generation))-round(min(tempdata2$generation)))*2+1))	
        tempgeneration <- rep(seq(round(min(tempdata2$generation)),round(max(tempdata2$generation)),by=0.5),length(unique(tempdata2$replicate)))
        newdata <- data.frame(replicate=tempreplicate,generation=tempgeneration)	
        tempreplicate <- as.numeric(tempreplicate)
        
        plot(tempdata2$generation,tempdata2$growth,xlab="generation",ylab="frequency",col=replicate,main=mer)
        if(checkpvalue>=0.05){
          for(j in unique(replicate)){		
            lines(tempdata2$generation[replicate==j],predict(m1)[replicate==j],col=j)
          }
        }	
        else if(checkpvalue<0.05){
          for(j in unique(replicate)){		
            lines(tempdata2$generation[replicate==j],predict(m2)[replicate==j],col=j)
          }						
        }	
        
      }	   
    }
    i <- i+1
    
  }
  
  if(plot==TRUE) dev.off()
  
  
  return(slope)
  
}

#summarizeSlope <- function(slopes){
#	if(length(slopes)<=15 ) return(mean(slopes))	
#	
#	if(length(slopes)>15){
#		temp <-hist(slopes,breaks=seq(min(slopes),max(slopes),by=(max(slopes)-min(slopes))/6),plot=FALSE)
#		maxcount <- which(temp$counts==max(temp$counts))
#		if(length(maxcount)>1) maxcount <- maxcount[1]
#		return((temp$breaks[maxcount]+temp$breaks[maxcount+1])/2)
#	}

#}

summarizeSlopedensity <- function(slopes){
  if(length(slopes)<=10 ) return(c(mean(slopes),sd(slopes)/sqrt(length(slopes)-1)))	
  
  if(length(slopes)>10){
    temp <- density(slopes)
    # calculate bootstrap standard error
    bootstrapest <- NULL
    for(b in 1:500){
      bootstrapslope <- sample(slopes,size=length(slopes),replace=TRUE)
      temp2 <- density(bootstrapslope)
      bootstrapest <- c(bootstrapest,temp2$x[which(temp2$y==max(temp2$y))])
    }
    
    return(c(temp$x[which(temp$y==max(temp$y))],sd(bootstrapest)))
  }
  
}
