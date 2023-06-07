#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# args[1] = scripts directory
# args[2] = working directory
# args[3] = condition1
# args[4] = condition2

library("Hmisc")
genome <- function (scriptdir, inputdir, inputfile, output){        
  
  ##genome-wide "screen" plot for fitted slope
  allgenes <- data.frame(genes = read.csv(paste(scriptdir,'/gene_locations_names_2',sep = ''), sep="\t", header = FALSE)[,1], number=read.csv(paste(scriptdir,'/gene_locations_names_2',sep = ''), sep="\t", header = FALSE)[,2])
  slopes <- data.frame(genes = read.csv(paste(inputdir, inputfile, sep=""), sep=",")[,1], avg = read.csv(paste(inputdir, inputfile, sep=""), sep=",")[,2], std = read.csv(paste(inputdir, inputfile, sep=""), sep=",")[,3])
  allslopes <- merge(allgenes, slopes, all.y=TRUE)
  allslopes <- allslopes[order(allslopes$number), ]
  
  pdf(paste(inputdir, output, ".pdf", sep=""), width=11, height=8.5)
  par(mfrow=c(4,1), mar=c(3,5,3,2)) 
  print(head(allslopes))
  
  #chromosome1
  plot(x=1, y = 1, type = 'n', xaxt="n", xlim = c(0, 1531933), ylim = c(-0.5, 0.5), xlab="", bty = 'n', main = paste(output, 'chr1', sep = ': '))
  segments(0,0,230218,0)
  chr1x <- vector()
  chr1y <- vector()
  
  for(i in 1:length(allslopes$genes)){
    chr = substring(allslopes[i,1],2,2)
    if(chr == 'A'){
      #		  if(allslopes[i,6] > 3){
      #			  if(allslopes[i,4] < 0.1){
      
      errbar(allslopes[i,2], allslopes[i,3], allslopes[i,3]+allslopes[i,4]*2, allslopes[i,3]-allslopes[i,4]*2, add = T, pch = '', xaxt = 'n', cap = 0)
      points(allslopes[i,2], allslopes[i,3], type = 'p', cex=0.7, col= "black", bg="lightskyblue", pch=21, xaxt="n", xlim = c(0, 1531933), ylim = c(-0.9, 0.5), xlab="", bty = 'n')
      points(x = c(151465), y = c(0), pch = 15, col = 'red')
      chr1x <- append(chr1x, allslopes[i,2])
      chr1y <- append(chr1y, allslopes[i,3])
      #	  }
      #	 }
    }
  }
  #lines(chr1x, chr1y, col = 'blue')
  
  
  #chromosome2	
  plot(x=1, y = 1, type = 'n', xaxt="n", xlim = c(0, 1531933), ylim = c(-0.5, 0.5), xlab="", bty = 'n', main = paste(output, 'chr2', sep = ': '))
  segments(0,0,813184,0)
  chr2x <- vector()
  chr2y <- vector()
  
  for(i in 1:length(allslopes$genes)){
    if(substring(allslopes[i,1],2,2) == 'B'){
      #			if(allslopes[i,6] > 3){
      #				if(allslopes[i,4] < 0.1){
      errbar(allslopes[i,2], allslopes[i,3], allslopes[i,3]+allslopes[i,4]*2, allslopes[i,3]-allslopes[i,4]*2, add = T, pch = '', xaxt = 'n', cap = 0)		
      points(allslopes[i,2], allslopes[i,3], type = 'p', cex=0.7, col= "black", bg="lightskyblue", pch=21, xaxt="n", xlim = c(0, 1531933), ylim = c(-0.9, 0.5), xlab="", bty = 'n')
      
      points(x = c(238207), y = c(0), pch = 15, col = 'red')
      chr2x <- append(chr2x, allslopes[i,2])
      chr2y <- append(chr2y, allslopes[i,3])
      #		}
      #	}
    }
  }
  #	segments(238207,0,513184,0, col = 'orange', lwd = 3)
  points(x = c(238207), y = c(0), pch = 15, col = 'red')
  #	lines(chr2x, chr2y, col = 'blue')
  
  #chromosome3	
  plot(x=1, y = 1, type = 'n', xaxt="n", xlim = c(0, 1531933), ylim = c(-0.5, 0.5), xlab="", bty = 'n', main = paste(output, 'chr3', sep = ': '))
  segments(0,0,316620,0)
  chr3x <- vector()
  chr3y <- vector()
  
  for(i in 1:length(allslopes$genes)){
    if(substring(allslopes[i,1],2,2) == 'C'){
      #			if(allslopes[i,6] > 3){
      #				if(allslopes[i,4] < 0.1){
      errbar(allslopes[i,2], allslopes[i,3], allslopes[i,3]+allslopes[i,4]*2, allslopes[i,3]-allslopes[i,4]*2, add = T, pch = '', xaxt = 'n', cap = 0)
      points(allslopes[i,2], allslopes[i,3], type = 'p', cex=0.7, col= "black", bg="lightskyblue", pch=21, xaxt="n", xlim = c(0, 1531933), ylim = c(-0.9, 0.5), xlab="", bty = 'n')
      
      points(x = c(114385), y = c(0), pch = 15, col = 'red')
      chr3x <- append(chr3x, allslopes[i,2])
      chr3y <- append(chr3y, allslopes[i,3])
      #		}
      #	}
    }
  }
  #	lines(chr3x, chr3y, col = 'blue')
  
  #chromosome4	
  plot(x=1, y = 1, type = 'n', xaxt="n", xlim = c(0, 1531933), ylim = c(-0.5, 0.5), xlab="", bty = 'n', main = paste(output, 'chr4', sep = ': '))
  segments(0,0,1531933,0)
  chr4x <- vector()
  chr4y <- vector()
  
  for(i in 1:length(allslopes$genes)){
    if(substring(allslopes[i,1],2,2) == 'D'){
      #			if(allslopes[i,6] > 3){
      #				if(allslopes[i,4] < 0.1){
      errbar(allslopes[i,2], allslopes[i,3], allslopes[i,3]+allslopes[i,4]*2, allslopes[i,3]-allslopes[i,4]*2, add = T, pch = '', xaxt = 'n', cap = 0)
      points(allslopes[i,2], allslopes[i,3], type = 'p', cex=0.7, col= "black", bg="lightskyblue", pch=21, xaxt="n", xlim = c(0, 1531933), ylim = c(-0.9, 0.5), xlab="", bty = 'n')
      
      chr4x <- append(chr4x, allslopes[i,2])
      chr4y <- append(chr4y, allslopes[i,3])
      #		}
      #	}
    }
  }
  #	lines(chr4x, chr4y, col = 'blue')
  #	segments(300000,0,1231933,0, col = 'orange', lwd = 3)
  #  segments(1154216,0,1161320,0, col = 'orange', lwd = 5)
  points(x = c(449711), y = c(0), pch = 15, col = 'red')
  #  segments(1443403,0,1441433,0,col = 'red', lwd = 5)
  
  #chromosome5
  plot(x=1, y = 1, type = 'n', xaxt="n", xlim = c(0, 1531933), ylim = c(-0.5, 0.5), xlab="", bty = 'n', main = paste(output, 'chr5', sep = ': '))
  segments(0,0,576874,0)
  chr5x <- vector()
  chr5y <- vector()
  
  for(i in 1:length(allslopes$genes)){
    chr = substring(allslopes[i,1],2,2)
    if(chr == 'E'){
      #		  if(allslopes[i,6] > 3){
      #			  if(allslopes[i,4] < 0.1){
      errbar(allslopes[i,2], allslopes[i,3], allslopes[i,3]+allslopes[i,4]*2, allslopes[i,3]-allslopes[i,4]*2, add = T, pch = '', xaxt = 'n', cap = 0)  
      points(allslopes[i,2], allslopes[i,3], type = 'p', cex=0.7, col= "black", bg="lightskyblue", pch=21, xaxt="n", xlim = c(0, 1531933), ylim = c(-0.9, 0.5), xlab="", bty = 'n')
      
      points(x = c(151987), y = c(0), pch = 15, col = 'red')
      chr5x <- append(chr5x, allslopes[i,2])
      chr5y <- append(chr5y, allslopes[i,3])
      #	  }
      #	 }
    }
  }
  # lines(chr5x, chr5y, col = 'blue')
  
  
  #chromosome6	
  plot(x=1, y = 1, type = 'n', xaxt="n", xlim = c(0, 1531933), ylim = c(-0.5, 0.5), xlab="", bty = 'n', main = paste(output, 'chr6', sep = ': '))
  segments(0,0,270161,0)
  chr6x <- vector()
  chr6y <- vector()
  
  for(i in 1:length(allslopes$genes)){
    if(substring(allslopes[i,1],2,2) == 'F'){
      #			if(allslopes[i,6] > 3){
      #				if(allslopes[i,4] < 0.1){
      errbar(allslopes[i,2], allslopes[i,3], allslopes[i,3]+allslopes[i,4]*2, allslopes[i,3]-allslopes[i,4]*2, add = T, pch = '', xaxt = 'n', cap = 0)		
      points(allslopes[i,2], allslopes[i,3], type = 'p', cex=0.7, col= "black", bg="lightskyblue", pch=21, xaxt="n", xlim = c(0, 1531933), ylim = c(-0.9, 0.5), xlab="", bty = 'n')
      
      points(x = c(148510), y = c(0), pch = 15, col = 'red')
      chr6x <- append(chr6x, allslopes[i,2])
      chr6y <- append(chr6y, allslopes[i,3])
      #		}
      #	}
    }
  }
  #	lines(chr6x, chr6y, col = 'blue')
  
  #chromosome7	
  plot(x=1, y = 1, type = 'n', xaxt="n", xlim = c(0, 1531933), ylim = c(-0.5, 0.5), xlab="", bty = 'n', main = paste(output, 'chr7', sep = ': '))
  segments(0,0,1090940,0)
  chr7x <- vector()
  chr7y <- vector()
  
  for(i in 1:length(allslopes$genes)){
    if(substring(allslopes[i,1],2,2) == 'G'){
      #			if(allslopes[i,6] > 3){
      #				if(allslopes[i,4] < 0.1){
      errbar(allslopes[i,2], allslopes[i,3], allslopes[i,3]+allslopes[i,4]*2, allslopes[i,3]-allslopes[i,4]*2, add = T, pch = '', xaxt = 'n', cap = 0)
      points(allslopes[i,2], allslopes[i,3], type = 'p', cex=0.7, col= "black", bg="lightskyblue", pch=21, xaxt="n", xlim = c(0, 1531933), ylim = c(-0.9, 0.5), xlab="", bty = 'n')
      
      points(x = c(497038), y = c(0), pch = 15, col = 'red')
      chr7x <- append(chr7x, allslopes[i,2])
      chr7y <- append(chr7y, allslopes[i,3])
      #		}
      #	}
    }
  }
  #	lines(chr7x, chr7y, col = 'blue')
  
  #chromosome8	
  plot(x=1, y = 1, type = 'n', xaxt="n", xlim = c(0, 1531933), ylim = c(-0.5, 0.5), xlab="", bty = 'n', main = paste(output, 'chr8', sep = ': '))
  segments(0,0,562643,0)
  chr8x <- vector()
  chr8y <- vector()
  
  for(i in 1:length(allslopes$genes)){
    if(substring(allslopes[i,1],2,2) == 'H'){
      #			if(allslopes[i,6] > 3){
      #				if(allslopes[i,4] < 0.1){
      errbar(allslopes[i,2], allslopes[i,3], allslopes[i,3]+allslopes[i,4]*2, allslopes[i,3]-allslopes[i,4]*2, add = T, pch = '', xaxt = 'n', cap = 0)
      points(allslopes[i,2], allslopes[i,3], type = 'p', cex=0.7, col= "black", bg="lightskyblue", pch=21, xaxt="n", xlim = c(0, 1531933), ylim = c(-0.9, 0.5), xlab="", bty = 'n')
      
      chr8x <- append(chr8x, allslopes[i,2])
      chr8y <- append(chr8y, allslopes[i,3])
      #		}
      #	}
    }
  }
  #	lines(chr8x, chr8y, col = 'blue')
  #	segments(300000,0,1231933,0, col = 'orange', lwd = 3)
  # segments(1154216,0,1161320,0, col = 'orange', lwd = 5)
  points(x = c(105703), y = c(0), pch = 15, col = 'red')
  # segments(1443403,0,1441433,0,col = 'red', lwd = 5)
  
  #chromosome9
  plot(x=1, y = 1, type = 'n', xaxt="n", xlim = c(0, 1531933), ylim = c(-0.5, 0.5), xlab="", bty = 'n', main = paste(output, 'chr9', sep = ': '))
  segments(0,0,439888,0)
  chr9x <- vector()
  chr9y <- vector()
  
  for(i in 1:length(allslopes$genes)){
    chr = substring(allslopes[i,1],2,2)
    if(chr == 'I'){
      #		  if(allslopes[i,6] > 3){
      #			  if(allslopes[i,4] < 0.1){
      errbar(allslopes[i,2], allslopes[i,3], allslopes[i,3]+allslopes[i,4]*2, allslopes[i,3]-allslopes[i,4]*2, add = T, pch = '', xaxt = 'n', cap = 0)  
      points(allslopes[i,2], allslopes[i,3], type = 'p', cex=0.7, col= "black", bg="lightskyblue", pch=21, xaxt="n", xlim = c(0, 1531933), ylim = c(-0.9, 0.5), xlab="", bty = 'n')
      
      points(x = c(355629), y = c(0), pch = 15, col = 'red')
      chr9x <- append(chr9x, allslopes[i,2])
      chr9y <- append(chr9y, allslopes[i,3])
      #	  }
      #	 }
    }
  }
  # lines(chr9x, chr9y, col = 'blue')
  
  
  #chromosome10	
  plot(x=1, y = 1, type = 'n', xaxt="n", xlim = c(0, 1531933), ylim = c(-0.5, 0.5), xlab="", bty = 'n', main = paste(output, 'chr10', sep = ': '))
  segments(0,0,745751,0)
  chr10x <- vector()
  chr10y <- vector()
  
  for(i in 1:length(allslopes$genes)){
    if(substring(allslopes[i,1],2,2) == 'J'){
      #			if(allslopes[i,6] > 3){
      #				if(allslopes[i,4] < 0.1){
      errbar(allslopes[i,2], allslopes[i,3], allslopes[i,3]+allslopes[i,4]*2, allslopes[i,3]-allslopes[i,4]*2, add = T, pch = '', xaxt = 'n', cap = 0)		
      points(allslopes[i,2], allslopes[i,3], type = 'p', cex=0.7, col= "black", bg="lightskyblue", pch=21, xaxt="n", xlim = c(0, 1531933), ylim = c(-0.9, 0.5), xlab="", bty = 'n')
      
      points(x = c(436425), y = c(0), pch = 15, col = 'red')
      chr10x <- append(chr10x, allslopes[i,2])
      chr10y <- append(chr10y, allslopes[i,3])
      #		}
      #	}
    }
  }
  #	lines(chr10x, chr10y, col = 'blue')
  
  #chromosome11	
  plot(x=1, y = 1, type = 'n', xaxt="n", xlim = c(0, 1531933), ylim = c(-0.5, 0.5), xlab="", bty = 'n', main = paste(output, 'chr11', sep = ': '))
  segments(0,0,666816,0)
  chr11x <- vector()
  chr11y <- vector()
  
  for(i in 1:length(allslopes$genes)){
    if(substring(allslopes[i,1],2,2) == 'K'){
      #			if(allslopes[i,6] > 3){
      #				if(allslopes[i,4] < 0.1){
      errbar(allslopes[i,2], allslopes[i,3], allslopes[i,3]+allslopes[i,4]*2, allslopes[i,3]-allslopes[i,4]*2, add = T, pch = '', xaxt = 'n', cap = 0)
      points(allslopes[i,2], allslopes[i,3], type = 'p', cex=0.7, col= "black", bg="lightskyblue", pch=21, xaxt="n", xlim = c(0, 1531933), ylim = c(-0.9, 0.5), xlab="", bty = 'n')
      
      points(x = c(440246), y = c(0), pch = 15, col = 'red')
      chr11x <- append(chr11x, allslopes[i,2])
      chr11y <- append(chr11y, allslopes[i,3])
      #		}
      #	}
    }
  }
  #	lines(chr11x, chr11y, col = 'blue')
  #	segments(294610,0,296193,0,col = 'orange', lwd = 5)
  
  #chromosome12	
  plot(x=1, y = 1, type = 'n', xaxt="n", xlim = c(0, 1531933), ylim = c(-0.5, 0.5), xlab="", bty = 'n', main = paste(output, 'chr12', sep = ': '))
  segments(0,0,1078177,0)
  chr12x <- vector()
  chr12y <- vector()
  
  for(i in 1:length(allslopes$genes)){
    if(substring(allslopes[i,1],2,2) == 'L'){
      #			if(allslopes[i,6] > 3){
      #				if(allslopes[i,4] < 0.1){
      errbar(allslopes[i,2], allslopes[i,3], allslopes[i,3]+allslopes[i,4]*2, allslopes[i,3]-allslopes[i,4]*2, add = T, pch = '', xaxt = 'n', cap = 0)
      points(allslopes[i,2], allslopes[i,3], type = 'p', cex=0.7, col= "black", bg="lightskyblue", pch=21, xaxt="n", xlim = c(0, 1531933), ylim = c(-0.9, 0.5), xlab="", bty = 'n')
      
      chr12x <- append(chr12x, allslopes[i,2])
      chr12y <- append(chr12y, allslopes[i,3])
      #		}
      #	}
    }
  }
  #	lines(chr12x, chr12y, col = 'blue')
  #	segments(300000,0,1231933,0, col = 'orange', lwd = 3)
  # segments(1154216,0,1161320,0, col = 'orange', lwd = 5)
  points(x = c(150947), y = c(0), pch = 15, col = 'red')
  # segments(1443403,0,1441433,0,col = 'red', lwd = 5)
  
  #chromosome13
  plot(x=1, y = 1, type = 'n', xaxt="n", xlim = c(0, 1531933), ylim = c(-0.5, 0.5), xlab="", bty = 'n', main = paste(output, 'chr13', sep = ': '))
  segments(0,0,924431,0)
  chr13x <- vector()
  chr13y <- vector()
  
  for(i in 1:length(allslopes$genes)){
    chr = substring(allslopes[i,1],2,2)
    if(chr == 'M'){
      #		  if(allslopes[i,6] > 3){
      #			  if(allslopes[i,4] < 0.1){
      errbar(allslopes[i,2], allslopes[i,3], allslopes[i,3]+allslopes[i,4]*2, allslopes[i,3]-allslopes[i,4]*2, add = T, pch = '', xaxt = 'n', cap = 0)  
      points(allslopes[i,2], allslopes[i,3], type = 'p', cex=0.7, col= "black", bg="lightskyblue", pch=21, xaxt="n", xlim = c(0, 1531933), ylim = c(-0.9, 0.5), xlab="", bty = 'n')
      
      points(x = c(268031), y = c(0), pch = 15, col = 'red')
      chr13x <- append(chr13x, allslopes[i,2])
      chr13y <- append(chr13y, allslopes[i,3])
      #	  }
      #	 }
    }
  }
  # lines(chr13x, chr13y, col = 'blue')
  
  
  #chromosome14	
  plot(x=1, y = 1, type = 'n', xaxt="n", xlim = c(0, 1531933), ylim = c(-0.5, 0.5), xlab="", bty = 'n', main = paste(output, 'chr14', sep = ': '))
  segments(0,0,784333,0)
  chr14x <- vector()
  chr14y <- vector()
  
  for(i in 1:length(allslopes$genes)){
    if(substring(allslopes[i,1],2,2) == 'N'){
      #			if(allslopes[i,6] > 3){
      #				if(allslopes[i,4] < 0.1){
      errbar(allslopes[i,2], allslopes[i,3], allslopes[i,3]+allslopes[i,4]*2, allslopes[i,3]-allslopes[i,4]*2, add = T, pch = '', xaxt = 'n', cap = 0)		
      points(allslopes[i,2], allslopes[i,3], type = 'p', cex=0.7, col= "black", bg="lightskyblue", pch=21, xaxt="n", xlim = c(0, 1531933), ylim = c(-0.9, 0.5), xlab="", bty = 'n')
      
      points(x = c(628758), y = c(0), pch = 15, col = 'red')
      chr14x <- append(chr14x, allslopes[i,2])
      chr14y <- append(chr14y, allslopes[i,3])
      #		}
      #	}
    }
  }
  #	lines(chr14x, chr14y, col = 'blue')
  
  #chromosome15	
  plot(x=1, y = 1, type = 'n', xaxt="n", xlim = c(0, 1531933), ylim = c(-0.5, 0.5), xlab="", bty = 'n', main = paste(output, 'chr15', sep = ': '))
  segments(0,0,1091291,0)
  chr15x <- vector()
  chr15y <- vector()
  
  for(i in 1:length(allslopes$genes)){
    if(substring(allslopes[i,1],2,2) == 'O'){
      #			if(allslopes[i,6] > 3){
      #				if(allslopes[i,4] < 0.1){
      errbar(allslopes[i,2], allslopes[i,3], allslopes[i,3]+allslopes[i,4]*2, allslopes[i,3]-allslopes[i,4]*2, add = T, pch = '', xaxt = 'n', cap = 0)
      points(allslopes[i,2], allslopes[i,3], type = 'p', cex=0.7, col= "black", bg="lightskyblue", pch=21, xaxt="n", xlim = c(0, 1531933), ylim = c(-0.9, 0.5), xlab="", bty = 'n')
      
      points(x = c(326702), y = c(0), pch = 15, col = 'red')
      chr15x <- append(chr15x, allslopes[i,2])
      chr15y <- append(chr15y, allslopes[i,3])
      #		}
      #	}
    }
  }
  #	lines(chr15x, chr15y, col = 'blue')
  
  #chromosome16	
  plot(x=1, y = 1, type = 'n', xaxt="n", xlim = c(0, 1531933), ylim = c(-0.5, 0.5), xlab="", bty = 'n', main = paste(output, 'chr16', sep = ': '))
  segments(0,0,948066,0)
  chr16x <- vector()
  chr16y <- vector()
  
  for(i in 1:length(allslopes$genes)){
    if(substring(allslopes[i,1],2,2) == 'P'){
      #			if(allslopes[i,6] > 3){
      #				if(allslopes[i,4] < 0.1){
      errbar(allslopes[i,2], allslopes[i,3], allslopes[i,3]+allslopes[i,4]*2, allslopes[i,3]-allslopes[i,4]*2, add = T, pch = '', xaxt = 'n', cap = 0)
      points(allslopes[i,2], allslopes[i,3], type = 'p', cex=0.7, col= "black", bg="lightskyblue", pch=21, xaxt="n", xlim = c(0, 1531933), ylim = c(-0.9, 0.5), xlab="", bty = 'n')
      
      chr16x <- append(chr16x, allslopes[i,2])
      chr16y <- append(chr16y, allslopes[i,3])
      #		}
      #	}
    }
  }
  #	lines(chr16x, chr16y, col = 'blue')
  #	segments(300000,0,1231933,0, col = 'orange', lwd = 3)
  # segments(1154216,0,1161320,0, col = 'orange', lwd = 5)
  points(x = c(555957), y = c(0), pch = 15, col = 'red')
  # segments(1443403,0,1441433,0,col = 'red', lwd = 5)
  
  #	dev.off()
}

plot1 <- genome(inputdir=args[2], inputfile=paste('/',args[3],"_slopes_mean_norm.csv", sep = ''), scriptdir = args[1],output=paste('/',args[3],"_plot_by_chr"))
plot1 <- genome(inputdir=args[2], inputfile=paste('/',args[4],"_slopes_mean_norm.csv", sep = ''), scriptdir = args[1],output=paste('/',args[4],"_plot_by_chr"))

