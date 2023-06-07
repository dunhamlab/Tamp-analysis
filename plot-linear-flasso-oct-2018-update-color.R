compare_linear_flasso_plot_v1 <- function(dat, dat.name, dat.varname, chr, centr, save=F) { 
  dat <- dat[dat$chr== chr, ]
  
  dat$Tamp.pos <- dat$start 
  dat$Tamp.pos[dat$arm == "R"] <- dat$stop[dat$arm == "R"]
  
  # Pull out data from left and right arms 
  dat.L <- dat[dat$arm == "L", ]
  dat.L <- dat.L[order(dat.L$Tamp.pos), ]
  dat.R <- dat[dat$arm == "R", ]
  dat.R <- dat.R[order(dat.R$Tamp.pos), ] 
  
  # If there are fewer than 7 data points on either side, 
  # then too few data to do a fused lasso analysis ... 
  if(nrow(dat.R) < 7 || nrow(dat.L) < 7) { 
    return()
  } 
  
  flasso.L <- fusedlasso1d(y=dat.L[, dat.varname], pos=dat.L$Tamp.pos)
  cv.flasso.L <-  cv.trendfilter(flasso.L)
  cv.flasso.L$i.2se <- which(cv.flasso.L$lambda == max(cv.flasso.L$lambda[cv.flasso.L$err <= min(cv.flasso.L$err) + 
                                                                            2*cv.flasso.L$se[which.min(cv.flasso.L$err)]]))
  flasso.fit.L.1se <- data.frame(Tamp.pos=dat.L$Tamp.pos, fit=flasso.L$fit[, cv.flasso.L$i.1se])
  flasso.fit.L.2se <- data.frame(Tamp.pos=dat.L$Tamp.pos, fit=flasso.L$fit[, cv.flasso.L$i.2se])
  
  flasso.R <- fusedlasso1d(y=dat.R[, dat.varname], pos=dat.R$Tamp.pos)
  cv.flasso.R <-  cv.trendfilter(flasso.R)
  cv.flasso.R$i.2se <- which(cv.flasso.R$lambda == max(cv.flasso.R$lambda[cv.flasso.R$err <= min(cv.flasso.R$err) + 
                                                                            2*cv.flasso.R$se[which.min(cv.flasso.R$err)]]))
  
  flasso.fit.R.1se <- data.frame(Tamp.pos=dat.R$Tamp.pos, fit=flasso.R$fit[, cv.flasso.R$i.1se])
  flasso.fit.R.2se <- data.frame(Tamp.pos=dat.R$Tamp.pos, fit=flasso.R$fit[, cv.flasso.R$i.2se])
  
  lmfit <- lm(as.formula(paste(dat.varname, "~", "Tamp.length", sep="")), data=dat)
  # Create data frame for use in plotting 
  lmfitplot <- data.frame(Tamp.pos= c(dat$Tamp.pos,centr[chr], centr[chr]+1e-16,dat$chrlength[1], 0), 
                            fit=predict(lmfit, data.frame(Tamp.length= c(dat$Tamp.length,  centr[chr], dat$chrlength[1]-centr[chr], 0, 0))))
  
  g <- ggplot() + geom_segment(aes(x=0,y=0,xend=chrom.lengths$end[chr],yend=0)) + geom_errorbar(aes_string(x="Tamp.pos", 
                                           ymin=paste(dat.varname, " - ", dat.varname, ".se", sep=""), 
                                           ymax=paste(dat.varname, " + ", dat.varname, ".se", sep="")), width=.1, 
                                data=dat[dat$chr == chr,],color='#C00000') +
    geom_point(aes_string(x="Tamp.pos", y=dat.varname), 
               data=dat[dat$chr == chr,], color='#C00000') + 
    geom_point(aes(x=centr[chr], y=0),size=4)  + 
    geom_line(data=lmfitplot,  aes(x=Tamp.pos, y=fit), linetype="dashed")+
    # geom_line(data=flasso.fit.L.1se, aes(x=Tamp.pos, y=fit), lty="dotted") +
    geom_line(data=flasso.fit.L.2se, aes(x=Tamp.pos, y=fit)) +
    #geom_line(data=flasso.fit.R.1se, aes(x=Tamp.pos, y=fit), lty="dotted") +
    geom_line(data=flasso.fit.R.2se, aes(x=Tamp.pos, y=fit)) +    
    ggtitle(paste(dat.name, ", Chr ", chr, sep="")) +
    labs(x= "starting nucleotide position of copied DNA", y = "relative fitness")+ 
    theme_minimal() +
    theme(axis.line.y = element_line(colour="black",size=0.4),axis.ticks.y= element_line(size=0.3))
  
  if(save) { 
    ggsave(paste(dat.varname, "-chr", chr, "-fit-v1-col.pdf", sep=""),  width =12, height=6)
  } else { 
    print(g)
  }
} 

####################################################################
compare_linear_flasso_plot_v2 <- function(dat, dat.name, dat.varname, chr, centr, col, save=F) { 
  dat <- dat[dat$chr== chr, ]

  dat$Tamp.pos <- dat$start 
  dat$Tamp.pos[dat$arm == "R"] <- dat$stop[dat$arm == "R"]
  dat$armL <- ifelse(dat$arm == "R", 1, 0)
  dat$Tamp.pos <- dat$Tamp.pos/1000
  dat$Tamp.length <- dat$Tamp.length/1000
  dat$chrlength <- dat$chrlength/1000
  centr = centr/1000
  chrom.lengths$end <- chrom.lengths$end/1000
  # Pull out data from left and right arms 
  dat.L <- dat[dat$arm == "L", ]
  dat.L <- dat.L[order(dat.L$Tamp.pos), ]
  dat.R <- dat[dat$arm == "R", ]
  dat.R <- dat.R[order(dat.R$Tamp.pos), ] 
  
  # If there are fewer than 7 data points on either side, 
  # then too few data to do a fused lasso analysis ... 
  if(nrow(dat.R) < 7 || nrow(dat.L) < 7) { 
    return()
  } 
  
  flasso.L <- fusedlasso1d(y=dat.L[, dat.varname], pos=dat.L$Tamp.pos)
  cv.flasso.L <-  cv.trendfilter(flasso.L)
  cv.flasso.L$i.2se <- which(cv.flasso.L$lambda == max(cv.flasso.L$lambda[cv.flasso.L$err <= min(cv.flasso.L$err) + 
                                                                            2*cv.flasso.L$se[which.min(cv.flasso.L$err)]]))
  flasso.fit.L.1se <- data.frame(Tamp.pos=dat.L$Tamp.pos, fit=flasso.L$fit[, cv.flasso.L$i.1se])
  flasso.fit.L.2se <- data.frame(Tamp.pos=dat.L$Tamp.pos, fit=flasso.L$fit[, cv.flasso.L$i.2se])
  
  flasso.R <- fusedlasso1d(y=dat.R[, dat.varname], pos=dat.R$Tamp.pos)
  cv.flasso.R <-  cv.trendfilter(flasso.R)
  cv.flasso.R$i.2se <- which(cv.flasso.R$lambda == max(cv.flasso.R$lambda[cv.flasso.R$err <= min(cv.flasso.R$err) + 
                                                                            2*cv.flasso.R$se[which.min(cv.flasso.R$err)]]))
  
  flasso.fit.R.1se <- data.frame(Tamp.pos=dat.R$Tamp.pos, fit=flasso.R$fit[, cv.flasso.R$i.1se])
  flasso.fit.R.2se <- data.frame(Tamp.pos=dat.R$Tamp.pos, fit=flasso.R$fit[, cv.flasso.R$i.2se])
  
  lmfit <- lm(as.formula(paste(dat.varname, "~", "Tamp.pos*armL", sep="")), data=dat)
  # Create data frame for use in plotting 
  lmfitplot <- data.frame(Tamp.pos= c(dat$Tamp.pos,centr[chr], centr[chr]+1e-16,
                                        dat$chrlength[1], 0), 
                            fit=predict(lmfit, 
                                        data.frame(Tamp.pos= c(dat$Tamp.pos,centr[chr], centr[chr]+1e-16,
                                                               dat$chrlength[1], 0), 
                                                   armL=c(dat$armL, 0, 1, 1, 0))))
  g <- ggplot() + geom_segment(aes(x=0,y=0,xend=chrom.lengths$end[chr],yend=0)) + geom_errorbar(aes_string(x="Tamp.pos", 
                                           ymin=paste(dat.varname, " - ", dat.varname, ".se", sep=""), 
                                           ymax=paste(dat.varname, " + ", dat.varname, ".se", sep="")), width=.1, 
                                data=dat[dat$chr == chr, ], color=col) +
    geom_point(aes_string(x="Tamp.pos", y=dat.varname), 
               data=dat[dat$chr == chr,],color=col) + 
    geom_point(aes(x=centr[chr], y=0),size=4)  + 
    geom_line(data=lmfitplot,  aes(x=Tamp.pos, y=fit), linetype="dashed")+
    # geom_line(data=flasso.fit.L.1se, aes(x=Tamp.pos, y=fit), lty="dotted") +
    geom_line(data=flasso.fit.L.2se, aes(x=Tamp.pos, y=fit), size=1.4) +
    #geom_line(data=flasso.fit.R.1se, aes(x=Tamp.pos, y=fit), lty="dotted") +
    geom_line(data=flasso.fit.R.2se, aes(x=Tamp.pos, y=fit), size=1.4) +    
    ggtitle(paste("Chromosome ", chr, sep="")) +
    labs(x= "starting nucleotide position of copied DNA (kb)", y = "Relative fitness")+ 
    theme_minimal() +
    theme(axis.line.y = element_line(colour="black",size=0.4),axis.ticks.y= element_line(size=0.3),text = element_text(size=20))
    # coord_cartesian(ylim=c(-0.5,0.5))
  
  if(save) { 
    ggsave(paste(dat.varname, "-chr", chr, "-fit-v2-col.pdf", sep=""),  width =12, height=6)
  } else { 
    print(g)
  }
} 

compare_linear_flasso_plot_v3 <- function(dat, dat.name, dat.varname, chr, centr, save=F) { 
  dat <- dat[dat$chr== chr, ]
  
  dat$Tamp.pos <- dat$start 
  dat$Tamp.pos[dat$arm == "R"] <- dat$stop[dat$arm == "R"]
  dat$armL <- ifelse(dat$arm == "R", 1, 0)
  dat$Tamp.posfit <- dat$Tamp.pos
  dat$Tamp.posfit[dat$arm == "R"] <- -dat$Tamp.posfit[dat$arm == "R"]
  
  # Pull out data from left and right arms 
  dat.L <- dat[dat$arm == "L", ]
  dat.L <- dat.L[order(dat.L$Tamp.pos), ]
  dat.R <- dat[dat$arm == "R", ]
  dat.R <- dat.R[order(dat.R$Tamp.pos), ] 
  
  # If there are fewer than 7 data points on either side, 
  # then too few data to do a fused lasso analysis ... 
  if(nrow(dat.R) < 7 || nrow(dat.L) < 7) { 
    return()
  } 
  
  flasso.L <- fusedlasso1d(y=dat.L[, dat.varname], pos=dat.L$Tamp.pos)
  cv.flasso.L <-  cv.trendfilter(flasso.L)
  cv.flasso.L$i.2se <- which(cv.flasso.L$lambda == max(cv.flasso.L$lambda[cv.flasso.L$err <= min(cv.flasso.L$err) + 
                                                                            2*cv.flasso.L$se[which.min(cv.flasso.L$err)]]))
  flasso.fit.L.1se <- data.frame(Tamp.pos=dat.L$Tamp.pos, fit=flasso.L$fit[, cv.flasso.L$i.1se])
  flasso.fit.L.2se <- data.frame(Tamp.pos=dat.L$Tamp.pos, fit=flasso.L$fit[, cv.flasso.L$i.2se])
  
  flasso.R <- fusedlasso1d(y=dat.R[, dat.varname], pos=dat.R$Tamp.pos)
  cv.flasso.R <-  cv.trendfilter(flasso.R)
  cv.flasso.R$i.2se <- which(cv.flasso.R$lambda == max(cv.flasso.R$lambda[cv.flasso.R$err <= min(cv.flasso.R$err) + 
                                                                            2*cv.flasso.R$se[which.min(cv.flasso.R$err)]]))
  
  flasso.fit.R.1se <- data.frame(Tamp.pos=dat.R$Tamp.pos, fit=flasso.R$fit[, cv.flasso.R$i.1se])
  flasso.fit.R.2se <- data.frame(Tamp.pos=dat.R$Tamp.pos, fit=flasso.R$fit[, cv.flasso.R$i.2se])
  
  # Fit the linear model 
  lmfit <- lm(as.formula(paste(dat.varname, "~", "Tamp.posfit + armL", sep="")), data=dat)
  # Create data frame for use in plotting 
  lmfitplot <- data.frame(Tamp.pos= c(dat$Tamp.pos,centr[chr], centr[chr]+1e-16, 
                                        0, dat$chrlength[1]), 
                            fit=predict(lmfit, 
                                        data.frame(Tamp.posfit= c(dat$Tamp.posfit,centr[chr], -(centr[chr]+1e-16), 0, 
                                                                  -dat$chrlength[1]), 
                                                   armL=c(dat$armL, 0, 1, 0, 1))))
  
  g <- ggplot() + geom_segment(aes(x=0,y=0,xend=chrom.lengths$end[chr],yend=0)) + geom_errorbar(aes_string(x="Tamp.pos", 
                                           ymin=paste(dat.varname, " - ", dat.varname, ".se", sep=""), 
                                           ymax=paste(dat.varname, " + ", dat.varname, ".se", sep="")), width=.1, 
                                data=dat[dat$chr == chr, ],color='#C00000') +
    geom_point(aes_string(x="Tamp.pos", y=dat.varname), 
               data=dat[dat$chr == chr,],color='#C00000') + 
    geom_point(aes(x=centr[chr], y=0))  + 
    geom_line(data=lmfitplot,  aes(x=Tamp.pos, y=fit), linetype="dashed")+
    # geom_line(data=flasso.fit.L.1se, aes(x=Tamp.pos, y=fit), lty="dotted") +
    geom_line(data=flasso.fit.L.2se, aes(x=Tamp.pos, y=fit)) +
    #geom_line(data=flasso.fit.R.1se, aes(x=Tamp.pos, y=fit), lty="dotted") +
    geom_line(data=flasso.fit.R.2se, aes(x=Tamp.pos, y=fit)) +    
    ggtitle(paste(dat.name, ", Chr ", chr, sep="")) +
    labs(x= "starting nucleotide position of copied DNA", y = "relative fitness")+ 
    theme_minimal() +
    theme(axis.line.y = element_line(colour="black",size=0.4),axis.ticks.y= element_line(size=0.3))
  
  if(save) { 
    ggsave(paste(dat.varname, "-chr", chr, "-fit-v3-col.pdf", sep=""),  width =12, height=6)
  } else { 
    print(g)
  }
} 


compare_linear_flasso_plot_v4 <- function(dat, dat.name, dat.varname, chr, centr, save=F) { 
  dat <- dat[dat$chr== chr, ]
  
  dat$Tamp.pos <- dat$start 
  dat$Tamp.pos[dat$arm == "R"] <- dat$stop[dat$arm == "R"]
  
  # Pull out data from left and right arms 
  dat.L <- dat[dat$arm == "L", ]
  dat.L <- dat.L[order(dat.L$Tamp.pos), ]
  dat.R <- dat[dat$arm == "R", ]
  dat.R <- dat.R[order(dat.R$Tamp.pos), ] 
  
  # If there are fewer than 7 data points on either side, 
  # then too few data to do a fused lasso analysis ... 
  if(nrow(dat.R) < 7 || nrow(dat.L) < 7) { 
    return()
  } 
  
  flasso.L <- fusedlasso1d(y=dat.L[, dat.varname], pos=dat.L$Tamp.pos)
  cv.flasso.L <-  cv.trendfilter(flasso.L)
  cv.flasso.L$i.2se <- which(cv.flasso.L$lambda == max(cv.flasso.L$lambda[cv.flasso.L$err <= min(cv.flasso.L$err) + 
                                                                            2*cv.flasso.L$se[which.min(cv.flasso.L$err)]]))
  flasso.fit.L.1se <- data.frame(Tamp.pos=dat.L$Tamp.pos, fit=flasso.L$fit[, cv.flasso.L$i.1se])
  flasso.fit.L.2se <- data.frame(Tamp.pos=dat.L$Tamp.pos, fit=flasso.L$fit[, cv.flasso.L$i.2se])
  
  flasso.R <- fusedlasso1d(y=dat.R[, dat.varname], pos=dat.R$Tamp.pos)
  cv.flasso.R <-  cv.trendfilter(flasso.R)
  cv.flasso.R$i.2se <- which(cv.flasso.R$lambda == max(cv.flasso.R$lambda[cv.flasso.R$err <= min(cv.flasso.R$err) + 
                                                                            2*cv.flasso.R$se[which.min(cv.flasso.R$err)]]))
  
  flasso.fit.R.1se <- data.frame(Tamp.pos=dat.R$Tamp.pos, 
                                 fit=flasso.R$fit[, cv.flasso.R$i.1se])
  flasso.fit.R.2se <- data.frame(Tamp.pos=dat.R$Tamp.pos, 
                                 fit=flasso.R$fit[, cv.flasso.R$i.2se])
  
  # Fit the linear model 
  lmfit <- lm(as.formula(paste(dat.varname, "~", "Tamp.length:(arm == \"R\") + Tamp.length", sep="")), 
              data=dat)
  # Create data frame for use in plotting 
  lmfitplot <- data.frame(Tamp.pos=c(0, centr[chr], 
                                     centr[chr]+1e-16, dat$chrlength[1]),
                          fit=predict(lmfit, data.frame(Tamp.length= c(0, centr[chr], 
                                                                       dat$chrlength[1] - (centr[chr] +1e-16), 0), 
                                                        arm=c("L", "L", "R", "R"))))
  
  g <- ggplot() + geom_segment(aes(x=0,y=0,xend=chrom.lengths$end[chr],yend=0)) + geom_errorbar(aes_string(x="Tamp.pos", 
                                           ymin=paste(dat.varname, " - ", dat.varname, ".se", sep=""), 
                                           ymax=paste(dat.varname, " + ", dat.varname, ".se", sep="")), width=.1, 
                                data=dat[dat$chr == chr, ],color='#C00000') +
    geom_point(aes_string(x="Tamp.pos", y=dat.varname), 
               data=dat[dat$chr == chr,],color='#C00000') + 
    geom_point(aes(x=centr[chr], y=0))  + 
    geom_line(data=lmfitplot,  aes(x=Tamp.pos, y=fit), linetype="dashed")+
   # geom_line(data=flasso.fit.L.1se, aes(x=Tamp.pos, y=fit), lty="dotted") +
    geom_line(data=flasso.fit.L.2se, aes(x=Tamp.pos, y=fit)) +
    #geom_line(data=flasso.fit.R.1se, aes(x=Tamp.pos, y=fit), lty="dotted") +
    geom_line(data=flasso.fit.R.2se, aes(x=Tamp.pos, y=fit)) +    
    ggtitle(paste(dat.name, ", Chr ", chr, sep="")) +
    labs(x= "starting nucleotide position of copied DNA", y = "relative fitness")+ 
    theme_minimal() +
    theme(axis.line.y = element_line(colour="black",size=0.4),axis.ticks.y= element_line(size=0.3))
  
  if(save) { 
    ggsave(paste(dat.varname, "-chr", chr, "-fit-v4-col.pdf", sep=""),  width =12, height=6)
    write.table(flasso.fit.R.2se, file = paste(dat.varname, "-chr", chr, "-fit-v4-R-arm.csv", sep=""), quote=FALSE,sep = ',',row.names = FALSE)
    write.table(flasso.fit.L.2se, file = paste(dat.varname, "-chr", chr, "-fit-v4-L-arm.csv", sep=""), quote=FALSE,sep = ',',row.names = FALSE)
  } else { 
    print(g)
  }
} 


