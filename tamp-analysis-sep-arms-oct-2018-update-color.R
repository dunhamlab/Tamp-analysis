library(genlasso)
library(flam)
library(ggplot2)

##### I changed this version to save breakpoints regardless of best fit #######


#### Functions ####
# Function for plotting the data from a setting specified by dat.name 
# where the variable name for rel. fitness is dat.varname 
# on a chromosome specified by chr, 
# where the positions of the centromeres are passed in by centr
plot_data <- function(dat, dat.name, dat.varname, chr, centr, save=F) {
  # Create variable indicating the nucleotide position of the 
  # non-telomeric end of the chromosome fragment copy 
  # = start on the left and = chromosome length - stop on the right 
  dat$Tamp.pos <- dat$start 
  dat$Tamp.pos[dat$arm == "R"] <- dat$stop[dat$arm == "R"]
  
  g <- ggplot() + geom_segment(aes(x=0,y=0,xend=chrom.lengths$end[chr],yend=0)) +
                                          geom_errorbar(aes_string(x="Tamp.pos", 
                                           ymin=paste(dat.varname, " - ", dat.varname, ".se", sep=""), 
                                           ymax=paste(dat.varname, " + ", dat.varname, ".se", sep="")), width=.1, 
                                data=dat[dat$chr == chr,], color='#960200') +
    geom_point(aes_string(x="Tamp.pos", y=dat.varname), 
               data=dat[dat$chr == chr,], color='#960200') + 
    geom_point(aes(x=centr[chr], y=0),size=4) + 
    ggtitle(paste(dat.name, ", Chr ", chr, sep="")) +
    labs(x= "starting nucleotide position of copied DNA", y = "relative fitness")+ 
    theme_minimal() +
    theme(axis.line.y = element_line(colour="black",size=0.4),axis.ticks.y= element_line(size=0.3))
  print(g)
  
  if(save) { 
    ggsave(paste(dat.varname, "-chr", chr, ".pdf", sep=""), width =12, height=6)
  }
}
setwd("~/Desktop/March2020_WFH/Witten_model/R-code")
source("plot-linear-flasso-oct-2018-update-color.R")
source("compare-linear-flasso-sep-arms-oct-2018-update.R")
# source("plot-linear-flasso-oct-2018-update-color-1se.R")
# source("compare-linear-flasso-sep-arms-oct-2018-update-1se.R")

#### Data ####
setwd("~/Desktop/March2020_WFH/Witten_model/data")
# anna's data
anna.gluc <- read.csv('171103_AS_glucose_slopes_locations.csv')
anna.sul <- read.csv('171103_AS_sulfate_slopes_locations.csv')
anna.phos <- read.csv('171103_AS_phosphate_slopes_locations.csv')

# abby's data
abby.D <- read.csv("MN_D_slopes_se_normF.csv")
abby.R <- read.csv("MN_R_slopes_se_normF.csv")

# abby.T30 <- read.csv('MN_T30_slopes_se_norm9.csv')
# abby.T37 <- read.csv('MN_T37_slopes_se_norm9.csv')

abby.low <- read.csv('MN_low_slopes_se_normF.csv')
abby.high <- read.csv('MN_high_slopes_se_normF.csv')
abby.stat <- read.csv("MN_stat_slopes_se_normF.csv")

# chromosome lengths and centromeres
chrom.lengths <- read.csv("chromosome_lengths.csv")

centr <- chrom.lengths$centromere

#### Plot data #### 
setwd("~/Desktop/March2020_WFH/Witten_model/data-plots")
for(chr in 1:16) {
  plot_data(anna.gluc, "glucose", "glucose", chr, centr, T)
  plot_data(anna.phos, "phosphate", "phosphate", chr, centr, T)
  plot_data(anna.sul, "sulfate", "sulfate", chr, centr, T)
  # plot_data(abby.low, "Exp1_T30", "Exp1_T30", chr, centr, T)
  # plot_data(abby.high, "Exp1_T37", "Exp1_T37", chr, centr, T)
  # plot_data(abby.stat,"MN_stat","MN_stat",chr,centr,T)
  # plot_data(abby.T30, "Exp2_T30","Exp2_T30",chr,centr,T)
  # plot_data(abby.T37, "Exp2_T37","Exp2_T37",chr,centr,T)
  # plot_data(abby.D, "MN_D", "MN_D", chr, centr, T)
  # plot_data(abby.R, "MN_R", "MN_R", chr, centr, T)
} 

#### Plot fits #### 
# for(chr in 1:16) {
  setwd("~/Desktop/March2020_WFH/Witten_model/fit-plots-v1-SameSlopesSameEnds-oct-2018-update/")
  # compare_linear_flasso_plot_v1(anna.gluc, "glucose", "glucose", chr, centr, T)
  # compare_linear_flasso_plot_v1(anna.phos, "phosphate", "phosphate", chr, centr, T)
  # compare_linear_flasso_plot_v1(anna.sul, "sulfate", "sulfate", chr, centr, T)
  # compare_linear_flasso_plot_v1(abby.low, "Exp1_T30", "Exp1_T30", chr, centr, T)
  # compare_linear_flasso_plot_v1(abby.high, "Exp1_T37", "Exp1_T37", chr, centr, T)
  # compare_linear_flasso_plot_v1(abby.stat,'MN_stat','MN_stat',chr,centr,save=T)
  # compare_linear_flasso_plot_v1(abby.T30,'Exp2_T30','Exp2_T30',chr,centr,save=T)
  # compare_linear_flasso_plot_v1(abby.T37,'Exp2_T37','Exp2_T37',chr,centr,save=T)
  # compare_linear_flasso_plot_v1(abby.D, "MN_D", "MN_D", chr, centr, save=T)
  # compare_linear_flasso_plot_v1(abby.R, "MN_R", "MN_R", chr, centr, save=T)
# }

  # compare_linear_flasso_plot_v2 <- function(dat, dat.name, dat.varname, chr, centr, col, save=F) { 
for(chr in 1:16) {
  setwd("~/Desktop/March2020_WFH/Witten_model/fit-plots-v2-DiffSlopesDiffEnds-oct-2018-update/")
  compare_linear_flasso_plot_v2(anna.gluc, "glucose", "glucose", chr, centr, T)
  compare_linear_flasso_plot_v2(anna.phos, "phosphate", "phosphate", chr, centr, T)
  compare_linear_flasso_plot_v2(anna.sul, "sulfate", "sulfate", chr, centr, T)
  compare_linear_flasso_plot_v2(dat = abby.low, dat.name="Exp1_T30", dat.varname="Exp1_T30", chr, centr=centr, col='#000066', save=T)
  compare_linear_flasso_plot_v2(abby.high, "Exp1_T37", "Exp1_T37", chr, centr, col='#ff9300', T)
  compare_linear_flasso_plot_v2(abby.stat,'MN_stat','MN_stat',chr,centr, col='#9a0794',save=T)
  compare_linear_flasso_plot_v2(abby.T30,'Exp2_T30','Exp2_T30',chr,centr,save=T)
  compare_linear_flasso_plot_v2(abby.T37,'Exp2_T37','Exp2_T37',chr,centr,save=T)
  compare_linear_flasso_plot_v2(abby.D, "MN_D", "MN_D", chr, centr, col='#107e7d', save=T)
  compare_linear_flasso_plot_v2(abby.R, "MN_R", "MN_R", chr, centr, col='#960200', save=T)
}
# for(chr in 1:16) {
#   setwd("~/Desktop/March2020_WFH/Witten_model/fit-plots-v3-SameSlopesDiffEnds-oct-2018-update/")
  # compare_linear_flasso_plot_v3(abby.low, "Exp1_T30", "Exp1_T30", chr, centr, T)
  # compare_linear_flasso_plot_v3(abby.high, "Exp1_T37", "Exp1_T37", chr, centr, T)
  # compare_linear_flasso_plot_v3(abby.stat,'MN_stat','MN_stat',chr,centr,save=T)
  # compare_linear_flasso_plot_v3(abby.T30,'Exp2_T30','Exp2_T30',chr,centr,save=T)
  # compare_linear_flasso_plot_v3(abby.T37,'Exp2_T37','Exp2_T37',chr,centr,save=T)
  # compare_linear_flasso_plot_v3(abby.D, "MN_D", "MN_D", chr, centr, save=T)
  # compare_linear_flasso_plot_v3(abby.R, "MN_R", "MN_R", chr, centr, save=T)
# }
# for(chr in 1:16) { 
#   setwd("~/Desktop/March2020_WFH/Witten_model/fit-plots-v4-DiffSlopesSameEnds-oct-2018-update/")
  # compare_linear_flasso_plot_v4(abby.low, "Exp1_T30", "Exp1_T30", chr, centr, T)
  # compare_linear_flasso_plot_v4(abby.high, "Exp1_T37", "Exp1_T37", chr, centr, T)
  # compare_linear_flasso_plot_v4(abby.stat,'MN_stat','MN_stat',chr,centr,save=T)
  # compare_linear_flasso_plot_v4(abby.T30,'Exp2_T30','Exp2_T30',chr,centr,save=T)
  # compare_linear_flasso_plot_v4(abby.T37,'Exp2_T37','Exp2_T37',chr,centr,save=T)
  # compare_linear_flasso_plot_v4(abby.D, "MN_D", "MN_D", chr, centr, save=T)
  # compare_linear_flasso_plot_v4(abby.R, "MN_R", "MN_R", chr, centr, save=T)
# }

#### Compare linear and piecewise constant on each arm separately ####
anna.sul.cvmse.L <- matrix(0, 10, 16)
anna.gluc.cvmse.L <- matrix(0, 10, 16)
anna.phos.cvmse.L <- matrix(0, 10, 16)
# abby.low.cvmse.L <- matrix(0, 10, 16)
# abby.high.cvmse.L <- matrix(0, 10, 16)
# abby.stat.cvmse.L <- matrix(0,10,16)
# abby.T30.cvmse.L <- matrix(0,10,16)
# abby.T37.cvmse.L <- matrix(0,10,16)
# abby.D.cvmse.L <- matrix(0, 10, 16)
# abby.R.cvmse.L <- matrix(0, 10, 16)

anna.sul.cvmse.R <- matrix(0, 10, 16)
anna.gluc.cvmse.R <- matrix(0, 10, 16)
anna.phos.cvmse.R <- matrix(0, 10, 16)
# abby.low.cvmse.R <- matrix(0, 10, 16)
# abby.high.cvmse.R <- matrix(0, 10, 16)
# abby.stat.cvmse.R <- matrix(0,10,16)
# abby.T30.cvmse.R <- matrix(0,10,16)
# abby.T37.cvmse.R <- matrix(0,10,16)
# abby.D.cvmse.R <- matrix(0, 10, 16)
# abby.R.cvmse.R <- matrix(0, 10, 16)

anna.sul.best.L <- rep("",16)
anna.gluc.best.L <- rep("",16)
anna.phos.best.L <- rep("",16)
# abby.low.best.L <- rep("",16)
# abby.high.best.L <- rep("",16)
# abby.stat.best.L <- rep("",16)
# abby.T30.best.L <- rep("",16)
# abby.T37.best.L <- rep("",16)
# abby.D.best.L <- rep("", 16)
# abby.R.best.L <- rep("", 16)

anna.sul.best.R <- rep("", 16)
anna.phos.best.R <- rep("",16)
anna.gluc.best.R <- rep("",16)
# abby.low.best.R <- rep("",16)
# abby.high.best.R <- rep("",16)
# abby.stat.best.R <- rep("",16)
# abby.T30.best.R <- rep("",16)
# abby.T37.best.R <- rep("",16)
# abby.D.best.R <- rep("", 16)
# abby.R.best.R <- rep("", 16)

anna.gluc.breaks <- data.frame(chr=character(), coord=character(), magnitude=character())
anna.phos.breaks <- data.frame(chr=character(), coord=character(), magnitude=character())
anna.sul.breaks <- data.frame(chr=character(), coord=character(), magnitude=character())
# abby.low.breaks <- data.frame(chr=character(), coord=character(), magnitude=character())
# abby.high.breaks <- data.frame(chr=character(), coord=character(), magnitude=character())
# abby.stat.breaks <- data.frame(chr=character(), coord=character(), magnitude=character())
# abby.T30.breaks <- data.frame(chr=character(), coord=character(), magnitude=character())
# abby.T37.breaks <- data.frame(chr=character(), coord=character(), magnitude=character())
# abby.D.breaks <- data.frame(chr=character(), coord=character(), magnitude=character())
# abby.R.breaks <- data.frame(chr=character(), coord=character(), magnitude=character())

all.models <- c("Piecewise Constant", "Linear Model 1", "Linear Model 2", 
                "Linear Model 3", "Linear Model 4")

#### CHANGE THIS LINE OF CODE IF YOU WANT TO CHANGE WHICH MODELS COMPARED ####
#### TO COMPARE ALL MODELS, SET IT TO 1:5 ####
#### TO COMPARE PIECEWISE CONSTANT and LM 4 FOR EXAMPLE, SET TO c(1, 5) #### 
models.considered.index <- c(1,3)

##### YOU WILL WANT TO REMEMBER TO CHANGE YOUR DIRECTORY SO YOU DON'T 
##### OVERWRITE RESULTS! #### 
setwd("~/Desktop/March2020_WFH/Witten_model/analysis-results/sep-arms-compare-PC-all-oct-2018-update-2se")

for(chr in 1:16) { 

  compare <- compare_linear_flasso_sep_arms(anna.sul,  "sulfate", chr, centr,
                                            all.models[models.considered.index])
  anna.sul.cvmse.L[,chr] <- compare$cv.mse.L
  anna.sul.cvmse.R[,chr] <- compare$cv.mse.R
  anna.sul.best.L[chr] <- compare$fits.best.L
  anna.sul.best.R[chr] <- compare$fits.best.R

  breakpts.L <- compare$breakpts.L
  breakpts.R <- compare$breakpts.R

  # if best is piecewise constant, save breakpts
  if(!is.na(anna.sul.best.L[chr])) {
  # if(anna.sul.best.L[chr] == "Piecewise Constant" & !is.na(anna.sul.best.L[chr])) {
    if(length(breakpts.L) != 0) {
      row.L <- as.integer(rownames(compare$lasso.fits.L[compare$lasso.fits.L$Tamp.pos %in% compare$breakpts.L,]))
      magnitude.L <- compare$lasso.fits.L$fit[row.L] - compare$lasso.fits.L$fit[(row.L-1)]
      anna.sul.break.L <- data.frame(chr=rep(chr,length(compare$breakpts.L)), coord=compare$breakpts.L,magnitude=magnitude.L)
      anna.sul.breaks <- rbind(anna.sul.breaks, anna.sul.break.L)
      # write.table(matrix(compare$breakpts.L, nrow=1), file=paste("anna.sul-chr", chr,
      #                                       "-left-arm-breakpts.txt",   sep=""),
      #           row.names=F, col.names=F)
    }
  }
  if(!is.na(anna.sul.best.R[chr])) {
    # if(anna.sul.best.R[chr] == "Piecewise Constant" & !is.na(anna.sul.best.R[chr])) {
    if(length(breakpts.R) != 0) {
      row.R <- as.integer(rownames(compare$lasso.fits.R[compare$lasso.fits.R$Tamp.pos %in% compare$breakpts.R,]))
      magnitude.R <- compare$lasso.fits.R$fit[row.R] - compare$lasso.fits.R$fit[(row.R-1)]
      anna.sul.break.R <- data.frame(chr=rep(chr,length(compare$breakpts.R)), coord=compare$breakpts.R,magnitude=magnitude.R)
      anna.sul.breaks <- rbind(anna.sul.breaks, anna.sul.break.R)
      #    #write.table(matrix(compare$breakpts.R, nrow=1), file=paste("anna.sul-chr", chr,
      #    #                                                          "-right-arm-breakpts.txt",   sep=""),
      #     #          row.names=F, col.names=F)
    }
  }

  compare <- compare_linear_flasso_sep_arms(anna.phos,  "phosphate", chr, centr,
                                            all.models[models.considered.index])
  anna.phos.cvmse.L[,chr] <- compare$cv.mse.L
  anna.phos.cvmse.R[,chr] <- compare$cv.mse.R
  anna.phos.best.L[chr] <- compare$fits.best.L
  anna.phos.best.R[chr] <- compare$fits.best.R

  breakpts.L <- compare$breakpts.L
  breakpts.R <- compare$breakpts.R

  # if best is piecewise constant, save breakpts
  if(!is.na(anna.phos.best.L[chr])) {
    # if(anna.phos.best.L[chr] == "Piecewise Constant" & !is.na(anna.phos.best.L[chr])) {
    if(length(breakpts.L) != 0) {
      row.L <- as.integer(rownames(compare$lasso.fits.L[compare$lasso.fits.L$Tamp.pos %in% compare$breakpts.L,]))
      magnitude.L <- compare$lasso.fits.L$fit[row.L] - compare$lasso.fits.L$fit[(row.L-1)]
      anna.phos.break.L <- data.frame(chr=rep(chr,length(compare$breakpts.L)), coord=compare$breakpts.L,magnitude=magnitude.L)
      anna.phos.breaks <- rbind(anna.phos.breaks, anna.phos.break.L)
      #    #write.table(matrix(compare$breakpts.L, nrow=1), file=paste("anna.phos-chr", chr,
      #    #                                                          "-left-arm-breakpts.txt",   sep=""),
      #     #          row.names=F, col.names=F)
    }
  }
  if(!is.na(anna.phos.best.R[chr])) {
  # if(anna.phos.best.R[chr] == "Piecewise Constant" & !is.na(anna.phos.best.R[chr])) {
    if(length(breakpts.R) != 0) {
      row.R <- as.integer(rownames(compare$lasso.fits.R[compare$lasso.fits.R$Tamp.pos %in% compare$breakpts.R,]))
      magnitude.R <- compare$lasso.fits.R$fit[row.R] - compare$lasso.fits.R$fit[(row.R-1)]
      anna.phos.break.R <- data.frame(chr=rep(chr,length(compare$breakpts.R)), coord=compare$breakpts.R,magnitude=magnitude.R)
      anna.phos.breaks <- rbind(anna.phos.breaks, anna.phos.break.R)
      #    #write.table(matrix(compare$breakpts.R, nrow=1), file=paste("anna.phos-chr", chr,
      #    #                                                          "-right-arm-breakpts.txt",   sep=""),
      #     #          row.names=F, col.names=F)
    }
  }

  compare <- compare_linear_flasso_sep_arms(anna.gluc,  "glucose", chr, centr,
                                            all.models[models.considered.index])
  anna.gluc.cvmse.L[,chr] <- compare$cv.mse.L
  anna.gluc.cvmse.R[,chr] <- compare$cv.mse.R
  anna.gluc.best.L[chr] <- compare$fits.best.L
  anna.gluc.best.R[chr] <- compare$fits.best.R

  breakpts.L <- compare$breakpts.L
  breakpts.R <- compare$breakpts.R

  # if best is piecewise constant, save breakpts
  if(!is.na(anna.gluc.best.L[chr])) {
  # if(anna.gluc.best.L[chr] == "Piecewise Constant" & !is.na(anna.gluc.best.L[chr])) {
    if(length(breakpts.L) != 0) {
      row.L <- as.integer(rownames(compare$lasso.fits.L[compare$lasso.fits.L$Tamp.pos %in% compare$breakpts.L,]))
      magnitude.L <- compare$lasso.fits.L$fit[row.L] - compare$lasso.fits.L$fit[(row.L-1)]
      anna.gluc.break.L <- data.frame(chr=rep(chr,length(compare$breakpts.L)), coord=compare$breakpts.L,magnitude=magnitude.L)
      anna.gluc.breaks <- rbind(anna.gluc.breaks, anna.gluc.break.L)
      # write.table(matrix(compare$breakpts.L, nrow=1), file=paste("anna.gluc-chr", chr,
      #                                       "-left-arm-breakpts.txt",   sep=""),
      #           row.names=F, col.names=F)
    }
  }
  if(!is.na(anna.gluc.best.R[chr])) {
  # if(anna.gluc.best.R[chr] == "Piecewise Constant" & !is.na(anna.gluc.best.R[chr])) {
    if(length(breakpts.R) != 0) {
      row.R <- as.integer(rownames(compare$lasso.fits.R[compare$lasso.fits.R$Tamp.pos %in% compare$breakpts.R,]))
      magnitude.R <- compare$lasso.fits.R$fit[row.R] - compare$lasso.fits.R$fit[(row.R-1)]
      anna.gluc.break.R <- data.frame(chr=rep(chr,length(compare$breakpts.R)), coord=compare$breakpts.R,magnitude=magnitude.R)
      anna.gluc.breaks <- rbind(anna.gluc.breaks, anna.gluc.break.R)
      # write.table(matrix(compare$breakpts.R, nrow=1), file=paste("anna.gluc-chr", chr,
      #                                                           "-right-arm-breakpts.txt",   sep=""),
      #           row.names=F, col.names=F)
    }
  }
  
  
  # compare <- compare_linear_flasso_sep_arms(abby.D,  "MN_D", chr, centr,
  #                                           all.models[models.considered.index])
  # abby.D.cvmse.L[, chr] <- compare$cv.mse.L
  # abby.D.cvmse.R[, chr] <- compare$cv.mse.R
  # abby.D.best.L[chr] <- compare$fits.best.L
  # abby.D.best.R[chr] <- compare$fits.best.R
  # 
  # breakpts.L <- compare$breakpts.L
  # breakpts.R <- compare$breakpts.R
  # 
  # # if best is piecewise constant, save breakpts
  # if(!is.na(abby.D.best.L[chr])) {
  # # if(abby.D.best.L[chr] == "Piecewise Constant" & !is.na(abby.D.best.L[chr])) {
  #   if(length(breakpts.L) != 0) {
  #     row.L <- as.integer(rownames(compare$lasso.fits.L[compare$lasso.fits.L$Tamp.pos %in% compare$breakpts.L,]))
  #     magnitude.L <- compare$lasso.fits.L$fit[row.L] - compare$lasso.fits.L$fit[(row.L-1)]
  #     abby.D.break.L <- data.frame(chr=rep(chr,length(compare$breakpts.L)), coord=compare$breakpts.L,magnitude=magnitude.L)
  #     abby.D.breaks <- rbind(abby.D.breaks, abby.D.break.L)
  #     #    #write.table(matrix(compare$breakpts.L, nrow=1), file=paste("abby.D-chr", chr, 
  #     #    #                                                          "-left-arm-breakpts.txt",   sep=""), 
  #     #     #          row.names=F, col.names=F)
  #   }
  # }
  # if(!is.na(abby.D.best.R[chr])) {
  # # if(abby.D.best.R[chr] == "Piecewise Constant" & !is.na(abby.D.best.R[chr])) {
  #   if(length(breakpts.R) != 0) {
  #     row.R <- as.integer(rownames(compare$lasso.fits.R[compare$lasso.fits.R$Tamp.pos %in% compare$breakpts.R,]))
  #     magnitude.R <- compare$lasso.fits.R$fit[row.R] - compare$lasso.fits.R$fit[(row.R-1)]
  #     abby.D.break.R <- data.frame(chr=rep(chr,length(compare$breakpts.R)), coord=compare$breakpts.R,magnitude=magnitude.R)
  #     abby.D.breaks <- rbind(abby.D.breaks, abby.D.break.R)
  #     #    # write.table(matrix(compare$breakpts.R, nrow=1), file=paste("abby.D-chr", chr, 
  #     #   #                                                           "-right-arm-breakpts.txt",   sep=""), 
  #     #    #           row.names=F, col.names=F)
  #   }
  # }
  # 
  # compare <- compare_linear_flasso_sep_arms(abby.R,  "MN_R", chr, centr,
  #                                           all.models[models.considered.index])
  # abby.R.cvmse.L[, chr] <- compare$cv.mse.L
  # abby.R.cvmse.R[, chr] <- compare$cv.mse.R
  # abby.R.best.L[chr] <- compare$fits.best.L
  # abby.R.best.R[chr] <- compare$fits.best.R
  # 
  # breakpts.L <- compare$breakpts.L
  # breakpts.R <- compare$breakpts.R
  # 
  # # if best is piecewise constant, save breakpts
  # if(!is.na(abby.R.best.L[chr])) {
  # # if(abby.R.best.L[chr] == "Piecewise Constant" & !is.na(abby.R.best.L[chr])) {
  #   if(length(breakpts.L) != 0) {
  #     row.L <- as.integer(rownames(compare$lasso.fits.L[compare$lasso.fits.L$Tamp.pos %in% compare$breakpts.L,]))
  #     magnitude.L <- compare$lasso.fits.L$fit[row.L] - compare$lasso.fits.L$fit[(row.L-1)]
  #     abby.R.break.L <- data.frame(chr=rep(chr,length(compare$breakpts.L)), coord=compare$breakpts.L,magnitude=magnitude.L)
  #     abby.R.breaks <- rbind(abby.R.breaks, abby.R.break.L)
  #     #    #write.table(matrix(compare$breakpts.L, nrow=1), file=paste("abby.R-chr", chr, 
  #     #    #                                                          "-left-arm-breakpts.txt",   sep=""), 
  #     #     #          row.names=F, col.names=F)
  #   }
  # }
  # if(!is.na(abby.R.best.R[chr])) {
  # # if(abby.R.best.R[chr] == "Piecewise Constant" & !is.na(abby.R.best.R[chr])) {
  #   if(length(breakpts.R) != 0) {
  #     row.R <- as.integer(rownames(compare$lasso.fits.R[compare$lasso.fits.R$Tamp.pos %in% compare$breakpts.R,]))
  #     magnitude.R <- compare$lasso.fits.R$fit[row.R] - compare$lasso.fits.R$fit[(row.R-1)]
  #     abby.R.break.R <- data.frame(chr=rep(chr,length(compare$breakpts.R)), coord=compare$breakpts.R,magnitude=magnitude.R)
  #     abby.R.breaks <- rbind(abby.R.breaks, abby.R.break.R)
  #     #    # write.table(matrix(compare$breakpts.R, nrow=1), file=paste("abby.R-chr", chr, 
  #     #   #                                                           "-right-arm-breakpts.txt",   sep=""), 
  #     #    #           row.names=F, col.names=F)
  #   }
  # }
  
  # compare <- compare_linear_flasso_sep_arms(abby.T30,  "Exp2_T30", chr, centr,
  #                                           all.models[models.considered.index])
  # abby.T30.cvmse.L[,chr] <- compare$cv.mse.L
  # abby.T30.cvmse.R[,chr] <- compare$cv.mse.R
  # abby.T30.best.L[chr] <- compare$fits.best.L
  # abby.T30.best.R[chr] <- compare$fits.best.R
  # 
  # breakpts.L <- compare$breakpts.L
  # breakpts.R <- compare$breakpts.R
  # 
  # # if best is piecewise constant, save breakpts
  # if(!is.na(abby.T30.best.L[chr])) {
  # # if(abby.T30.best.L[chr] == "Piecewise Constant" & !is.na(abby.T30.best.L[chr])) {
  #   if(length(breakpts.L) != 0) {
  #     row.L <- as.integer(rownames(compare$lasso.fits.L[compare$lasso.fits.L$Tamp.pos %in% compare$breakpts.L,]))
  #     magnitude.L <- compare$lasso.fits.L$fit[row.L] - compare$lasso.fits.L$fit[(row.L-1)]
  #     abby.T30.break.L <- data.frame(chr=rep(chr,length(compare$breakpts.L)), coord=compare$breakpts.L,magnitude=magnitude.L)
  #     abby.T30.breaks <- rbind(abby.T30.breaks, abby.T30.break.L)
  #     #    #write.table(matrix(compare$breakpts.L, nrow=1), file=paste("abby.T30-chr", chr,
  #     #    #                                                          "-left-arm-breakpts.txt",   sep=""),
  #     #     #          row.names=F, col.names=F)
  #   }
  # }
  # if(!is.na(abby.T30.best.R[chr])) {
  # # if(abby.T30.best.R[chr] == "Piecewise Constant" & !is.na(abby.T30.best.R[chr])) {
  #   if(length(breakpts.R) != 0) {
  #     row.R <- as.integer(rownames(compare$lasso.fits.R[compare$lasso.fits.R$Tamp.pos %in% compare$breakpts.R,]))
  #     magnitude.R <- compare$lasso.fits.R$fit[row.R] - compare$lasso.fits.R$fit[(row.R-1)]
  #     abby.T30.break.R <- data.frame(chr=rep(chr,length(compare$breakpts.R)), coord=compare$breakpts.R,magnitude=magnitude.R)
  #     abby.T30.breaks <- rbind(abby.T30.breaks, abby.T30.break.R)
  #     #    # write.table(matrix(compare$breakpts.R, nrow=1), file=paste("abby.T30-chr", chr,
  #     #   #                                                           "-right-arm-breakpts.txt",   sep=""),
  #     #    #           row.names=F, col.names=F)
  #   }
  # }
  # 
  # compare <- compare_linear_flasso_sep_arms(abby.T37,  "Exp2_T37", chr, centr,
  #                                           all.models[models.considered.index])
  # abby.T37.cvmse.L[,chr] <- compare$cv.mse.L
  # abby.T37.cvmse.R[,chr] <- compare$cv.mse.R
  # abby.T37.best.L[chr] <- compare$fits.best.L
  # abby.T37.best.R[chr] <- compare$fits.best.R
  # 
  # breakpts.L <- compare$breakpts.L
  # breakpts.R <- compare$breakpts.R
  # 
  # # if best is piecewise constant, save breakpts
  # if(!is.na(abby.T37.best.L[chr])) {
  # # if(abby.T37.best.L[chr] == "Piecewise Constant" & !is.na(abby.T37.best.L[chr])) {
  #   if(length(breakpts.L) != 0) {
  #     row.L <- as.integer(rownames(compare$lasso.fits.L[compare$lasso.fits.L$Tamp.pos %in% compare$breakpts.L,]))
  #     magnitude.L <- compare$lasso.fits.L$fit[row.L] - compare$lasso.fits.L$fit[(row.L-1)]
  #     abby.T37.break.L <- data.frame(chr=rep(chr,length(compare$breakpts.L)), coord=compare$breakpts.L,magnitude=magnitude.L)
  #     abby.T37.breaks <- rbind(abby.T37.breaks, abby.T37.break.L)
  #     #    #write.table(matrix(compare$breakpts.L, nrow=1), file=paste("abby.T37-chr", chr,
  #     #    #                                                          "-left-arm-breakpts.txt",   sep=""),
  #     #     #          row.names=F, col.names=F)
  #   }
  # }
  # if(!is.na(abby.T37.best.R[chr])) {
  # # if(abby.T37.best.R[chr] == "Piecewise Constant" & !is.na(abby.T37.best.R[chr])) {
  #   if(length(breakpts.R) != 0) {
  #     row.R <- as.integer(rownames(compare$lasso.fits.R[compare$lasso.fits.R$Tamp.pos %in% compare$breakpts.R,]))
  #     magnitude.R <- compare$lasso.fits.R$fit[row.R] - compare$lasso.fits.R$fit[(row.R-1)]
  #     abby.T37.break.R <- data.frame(chr=rep(chr,length(compare$breakpts.R)), coord=compare$breakpts.R,magnitude=magnitude.R)
  #     abby.T37.breaks <- rbind(abby.T37.breaks, abby.T37.break.R)
  #     #    #write.table(matrix(compare$breakpts.R, nrow=1), file=paste("abby.T37-chr", chr,
  #     #    #                                                          "-right-arm-breakpts.txt",   sep=""),
  #     #     #          row.names=F, col.names=F)
  #   }
  # }
  # 
  # compare <- compare_linear_flasso_sep_arms(abby.low,  "Exp1_T30", chr, centr,
  #                                           all.models[models.considered.index])
  # abby.low.cvmse.L[,chr] <- compare$cv.mse.L
  # abby.low.cvmse.R[,chr] <- compare$cv.mse.R
  # abby.low.best.L[chr] <- compare$fits.best.L
  # abby.low.best.R[chr] <- compare$fits.best.R
  # 
  # breakpts.L <- compare$breakpts.L
  # breakpts.R <- compare$breakpts.R
  # 
  # # if best is piecewise constant, save breakpts
  # if(!is.na(abby.low.best.L[chr])) {
  #   # if(abby.low.best.L[chr] == "Piecewise Constant" & !is.na(abby.low.best.L[chr])) {
  #   if(length(breakpts.L) != 0) {
  #     row.L <- as.integer(rownames(compare$lasso.fits.L[compare$lasso.fits.L$Tamp.pos %in% compare$breakpts.L,]))
  #     magnitude.L <- compare$lasso.fits.L$fit[row.L] - compare$lasso.fits.L$fit[(row.L-1)]
  #     abby.low.break.L <- data.frame(chr=rep(chr,length(compare$breakpts.L)), coord=compare$breakpts.L,magnitude=magnitude.L)
  #     abby.low.breaks <- rbind(abby.low.breaks, abby.low.break.L)
  #     #    #write.table(matrix(compare$breakpts.L, nrow=1), file=paste("abby.low-chr", chr,
  #     #    #                                                          "-left-arm-breakpts.txt",   sep=""),
  #     #     #          row.names=F, col.names=F)
  #   }
  # }
  # if(!is.na(abby.low.best.R[chr])) {
  #   # if(abby.low.best.R[chr] == "Piecewise Constant" & !is.na(abby.low.best.R[chr])) {
  #   if(length(breakpts.R) != 0) {
  #     row.R <- as.integer(rownames(compare$lasso.fits.R[compare$lasso.fits.R$Tamp.pos %in% compare$breakpts.R,]))
  #     magnitude.R <- compare$lasso.fits.R$fit[row.R] - compare$lasso.fits.R$fit[(row.R-1)]
  #     abby.low.break.R <- data.frame(chr=rep(chr,length(compare$breakpts.R)), coord=compare$breakpts.R,magnitude=magnitude.R)
  #     abby.low.breaks <- rbind(abby.low.breaks, abby.low.break.R)
  #     #    # write.table(matrix(compare$breakpts.R, nrow=1), file=paste("abby.low-chr", chr,
  #     #   #                                                           "-right-arm-breakpts.txt",   sep=""),
  #     #    #           row.names=F, col.names=F)
  #   }
  # }
  # 
  # compare <- compare_linear_flasso_sep_arms(abby.high,  "Exp1_T37", chr, centr,
  #                                           all.models[models.considered.index])
  # abby.high.cvmse.L[,chr] <- compare$cv.mse.L
  # abby.high.cvmse.R[,chr] <- compare$cv.mse.R
  # abby.high.best.L[chr] <- compare$fits.best.L
  # abby.high.best.R[chr] <- compare$fits.best.R
  # 
  # breakpts.L <- compare$breakpts.L
  # breakpts.R <- compare$breakpts.R
  # 
  # # if best is piecewise constant, save breakpts
  # if(!is.na(abby.high.best.L[chr])) {
  #   # if(abby.high.best.L[chr] == "Piecewise Constant" & !is.na(abby.high.best.L[chr])) {
  #   if(length(breakpts.L) != 0) {
  #     row.L <- as.integer(rownames(compare$lasso.fits.L[compare$lasso.fits.L$Tamp.pos %in% compare$breakpts.L,]))
  #     magnitude.L <- compare$lasso.fits.L$fit[row.L] - compare$lasso.fits.L$fit[(row.L-1)]
  #     abby.high.break.L <- data.frame(chr=rep(chr,length(compare$breakpts.L)), coord=compare$breakpts.L,magnitude=magnitude.L)
  #     abby.high.breaks <- rbind(abby.high.breaks, abby.high.break.L)
  #     #    #write.table(matrix(compare$breakpts.L, nrow=1), file=paste("abby.high-chr", chr,
  #     #    #                                                          "-left-arm-breakpts.txt",   sep=""),
  #     #     #          row.names=F, col.names=F)
  #   }
  # }
  # if(!is.na(abby.high.best.R[chr])) {
  #   # if(abby.high.best.R[chr] == "Piecewise Constant" & !is.na(abby.high.best.R[chr])) {
  #   if(length(breakpts.R) != 0) {
  #     row.R <- as.integer(rownames(compare$lasso.fits.R[compare$lasso.fits.R$Tamp.pos %in% compare$breakpts.R,]))
  #     magnitude.R <- compare$lasso.fits.R$fit[row.R] - compare$lasso.fits.R$fit[(row.R-1)]
  #     abby.high.break.R <- data.frame(chr=rep(chr,length(compare$breakpts.R)), coord=compare$breakpts.R,magnitude=magnitude.R)
  #     abby.high.breaks <- rbind(abby.high.breaks, abby.high.break.R)
  #     #    #write.table(matrix(compare$breakpts.R, nrow=1), file=paste("abby.high-chr", chr,
  #     #    #                                                          "-right-arm-breakpts.txt",   sep=""),
  #     #     #          row.names=F, col.names=F)
  #   }
  # }
  # ####
  # compare <- compare_linear_flasso_sep_arms(abby.stat,  "MN_stat", chr, centr,
  #                                           all.models[models.considered.index])
  # abby.stat.cvmse.L[,chr] <- compare$cv.mse.L
  # abby.stat.cvmse.R[,chr] <- compare$cv.mse.R
  # abby.stat.best.L[chr] <- compare$fits.best.L
  # abby.stat.best.R[chr] <- compare$fits.best.R
  # 
  # breakpts.L <- compare$breakpts.L
  # breakpts.R <- compare$breakpts.R
  # 
  # # if best is piecewise constant, save breakpts
  # if(!is.na(abby.stat.best.L[chr])) {
  # # if(abby.stat.best.L[chr] == "Piecewise Constant" & !is.na(abby.stat.best.L[chr])) {
  #   if(length(breakpts.L) != 0) {
  #     row.L <- as.integer(rownames(compare$lasso.fits.L[compare$lasso.fits.L$Tamp.pos %in% compare$breakpts.L,]))
  #     magnitude.L <- compare$lasso.fits.L$fit[row.L] - compare$lasso.fits.L$fit[(row.L-1)]
  #     abby.stat.break.L <- data.frame(chr=rep(chr,length(compare$breakpts.L)), coord=compare$breakpts.L,magnitude=magnitude.L)
  #     abby.stat.breaks <- rbind(abby.stat.breaks, abby.stat.break.L)
  #     #    #write.table(matrix(compare$breakpts.L, nrow=1), file=paste("abby.stat-chr", chr,
  #     #    #                                                          "-left-arm-breakpts.txt",   sep=""),
  #     #     #          row.names=F, col.names=F)
  #   }
  # }
  # if(!is.na(abby.stat.best.R[chr])) {
  # # if(abby.stat.best.R[chr] == "Piecewise Constant" & !is.na(abby.stat.best.R[chr])) {
  #   if(length(breakpts.R) != 0) {
  #     row.R <- as.integer(rownames(compare$lasso.fits.R[compare$lasso.fits.R$Tamp.pos %in% compare$breakpts.R,]))
  #     magnitude.R <- compare$lasso.fits.R$fit[row.R] - compare$lasso.fits.R$fit[(row.R-1)]
  #     abby.stat.break.R <- data.frame(chr=rep(chr,length(compare$breakpts.R)), coord=compare$breakpts.R,magnitude=magnitude.R)
  #     abby.stat.breaks <- rbind(abby.stat.breaks, abby.stat.break.R)
  #     #    #write.table(matrix(compare$breakpts.R, nrow=1), file=paste("abby.stat-chr", chr,
  #     #    #                                                          "-right-arm-breakpts.txt",   sep=""),
  #     #     #          row.names=F, col.names=F)
  #   }
  # }
  
}

write.csv(anna.sul.breaks, file='anna.sul.breaks.csv')
write.csv(anna.phos.breaks, file='anna.phos.breaks.csv')
write.csv(anna.gluc.breaks, file='anna.gluc.breaks.csv')
# write.csv(abby.low.breaks, file='abby.low.breaks.csv')
# write.csv(abby.high.breaks, file='abby.high.breaks.csv')
# write.csv(abby.T30.breaks, file='abby.T30.breaks.csv')
# write.csv(abby.T37.breaks, file='abby.T37.breaks.csv')
# write.csv(abby.D.breaks, file='abby.D.breaks.csv')
# write.csv(abby.R.breaks, file='abby.R.breaks.csv')
# write.csv(abby.stat.breaks, file='abby.stat.breaks.csv')


row.labels <- c("Mean CV-MSE, Piecewise Constant Model", 
               "S.E. CV-MSE, Piecewise Constant Model",
               "Mean CV-MSE, Linear Model 1", 
               "S.E. CV-MSE, Linear Model 1", 
               "Mean CV-MSE, Linear Model 2", 
               "S.E. CV-MSE, Linear Model 2", 
               "Mean CV-MSE, Linear Model 3", 
               "S.E. CV-MSE, Linear Model 3",
               "Mean CV-MSE, Linear Model 4", 
               "S.E. CV-MSE, Linear Model 4", 
               "Which Model Fits Best")


anna.sul.results.L <- rbind(data.frame(anna.sul.cvmse.L), anna.sul.best.L)
names(anna.sul.results.L) <- paste("Chr", 1:16)
# Only show the models we care about
anna.sul.results.L <- anna.sul.results.L[c(sort(c(models.considered.index*2,
                                                  models.considered.index*2-1)), 11), ]
rownames(anna.sul.results.L) <- row.labels[c(sort(c(models.considered.index*2,
                                                    models.considered.index*2-1)), 11)]
write.csv(anna.sul.results.L, file="sulfate-results-left-arm.csv")

anna.sul.results.R <- rbind(data.frame(anna.sul.cvmse.R), anna.sul.best.R)
names(anna.sul.results.R) <- paste("Chr", 1:16)
# Only show the models we care about
anna.sul.results.R <- anna.sul.results.R[c(sort(c(models.considered.index*2,
                                                  models.considered.index*2-1)), 11), ]
rownames(anna.sul.results.R) <- row.labels[c(sort(c(models.considered.index*2,
                                                    models.considered.index*2-1)), 11)]
write.csv(anna.sul.results.R, file="sulfate-results-right-arm.csv")


anna.phos.results.L <- rbind(data.frame(anna.phos.cvmse.L), anna.phos.best.L)
names(anna.phos.results.L) <- paste("Chr", 1:16)
# Only show the models we care about
anna.phos.results.L <- anna.phos.results.L[c(sort(c(models.considered.index*2,
                                                    models.considered.index*2-1)), 11), ]
rownames(anna.phos.results.L) <- row.labels[c(sort(c(models.considered.index*2,
                                                     models.considered.index*2-1)), 11)]
write.csv(anna.phos.results.L, file="phosphate-results-left-arm.csv")

anna.phos.results.R <- rbind(data.frame(anna.phos.cvmse.R), anna.phos.best.R)
names(anna.phos.results.R) <- paste("Chr", 1:16)
# Only show the models we care about
anna.phos.results.R <- anna.phos.results.R[c(sort(c(models.considered.index*2,
                                                    models.considered.index*2-1)), 11), ]
rownames(anna.phos.results.R) <- row.labels[c(sort(c(models.considered.index*2,
                                                     models.considered.index*2-1)), 11)]
write.csv(anna.phos.results.R, file="phosphate-results-right-arm.csv")


anna.gluc.results.L <- rbind(data.frame(anna.gluc.cvmse.L), anna.gluc.best.L)
names(anna.gluc.results.L) <- paste("Chr", 1:16)
# Only show the models we care about
anna.gluc.results.L <- anna.gluc.results.L[c(sort(c(models.considered.index*2,
                                                    models.considered.index*2-1)), 11), ]
rownames(anna.gluc.results.L) <- row.labels[c(sort(c(models.considered.index*2,
                                                     models.considered.index*2-1)), 11)]
write.csv(anna.gluc.results.L, file="glucose-results-left-arm.csv")

anna.gluc.results.R <- rbind(data.frame(anna.gluc.cvmse.R), anna.gluc.best.R)
names(anna.gluc.results.R) <- paste("Chr", 1:16)
# Only show the models we care about
anna.gluc.results.R <- anna.gluc.results.R[c(sort(c(models.considered.index*2,
                                                    models.considered.index*2-1)), 11), ]
rownames(anna.gluc.results.R) <- row.labels[c(sort(c(models.considered.index*2,
                                                     models.considered.index*2-1)), 11)]
write.csv(anna.gluc.results.R, file="glucose-results-right-arm.csv")


# abby.D.results.L <- rbind(data.frame(abby.D.cvmse.L), abby.D.best.L)
# names(abby.D.results.L) <- paste("Chr", 1:16)
# # Only show the models we care about
# abby.D.results.L <- abby.D.results.L[c(sort(c(models.considered.index*2,
#                                                   models.considered.index*2-1)), 11), ]
# rownames(abby.D.results.L) <- row.labels[c(sort(c(models.considered.index*2,
#                                                     models.considered.index*2-1)), 11)]
# write.csv(abby.D.results.L, file="MN_D-results-left-arm.csv")
# 
# abby.D.results.R <- rbind(data.frame(abby.D.cvmse.R), abby.D.best.R)
# names(abby.D.results.R) <- paste("Chr", 1:16)
# # Only show the models we care about
# abby.D.results.R <- abby.D.results.R[c(sort(c(models.considered.index*2,
#                                                   models.considered.index*2-1)), 11), ]
# rownames(abby.D.results.R) <- row.labels[c(sort(c(models.considered.index*2,
#                                                     models.considered.index*2-1)), 11)]
# write.csv(abby.D.results.R, file="MN_D-results-right-arm.csv")
# 
# abby.R.results.L <- rbind(data.frame(abby.R.cvmse.L), abby.R.best.L)
# names(abby.R.results.L) <- paste("Chr", 1:16)
# # Only show the models we care about
# abby.R.results.L <- abby.R.results.L[c(sort(c(models.considered.index*2,
#                                                   models.considered.index*2-1)), 11), ]
# rownames(abby.R.results.L) <- row.labels[c(sort(c(models.considered.index*2,
#                                                     models.considered.index*2-1)), 11)]
# write.csv(abby.R.results.L, file="MN_R-results-left-arm.csv")
# 
# abby.R.results.R <- rbind(data.frame(abby.R.cvmse.R), abby.R.best.R)
# names(abby.R.results.R) <- paste("Chr", 1:16)
# # Only show the models we care about
# abby.R.results.R <- abby.R.results.R[c(sort(c(models.considered.index*2,
#                                                   models.considered.index*2-1)), 11), ]
# rownames(abby.R.results.R) <- row.labels[c(sort(c(models.considered.index*2,
#                                                     models.considered.index*2-1)), 11)]
# write.csv(abby.R.results.R, file="MN_R-results-right-arm.csv")
# 
# # abby.T30.results.L <- rbind(data.frame(abby.T30.cvmse.L), abby.T30.best.L)
# # names(abby.T30.results.L) <- paste("Chr", 1:16)
# # # Only show the models we care about
# # abby.T30.results.L <- abby.T30.results.L[c(sort(c(models.considered.index*2,
# #                                               models.considered.index*2-1)), 11), ]
# # rownames(abby.T30.results.L) <- row.labels[c(sort(c(models.considered.index*2,
# #                                                   models.considered.index*2-1)), 11)]
# # write.csv(abby.T30.results.L, file="Exp2_T30-results-left-arm.csv")
# # 
# # abby.T30.results.R <- rbind(data.frame(abby.T30.cvmse.R), abby.T30.best.R)
# # names(abby.T30.results.R) <- paste("Chr", 1:16)
# # # Only show the models we care about
# # abby.T30.results.R <- abby.T30.results.R[c(sort(c(models.considered.index*2,
# #                                               models.considered.index*2-1)), 11), ]
# # rownames(abby.T30.results.R) <- row.labels[c(sort(c(models.considered.index*2,
# #                                                   models.considered.index*2-1)), 11)]
# # write.csv(abby.T30.results.R, file="Exp2_T30-results-right-arm.csv")
# # 
# # 
# # abby.T37.results.L <- rbind(data.frame(abby.T37.cvmse.L), abby.T37.best.L)
# # names(abby.T37.results.L) <- paste("Chr", 1:16)
# # # Only show the models we care about
# # abby.T37.results.L <- abby.T37.results.L[c(sort(c(models.considered.index*2,
# #                                                   models.considered.index*2-1)), 11), ]
# # rownames(abby.T37.results.L) <- row.labels[c(sort(c(models.considered.index*2,
# #                                                     models.considered.index*2-1)), 11)]
# # write.csv(abby.T37.results.L, file="Exp2_T37-results-left-arm.csv")
# # 
# # abby.T37.results.R <- rbind(data.frame(abby.T37.cvmse.R), abby.T37.best.R)
# # names(abby.T37.results.R) <- paste("Chr", 1:16)
# # # Only show the models we care about
# # abby.T37.results.R <- abby.T37.results.R[c(sort(c(models.considered.index*2,
# #                                                   models.considered.index*2-1)), 11), ]
# # rownames(abby.T37.results.R) <- row.labels[c(sort(c(models.considered.index*2,
# #                                                     models.considered.index*2-1)), 11)]
# # write.csv(abby.T37.results.R, file="Exp2_T37-results-right-arm.csv")
# 
# 
# abby.stat.results.L <- rbind(data.frame(abby.stat.cvmse.L), abby.stat.best.L)
# names(abby.stat.results.L) <- paste("Chr", 1:16)
# # Only show the models we care about
# abby.stat.results.L <- abby.stat.results.L[c(sort(c(models.considered.index*2,
#                                                   models.considered.index*2-1)), 11), ]
# rownames(abby.stat.results.L) <- row.labels[c(sort(c(models.considered.index*2,
#                                                     models.considered.index*2-1)), 11)]
# write.csv(abby.stat.results.L, file="MN_stat-results-left-arm.csv")
# 
# abby.stat.results.R <- rbind(data.frame(abby.stat.cvmse.R), abby.stat.best.R)
# names(abby.stat.results.R) <- paste("Chr", 1:16)
# # Only show the models we care about
# abby.stat.results.R <- abby.stat.results.R[c(sort(c(models.considered.index*2,
#                                                   models.considered.index*2-1)), 11), ]
# rownames(abby.stat.results.R) <- row.labels[c(sort(c(models.considered.index*2,
#                                                     models.considered.index*2-1)), 11)]
# write.csv(abby.stat.results.R, file="MN_stat-results-right-arm.csv")
# 
# abby.low.results.L <- rbind(data.frame(abby.low.cvmse.L), abby.low.best.L)
# names(abby.low.results.L) <- paste("Chr", 1:16)
# # Only show the models we care about
# abby.low.results.L <- abby.low.results.L[c(sort(c(models.considered.index*2,
#                                                   models.considered.index*2-1)), 11), ]
# rownames(abby.low.results.L) <- row.labels[c(sort(c(models.considered.index*2,
#                                                     models.considered.index*2-1)), 11)]
# write.csv(abby.low.results.L, file="Exp1_T30-results-left-arm.csv")
# 
# abby.low.results.R <- rbind(data.frame(abby.low.cvmse.R), abby.low.best.R)
# names(abby.low.results.R) <- paste("Chr", 1:16)
# # Only show the models we care about
# abby.low.results.R <- abby.low.results.R[c(sort(c(models.considered.index*2,
#                                                   models.considered.index*2-1)), 11), ]
# rownames(abby.low.results.R) <- row.labels[c(sort(c(models.considered.index*2,
#                                                     models.considered.index*2-1)), 11)]
# write.csv(abby.low.results.R, file="Exp1_T30-results-right-arm.csv")
# 
# 
# abby.high.results.L <- rbind(data.frame(abby.high.cvmse.L), abby.high.best.L)
# names(abby.high.results.L) <- paste("Chr", 1:16)
# # Only show the models we care about
# abby.high.results.L <- abby.high.results.L[c(sort(c(models.considered.index*2,
#                                                   models.considered.index*2-1)), 11), ]
# rownames(abby.high.results.L) <- row.labels[c(sort(c(models.considered.index*2,
#                                                     models.considered.index*2-1)), 11)]
# write.csv(abby.high.results.L, file="Exp1_T37-results-left-arm.csv")
# 
# abby.high.results.R <- rbind(data.frame(abby.high.cvmse.R), abby.high.best.R)
# names(abby.high.results.R) <- paste("Chr", 1:16)
# # Only show the models we care about
# abby.high.results.R <- abby.high.results.R[c(sort(c(models.considered.index*2,
#                                                   models.considered.index*2-1)), 11), ]
# rownames(abby.high.results.R) <- row.labels[c(sort(c(models.considered.index*2,
#                                                     models.considered.index*2-1)), 11)]
# write.csv(abby.high.results.R, file="Exp1_T37-results-right-arm.csv")


save(abby.D.results.L, abby.D.results.R,
     abby.R.results.L, abby.R.results.R,
     # abby.T30.results.L, abby.T30.results.R,
     # abby.T37.results.L, abby.T37.results.R,
     abby.low.results.L, abby.low.results.R,
     abby.high.results.L, abby.high.results.R,
     abby.stat.results.L, abby.stat.results.R,
     file="results.Rdata")

