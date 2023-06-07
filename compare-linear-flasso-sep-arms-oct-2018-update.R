compare_linear_flasso_sep_arms_v1 <- function(dat, dat.varname, chr, centr) { 
  # Create variable indicating the nucleotide position of the 
  # non-telomeric end of the chromosome fragment copy 
  # = start on the left and = chromosome length - stop on the right 
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
    return(list(compare.L=rep(NA, 4), 
                compare.R=rep(NA, 4), 
                flasso.fit.L=NA, flasso.fit.R=NA, lm.fit=NA, lm.fit.plot=NA) )
  } 
  
  dat$armL <- ifelse(dat$arm == "R", 1, 0)
  dat$Tamp.posfit <- dat$Tamp.pos
  dat$Tamp.posfit[dat$arm == "R"] <- dat$Tamp.posfit[dat$arm == "R"] - dat$chrlength[1]
  
  # Fit the linear model 
  lmfit <- lm(as.formula(paste(dat.varname, "~", "Tamp.length", sep="")), data=dat)
  # Create data frame for use in plotting 
  lmfitplot <- data.frame(Tamp.pos= c(dat$Tamp.pos,centr[chr], centr[chr]+1e-16,dat$chrlength[1], 0), 
                          fit=predict(lmfit, data.frame(Tamp.length= c(dat$Tamp.length,  centr[chr], dat$chrlength[1]-centr[chr], 0, 0))))
  
  # Fit the fused lasso using 5-fold CV-MSE.
  # Tuning parameter chosen by 1-se rule. 
  foldid.L <- c(0, rep(seq(1, 5), nrow(dat.L) - 2)[seq(1, nrow(dat.L) - 2)], 0)
  foldid.R <- c(0, rep(seq(1, 5), nrow(dat.R) - 2)[seq(1, nrow(dat.R) - 2)], 0)
  flasso.L <- fusedlasso1d(y=dat.L[, dat.varname],  pos=dat.L$Tamp.pos)
  flasso.R <- fusedlasso1d(y=dat.R[, dat.varname],  pos=dat.R$Tamp.pos)
  cv.flasso.L <-  cv.trendfilter(flasso.L)
  cv.flasso.L$i.2se <- which(cv.flasso.L$lambda == max(cv.flasso.L$lambda[cv.flasso.L$err <= min(cv.flasso.L$err) + 
                                                                            2*cv.flasso.L$se[which.min(cv.flasso.L$err)]]))
  flasso.fit.L <- data.frame(Tamp.pos=dat.L$Tamp.pos, fit=flasso.L$fit[, cv.flasso.L$i.2se])
  cv.flasso.R <-  cv.trendfilter(flasso.R)
  cv.flasso.R$i.2se <- which(cv.flasso.R$lambda == max(cv.flasso.R$lambda[cv.flasso.R$err <= min(cv.flasso.R$err) + 
                                                                            2*cv.flasso.R$se[which.min(cv.flasso.R$err)]]))
  
  flasso.fit.R <- data.frame(Tamp.pos=dat.R$Tamp.pos,  fit=flasso.R$fit[, cv.flasso.R$i.2se]) 
  
  breakpt.L <- flasso.fit.L$Tamp.pos[which(abs(diff(flasso.fit.L$fit,1)) > 1e-15)+1]
  breakpt.R <- flasso.fit.R$Tamp.pos[which(abs(diff(flasso.fit.R$fit,1)) > 1e-15)+1]
  
  # Compute mean and standard error of CV-MSE for linear regression over 5 folds
  mse.lm.L <- rep(0, 5)
  mse.lm.R <- rep(0, 5)
  for(i in 1:5) { 
    cv.lm.fit <- lm(as.formula(paste(dat.varname, "~Tamp.length", sep="")),
                    data = dat[c(which(foldid.L !=i), which(foldid.R !=i) + length(foldid.L)), ])
    cv.lm.fitted.L <- predict(cv.lm.fit,  newdata=dat[foldid.L ==i, ])
    cv.lm.fitted.R <- predict(cv.lm.fit,  newdata=dat[which(foldid.R ==i) + length(foldid.L), ])
    
    mse.lm.L[i] <- mean((cv.lm.fitted.L -  dat[foldid.L ==i, dat.varname])^2) 
    mse.lm.R[i] <- mean((cv.lm.fitted.R - dat[which(foldid.R ==i) + length(foldid.L), dat.varname])^2)
  }  
  
  compare.L <- c( mean(mse.lm.L), 
                  sd(mse.lm.L)/sqrt(5), 
                  cv.flasso.L$err[cv.flasso.L$i.2se], 
                  cv.flasso.L$se[cv.flasso.L$i.2se])
  
  compare.R <- c(mean(mse.lm.R), 
                 sd(mse.lm.R)/sqrt(5), 
                 cv.flasso.R$err[cv.flasso.R$i.2se], 
                 cv.flasso.R$se[cv.flasso.R$i.2se])
  
  # Return results 
  return(list(compare.L=compare.L, compare.R=compare.R, 
              flasso.fit.L=flasso.fit.L,  flasso.fit.R=flasso.fit.R, 
              breakpt.L=breakpt.L, breakpt.R=breakpt.R,
              lm.fit=lmfit, lm.fit.plot=lmfitplot))
}

compare_linear_flasso_sep_arms_v2 <- function(dat, dat.varname, chr, centr) { 
  # Create variable indicating the nucleotide position of the 
  # non-telomeric end of the chromosome fragment copy 
  # = start on the left and = chromosome length - stop on the right 
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
    return(list(compare.L=rep(NA, 4), 
                compare.R=rep(NA, 4), 
                flasso.fit.L=NA, flasso.fit.R=NA, lm.fit=NA, lm.fit.plot=NA) )
  } 
  
  dat$armL <- ifelse(dat$arm == "R", 1, 0)
  dat$Tamp.posfit <- dat$Tamp.pos
  dat$Tamp.posfit[dat$arm == "R"] <- dat$Tamp.posfit[dat$arm == "R"] - dat$chrlength[1]

    # Fit the linear model 
  lmfit <- lm(as.formula(paste(dat.varname, "~", "Tamp.pos*armL", sep="")), data=dat)
  # Create data frame for use in plotting 
  lmfitplot <- data.frame(Tamp.pos= c(dat$Tamp.pos,centr[chr], centr[chr]+1e-16,
                                      dat$chrlength[1], 0), 
                          fit=predict(lmfit, 
                                      data.frame(Tamp.pos= c(dat$Tamp.pos,centr[chr], centr[chr]+1e-16,
                                                             dat$chrlength[1], 0), 
                                                 armL=c(dat$armL, 0, 1, 1, 0))))
  
  # Fit the fused lasso using 5-fold CV-MSE.
  # Tuning parameter chosen by 1-se rule. 
  foldid.L <- c(0, rep(seq(1, 5), nrow(dat.L) - 2)[seq(1, nrow(dat.L) - 2)], 0)
  foldid.R <- c(0, rep(seq(1, 5), nrow(dat.R) - 2)[seq(1, nrow(dat.R) - 2)], 0)
  flasso.L <- fusedlasso1d(y=dat.L[, dat.varname],  pos=dat.L$Tamp.pos)
  flasso.R <- fusedlasso1d(y=dat.R[, dat.varname],  pos=dat.R$Tamp.pos)
  cv.flasso.L <-  cv.trendfilter(flasso.L)
  cv.flasso.L$i.2se <- which(cv.flasso.L$lambda == max(cv.flasso.L$lambda[cv.flasso.L$err <= min(cv.flasso.L$err) + 
                                                                            2*cv.flasso.L$se[which.min(cv.flasso.L$err)]]))
  flasso.fit.L <- data.frame(Tamp.pos=dat.L$Tamp.pos, fit=flasso.L$fit[, cv.flasso.L$i.2se])
  cv.flasso.R <-  cv.trendfilter(flasso.R)
  cv.flasso.R$i.2se <- which(cv.flasso.R$lambda == max(cv.flasso.R$lambda[cv.flasso.R$err <= min(cv.flasso.R$err) + 
                                                                            2*cv.flasso.R$se[which.min(cv.flasso.R$err)]]))
  flasso.fit.R <- data.frame(Tamp.pos=dat.R$Tamp.pos,  fit=flasso.R$fit[, cv.flasso.R$i.2se]) 
  
  breakpt.L <- flasso.fit.L$Tamp.pos[which(abs(diff(flasso.fit.L$fit,1)) > 1e-15)+1]
  breakpt.R <- flasso.fit.R$Tamp.pos[which(abs(diff(flasso.fit.R$fit,1)) > 1e-15)+1]
  
  # Compute mean and standard error of CV-MSE for linear regression over 5 folds
  mse.lm.L <- rep(0, 5)
  mse.lm.R <- rep(0, 5)
  for(i in 1:5) { 
    cv.lm.fit <- lm(as.formula(paste(dat.varname, "~", "Tamp.pos*armL", sep="")),
                    data = dat[c(which(foldid.L !=i), which(foldid.R !=i) + length(foldid.L)), ])
    cv.lm.fitted.L <- predict(cv.lm.fit,  newdata=dat[foldid.L ==i, ])
    cv.lm.fitted.R <- predict(cv.lm.fit,  newdata=dat[which(foldid.R ==i) + length(foldid.L), ])
    
    mse.lm.L[i] <- mean((cv.lm.fitted.L -  dat[foldid.L ==i, dat.varname])^2) 
    mse.lm.R[i] <- mean((cv.lm.fitted.R - dat[which(foldid.R ==i) + length(foldid.L), dat.varname])^2)
  }  
  
  compare.L <- c( mean(mse.lm.L), 
                  sd(mse.lm.L)/sqrt(5), 
                  cv.flasso.L$err[cv.flasso.L$i.2se], 
                  cv.flasso.L$se[cv.flasso.L$i.2se])
  
  compare.R <- c(mean(mse.lm.R), 
                 sd(mse.lm.R)/sqrt(5), 
                 cv.flasso.R$err[cv.flasso.R$i.2se], 
                 cv.flasso.R$se[cv.flasso.R$i.2se])
  
  # Return results 
  return(list(compare.L=compare.L, compare.R=compare.R, 
              flasso.fit.L=flasso.fit.L,  flasso.fit.R=flasso.fit.R, 
              breakpt.L=breakpt.L, breakpt.R=breakpt.R,
              lm.fit=lmfit, lm.fit.plot=lmfitplot))
}

compare_linear_flasso_sep_arms_v3 <- function(dat, dat.varname, chr, centr) { 
  # Create variable indicating the nucleotide position of the 
  # non-telomeric end of the chromosome fragment copy 
  # = start on the left and = chromosome length - stop on the right 
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
    return(list(compare.L=rep(NA, 4), 
                compare.R=rep(NA, 4), 
                flasso.fit.L=NA, flasso.fit.R=NA, lm.fit=NA, lm.fit.plot=NA) )
  } 
  
  dat$armL <- ifelse(dat$arm == "R", 1, 0)
  dat$Tamp.posfit <- dat$Tamp.pos
  dat$Tamp.posfit[dat$arm == "R"] <- -dat$Tamp.posfit[dat$arm == "R"]
  # Fit the linear model 
  lmfit <- lm(as.formula(paste(dat.varname, "~", "Tamp.posfit + armL", sep="")), data=dat)
  # Create data frame for use in plotting 
  lmfitplot <- data.frame(Tamp.pos= c(dat$Tamp.pos,centr[chr], centr[chr]+1e-16, 
                                      0, dat$chrlength[1]), 
                          fit=predict(lmfit, 
                                      data.frame(Tamp.posfit= c(dat$Tamp.posfit,centr[chr], -(centr[chr]+1e-16), 0, 
                                                                -dat$chrlength[1]), 
                                                 armL=c(dat$armL, 0, 1, 0, 1))))
  
  # Fit the fused lasso using 5-fold CV-MSE.
  # Tuning parameter chosen by 1-se rule. 
  foldid.L <- c(0, rep(seq(1, 5), nrow(dat.L) - 2)[seq(1, nrow(dat.L) - 2)], 0)
  foldid.R <- c(0, rep(seq(1, 5), nrow(dat.R) - 2)[seq(1, nrow(dat.R) - 2)], 0)
  flasso.L <- fusedlasso1d(y=dat.L[, dat.varname],  pos=dat.L$Tamp.pos)
  flasso.R <- fusedlasso1d(y=dat.R[, dat.varname],  pos=dat.R$Tamp.pos)
  cv.flasso.L <-  cv.trendfilter(flasso.L)
  cv.flasso.L$i.2se <- which(cv.flasso.L$lambda == max(cv.flasso.L$lambda[cv.flasso.L$err <= min(cv.flasso.L$err) + 
                                                                            2*cv.flasso.L$se[which.min(cv.flasso.L$err)]]))
  flasso.fit.L <- data.frame(Tamp.pos=dat.L$Tamp.pos, fit=flasso.L$fit[, cv.flasso.L$i.2se])
  cv.flasso.R <-  cv.trendfilter(flasso.R)
  
  cv.flasso.R$i.2se <- which(cv.flasso.R$lambda == max(cv.flasso.R$lambda[cv.flasso.R$err <= min(cv.flasso.R$err) + 
                                                                            2*cv.flasso.R$se[which.min(cv.flasso.R$err)]]))
  
  flasso.fit.R <- data.frame(Tamp.pos=dat.R$Tamp.pos,  fit=flasso.R$fit[, cv.flasso.R$i.2se]) 
  
  breakpt.L <- flasso.fit.L$Tamp.pos[which(abs(diff(flasso.fit.L$fit,1)) > 1e-15)+1]
  breakpt.R <- flasso.fit.R$Tamp.pos[which(abs(diff(flasso.fit.R$fit,1)) > 1e-15)+1]
  
  # Compute mean and standard error of CV-MSE for linear regression over 5 folds
  mse.lm.L <- rep(0, 5)
  mse.lm.R <- rep(0, 5)
  for(i in 1:5) { 
    cv.lm.fit <- lm(as.formula(paste(dat.varname, "~", "Tamp.posfit + armL", sep="")),
                    data = dat[c(which(foldid.L !=i), which(foldid.R !=i) + length(foldid.L)), ])
    cv.lm.fitted.L <- predict(cv.lm.fit,  newdata=dat[foldid.L ==i, ])
    cv.lm.fitted.R <- predict(cv.lm.fit,  newdata=dat[which(foldid.R ==i) + length(foldid.L), ])
    
    mse.lm.L[i] <- mean((cv.lm.fitted.L -  dat[foldid.L ==i, dat.varname])^2) 
    mse.lm.R[i] <- mean((cv.lm.fitted.R - dat[which(foldid.R ==i) + length(foldid.L), dat.varname])^2)
  }  
  
  compare.L <- c( mean(mse.lm.L), 
                  sd(mse.lm.L)/sqrt(5), 
                  cv.flasso.L$err[cv.flasso.L$i.2se], 
                  cv.flasso.L$se[cv.flasso.L$i.2se])
  
  compare.R <- c(mean(mse.lm.R), 
                 sd(mse.lm.R)/sqrt(5), 
                 cv.flasso.R$err[cv.flasso.R$i.2se], 
                 cv.flasso.R$se[cv.flasso.R$i.2se])
  
  # Return results 
  return(list(compare.L=compare.L, compare.R=compare.R, 
              flasso.fit.L=flasso.fit.L,  flasso.fit.R=flasso.fit.R, 
              breakpt.L=breakpt.L, breakpt.R=breakpt.R,
              lm.fit=lmfit, lm.fit.plot=lmfitplot))
}

#####################################################

compare_linear_flasso_sep_arms_v4 <- function(dat, dat.varname, chr, centr) { 
  # Create variable indicating the nucleotide position of the 
  # non-telomeric end of the chromosome fragment copy 
  # = start on the left and = chromosome length - stop on the right 
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
    return(list(compare.L=rep(NA, 4), 
                compare.R=rep(NA, 4), 
                flasso.fit.L=NA, flasso.fit.R=NA, lm.fit=NA, lm.fit.plot=NA) )
  } 
  
  
  dat$armL <- ifelse(dat$arm == "R", 1, 0)
  dat$Tamp.posfit <- dat$Tamp.pos
  dat$Tamp.posfit[dat$arm == "R"] <- dat$Tamp.posfit[dat$arm == "R"] - dat$chrlength[1]
  # Fit the linear model 
  lmfit <- lm(as.formula(paste(dat.varname, "~", "Tamp.posfit:armL + Tamp.posfit", sep="")), data=dat)
  # Create data frame for use in plotting 
  lmfitplot <- data.frame(Tamp.pos= c(dat$Tamp.pos,centr[chr], centr[chr]+1e-16, dat$chrlength[1], 0), 
                          fit=predict(lmfit, data.frame(Tamp.posfit= c(dat$Tamp.posfit,centr[chr], centr[chr]+1e-16 - dat$chrlength[1], 0, 0), armL=c(dat$armL, 0, 1, 1, 0))))
  
  # Fit the fused lasso using 5-fold CV-MSE.
  # Tuning parameter chosen by 2-se rule. 
  foldid.L <- c(0, rep(seq(1, 5), nrow(dat.L) - 2)[seq(1, nrow(dat.L) - 2)], 0)
  foldid.R <- c(0, rep(seq(1, 5), nrow(dat.R) - 2)[seq(1, nrow(dat.R) - 2)], 0)
  flasso.L <- fusedlasso1d(y=dat.L[, dat.varname],  pos=dat.L$Tamp.pos)
  flasso.R <- fusedlasso1d(y=dat.R[, dat.varname],  pos=dat.R$Tamp.pos)
  cv.flasso.L <-  cv.trendfilter(flasso.L)
  cv.flasso.L$i.2se <- which(cv.flasso.L$lambda == max(cv.flasso.L$lambda[cv.flasso.L$err <= min(cv.flasso.L$err) + 
                                                                            2*cv.flasso.L$se[which.min(cv.flasso.L$err)]]))
  
  flasso.fit.L <- data.frame(Tamp.pos=dat.L$Tamp.pos, fit=flasso.L$fit[, cv.flasso.L$i.2se])
  cv.flasso.R <-  cv.trendfilter(flasso.R)
  cv.flasso.R$i.2se <- which(cv.flasso.R$lambda == max(cv.flasso.R$lambda[cv.flasso.R$err <= min(cv.flasso.R$err) + 
                                                                            2*cv.flasso.R$se[which.min(cv.flasso.R$err)]]))
  flasso.fit.R <- data.frame(Tamp.pos=dat.R$Tamp.pos,  fit=flasso.R$fit[, cv.flasso.R$i.2se]) 
  
  breakpt.L <- flasso.fit.L$Tamp.pos[which(abs(diff(flasso.fit.L$fit,1)) > 1e-15)+1]
  breakpt.R <- flasso.fit.R$Tamp.pos[which(abs(diff(flasso.fit.R$fit,1)) > 1e-15)+1]
  
  # Compute mean and standard error of CV-MSE for linear regression over 5 folds
  mse.lm.L <- rep(0, 5)
  mse.lm.R <- rep(0, 5)
  for(i in 1:5) { 
    cv.lm.fit <- lm(as.formula(paste(dat.varname, "~", "Tamp.posfit:armL + Tamp.posfit", sep="")),
                    data = dat[c(which(foldid.L !=i), which(foldid.R !=i) + length(foldid.L)), ])
    cv.lm.fitted.L <- predict(cv.lm.fit,  newdata=dat[foldid.L ==i, ])
    cv.lm.fitted.R <- predict(cv.lm.fit,  newdata=dat[which(foldid.R ==i) + length(foldid.L), ])
    
    mse.lm.L[i] <- mean((cv.lm.fitted.L -  dat[foldid.L ==i, dat.varname])^2) 
    mse.lm.R[i] <- mean((cv.lm.fitted.R - dat[which(foldid.R ==i) + length(foldid.L), dat.varname])^2)
  }  
  
  compare.L <- c( mean(mse.lm.L), 
                  sd(mse.lm.L)/sqrt(5), 
                  cv.flasso.L$err[cv.flasso.L$i.2se], 
                  cv.flasso.L$se[cv.flasso.L$i.2se])
  
  compare.R <- c(mean(mse.lm.R), 
                 sd(mse.lm.R)/sqrt(5), 
                 cv.flasso.R$err[cv.flasso.R$i.2se], 
                 cv.flasso.R$se[cv.flasso.R$i.2se])
  
  # Return results 
  return(list(compare.L=compare.L, compare.R=compare.R, 
              flasso.fit.L=flasso.fit.L,  flasso.fit.R=flasso.fit.R, 
              breakpt.L=breakpt.L, breakpt.R=breakpt.R,
              lm.fit=lmfit, lm.fit.plot=lmfitplot))
}

compare_linear_flasso_sep_arms <- function(dat, dat.varname, chr, centr, 
                                  models=c("Piecewise Constant", "Linear Model 1", "Linear Model 2", 
                                           "Linear Model 3", "Linear Model 4")) { 
  cv.errs1 <- compare_linear_flasso_sep_arms_v1(dat, dat.varname, chr, centr)
  cv.errs2 <- compare_linear_flasso_sep_arms_v2(dat, dat.varname, chr, centr)
  cv.errs3 <- compare_linear_flasso_sep_arms_v3(dat, dat.varname, chr, centr)
  cv.errs4 <- compare_linear_flasso_sep_arms_v4(dat, dat.varname, chr, centr)
  
  
  if(all(c(!is.na(cv.errs1$compare.L), !is.na(cv.errs1$compare.R)))) { 
    # mean cv err for all models 
    mean.L <- c(cv.errs1$compare.L[3], cv.errs1$compare.L[1], cv.errs2$compare.L[1], 
              cv.errs3$compare.L[1], cv.errs4$compare.L[1])
    # mean cv err for the models you're interested in
    all.models <- c("Piecewise Constant", "Linear Model 1", "Linear Model 2", 
                    "Linear Model 3", "Linear Model 4")
    mean.L <- mean.L[all.models %in% models]
    
    best.L <- models[which.min(mean.L)]
    
    compare.L <- c(cv.errs1$compare.L[3:4], cv.errs1$compare.L[1:2], cv.errs2$compare.L[1:2], 
                 cv.errs3$compare.L[1:2], cv.errs4$compare.L[1:2]) 
    
    # having it return breakpoints regardless of fit
    breakpts.L <- cv.errs1$breakpt.L
    lasso.fits.L <- cv.errs1$flasso.fit.L
    # if(best.L == "Piecewise Constant") { 
    #   breakpts.L <- cv.errs1$breakpt.L
    #   lasso.fits.L <- cv.errs1$flasso.fit.L
    # } else { 
    #   breakpts.L <- NA
    #   lasso.fits.L <- NA   
    # }
    
    names(compare.L) <- c("Mean CV-MSE, Piecewise Constant Model", 
                        "S.E. CV-MSE, Piecewise Constant Model",
                        "Mean CV-MSE, Linear Model 1", 
                        "S.E. CV-MSE, Linear Model 1", 
                        "Mean CV-MSE, Linear Model 2", 
                        "S.E. CV-MSE, Linear Model 2", 
                        "Mean CV-MSE, Linear Model 3", 
                        "S.E. CV-MSE, Linear Model 3",
                        "Mean CV-MSE, Linear Model 4", 
                        "S.E. CV-MSE, Linear Model 4")
    
    # mean cv err for all models 
    mean.R <- c(cv.errs1$compare.R[3], cv.errs1$compare.R[1], cv.errs2$compare.R[1], 
                cv.errs3$compare.R[1], cv.errs4$compare.R[1])
    # mean cv err for the models you're interested in
    all.models <- c("Piecewise Constant", "Linear Model 1", "Linear Model 2", 
                    "Linear Model 3", "Linear Model 4")
    mean.R <- mean.R[all.models %in% models]
    
    best.R <- models[which.min(mean.R)]
    
    compare.R <- c(cv.errs1$compare.R[3:4], cv.errs1$compare.R[1:2], cv.errs2$compare.R[1:2], 
                   cv.errs3$compare.R[1:2], cv.errs4$compare.R[1:2]) 
    
    # having it return breakpoints regardless of fit
    breakpts.R <- cv.errs1$breakpt.R
    lasso.fits.R <- cv.errs1$flasso.fit.R
    # if(best.R == "Piecewise Constant") { 
    #   breakpts.R <- cv.errs1$breakpt.R
    #   lasso.fits.R <- cv.errs1$flasso.fit.R
    # } else { 
    #   breakpts.R <- NA
    #   lasso.fits.R <- NA   
    # }
    
    names(compare.R) <- c("Mean CV-MSE, Piecewise Constant Model", 
                          "S.E. CV-MSE, Piecewise Constant Model",
                          "Mean CV-MSE, Linear Model 1", 
                          "S.E. CV-MSE, Linear Model 1", 
                          "Mean CV-MSE, Linear Model 2", 
                          "S.E. CV-MSE, Linear Model 2", 
                          "Mean CV-MSE, Linear Model 3", 
                          "S.E. CV-MSE, Linear Model 3",
                          "Mean CV-MSE, Linear Model 4", 
                          "S.E. CV-MSE, Linear Model 4")
    return(list(cv.mse.L=compare.L, fits.best.L=best.L, models.considered=models, 
                breakpts.L=breakpts.L, 
                cv.mse.R=compare.R, fits.best.R=best.R, 
                breakpts.R=breakpts.R,lasso.fits.L=lasso.fits.L,lasso.fits.R=lasso.fits.R))
  } else { 
    return(list(cv.mse.L=rep(NA, 10), 
                cv.mse.R=rep(NA, 10), fits.best.L=NA,
                fits.best.R=NA, models.considered=NA, 
                breakpts.L=NA, 
                breakpts.R=NA,lasso.fits.L=NA,lasso.fits.R=NA))  
  }
}

