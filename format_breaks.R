# open data files for gene coordinate info
setwd('/Users/AbbyKeller/Desktop/March2020_WFH/Witten_model/data/')

glucose <- read.csv('171103_AS_glucose_slopes_locations.csv')
phosphate <- read.csv('171103_AS_phosphate_slopes_locations.csv')
sulfate <- read.csv('171103_AS_sulfate_slopes_locations.csv')
DMSO <- read.csv('MN_D_slopes_se_normF.csv')
Rad <- read.csv('MN_R_slopes_se_normF.csv')
low <- read.csv('MN_low_slopes_se_normF.csv')
high <- read.csv('MN_high_slopes_se_normF.csv')
stat <- read.csv('MN_stat_slopes_se_normF.csv')

# order data to read left to right (increasing coordinate) along chromosome
glucose <- glucose[order(glucose$chr,glucose$number),]
phosphate <- phosphate[order(phosphate$chr,phosphate$number),]
sulfate <- sulfate[order(sulfate$chr,sulfate$number),]
DMSO <- DMSO[order(DMSO$chr,DMSO$number),]
Rad <- Rad[order(Rad$chr,Rad$number),]
low <- low[order(low$chr,low$number),]
high <- high[order(high$chr,high$number),]
stat <- stat[order(stat$chr,stat$number),]

# read in breakpoint data
setwd('/Users/AbbyKeller/Desktop/March2020_WFH/Witten_model/analysis-results/sep-arms-compare-PC-all-oct-2018-update-2se/')

gluc.breaks <- read.csv('anna.gluc.breaks.csv')
phos.breaks <- read.csv('anna.phos.breaks.csv')
sulf.breaks <- read.csv('anna.sul.breaks.csv')
D.breaks <- read.csv('abby.D.breaks.csv')
R.breaks <- read.csv('abby.R.breaks.csv')
low.breaks <- read.csv('abby.low.breaks.csv')
high.breaks <- read.csv('abby.high.breaks.csv')
stat.breaks <- read.csv('abby.stat.breaks.csv')

glu.frame <- data.frame(chromosome=character(),arm=character(),break.start=character(),break.stop=character(),magnitude=character(),step.type=character())
j <- 1
for (row in 1:nrow(glucose)) {
  if (j <= nrow(gluc.breaks)) {
    arm <- substr(as.character(glucose$genes[row]),3,3)
    if (arm == 'L') {
     match <- glucose$start[row]
    } else {
      match <- glucose$stop[row]
    }
    print(match)
    if (match == gluc.breaks$coord[j]) {
      print('hey!')
      print(glucose$number[row])
      chr <- glucose$chr[row]
      break.start <- glucose$start[(row-1)]
      break.end <- glucose$start[row]
      magnitude <- gluc.breaks$magnitude[j]
      if (arm == 'L') {
        if (magnitude > 0) {
          break.type <- 'upstep'
        } else {
          break.type <- 'downstep'
        }
      } else {
        if (magnitude < 0) {
          break.type <- 'upstep'
        } else {
          break.type <- 'downstep'
        }
      }
      new.frame <- data.frame(chromosome=chr,arm=arm,break.start=break.start,break.stop=break.end,magnitude=magnitude,step.type=break.type)
      glu.frame <- rbind(glu.frame, new.frame)
      j <- j+1
    } else {
      print('nope')
    }
  }
}

pho.frame <- data.frame(chromosome=character(),arm=character(),break.start=character(),break.stop=character(),magnitude=character(),step.type=character())
j <- 1
for (row in 1:nrow(phosphate)) {
  if (j <= nrow(phos.breaks)) {
    arm <- substr(as.character(phosphate$genes[row]),3,3)
    if (arm == 'L') {
      match <- phosphate$start[row]
    } else {
      match <- phosphate$stop[row]
    }
    print(match)
    if (match == phos.breaks$coord[j]) {
      print('hey!')
      print(phosphate$number[row])
      chr <- phosphate$chr[row]
      break.start <- phosphate$start[(row-1)]
      break.end <- phosphate$start[row]
      magnitude <- phos.breaks$magnitude[j]
      if (arm == 'L') {
        if (magnitude > 0) {
          break.type <- 'upstep'
        } else {
          break.type <- 'downstep'
        }
      } else {
        if (magnitude < 0) {
          break.type <- 'upstep'
        } else {
          break.type <- 'downstep'
        }
      }
      new.frame <- data.frame(chromosome=chr,arm=arm,break.start=break.start,break.stop=break.end,magnitude=magnitude,step.type=break.type)
      pho.frame <- rbind(pho.frame, new.frame)
      j <- j+1
    } else {
      print('nope')
    }
  }
}

sul.frame <- data.frame(chromosome=character(),arm=character(),break.start=character(),break.stop=character(),magnitude=character(),step.type=character())
j <- 1
for (row in 1:nrow(sulfate)) {
  if (j <= nrow(sulf.breaks)) {
    arm <- substr(as.character(sulfate$genes[row]),3,3)
    if (arm == 'L') {
      match <- sulfate$start[row]
    } else {
      match <- sulfate$stop[row]
    }
    print(match)
    if (match == sulf.breaks$coord[j]) {
      print('hey!')
      print(sulfate$number[row])
      chr <- sulfate$chr[row]
      break.start <- sulfate$start[(row-1)]
      break.end <- sulfate$start[row]
      magnitude <- sulf.breaks$magnitude[j]
      if (arm == 'L') {
        if (magnitude > 0) {
          break.type <- 'upstep'
        } else {
          break.type <- 'downstep'
        }
      } else {
        if (magnitude < 0) {
          break.type <- 'upstep'
        } else {
          break.type <- 'downstep'
        }
      }
      new.frame <- data.frame(chromosome=chr,arm=arm,break.start=break.start,break.stop=break.end,magnitude=magnitude,step.type=break.type)
      sul.frame <- rbind(sul.frame, new.frame)
      j <- j+1
    } else {
      print('nope')
    }
  }
}

D.frame <- data.frame(chromosome=character(),arm=character(),break.start=character(),break.stop=character(),magnitude=character(),step.type=character())
j <- 1
for (row in 1:nrow(DMSO)) {
  if (j <= nrow(D.breaks)) {
    arm <- substr(as.character(DMSO$genes[row]),3,3)
    if (arm == 'L') {
      match <- DMSO$start[row]
    } else {
      match <- DMSO$stop[row]
    }
    print(match)
    if (match == D.breaks$coord[j]) {
      print('hey!')
      print(DMSO$number[row])
      chr <- DMSO$chr[row]
      break.start <- DMSO$start[(row-1)]
      break.end <- DMSO$start[row]
      magnitude <- D.breaks$magnitude[j]
      if (arm == 'L') {
        if (magnitude > 0) {
          break.type <- 'upstep'
        } else {
          break.type <- 'downstep'
        }
      } else {
        if (magnitude < 0) {
          break.type <- 'upstep'
        } else {
          break.type <- 'downstep'
        }
      }
      new.frame <- data.frame(chromosome=chr,arm=arm,break.start=break.start,break.stop=break.end,magnitude=magnitude,step.type=break.type)
      D.frame <- rbind(D.frame, new.frame)
      j <- j+1
    } else {
      print('nope')
    }
  }
}

R.frame <- data.frame(chromosome=character(),arm=character(),break.start=character(),break.stop=character(),magnitude=character(),step.type=character())
j <- 1
for (row in 1:nrow(Rad)) {
  if (j <= nrow(R.breaks)) {
    arm <- substr(as.character(Rad$genes[row]),3,3)
    if (arm == 'L') {
      match <- Rad$start[row]
    } else {
      match <- Rad$stop[row]
    }
    print(match)
    if (match == R.breaks$coord[j]) {
      print('hey!')
      print(Rad$number[row])
      chr <- Rad$chr[row]
      break.start <- Rad$start[(row-1)]
      break.end <- Rad$start[row]
      magnitude <- R.breaks$magnitude[j]
      if (arm == 'L') {
        if (magnitude > 0) {
          break.type <- 'upstep'
        } else {
          break.type <- 'downstep'
        }
      } else {
        if (magnitude < 0) {
          break.type <- 'upstep'
        } else {
          break.type <- 'downstep'
        }
      }
      new.frame <- data.frame(chromosome=chr,arm=arm,break.start=break.start,break.stop=break.end,magnitude=magnitude,step.type=break.type)
      R.frame <- rbind(R.frame, new.frame)
      j <- j+1
    } else {
      print('nope')
    }
  }
}

low.frame <- data.frame(chromosome=character(),arm=character(),break.start=character(),break.stop=character(),magnitude=character(),step.type=character())
j <- 1
for (row in 1:nrow(low)) {
  if (j <= nrow(low.breaks)) {
    arm <- substr(as.character(low$genes[row]),3,3)
    if (arm == 'L') {
      match <- low$start[row]
    } else {
      match <- low$stop[row]
    }
    print(match)
    if (match == low.breaks$coord[j]) {
      print('hey!')
      print(low$number[row])
      chr <- low$chr[row]
      break.start <- low$start[(row-1)]
      break.end <- low$start[row]
      magnitude <- low.breaks$magnitude[j]
      if (arm == 'L') {
        if (magnitude > 0) {
          break.type <- 'upstep'
        } else {
          break.type <- 'downstep'
        }
      } else {
        if (magnitude < 0) {
          break.type <- 'upstep'
        } else {
          break.type <- 'downstep'
        }
      }
      new.frame <- data.frame(chromosome=chr,arm=arm,break.start=break.start,break.stop=break.end,magnitude=magnitude,step.type=break.type)
      low.frame <- rbind(low.frame, new.frame)
      j <- j+1
    } else {
      print('nope')
    }
  }
}

high.frame <- data.frame(chromosome=character(),arm=character(),break.start=character(),break.stop=character(),magnitude=character(),step.type=character())
j <- 1
for (row in 1:nrow(high)) {
  if (j <= nrow(high.breaks)) {
    arm <- substr(as.character(high$genes[row]),3,3)
    if (arm == 'L') {
      match <- high$start[row]
    } else {
      match <- high$stop[row]
    }
    print(match)
    if (match == high.breaks$coord[j]) {
      print('hey!')
      print(high$number[row])
      chr <- high$chr[row]
      break.start <- high$start[(row-1)]
      break.end <- high$start[row]
      magnitude <- high.breaks$magnitude[j]
      if (arm == 'L') {
        if (magnitude > 0) {
          break.type <- 'upstep'
        } else {
          break.type <- 'downstep'
        }
      } else {
        if (magnitude < 0) {
          break.type <- 'upstep'
        } else {
          break.type <- 'downstep'
        }
      }
      new.frame <- data.frame(chromosome=chr,arm=arm,break.start=break.start,break.stop=break.end,magnitude=magnitude,step.type=break.type)
      high.frame <- rbind(high.frame, new.frame)
      j <- j+1
    } else {
      print('nope')
    }
  }
}

stat.frame <- data.frame(chromosome=character(),arm=character(),break.start=character(),break.stop=character(),magnitude=character(),step.type=character())
j <- 1
for (row in 1:nrow(stat)) {
  if (j <= nrow(stat.breaks)) {
    arm <- substr(as.character(stat$genes[row]),3,3)
    if (arm == 'L') {
      match <- stat$start[row]
    } else {
      match <- stat$stop[row]
    }
    print(match)
    if (match == stat.breaks$coord[j]) {
      print('hey!')
      print(stat$number[row])
      chr <- stat$chr[row]
      break.start <- stat$start[(row-1)]
      break.end <- stat$start[row]
      magnitude <- stat.breaks$magnitude[j]
      if (arm == 'L') {
        if (magnitude > 0) {
          break.type <- 'upstep'
        } else {
          break.type <- 'downstep'
        }
      } else {
        if (magnitude < 0) {
          break.type <- 'upstep'
        } else {
          break.type <- 'downstep'
        }
      }
      new.frame <- data.frame(chromosome=chr,arm=arm,break.start=break.start,break.stop=break.end,magnitude=magnitude,step.type=break.type)
      stat.frame <- rbind(stat.frame, new.frame)
      j <- j+1
    } else {
      print('nope')
    }
  }
}

write.csv(sul.frame, file = 'sulfate_breakMags.csv')
write.csv(pho.frame, file = 'phosphate_breakMags.csv')
write.csv(glu.frame, file = 'glucose_breakMags.csv')
# write.csv(D.frame,file = 'DMSO_breakMags.csv')
# write.csv(R.frame,file = 'Radicicol_breakMags.csv')
# write.csv(low.frame,file = 'low_breakMags.csv')
# write.csv(high.frame,file = 'high_breakMags.csv')
# write.csv(stat.frame,file = 'stat_breakMags.csv')


glu.magadj <- NULL
for (i in 1:length(glu.frame$magnitude)){
  if (glu.frame$step.type[i] == 'upstep'){
    glu.magadj = c(glu.magadj,abs(glu.frame$magnitude[i]))
  }else if (glu.frame$step.type[i] == 'downstep'){
    glu.magadj = c(glu.magadj,(0-abs(glu.frame$magnitude[i])))
  }
  i = i+1
}
hist(glu.magadj,breaks = 30, col = '#e41a1c', main = 'Glucose', xlim = c(-0.3,0.3),xlab = 'Step Magnitude')
abline(v=0.05)
abline(v=-0.05)

pho.magadj <- NULL
for (i in 1:length(pho.frame$magnitude)){
  if (pho.frame$step.type[i] == 'upstep'){
    pho.magadj = c(pho.magadj,abs(pho.frame$magnitude[i]))
  }else if (pho.frame$step.type[i] == 'downstep'){
    pho.magadj = c(pho.magadj,(0-abs(pho.frame$magnitude[i])))
  }
  i = i+1
}
hist(pho.magadj,breaks = 30, col = '#4daf4a', main = 'Phosphate', xlim = c(-0.3,0.3),xlab = 'Step Magnitude')
abline(v=0.05)
abline(v=-0.05)

sul.magadj <- NULL
for (i in 1:length(sul.frame$magnitude)){
  if (sul.frame$step.type[i] == 'upstep'){
    sul.magadj = c(sul.magadj,abs(sul.frame$magnitude[i]))
  }else if (sul.frame$step.type[i] == 'downstep'){
    sul.magadj = c(sul.magadj,(0-abs(pho.frame$magnitude[i])))
  }
  i = i+1
}
hist(sul.magadj,breaks = 30, col = '#377eb8', main = 'Sulfate', xlim = c(-0.3,0.3),xlab = 'Step Magnitude')
abline(v=0.05)
abline(v=-0.05)


D.magadj <- NULL
for (i in 1:length(D.frame$magnitude)){
  if (D.frame$step.type[i] == 'upstep'){
    D.magadj = c(D.magadj,abs(D.frame$magnitude[i]))
  }else if (D.frame$step.type[i] == 'downstep'){
    D.magadj = c(D.magadj,(0-abs(D.frame$magnitude[i])))
  }
  i = i+1
}
hist(D.magadj,breaks = 30, col = '#00cccc', main = 'DMSO', xlim = c(-0.3,0.3),xlab = 'Step Magnitude')

R.magadj <- NULL
for (i in 1:length(R.frame$magnitude)){
  if (R.frame$step.type[i] == 'upstep'){
    R.magadj = c(R.magadj,abs(R.frame$magnitude[i]))
  }else if (R.frame$step.type[i] == 'downstep'){
    R.magadj = c(R.magadj,(0-abs(R.frame$magnitude[i])))
  }
  i = i+1
}
hist(R.magadj,breaks = 30, col = '#990066', main = 'Radicicol', xlim = c(-0.3,0.3),xlab = 'Step Magnitude')


low.magadj <- NULL
for (i in 1:length(low.frame$magnitude)){
  if (low.frame$step.type[i] == 'upstep'){
    low.magadj = c(low.magadj,abs(low.frame$magnitude[i]))
  }else if (low.frame$step.type[i] == 'downstep'){
    low.magadj = c(low.magadj,(0-abs(low.frame$magnitude[i])))
  }
  i = i+1
}
hist(low.magadj,breaks = 30, col = '#000066', main = 'T30', xlim = c(-0.3,0.3),xlab = 'Step Magnitude')

high.magadj <- NULL
for (i in 1:length(high.frame$magnitude)){
  if (high.frame$step.type[i] == 'upstep'){
    high.magadj = c(high.magadj,abs(high.frame$magnitude[i]))
  }else if (high.frame$step.type[i] == 'downstep'){
    high.magadj = c(high.magadj,(0-abs(high.frame$magnitude[i])))
  }
  i = i+1
}
hist(high.magadj,breaks = 30, col = '#ff9300', main = 'T37', xlim = c(-0.3,0.3),xlab = 'Step Magnitude')

stat.magadj <- NULL
for (i in 1:length(stat.frame$magnitude)){
  if (stat.frame$step.type[i] == 'upstep'){
    stat.magadj = c(stat.magadj,abs(stat.frame$magnitude[i]))
  }else if (stat.frame$step.type[i] == 'downstep'){
    stat.magadj = c(stat.magadj,(0-abs(stat.frame$magnitude[i])))
  }
  i = i+1
}
hist(stat.magadj,breaks = 30, col = '#9A0794', main = 'Stationary', xlim = c(-0.3,0.3),xlab = 'Step Magnitude')

par(mfrow=c(1,3),mar=c(4,3,4,1))
hist(glu.magadj,breaks = 30, col = '#e41a1c', main = 'Glucose', xlim = c(-0.3,0.3),xlab = 'Step Magnitude')
abline(v=0.05)
abline(v=-0.05)
hist(pho.magadj,breaks = 30, col = '#4daf4a', main = 'Phosphate', xlim = c(-0.3,0.3),xlab = 'Step Magnitude')
abline(v=0.05)
abline(v=-0.05)
hist(sul.magadj,breaks = 30, col = '#377eb8', main = 'Sulfate', xlim = c(-0.3,0.3),xlab = 'Step Magnitude')
abline(v=0.05)
abline(v=-0.05)



par(mfrow=c(2,3),mar=c(5,3,5,3))
hist(stat.magadj,breaks = 30, col = '#9A0794', main = 'Stationary', xlim = c(-0.3,0.3),xlab = 'Step Magnitude')
hist(low.magadj,breaks = 30, col = '#000066', main = 'T30', xlim = c(-0.3,0.3),xlab = 'Step Magnitude')
hist(high.magadj,breaks = 30, col = '#ff9300', main = 'T37', xlim = c(-0.3,0.3),xlab = 'Step Magnitude')
hist(D.magadj,breaks = 30, col = '#00cccc', main = 'DMSO', xlim = c(-0.3,0.3),xlab = 'Step Magnitude')
hist(R.magadj,breaks = 30, col = '#990066', main = 'Radicicol', xlim = c(-0.3,0.3),xlab = 'Step Magnitude')

all = data.frame(condition = c(rep('stationary',length(stat.magadj)),rep('T30',length(low.magadj))),magnitude = c(stat.magadj,low.magadj))
p <- ggplot(all,aes(x=condition,y=magnitude,fill=condition)) + geom_violin() + geom_boxplot(width=.1,position = position_dodge(.9)) + theme_light()
p
