#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## args
# args[1] = scripts directory
# args[2] = working directory
# args[3] = group1 (MN_D)
# args[4] = group2 (MN_R)
# args[5] = exp name (DMSOvRadicicol)

conditions = unlist(strsplit(args[5],'v',fixed = TRUE))

library(ggplot2)

# plot across whole genome
pdf(paste(args[2],'/',args[5],'_fitness_overview.pdf',sep=''), width=20, height=8.5)
par(mfrow=c(2,1), mar=c(3,5,3,2))

con1 = read.csv(paste(args[2],'/',args[3],'_slopes_se_norm.csv',sep=''))
con2 = read.csv(paste(args[2],'/',args[4],'_slopes_se_norm.csv',sep=''))

chromosomes = c(0,230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066)
centromeres = c(151465, 238207, 114385, 449711, 151987, 148510, 497038, 105703, 355629, 436425, 440246, 150947, 268031, 628758, 326702, 555957)

chr1x <- vector()
chr1y <- vector()
chr2x <- vector()
chr2y <- vector()

for(i in 1:length(con1$genes)){
  n = as.integer(con1$chr[i])
  coordinate = as.integer(con1$number[i])+as.integer(sum(chromosomes[1:n]))
  chr1x = append(chr1x, coordinate)
  chr1y = append(chr1y, con1[i,4])
}

for(i in 1:length(con2$genes)){
  n = as.integer(con2$chr[i])
  coordinate = as.integer(con2$number[i])+as.integer(sum(chromosomes[1:n]))
  chr2x = append(chr2x, coordinate)
  chr2y = append(chr2y, con2[i,4])
}

plot(chr1x, chr1y, type = 'h', cex=0.7, col= "#00cccc", pch=15, xaxt="n", main = conditions[1], xlab = 'Genome coordinate', ylab = 'Relative Fitness', ylim = c(-0.45, 0.45), xaxs = 'i')
abline(h=0)
abline(v=230218.5, lty = 2)
abline(v=1043402.5, lty = 2)
abline(v=1360022.5, lty = 2)
abline(v=2891955.5, lty = 2)
abline(v=3468829.5, lty = 2)
abline(v=3738990.5, lty = 2)
abline(v=4829930.5, lty = 2)
abline(v=5392573.5, lty = 2)
abline(v=5832461.5, lty = 2)
abline(v=6578212.5, lty = 2)
abline(v=7245028.5, lty = 2)
abline(v=8323205.5, lty = 2)
abline(v=9247636.5, lty = 2)
abline(v=10031969.5, lty = 2)
abline(v=11123260.5, lty = 2) 
points(x = c(151465, 468425, 1157787, 1809733, 3043942, 3617339, 4235910, 4935516, 5748202, 6268768, 7018341, 7395856, 8591236, 9876394, 10358553, 11679217), y = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), pch = 15, col = 'black')

plot(chr2x, chr2y, type = 'h', cex=0.7, col= "#990066", pch=15, xaxt="n", main = conditions[2], xlab = 'Genome coordinate', ylab = 'Relative Fitness', ylim = c(-0.45, 0.45), xaxs = 'i')
abline(h=0)
abline(v=230218.5, lty = 2)
abline(v=1043402.5, lty = 2)
abline(v=1360022.5, lty = 2)
abline(v=2891955.5, lty = 2)
abline(v=3468829.5, lty = 2)
abline(v=3738990.5, lty = 2)
abline(v=4829930.5, lty = 2)
abline(v=5392573.5, lty = 2)
abline(v=5832461.5, lty = 2)
abline(v=6578212.5, lty = 2)
abline(v=7245028.5, lty = 2)
abline(v=8323205.5, lty = 2)
abline(v=9247636.5, lty = 2)
abline(v=10031969.5, lty = 2)
abline(v=11123260.5, lty = 2) 
points(x = c(151465, 468425, 1157787, 1809733, 3043942, 3617339, 4235910, 4935516, 5748202, 6268768, 7018341, 7395856, 8591236, 9876394, 10358553, 11679217), y = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), pch = 15, col = 'black')

dev.off()

#plot by length
pdf(file = paste(args[2],'/',args[3],'_byLength.pdf',sep=''), width=6, height=4)

data = Reduce(function(x, y) merge(x, y, all=TRUE), list(con1, con2))
p = ggplot(data)
pp = p + geom_point(aes(x = Tamp.length, y = data[,9], color = factor(chr)), na.rm = TRUE)+scale_x_continuous(limits=c(0,1100000))+scale_y_continuous(limits=c(-0.5,0.5)) + ylab(conditions[1])
pp
dev.off()

pdf(file = paste(args[2],'/',args[4],'_byLength.pdf',sep=''), width=6, height=4)

p = ggplot(data)
pp = p + geom_point(aes(x = Tamp.length, y = data[,11], color = factor(chr)), na.rm = TRUE)+scale_x_continuous(limits=c(0,1100000))+scale_y_continuous(limits=c(-0.5,0.5)) + ylab(conditions[2])
pp
dev.off()

#####################################################
refFile <- read.csv('/Users/ankeller/Dropbox/190404_Witten_model/fit-plots-v4-DiffSlopesSameEnds-oct-2018-update/190411_fit_reference.csv')
chrArms <- c('BR','BL','CR','CL','DR','DL','ER','EL','FR','FL','GR','GL','HR','HL','IR','IL','JR','JL','KR','KL','LR','LL','MR','ML','NR','NL','OR','OL','OR','PL')
PC.rad.List <- c()
i <- 1
for (call in refFile$Radicicol) {
  if (call == 'PC') {
    PC.rad.List <- c(PC.rad.List,chrArms[i])
  }
  i <- i + 1
}

PC.rad <- data[which(substr(data$genes,2,3) %in% PC.rad.List),]
LM.rad <- data[which(!(substr(data$genes,2,3) %in% PC.rad.List)),]

p = ggplot(PC.rad) + aes(x=Tamp.length, y=PC.rad[,11],color = factor(chr))
pp = p + geom_point(na.rm = TRUE) + facet_wrap(~chr) + scale_x_continuous(limits=c(0,1100000))+scale_y_continuous(limits=c(-0.5,0.5)) + ylab(conditions[2])
pp 
slope.PC <- coef(lm(PC.rad$MN_R~PC.rad$Tamp.length))[2]


p = ggplot(LM.rad)
pp = p + geom_point(aes(x=Tamp.length, y=LM.rad[,11], color = factor(chr)), na.rm = TRUE)+scale_x_continuous(limits=c(0,1100000))+scale_y_continuous(limits=c(-0.5,0.5)) + ylab(conditions[2])
pp
slope.LM <- coef(lm(LM.rad$MN_R~LM.rad$Tamp.length))[2]



#####################################################

# compare conditions

allgenes <- data.frame(genes = read.csv(paste(args[1],'/gene_locations_names_chr',sep=''), sep="\t", header = FALSE)[,1])
con1 = read.csv(paste(args[2],'/',args[3],'_slopes_se_norm.csv',sep=''))
con2 = read.csv(paste(args[2],'/',args[4],'_slopes_se_norm.csv',sep=''))

data = Reduce(function(x, y) merge(x, y, all=TRUE), list(allgenes, con1, con2))
write.csv(data, paste(args[2],'/',args[5],'_data_summary.csv',sep=''))

pdf(file = paste(args[2],'/',args[5],'_compare.pdf',sep=''), width=6, height=4)
q = ggplot(data)
qq = q + geom_point(aes(x = data[,9], y = data[,11], color = factor(chr)), na.rm = TRUE)+scale_x_continuous(limits=c(-0.4, 0.4))+scale_y_continuous(limits=c(-0.4,0.4))
qq
dev.off()

