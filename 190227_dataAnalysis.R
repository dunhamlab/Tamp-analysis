#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## args
# args[1] = scripts directory
# args[2] = working directory
# args[3] = group (MN_R)

# read in files
datafile = read.csv(paste(args[2],'/',args[3],'_slopes_remove_NAs.csv',sep = ''))

# make histograms (make plots new files rather than to the workspace)
pdf(file = paste(args[2],'/',args[3],'_hist.pdf',sep=''), width=6, height=4)
hist(datafile$standard.error,breaks = 100)
data.avg = mean(datafile$standard.error)
data.sd = sd(datafile$standard.error)
abline(v=(data.avg + data.sd))
data.cutoff = data.avg + data.sd
print(data.cutoff)
dev.off()

# filter out genes with sterr higher than cutoff
data.filt = datafile[which(datafile$standard.error < data.cutoff),]
data.filtout = datafile[which(datafile$standard.error > data.cutoff),]
filt.out = nrow(datafile) - nrow(data.filt)
print(filt.out)

# save intermediate files
write.csv(data.filtout, file = paste(args[2],'/',args[3],'_slopes_filtOut.csv',sep=''))
write.csv(data.filt, file = paste(args[2],'/',args[3],'_slopes_se_cutoff.csv',sep=''))

# get new mean of remaining genes
slope.avg = mean(data.filt$slopes)

# adjust slopes so mean is 0
norm.data = data.filt
norm.data$slopes = norm.data$slopes - slope.avg
print(mean(norm.data$slopes))

write.csv(norm.data, file = paste(args[2],'/',args[3],'_slopes_mean_norm.csv',sep=''),row.names = FALSE,quote = FALSE)






