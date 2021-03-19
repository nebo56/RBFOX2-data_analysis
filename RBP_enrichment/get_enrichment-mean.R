library("data.table")
args<-commandArgs(TRUE)

# Script will calculate crosslink enrichment of target vs control region from input BED format

# eCLIP - control crosslinks
######################################
eCLIP.1.control <- fread(args[1], sep='\t')
eCLIP.1.control.sum <- aggregate(eCLIP.1.control$V8, list(eCLIP.1.control$V7), FUN=sum) #because we are only interested in mean we summed both strand together even that the order of enrichment is different
eCLIP.1.control.mean <- mean(eCLIP.1.control.sum$x)
remove(eCLIP.1.control)

# eCLIP - peak crosslinks
######################################
eCLIP.1 <- fread(args[2], sep='\t')
eCLIP.1.sum <- aggregate(eCLIP.1$V8, list(eCLIP.1$V7), FUN=sum) #because we are only interested in mean we summed both strand together even that the order of enrichment is different
eCLIP.1.mean <- mean(eCLIP.1.sum$x)
remove(eCLIP.1)

# normalised peak
norm.peak.enr <- t(data.frame(c(args[3],(eCLIP.1.mean / eCLIP.1.control.mean))))  #normalise average of crosslink inside the peak with the average control from the control region
write.table(norm.peak.enr, args[4], sep='\t', quote = FALSE, append = TRUE, col.names = FALSE, row.names = FALSE)


