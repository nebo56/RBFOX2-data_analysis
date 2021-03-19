library("ggplot2")
library("smoother")
library("cowplot")

args<-commandArgs(TRUE)
smoothing_window <- 25

# replicate 1
motif.rep1.xl <- read.table("L1_S1_L001_R1_001.fastq-q255.rmDup-xl-sum-flanked300-motif_density.txt", header=FALSE, sep=",",stringsAsFactors=F,na.strings=c("NA", "-", "?"))
motif.rep1.xl <- t(motif.rep1.xl)
colnames(motif.rep1.xl)[1] <- "coverage"
motif.rep1.xl <- as.data.frame(motif.rep1.xl)
motif.rep1.xl$map <- c(-300:300)
motif.rep1.xl$smooth <- smth(motif.rep1.xl$coverage, window = smoothing_window, method = "gaussian")

# replicate 2
motif.rep2.xl <- read.table("L2_S2_L001_R1_001.fastq-q255.rmDup-xl-sum-flanked300-motif_density.txt", header=FALSE, sep=",",stringsAsFactors=F,na.strings=c("NA", "-", "?"))
motif.rep2.xl <- t(motif.rep2.xl)
colnames(motif.rep2.xl)[1] <- "coverage"
motif.rep2.xl <- as.data.frame(motif.rep2.xl)
motif.rep2.xl$map <- c(-300:300)
motif.rep2.xl$smooth <- smth(motif.rep2.xl$coverage, window = smoothing_window, method = "gaussian")

# replicate 3
motif.rep3.xl <- read.table("L3_S3_L001_R1_001.fastq-q255.rmDup-xl-sum-flanked300-motif_density.txt", header=FALSE, sep=",",stringsAsFactors=F,na.strings=c("NA", "-", "?"))
motif.rep3.xl <- t(motif.rep3.xl)
colnames(motif.rep3.xl)[1] <- "coverage"
motif.rep3.xl <- as.data.frame(motif.rep3.xl)
motif.rep3.xl$map <- c(-300:300)
motif.rep3.xl$smooth <- smth(motif.rep3.xl$coverage, window = smoothing_window, method = "gaussian")

# adjust sizes for the main figure
tans <- 1.0
g_size <- 0.4
adj <- 0.4

# save output as PDF
pdf(args[4])

# plot
ggRNAmap.smooth <- ggplot() + theme_cowplot() + theme(plot.title = element_text(size = 8, face = "plain")) + background_grid(major = "xy", minor = "none") +
  geom_line(aes(motif.rep1.xl$map, as.vector(motif.rep1.xl$smooth), colour="rep1"),size=g_size, alpha=tans) +
  geom_line(aes(motif.rep2.xl$map, as.vector(motif.rep2.xl$smooth), colour="rep2"),size=g_size, alpha=tans) + 
  geom_line(aes(motif.rep3.xl$map, as.vector(motif.rep3.xl$smooth), colour="rep3"),size=g_size, alpha=tans) + 
  ggtitle("RBFOX2 (AG)U)GCAUG motif enrichment around eiCLIP crosslinks") + 
  xlab("position relative to eiCLIP cDNA-start") + 
  ylab("normalised motif coverage") + 
  theme(text=element_text(size=8),axis.text=element_text(size=8), axis.title=element_text(size=8,face="plain")) + 
  scale_colour_manual(values=c("#C90000","#788285","#0952B0")) +
  scale_x_continuous(limits = c(-300, 300)) 
ggRNAmap.smooth


