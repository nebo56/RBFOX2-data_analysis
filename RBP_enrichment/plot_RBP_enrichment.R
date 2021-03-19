library("ggplot2")
args <- commandArgs(TRUE)

# Plot RBP enrichments as histograms

### SmB peaks ###
peak2 <- read.table(args[1], sep="\t", header=TRUE)

pdf(args[2], width = 12)

### peak2 ###
ggplot(peak2, aes(x = reorder(RBP, -norm.enrichment), y = norm.enrichment)) + theme_bw() + # sorted by enrichment
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("eCLIP: RBP enrichment around splice sites") +
  xlab("RBP") + 
  ylab("crosslink enrichment") +
  geom_bar(stat = "identity")


