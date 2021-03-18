library(readxl)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(tidyr)
library(gridExtra)
library(dplyr)

set.seed(357)

# setwd
setwd("/Users/zacc/USyd/UMI_repeatexpansions")

# read in data
dat <- read_xlsx("/Users/zacc/USyd/UMI_repeatexpansions/ngs_metrics.xlsx",sheet = 1)

# calculate library size
mean(dat$size_bp_bioanalyser)
sd(dat$size_bp_bioanalyser)

t.test(dat$size_bp_bioanalyser[dat$assay == "sWGS"],
       dat$size_bp_bioanalyser[dat$assay == "repeSeq"])

# reads v pooling
mycolours <- c("repeSeq" = "dodgerblue", "sWGS" = "grey")
p1 <- ggplot(dat, aes(x=assay, y=reads_r1/1000000, fill=assay)) + 
  geom_point(position = position_jitterdodge(),alpha=0.5, aes(colour = assay)) + 
  geom_boxplot(alpha=0.6,outlier.alpha=0,) + ylab("Reads (M)") + theme_bw() + xlab("") + 
  scale_fill_manual(values=mycolours) +
  scale_color_manual(values = mycolours) +
  theme(legend.position="none") 

dat$reads.pool <- (dat$reads_r1/sum(dat$reads_r1))/(dat$prc_seq_pool/100)

p2 <- ggplot(dat, aes(x=assay, y=reads.pool, fill=assay)) + 
  geom_point(position = position_jitterdodge(),alpha=0.5, aes(colour = assay)) + 
  geom_boxplot(alpha=0.6,outlier.alpha=0,) + ylab("Observed / Expected") + theme_bw() + xlab("") + 
  scale_fill_manual(values=mycolours) +
  scale_color_manual(values = mycolours) +
  theme(legend.position="none") 

pdf("boxplot_reads.pdf")
grid.arrange(p1,p2,nrow=2, ncol=4)
dev.off()

# correlation between sWGS and repeSeq
mean(dat$reads_r1[dat$assay == "sWGS"])
mean(dat$reads_r1[dat$assay == "repeSeq"])

cor.test(dat$reads_r1[dat$assay == "sWGS"],
         dat$reads_r1[dat$assay == "repeSeq"])

# reads processed
dat$reads <- dat$reads_r1 - dat$reads_collapse_r1
dat$unique_r1 <- dat$reads_collapse_r1 - dat$reads_collapse_trim_r1

mydata <- pivot_longer(dat,c("reads","reads_collapse_trim_r1"))
mydata$name <- factor(mydata$name)
levels(mydata$name) <- c("reads","reads_collapse_trim_r1")

mycolours <- c("reads" = "black", "reads_collapse_trim_r1" = "grey")
p3<-ggplot(mydata, aes(fill=name, y=value/1000000, x=sample_name)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=mycolours) +
  ggtitle("") + ylab("reads (M)") +
  theme_minimal() + theme(legend.position="top") +
  xlab("")

pdf("barplot_reads_retained.pdf")
grid.arrange(p3,nrow=2, ncol=2)
dev.off()

(1-mean(dat$reads_collapse_trim_r1[dat$assay == "sWGS"]/dat$reads_r1[dat$assay == "sWGS"])) * 100
(1-mean(dat$reads_collapse_trim_r1[dat$assay == "repeSeq"]/dat$reads_r1[dat$assay == "repeSeq"])) * 100

t.test(c(dat$reads_collapse_trim_r1[dat$assay == "sWGS"]/dat$reads_r1[dat$assay == "sWGS"]),
       c(dat$reads_collapse_trim_r1[dat$assay == "repeSeq"]/dat$reads_r1[dat$assay == "repeSeq"]))

t.test(dat$dup_prc_r1_collapse_trim[dat$assay == "sWGS"],
       dat$dup_prc_r1_collapse_trim[dat$assay == "repeSeq"])

# sequencing metrics
mycolours <- c("repeSeq" = "dodgerblue", "sWGS" = "grey")

dat$mean_gc_ct <- rowMeans(dat[,c("gc_prc_r1_collapse_trim","gc_prc_r2_collapse_trim")])
p1 <- ggplot(dat, aes(x=assay, y=mean_gc_ct*100, fill=assay)) + 
  geom_point(position = position_jitterdodge(),alpha=0.5, aes(colour = assay)) + 
  geom_boxplot(alpha=0.6,outlier.alpha=0,) + ylab("GC %") + theme_bw() + xlab("") + 
  scale_fill_manual(values=mycolours) +
  scale_color_manual(values = mycolours) +
  theme(legend.position="none") + ylim(40,43.5)

mycolours <- c("repeSeq" = "dodgerblue", "sWGS" = "grey")

p2 <- ggplot(dat, aes(x=assay, y=length_r1_collapse_trim, fill=assay)) + 
  geom_point(position = position_jitterdodge(),alpha=0.5, aes(colour = assay)) + 
  geom_boxplot(alpha=0.6,outlier.alpha=0,) + ylab("read length (bp)") + theme_bw() + xlab("") + 
  scale_fill_manual(values=mycolours) +
  scale_color_manual(values = mycolours) +
  theme(legend.position="none") 

pdf("barplot_gc_&_readlength.pdf")
grid.arrange(p1,p2,nrow=2, ncol=4)
dev.off()

# reads with "CCGGGG"
mycolours <- c("repeSeq" = "dodgerblue", "sWGS" = "grey")

p1 <- ggplot(dat, aes(x=assay, y=(dat$CCGGGG_r1_collapse_trim/dat$reads_collapse_trim_r1) * 100, fill=assay)) + 
  geom_point(position = position_jitterdodge(),alpha=0.5, aes(colour = assay)) + 
  geom_boxplot(alpha=0.6,outlier.alpha=0,) + ylab("Reads w 'CCGGGG' (%)") + theme_bw() + xlab("") + 
  scale_fill_manual(values=mycolours) +
  scale_color_manual(values = mycolours) +
  theme(legend.position="none") 

p2 <- ggplot(dat, aes(x=assay, y=(dat$CCGGGGCCGGGG_r1_collapse_trim/dat$reads_collapse_trim_r1) * 100, fill=assay)) + 
  geom_point(position = position_jitterdodge(),alpha=0.5, aes(colour = assay)) + 
  geom_boxplot(alpha=0.6,outlier.alpha=0,) + ylab("Reads w 'CCGGGGCCGGGG' (%)") + theme_bw() + xlab("") + 
  scale_fill_manual(values=mycolours) +
  scale_color_manual(values = mycolours) +
  theme(legend.position="none") 

t.test(c((dat$CCGGGG_r1_collapse_trim/dat$reads_collapse_trim_r1) * 100)[dat$assay == "sWGS"],
       c((dat$CCGGGG_r1_collapse_trim/dat$reads_collapse_trim_r1) * 100)[dat$assay == "repeSeq"])

t.test(c((dat$CCGGGGCCGGGG_r1_collapse_trim/dat$reads_collapse_trim_r1) * 100)[dat$assay == "sWGS"],
       c((dat$CCGGGGCCGGGG_r1_collapse_trim/dat$reads_collapse_trim_r1) * 100)[dat$assay == "repeSeq"])

# reads mapping capture target
mycolours <- c("repeSeq" = "dodgerblue", "sWGS" = "grey")

p3 <- ggplot(dat, aes(x=assay, y=(dat$reads_art_pulldown/dat$reads_collapse_trim_r1) * 100, fill=assay)) + 
  geom_point(position = position_jitterdodge(),alpha=0.5, aes(colour = assay)) + 
  geom_boxplot(alpha=0.6,outlier.alpha=0,) + ylab("Reads mapping capture (%)") + theme_bw() + xlab("") + 
  scale_fill_manual(values=mycolours) +
  scale_color_manual(values = mycolours) +
  theme(legend.position="none") 

p4 <- ggplot(dat, aes(x=assay, y=(dat$reads_art_c9region/dat$reads_collapse_trim_r1) * 100, fill=assay)) + 
  geom_point(position = position_jitterdodge(),alpha=0.5, aes(colour = assay)) + 
  geom_boxplot(alpha=0.6,outlier.alpha=0,) + ylab("Reads mapping c9orf72 (%)") + theme_bw() + xlab("") + 
  scale_fill_manual(values=mycolours) +
  scale_color_manual(values = mycolours) +
  theme(legend.position="none") 

t.test(c((dat$reads_art_pulldown/dat$reads_collapse_trim_r1) * 100)[dat$assay == "sWGS"],
       c((dat$reads_art_pulldown/dat$reads_collapse_trim_r1) * 100)[dat$assay == "repeSeq"])

t.test(c((dat$reads_art_c9region/dat$reads_collapse_trim_r1) * 100)[dat$assay == "sWGS"],
       c((dat$reads_art_c9region/dat$reads_collapse_trim_r1) * 100)[dat$assay == "repeSeq"])


pdf("boxplot_capture_enrich.pdf")
grid.arrange(p1,p2,p3,p4,nrow=2, ncol=4)
dev.off()

# enrichment between repeSeq vs intra-sample
dat2 <- dat[dat$assay == "repeSeq" & dat$source == "human",]
dat1 <- dat[dat$assay == "sWGS" & dat$source == "human",]

dat2$repeSeq_intersect <- as.numeric(dat2$intersect_all_repeSeq) / sum(dat2$bases_covered)
dat2$sample_intersect <- as.numeric(dat2$intersect_repewgs) / (dat2$bases_covered + dat1$bases_covered)

mydata <- pivot_longer(dat2,c("repeSeq_intersect","sample_intersect"))
mydata$name <- factor(mydata$name)
levels(mydata$name) <- c("repeSeq_intersect","sample_intersect")
mycolours <- c("repeSeq_intersect" = "dodgerblue", "sample_intersect" = "grey")

p3<-ggplot(mydata, aes(fill=name, y=value * 100, x=sample_name)) + 
  geom_bar(position='dodge', stat="identity") +
  scale_fill_manual(values=mycolours) +
  ggtitle("") + ylab("intersected bases (%)") +
  theme_minimal() + theme(legend.position="top") +
  xlab("")

pdf("barplot_intersect.pdf")
grid.arrange(p3,nrow=2, ncol=4)
dev.off()
