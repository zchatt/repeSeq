library(LOLA)
library(simpleCache)
library(GenomicRanges)
require(data.table)

set.seed(357)

setwd("/Users/zacc/USyd/UMI_repeatexpansions/lola")

# loading regions from LOLA db
regionDB = loadRegionDB("/Users/zacc/USyd/UMI_repeatexpansions/lola/nm/t1/resources/regions/LOLACore/hg19")

# load in-house regions
regionSetA = readBed("inter_all_UMB1-2_merged.bed")
regionSetB = readBed("inter_all_UMB2-2_merged.bed")
regionSetC = readBed("inter_all_UMB3-2_merged.bed")
regionSetD = readBed("inter_all_UMB4-2_merged.bed")

userSets = GRangesList(regionSetA, regionSetB, regionSetC, regionSetD)

# set universe
universe_swgs = readBed("scwgs_merged.bed") # this universe is from 4 x swgs experiments merged
universe_repe = readBed("repe_merged.bed") # this universe is from 4 x repeseq experiments merged

# run enrichment
res = runLOLA(userSets, universe_swgs, regionDB, cores=1,redefineUserSets=TRUE)
resr = runLOLA(userSets, universe_repe, regionDB, cores=1,redefineUserSets=TRUE)

# analye results
res[res$userSet == 1,]
res[res$userSet == 2,]
res[res$userSet == 3,]
res[res$userSet == 4,]

resr[resr$userSet == 1 & resr$collection == "ucsc_features",]
resr[resr$userSet == 2 & resr$collection == "ucsc_features",]
resr[resr$userSet == 3 & resr$collection == "ucsc_features",]
resr[resr$userSet == 4 & resr$collection == "ucsc_features",]

# plotting
# ucsc features
dplot <- rbind(res[res$collection == "ucsc_features",],
               resr[resr$collection == "ucsc_features",])
dplot$universe <- c(rep("sWGS", nrow(res[res$collection == "ucsc_features",])),
                    rep("repeSeq",nrow(resr[resr$collection == "ucsc_features",])))
dplot$description <- gsub("UCSC ","",dplot$description)

dplot$userSet1 <- as.character(paste("repeSeq",dplot$userSet))
mycolours <- c("repeSeq" = "grey", "sWGS" = "black")
dplot$assay <- rep("repeSeq",nrow(dplot))

p3 <- ggplot(dplot, aes(x=reorder(description,-oddsRatio,na.rm = TRUE), y=oddsRatio, fill=factor(universe))) + 
  geom_point(position = position_jitterdodge(),alpha=0.5, aes(colour = assay)) + 
  geom_boxplot(alpha=0.6,outlier.alpha=0,) + ylab("Odds-Ratio") + theme_bw() + xlab("") + 
  scale_fill_manual(values=mycolours) +
  scale_color_manual(values = mycolours) +
  theme(legend.position="none")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("boxplot_ucsc_feature_enrich.pdf")
grid.arrange(p3,nrow=2, ncol=2)
dev.off()

# encode tfbs
dplot <- rbind(res[res$collection == "encode_tfbs",],
               resr[resr$collection == "encode_tfbs",])
dplot$universe <- c(rep("sWGS", nrow(res[res$collection == "encode_tfbs",])),
                    rep("repeSeq",nrow(resr[resr$collection == "encode_tfbs",])))
dplot$description <- gsub("ChIP ","",dplot$description)

dplot$userSet1 <- as.character(paste("repeSeq",dplot$userSet))
mycolours <- c("repeSeq" = "grey", "sWGS" = "black")
dplot$assay <- rep("repeSeq",nrow(dplot))

# select top 10 for each sample
dplot$u_a <- paste(dplot$userSet,dplot$antibody)
d <- data.table(dplot, key="u_a")
dplot <- d[, head(.SD,1), by=u_a]

dplot <- dplot[1:40,]

p4 <- ggplot(dplot, aes(x=reorder(antibody,pValueLog,na.rm = TRUE), y=pValueLog, fill=factor(universe))) + 
  geom_point(position = position_jitterdodge(),alpha=0.5, aes(colour = assay)) + 
  geom_boxplot(alpha=0.6,outlier.alpha=0,) + ylab("Odds-Ratio") + theme_bw() + xlab("") + 
  scale_fill_manual(values=mycolours) +
  scale_color_manual(values = mycolours) +
  theme(legend.position="none")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

