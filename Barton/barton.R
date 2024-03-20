## Analysis of the Barton data

## KO gene: YOR290C

# Loading libraries, data, and other prep---------------------------------------

library(ALDEx2)
library(tidyverse)
library(ggrepel)

set.seed(2022)
##Loading in the RDA data
load(file.path("Barton", "data", "barton.d.agg.t.Rda"))
dim(d.agg.t)
##6328 taxa over 96 samples

##The only metadata is "SNF" versus "WT"
##Samples to discard. Barton did not include these samples because of data quality issues.
sampsToDiscard <- c("WT.21", "WT.22", "WT.25", "WT.28", "WT.34", "WT.36", "SNF2.6", "SNF2.13", "SNF2.25", "SNF2.35")

##Discarding those samples in the data set
d.agg.t = d.agg.t[,!(colnames(d.agg.t) %in% sampsToDiscard)]

## strip out zero count genes in all samples
geneCounts = d.agg.t[rowSums(d.agg.t)>0,]

#-------------------------------------------------------------------------------

# Calculating the estimate of \theta^\perp--------------------------------------

## Calculating the geometric means
normalized_GC <- apply(geneCounts+.5, 2, FUN = function(x){x/sum(x)})
gm_GC <- apply(normalized_GC, 2, FUN = function(x){mean(log2(x))})
graph_df <- data.frame("gm" = gm_GC, "group" = c(rep("SNF2", 44), rep("WT", 42)))

#LFC = WT - SNF2
##Implied by the CLR on the observed data
mean(-1*gm_GC[45:86]) - mean(-1*gm_GC[1:44])

##Implied by the CLR on the Dirichlet data
conds = c(rep("SNF", 44), rep("WT", 42))
clr.samps <- aldex.clr(geneCounts, conds, mc.samples = 1000, denom = "all")
dir.sams <- clr.samps@dirichletData

gm_imp <- rep(NA,1000)

for(i in 1:1000){
  tmp <- matrix(NA,nrow = nrow(geneCounts), ncol = ncol(geneCounts))
  for(k in 1:86){
    tmp[,k] <- dir.sams[[k]][,i]
  }
  tmp_gm <- apply(tmp, 2, FUN = function(x){mean(log2(x))})
  gm_imp[i] <- mean(-1*mean(tmp_gm[45:86])) - mean(-1*mean(tmp_gm[1:44]))
}

mean(gm_imp)

#-------------------------------------------------------------------------------


## Comparing the different total models-----------------------------------------
gamma <- c(0, .25, .5, .75, 1, 1.25, 1.5, 1.75, 2)
sig.mod <- list()
sig.mod.inc <- list()
sign.mod <- list()
sign.mod.inc <- list()

for(i in 1:length(gamma)){
  mu <- c(rep(0.9,44), rep(1, 42))
  gamma.inc <- ALDEx2::aldex.makeScaleMatrix(gamma[i], mu, conds, log = FALSE, mc.samples = 500)
  
  print("Running clr model...")
  mod <- aldex(geneCounts, conds, mc.samples = 500, gamma = gamma[i], effect = TRUE)
  
  print("Running increase model...")
  mod.inc <- aldex(geneCounts, conds, gamma = gamma.inc, effect = TRUE, mc.samples = 500)
  
  sig.mod[[i]] <- rownames(mod[mod$we.eBH <= 0.05,])
  sign.mod[[i]] <- mod[, "effect"]
  sig.mod.inc[[i]] <- rownames(mod.inc[mod.inc$we.eBH <= 0.05,])
  sign.mod.inc[[i]] <- mod.inc[, "effect"]
  print(i)
}

##Making the graph
allTaxa <- rownames(geneCounts)
tax.nonSig <- c()
tax.inBoth <- c()
tax.inCLR <- c()
tax.inINC <- c()
sign.same <- c()
sign.diff <- c()

for(i in 1:length(gamma)){
  inEither <- unique(c(sig.mod[[i]], sig.mod.inc[[i]]))
  inBoth <- intersect(sig.mod[[i]], sig.mod.inc[[i]])
  tmp.nonSig <- sum(!(allTaxa %in% inEither))
  tmp.inBoth <- length(inBoth)
  tmp.inCLR <- length(setdiff(sig.mod[[i]], sig.mod.inc[[i]]))
  tmp.inINC <- length(setdiff(sig.mod.inc[[i]], sig.mod[[i]]))
  
  ##Now calculating the sign of those in both
  inds <- which(allTaxa %in% inBoth)
  sign.same.tmp <- sum(sign(sign.mod[[i]][inds]) == sign(sign.mod.inc[[i]][inds]))
  sign.diff.tmp <- sum(sign(sign.mod[[i]][inds]) != sign(sign.mod.inc[[i]][inds]))
  
  ##Concatenating
  tax.nonSig <- c(tax.nonSig, tmp.nonSig)
  tax.inBoth <- c(tax.inBoth, tmp.inBoth)
  tax.inCLR <- c(tax.inCLR, tmp.inCLR)
  tax.inINC <- c(tax.inINC, tmp.inINC)
  sign.same <- c(sign.same, sign.same.tmp)
  sign.diff <- c(sign.diff, sign.diff.tmp)
}


graph.df <- data.frame("Gamma" = rep(gamma, times = 4), "Count" = c(tax.nonSig, tax.inBoth, tax.inCLR, tax.inINC), "Group" = rep(c("Neither", "Both", "CLR Only", "Informed Only"), each = length(gamma)))

graph.df <- graph.df %>%
  mutate(gamma = Gamma) %>%
  mutate(Gamma = factor(Gamma)) %>%
  mutate(Group = factor(Group)) %>%
  mutate(Percent = (Count/nrow(geneCounts)) * 100) %>%
  filter(Group != "Neither")

graph.df$Group <- ordered(graph.df$Group, levels = c( "Informed Only", "CLR Only", "Both"))


ggplot(graph.df, aes(x= Gamma, y = Percent, group = Group, fill = Group)) +
  geom_bar(position="stack", stat="identity", width=.975) +
  scale_fill_manual(values = c("royalblue3", "red3", "orange")) +
  xlab(expression(gamma)) +
  ylab("Percent in Category") +
  theme_classic() +
  labs(fill='Significant In:') +
  geom_hline(yintercept=c(0,12.5,25,37.5,50,62.5,75,87.5,100), col = "lightgrey",linewidth=.5) +
  geom_segment(aes(x = 0.5, xend  = 0.5, y = 0, yend = 100), color = "lightgrey", lwd = 0.5)+
  geom_segment(aes(x = 1.5, xend  = 1.5, y = 0, yend = 100), color = "lightgrey", lwd = 0.5)+
  geom_segment(aes(x = 2.5, xend  = 2.5, y = 0, yend = 100), color = "lightgrey", lwd = 0.5)+
  geom_segment(aes(x = 3.5, xend  = 3.5, y = 0, yend = 100), color = "lightgrey", lwd = 0.5)+
  geom_segment(aes(x = 4.5, xend  = 4.5, y = 0, yend = 100), color = "lightgrey", lwd = 0.5)+
  geom_segment(aes(x = 5.5, xend  = 5.5, y = 0, yend = 100), color = "lightgrey", lwd = 0.5)+
  geom_segment(aes(x = 6.5, xend  = 6.5, y = 0, yend = 100), color = "lightgrey", lwd = 0.5)+
  geom_segment(aes(x = 7.5, xend  = 7.5, y = 0, yend = 100), color = "lightgrey", lwd = 0.5)+
  geom_segment(aes(x = 8.5, xend  = 8.5, y = 0, yend = 100), color = "lightgrey", lwd = 0.5)+
  scale_y_continuous(sec.axis = sec_axis(~.*63.27, name = "Number in Category")) + 
  theme(text=element_text(size=21)) 


ggsave(file.path("Barton", "results", "Barton_expandedGammaDiagramnoLegend.pdf"), width = 8.5, height = 5)

graph.df <- graph.df %>%
  filter(gamma >=1) %>%
  group_by(gamma) %>%
  ungroup()


##Creating inlay for smaller values
ggplot(graph.df, aes(x= Gamma, y = Percent, group = Group, fill = Group)) +
  geom_bar(position="stack", stat="identity", width=.975) +
  scale_fill_manual(values = c("royalblue3", "red3", "orange")) +
  xlab(expression(gamma)) +
  ylab("Percent in Category") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  labs(fill='Significant In:') +
  coord_cartesian(ylim=c(0, 5)) +
  geom_hline(yintercept=c(0,1,2,3,4,5), col = "lightgrey",linewidth=.5) +
  geom_segment(aes(x = 0.5, xend  = 0.5, y = 0, yend = 100), color = "lightgrey", lwd = 0.5)+
  geom_segment(aes(x = 1.5, xend  = 1.5, y = 0, yend = 100), color = "lightgrey", lwd = 0.5)+
  geom_segment(aes(x = 2.5, xend  = 2.5, y = 0, yend = 100), color = "lightgrey", lwd = 0.5)+
  geom_segment(aes(x = 3.5, xend  = 3.5, y = 0, yend = 100), color = "lightgrey", lwd = 0.5)+
  geom_segment(aes(x = 4.5, xend  = 4.5, y = 0, yend = 100), color = "lightgrey", lwd = 0.5)+
  geom_segment(aes(x = 5.5, xend  = 5.5, y = 0, yend = 100), color = "lightgrey", lwd = 0.5)+
  scale_y_continuous(sec.axis = sec_axis(~.*63.27, name = "Number in Category")) +
  theme(text=element_text(size=21)) 


ggsave(file.path("Barton", "results", "Barton_expandedGammaDiagramZoomnoLegend.pdf"), width = 7, height = 5)

#-------------------------------------------------------------------------------
