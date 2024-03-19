## Code to re-create the simulations

## Specifying some helper functions

## Functions for sampling the data ---------------------------------------------
create_true_abundances <- function(d, n){
  dd <- length(d)/2
  dat <- d %>%
    sapply(function(x) rpois(n, lambda=x)) %>% 
    t() %>%
    as.data.frame() %>%
    split(rep(1:2, each=dd)) %>%
    purrr::map(~`rownames<-`(.x, paste0("Taxa", 1:dd))) %>%
    purrr::map(t) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    cbind(Condition=factor(rep(c(1, 2), each=n)), .) %>%
    #cbind(Condition=factor(rep(c("Pre", "Post"), each=n), levels = c("Pre", "Post")), .) %>%
    `rownames<-`(., NULL)
  return(dat)
}

resample_data <- function(dat, seq.depth){
  ddat <- driver::miniclo(as.matrix(dat[,-1]))
  for (i in 1:nrow(dat)){
    dat[i,-1] <- rmultinom(1, size=seq.depth, prob=ddat[i,])
  }
  return(dat)
}
#-------------------------------------------------------------------------------

# Functions for assessing significance------------------------------------------
sig_code <- function(sig, Taxa, truth){
  out <- rep("TN", length(Taxa))
  out[sig &(Taxa %in% truth)] <- "TP" # True Positives
  out[sig & (out!="TP")] <- "FP" # False Positives
  out[!sig & (Taxa %in% truth)] <- "FN" # False Negatives
  return(out)
}

append_sig <- function(fit, f){
  s <- f(fit)$category
  mutate(fit, sig=ifelse(category%in%s, TRUE, FALSE))
}
#-------------------------------------------------------------------------------

# Functions for plotting--------------------------------------------------------
strip_end_num <- function(x){
  x <- as.character(x)
  as.numeric(str_extract(x, "[[:digit:]]+"))
}


plot_count <- function(dat){
  gather(dat, Taxa, Count, -Condition) %>% 
    mutate(Taxa=strip_end_num(Taxa)) %>% 
    mutate(Taxa=factor(Taxa)) %>% 
    mutate(Condition = ifelse(Condition == 1, "Pre", "Post")) %>%
    ggplot(aes(x=Taxa, y=Count)) +
    geom_boxplot(aes(fill = Condition, color = Condition), position=position_dodge(width=1), 
                 linewidth=1)+
    scale_y_log10() +
    theme(text = element_text(size=16)) +
    theme_bw() +
    scale_color_npg() +
    scale_fill_npg()+
    #scale_fill_manual(values = c("#fdae61", "#2b83ba")) + 
    #scale_color_manual(values = c("#fdae61", "#2b83ba")) +
    labs(color='Antibiotic\nTreatment') +
    labs(fill='Antibiotic\nTreatment') 
}

plot_sig2 <- function(rrs, truth, ...){
  rrs$dat <- NULL
  names(rrs) <- model.names[names(rrs)]
  rrs.aldex <- rrs[5:10]
  rrs.others <- rrs[1:4]
  p.aldex <- bind_rows(rrs.aldex, .id="Model") %>% 
    dplyr::select(Model, category, sig) %>% 
    mutate(Taxa = category) %>% 
    mutate(Taxa=strip_end_num(Taxa)) %>% 
    mutate(sigcode = sig_code(sig, Taxa, truth)) %>% 
    mutate(Taxa=factor(Taxa), sigcode=factor(sigcode, 
                                             levels=c("TP", "TN", 
                                                      "FP", "FN"))) %>% 
    mutate(Model=factor(Model, levels=model.name.levels)) %>% 
    ggplot(aes(x=Taxa, y=Model)) +
    geom_tile_pattern(aes(fill=sigcode, pattern = sigcode), color="darkgrey",pattern_fill = 'grey',pattern_colour  = 'grey', pattern_density = 0.015) +
    theme_minimal() +
    theme(panel.grid = element_blank(), 
          legend.title=element_blank(),
          text = element_text(size=16),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("ALDEx2")+
    scale_pattern_manual(values = c(TP = "none", TN = "none", FP = "none", FN = "stripe")) +
    scale_fill_manual(values= c("black", "white", "grey", "white"))
  
  p.others <- bind_rows(rrs.others, .id="Model") %>% 
    dplyr::select(Model, category, sig) %>% 
    mutate(Taxa = category) %>% 
    mutate(Taxa=strip_end_num(Taxa)) %>% 
    mutate(sigcode = sig_code(sig, Taxa, truth)) %>% 
    mutate(Taxa=factor(Taxa), sigcode=factor(sigcode, 
                                             levels=c("TP", "TN", 
                                                      "FP", "FN"))) %>% 
    mutate(Model=factor(Model, levels=model.name.levels)) %>% 
    ggplot(aes(x=Taxa, y=Model)) +
    geom_tile_pattern(aes(fill=sigcode, pattern = sigcode), color="darkgrey",pattern_fill = 'grey',pattern_colour  = 'grey', pattern_density = 0.015) +
    theme_minimal() +
    theme(panel.grid = element_blank(), 
          legend.title=element_blank(),
          text = element_text(size=16),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("Others") +
    scale_pattern_manual(values = c(TP = "none", TN = "none", FP = "none", FN = "stripe")) +
    scale_fill_manual(values= c("black", "white", "grey", "white"))
  
  p.comb <- list(p.aldex, p.others)
  return(p.comb)
}

#-------------------------------------------------------------------------------

# Functions for running the tested methods--------------------------------------
run_aldex2 <- function(dat, denom="all", gamma= NULL, mc.samples = 500){
  countdata <- t(dat[,-1,drop=F])
  colnames(countdata) <- paste0("n", 1:ncol(countdata))
  aldex.fit <- aldex(countdata, as.character(dat$Condition), denom=denom, mc.samples = mc.samples, gamma = gamma) %>% 
    as.data.frame() %>% 
    rownames_to_column("category") %>% 
    dplyr::select(category, effect, we.eBH) %>% 
    mutate(padj=we.eBH) %>% 
    mutate(mean=effect) %>% 
    mutate(low=NA, high=NA) %>%
    mutate(sig = ifelse(padj < 0.05, TRUE, FALSE))
  return(aldex.fit)
}

run_deseq2 <- function(dat){
  coldata <- dat[,"Condition",drop=F]
  countdata <- t(dat[,-1,drop=F])
  rownames(coldata) <- colnames(countdata) <- paste0("n", 1:ncol(countdata))
  dds <- DESeqDataSetFromMatrix(countData=countdata,
                                colData=coldata,
                                design = ~Condition)
  dds <- DESeq(dds)
  res <- results(dds)
  fit <- res %>% 
    as.data.frame() %>% 
    rownames_to_column("category") %>% 
    dplyr::select(category, log2FoldChange, padj, lfcSE) %>% 
    mutate(low = log2FoldChange - 1.96*lfcSE, 
           high = log2FoldChange + 1.96*lfcSE) %>% 
    mutate(mean=log2FoldChange) %>%
    mutate(sig = ifelse(padj < 0.05, TRUE, FALSE))
  return(fit)
}

run_limma <- function(dat){
  coldata <- dat[,"Condition",drop=F]
  countdata <- t(dat[,-1,drop=F])
  rownames(coldata) <- colnames(countdata) <- paste0("n", 1:ncol(countdata))
  rownames(countdata) <- paste0("Taxa", 1:nrow(countdata))
  design <- model.matrix(~coldata$Condition)
  
  y <- DGEList(counts=countdata)
  y <- calcNormFactors(y)
  design <- model.matrix(~0+Condition, data = coldata)
  v <- voom(y,design,plot=FALSE)
  fit <- lmFit(v,design)
  contr <- makeContrasts(Condition1 - Condition2, levels = colnames(coef(fit)))
  contr.fit <- contrasts.fit(fit, contr)
  contr.fit <- eBayes(contr.fit)
  et <- topTable(contr.fit,number =20, p.value = 1) %>%
    rownames_to_column("category") %>%
    mutate(sig = ifelse(adj.P.Val < 0.05, TRUE, FALSE))
  
  return(et)
}

run_edger <- function(dat){
  coldata <- dat[,"Condition",drop=F]
  countdata <- t(dat[,-1,drop=F])
  rownames(coldata) <- colnames(countdata) <- paste0("n", 1:ncol(countdata))
  rownames(countdata) <- paste0("Taxa", 1:nrow(countdata))
  design <- model.matrix(~coldata$Condition)
  y <- DGEList(counts=countdata,group=c(coldata$Condition))
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  et <- exactTest(y)
  et <- topTags(et, n=20, p.value = 1) %>%
    as.data.frame() %>%
    rownames_to_column("category") %>%
    mutate(sig = ifelse(FDR < 0.05, TRUE, FALSE))
  return(et)
}

run_baySeq <- function(dat){
  coldata <- dat[,"Condition",drop=F]
  countdata <- t(dat[,-1,drop=F])
  rownames(coldata) <- colnames(countdata) <- paste0("n", 1:ncol(countdata))
  rownames(countdata) <- paste0("Taxa", 1:nrow(countdata))
  groups <- list(NDE=rep(1, nrow(coldata)), 
                 DE=as.numeric(coldata$Condition))
  dds <- new("countData", data = countdata, groups = groups, replicates = as.character(coldata$Condition))
  libsizes(dds) <- getLibsizes(dds)
  dds <- getPriors.NB(dds, cl = NULL)
  dds <- getLikelihoods(dds,cl=NULL,bootStraps=3,verbose=FALSE)
  res <- topCounts(dds,group="DE", number = 20)
  
  ##Merge to the data to see which taxa because the names are uninformative
  res <- res[,c("n1","n2","n3","n4","n5","likes","DE","FDR.DE","FWER.DE")]
  taxa_names <- countdata[,c("n1","n2","n3","n4","n5")] %>%
    as.data.frame() %>%
    rownames_to_column("category") %>%
    join(res, by = c("n1","n2","n3","n4","n5"), type = "right") %>%
    mutate(sig = ifelse(FWER.DE < 0.05, TRUE, FALSE))

  return(taxa_names)
  
}

#-------------------------------------------------------------------------------

# Analysis code-----------------------------------------------------------------

# The directory change assumes you are using the R project.

setwd("simulation")
###Loading libraries
library(plyr)
library(dplyr)
library(phyloseq)
library(driver)
library(tidyverse)
library(ALDEx2)
library(DESeq2)
library(limma)
library(edgeR)
library(baySeq)
library(gghighlight)
library(cowplot)
library(ggplot2)
library(magrittr)
library(stringi)
library(directlabels)
library(ggpattern)
library(latex2exp)
library(ggrepel)
library(ggsci)
set.seed(12345)

###Scaled ALDEx2 simulation

# Code to benchmark performance for one sample size (n = 50)--------------------
runBenchmark <- function(d, n, seq.depth, pval = 0.05, mc.samples = 500){
  dd <- length(d)/2
  truth1 <- !(d[1:dd]==d[(dd+1):(2*dd)])##testing if the mean is different
  truth2 <- (1:dd)[truth1]##Locations of the differences
  dat <- create_true_abundances(d, n=n)
  log.totals=log2(rowSums(dat[,-1]))
  print("Theta^perp for the true...")
  print(mean(log.totals[(n+1):(2*n)]) - mean(log.totals[1:n]))
  
  rdat <- resample_data(dat, seq.depth=seq.depth)
  
  norm_rdat <- apply(rdat[,-1]+.5, 1, FUN = function(x){x/sum(x)})
  gm_rdat <- apply(norm_rdat, 2, FUN = function(x){mean(log2(x))})
  
  print("Scale implied by the CLR normalization...")
  clr.samps <- aldex.clr(t(rdat[,-1]), as.character(rdat$Condition), mc.samples = 500, denom = "all")
  dir.sams <- clr.samps@dirichletData
  
  gm_imp <- rep(NA,500)
  for(i in 1:500){
    tmp <- matrix(NA,nrow = 20, ncol = 100)
    for(k in 1:100){
      tmp[,k] <- dir.sams[[k]][,i]
    }
    tmp_gm <- apply(tmp, 2, FUN = function(x){mean(log2(x))})
    gm_imp[i] <- mean(-1*mean(tmp_gm[(n+1):(2*n)])) - mean(-1*mean(tmp_gm[1:n]))
  }
  print("Mean of the CLR norm across Dirichlet samples...")
  print(mean(gm_imp))
  print("Mean of the CLR norm on observed data...")
  print(mean(-1*gm_rdat[(n+1):(2*n)]) - mean(-1*gm_rdat[1:n]))
  
  # Running the other models
  afit <- run_aldex2(rdat, mc.samples = mc.samples)

  dfit <- run_deseq2(rdat)
  
  lfit <- run_limma(rdat)
  
  bsfit <- run_baySeq(rdat)
  
  efit <- run_edger(rdat)
  
  # Running new scale samples
  tfit.mini <- run_aldex2(rdat, mc.samples = mc.samples, gamma = 1e-3)
  tfit.small <- run_aldex2(rdat, mc.samples = mc.samples, gamma = .0725)
  tfit.mid <- run_aldex2(rdat, mc.samples = mc.samples, gamma =  .25)
  tfit.large <- run_aldex2(rdat, mc.samples = mc.samples, gamma = .5)
  scale.samps <- aldex.makeScaleMatrix(gamma = .25, mu = c(rep(1,n),rep(.9,n)), conditions = as.character(rdat[,1]), log = FALSE, mc.samples = mc.samples)
  tfit.informed <- run_aldex2(rdat, mc.samples = mc.samples, gamma = scale.samps)
  rrs <- list(dat=rdat, bsfit = bsfit, dfit = dfit, efit = efit, lfit = lfit, afit=afit, tfit.mini = tfit.mini, tfit.small = tfit.small, tfit.mid = tfit.mid, tfit.large = tfit.large, tfit.informed = tfit.informed)

  # Plotting
  p1 <- plot_count(dat)
  p2 <- plot_sig2(rrs, truth=truth2)
  p1 <- p1+theme(axis.title.x = element_blank(), 
                 axis.text.x = element_blank(), 
                 axis.ticks.x=element_blank(),
                 text = element_text(size=18))
  p <- plot_grid(p1, p2[[1]], p2[[2]], nrow=3, align="v", rel_heights=c(1.7, 1,1))
  p
  #ggsave(p1, file.path("Results","aldex_deseq_failures.pdf"))
}# end of function

###Setting the data parameters for all simulations
d <- c(4000, 4000, 4000, 4000, 4000, 400,400,400,400,4000,400,500,500,500,400,400,400,400,400,400, # Pre
       4000, 4000, 3000, 2000, 4000, 400,400,400,400,4000,400,500,500,500,200,400,400,400,400,100) # Post

# Specifying models and model names
model.names <- c("afit"="Original",
                 "dfit" = "DESeq2",
                 "bsfit" = "BaySeq",
                 "lfit" = "limma",
                 "efit" = "edgeR",
                 "tfit.mini" = "Gamma = 0",
                 "tfit.small" = "Gamma = 0.07",
                 "tfit.mid"= "Gamma = 0.25",
                 "tfit.large" = "Gamma = 0.5",
                 "tfit.informed" = "Informed")
model.name.levels <- c("Informed", "Gamma = 0.5", "Gamma = 0.25","Gamma = 0.07", "Gamma = 0", "Original","limma","edgeR", "DESeq2", "BaySeq")

##If this throws an 'Error in seq.default(from, to by)' error, increase size of your plot viewer.
runBenchmark(d, n = 50, seq.depth = 5000, mc.samples = 2000)
ggsave(file.path("results", "sim-res-by-method.pdf"), height = 7, units = "in", width = 7)

#-------------------------------------------------------------------------------

#Code to create plot for FDR by method over sample size-------------------------

##FDR by method over sample sizes
##To create the line graph

## This function will be used to run over a single sample size
runBenchmark_DF <- function(d, n, seq.depth, pval = 0.05, mc.samples = 500){
  dd <- length(d)/2
  truth1 <- !(d[1:dd]==d[(dd+1):(2*dd)])##testing if the mean is different
  truth2 <- (1:dd)[truth1]##Locations of the differences
  dat <- create_true_abundances(d, n=n)
  log.totals=log2(rowSums(dat[,-1]))
  
  taxa_truth <- paste0("Taxa", truth2)
  rdat <- resample_data(dat, seq.depth=seq.depth)
  sample.totals=rowSums(dat[,-1])
  
  #Storing results in order
  tp <- rep(NA, 7)
  fp <- rep(NA, 7)
  afit <- run_aldex2(rdat, mc.samples = mc.samples) %>%
    filter(padj < 0.05)
  tp[1] <- sum(afit$category %in% taxa_truth)
  fp[1] <- sum(!(afit$category %in% taxa_truth))
  
  dfit <- run_deseq2(rdat) %>%
    filter(padj < 0.05)
  tp[2] <- sum(dfit$category %in% taxa_truth)
  fp[2] <- sum(!(dfit$category %in% taxa_truth))
  
  lfit <- run_limma(rdat) %>%
    filter(adj.P.Val < 0.05)
  tp[3] <- sum(lfit$category %in% taxa_truth)
  fp[3] <- sum(!(lfit$category %in% taxa_truth))
  
  bsfit <- run_baySeq(rdat) %>%
    filter(FWER.DE < 0.05)
  tp[4] <- sum(bsfit$category %in% taxa_truth)
  fp[4] <- sum(!(bsfit$category %in% taxa_truth))
  
  efit <- run_edger(rdat) %>%
    filter(FDR < 0.05)
  tp[5] <- sum(efit$category %in% taxa_truth)
  fp[5] <- sum(!(efit$category %in% taxa_truth))

  tfit.mid <- run_aldex2(rdat, mc.samples = mc.samples, gamma =  .5)%>%
    filter(padj < 0.05)
  tp[6] <- sum(tfit.mid$category %in% taxa_truth)
  fp[6] <- sum(!(tfit.mid$category %in% taxa_truth))
  
  
  scale.samps <- aldex.makeScaleMatrix(gamma = .25, mu = c(rep(1,n),rep(.9,n)), conditions = as.character(rdat[,1]), log = FALSE, mc.samples = mc.samples)
  tfit.informed <- run_aldex2(rdat, mc.samples = mc.samples, gamma = scale.samps)%>%
    filter(padj < 0.05)
  tp[7] <- sum(tfit.informed$category %in% taxa_truth)
  fp[7] <- sum(!(tfit.informed$category %in% taxa_truth))
  
  
  
  res <- data.frame("tp" = tp, "fp" = fp, "method" = c("ALDEx2", "DESeq2", "limma", "baySeq", "edgeR", "ALDEx2 (gamma = 0.5)", "Informed"))
  return(res)
}# end of function

##Running over sample sizes
n_to_test <- c(5,10,25,50,75,100,125,150,175,200,225,250,275,300)
# Empty data frame to store results
benchmark_df <- data.frame()
k <- 3 #This is the number of replicates
for(j in 1:k){
  for(i in 1:length(n_to_test)){
    tmp <- runBenchmark_DF(d, n = n_to_test[i], seq.depth = 5000, mc.samples = 2000)
    tmp$n <- rep(n_to_test[i], nrow(tmp))
    tmp$k <- rep(j, nrow(tmp))
    benchmark_df <- rbind(benchmark_df, tmp)
    print(i)
    
  }
}

## Plotting
benchmark_df$fdr <- benchmark_df$fp/(benchmark_df$tp + benchmark_df$fp)
benchmark_df$fdr <- ifelse(is.na(benchmark_df$fdr), 0 , benchmark_df$fdr)
benchmark_df$label <- ifelse(benchmark_df$n == 300 & benchmark_df$k == 5, benchmark_df$method, NA)
benchmark_df$method <- ifelse(benchmark_df$method == "Informed", "ALDEx2 (Informed)", benchmark_df$method)
ggplot(benchmark_df, aes(x=n, y = fdr, group = method, color = method)) +
  #geom_line(aes(linetype = method), lwd = 0.5) +
  geom_point(aes(shape = method), alpha = 0.9) +
  geom_smooth(method = "loess", alpha = 0.05, linewidth = 1, span = 1) +
  scale_color_npg() +
  scale_shape_manual(values = c(1, 15, 3, 4, 5, 6, 7)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        text = element_text(size=16)) +
  xlab("Sample Size") +
  ylim(c(0,1)) +
  ylab("False Discovery Rate")
  #geom_label_repel(aes(label = label),
  #                 na.rm = TRUE) 

ggsave(file.path("results", "sim-fdr-by-method-by-size-shape.pdf"), height = 7, units = "in", width = 10)


#-------------------------------------------------------------------------------

#Plot to create gamma plot (supplement)-----------------------------------------

## Gamma plot

# Sampling one data set
dd <- length(d)/2
truth1 <- !(d[1:dd]==d[(dd+1):(2*dd)])##testing if the mean is different
truth2 <- (1:dd)[truth1]##Locations of the differences
dat <- create_true_abundances(d, n=50)
rdat <- resample_data(dat, seq.depth=5000)

# Running aldex2 on different values of gamma
countdata <- t(rdat[,-1,drop=F])
colnames(countdata) <- paste0("n", 1:ncol(countdata))
clr <- aldex.clr(countdata, as.character(rdat$Condition), mc.samples = 1000, gamma = 1e-3)
sen_res <- aldex.senAnalysis(clr, gamma = c(1e-3,.01, 0.025, 0.05, 0.075, .1,.2,.3,.4,.5,.6,.7,.8,.9,1))

# Transforming from list to data frame

gamma <- as.numeric(sub("gamma_", "", names(sen_res)))
B <- matrix(NA, nrow = length(sen_res), ncol = dim(sen_res[[1]])[1])
pvals <- matrix(NA, nrow = length(sen_res), ncol = dim(sen_res[[1]])[1])

for(i in 1:length(sen_res)){
  pvals[i,] <- sen_res[[i]]$we.ep
  B[i,] <- sen_res[[i]]$effect
}

# Prepping for plotting and plottind data

P <- as.data.frame(pvals) 
P$gamma <- gamma
P <- P[,c(ncol(P), 1:(ncol(P)-1))]
P <- reshape(P, idvar = "gamma",
             varying = list(2:ncol(P)),
             timevar = "Sequence",
             v.names = "pval", direction = "long")
P$Sequence = paste0("V", P$Sequence)

P.toLabel <- P[(P$pval < 0.1),]
P.toLabel <- P.toLabel[order(-P.toLabel$gamma),]
seq_to_label <- as.numeric(sub("V", "", unique(P.toLabel$Sequence)))

taxa_to_label = as.vector(na.omit(seq_to_label))

B.graph <- as.data.frame(B)
B.graph$gamma <- gamma
B.graph <- B.graph[,c(ncol(B.graph), 1:(ncol(B.graph)-1))]
B.graph <- reshape(B.graph, idvar = "gamma",
                   varying = list(2:ncol(B.graph)),
                   timevar = "Sequence",
                   v.names = "Effect", direction = "long")
B.graph$Sequence <- paste0("V", B.graph$Sequence)
B.graph <- base::merge(B.graph, P, by = c("gamma", "Sequence"))
B.graph$Sequence <- sub("V", "", B.graph$Sequence)
B.graph$labl = ifelse(B.graph$Sequence %in% taxa_to_label, B.graph$Sequence, NA)



library(gghighlight)
library(ggrepel)
ggplot(B.graph, aes(x=gamma, y=Effect, group = Sequence)) +
  geom_line() +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  gghighlight(pval < 0.05) +
  xlab(expression(gamma)) +
  ylab("Effect Size") +
  theme(text = element_text(size = 16))

ggsave(file.path("results", "sim-gamma.pdf"), height = 4, units = "in", width = 7)

#-------------------------------------------------------------------------------

# Testing a geo. mean scale model (noise added directly) -----------------------

## Scale model for adding noise to the geometric mean directly
scale.model <- function(gamma, conds, p, mc.samples){
  scale_samples <- matrix(NA, length(p), mc.samples) ## empty container
  
  for(i in 1:length(p)){
    geo_means <- log(apply(p[[i]],2,ALDEx2:::gm))
    noise <- stats::rnorm(mc.samples, 0, gamma)

    scale_samples[i,] <- geo_means + noise
  }
  scale_samples <- log2(exp(scale_samples))
  return(scale_samples)
}

## A (slightly) updated version of the aldex.clr function that uses this scale model
new.clr <- function( reads, conds, mc.samples=128, denom="all", 
                     verbose=FALSE, useMC=FALSE, gamma = NULL) {
  #  invocation:
  #  use selex dataset from ALDEx2 library
  #  x <- aldex.clr( reads, conds, mc.samples=128, denom="all", verbose=FALSE, useMC=FALSE )
  #  this function generates the centre log-ratio transform of Monte-Carlo instances
  #  drawn from the Dirichlet distribution.
  
  # INPUT
  # The 'reads' data.frame MUST have row
  # and column names that are unique, and
  # looks like the following:
  #
  #              T1a T1b  T2  T3  N1  N2
  #   Gene_00001   0   0   2   0   0   1
  #   Gene_00002  20   8  12   5  19  26
  #   Gene_00003   3   0   2   0   0   0
  #       ... many more rows ...
  #
  # ---------------------------------------------------------------------
  
  # OUTPUT
  # The output returned is a list (x) that contains Monte-Carlo instances of
  # the centre log-ratio transformed values for each sample
  # Access to values
  # sample IDs: names(x)
  # number of features (genes, OTUs): length(x[[1]][,1])
  # number of Monte-Carlo Dirichlet instances: length(x[[1]][1,])
  # feature names: rownames(x[[1]])
  
  # Fully validate and coerce the data into required formats
  # coerce SummarizedExperiment reads into data.frame
  summarizedExperiment=FALSE
  if (summarizedExperiment) {
    reads <- data.frame(as.list(assays(reads,withDimnames=TRUE)))
    if (verbose) {
      message("converted SummarizedExperiment read count object into data frame")
    }
  }
  # make sure the conditions vector or matrix is reasonable
  if(missing(conds)){
    
    message("no conditions provided: forcing denom = 'all'")
    message("no conditions provided: forcing conds = 'NA'")
    denom <- "all"
    conds <- rep("NA", ncol(reads))
    
  }
  
  # if a model matrix is supplied, then aldex.effect is not valid
  # force the use of either all for the denominator
  # or
  # the use of a user-supplied denominator
  if(is(conds, "vector")){
    message("conditions vector supplied")
  }    
  else if(is(conds, "matrix")  & all(round(conds) == conds)){
    message("integer matrix provided")
    if(is.vector(denom, mode="numeric")){
      message("user-defined denominator used")
    } else if (denom == "all"){
      message("using all features for denominator")
    } else {
      stop("please supply the desired vector of indices for the denominator")
    }
    #     if(conds.col == 0){
    #       message("conditions provided as matrix: selecting first column for aldex.clr")
    #       conds <- as.character(conds[,1])
    #     }else{
    #       message("conditions provided as matrix: user selected column for aldex.clr")
    #       if(is.numeric(conds.col)){
    #         print(conds[,conds.col])
    #         conds <- as.vector(conds[,conds.col])
    #         prints(conds)
    #         print(length(conds))
    #       }
    #     }
  }
  
  if(ncol(reads) != length(conds) & !is(conds, "matrix")){
    print(length(conds))
    print(ncol(reads))
    stop("mismatch between number of samples and condition vector")
  }
  
  # make sure that the multicore package is in scope and return if available
  has.BiocParallel <- FALSE
  if ("BiocParallel" %in% rownames(installed.packages()) & useMC){
    message("multicore environment is is OK -- using the BiocParallel package")
    #require(BiocParallel)
    has.BiocParallel <- TRUE
  }
  else {
    message("operating in serial mode")
  }
  
  # make sure that mc.samples is an integer, despite it being a numeric type value
  mc.samples <- as.numeric(as.integer(mc.samples))
  
  #  remove all rows with reads less than the minimum set by minsum
  minsum <- 0
  
  # remove any row in which the sum of the row is 0
  z <- as.numeric(apply(reads, 1, sum))
  reads <- as.data.frame( reads[(which(z > minsum)),]  )
  
  if (verbose) message("removed rows with sums equal to zero")
  
  
  #  SANITY CHECKS ON THE DATA INPUT
  if ( any( round(reads) != reads ) ) stop("not all reads are integers")
  if ( any( reads < 0 ) )             stop("one or more reads are negative")
  
  for ( col in names(reads) ) {
    if ( any( ! is.finite( reads[[col]] ) ) )  stop("one or more reads are not finite")
  }
  
  if ( length(rownames(reads)) == 0 ) stop("rownames(reads) cannot be empty")
  if ( length(colnames(reads)) == 0 ) stop("colnames(reads) cannot be empty")
  
  if ( length(rownames(reads)) != length(unique(rownames(reads))) ) stop ("row names are not unique")
  if ( length(colnames(reads)) != length(unique(colnames(reads))) ) stop ("col names are not unique")
  if ( mc.samples < 128 ) warning("values are unreliable when estimated with so few MC smps")
  
  # add a prior expection to all remaining reads that are 0
  # this should be by a Count Zero Multiplicative approach, but in practice
  # this is not necessary because of the large number of features
  prior <- 0.5
  
  # This extracts the set of features to be used in the geometric mean computation
  # returns a list of features
  if(is.null(gamma)){
    feature.subset <- aldex.set.mode(reads, conds, denom)
    if ( length(feature.subset[[1]]) == 0 ) stop("No low variance, high abundance features in common between conditions\nPlease choose another denomiator.")
  } else{
    feature.subset <- vector()
  }
  
  reads <- reads + prior
  
  if (verbose == TRUE) message("data format is OK")
  
  # ---------------------------------------------------------------------
  # Generate a Monte Carlo instance of the frequencies of each sample via the Dirichlet distribution,
  # returns frequencies for each feature in each sample that are consistent with the
  # feature count observed as a proportion of the total counts per sample given
  # technical variation (i.e. proportions consistent with error observed when resequencing the same library)
  
  nr <- nrow( reads )
  rn <- rownames( reads )
  
  #this returns a list of proportions that are consistent with the number of reads per feature and the
  #total number of reads per sample
  
  # environment test, runs in multicore if possible
  if (has.BiocParallel){
    p <- bplapply( reads ,
                   function(col) {
                     q <- t( ALDEx2:::rdirichlet( mc.samples, col ) ) ;
                     rownames(q) <- rn ;
                     q })
    names(p) <- names(reads)
  }
  else{
    p <- lapply( reads ,
                 function(col) {
                   q <- t( ALDEx2:::rdirichlet( mc.samples, col ) ) ;
                   rownames(q) <- rn ; q } )
  }
  
  # sanity check on the data, should never fail
  for ( i in 1:length(p) ) {
    if ( any( ! is.finite( p[[i]] ) ) ) stop("non-finite frequencies estimated")
  }
  
  if (verbose == TRUE) message("dirichlet samples complete")
  
  # ---------------------------------------------------------------------
  # Add scale samples (if desired)
  # Checking the size of the scale samples
  
  if(!is.null(gamma)){
    message("aldex.scaleSim: adjusting samples to reflect scale uncertainty.")
    l2p <- list()
    if(length(gamma) == 1){ ##Add uncertainty around the scale samples
      
      ## grabbing samples from the default scale model
      if(verbose) message("sampling from the default scale model.")
      scale_samples <- scale.model(gamma, conds, p, mc.samples)
      
      for(i in 1:length(p)){
        l2p[[i]] <- sweep(log2(p[[i]]), 2,  scale_samples[i,], "-")
      }
    }
    names(l2p) <- names(p)
  }
  
  # ---------------------------------------------------------------------
  # Take the log2 of the frequency and subtract the geometric mean log2 frequency per sample
  # i.e., do a centered logratio transformation as per Aitchison
  
  # apply the function over elements in a list, that contains an array
  if(is.null(gamma)){
    scale_samples <- NULL
    # DEFAULT
    if (is.list(feature.subset)) {
      # ZERO only
      feat.result <- vector("list", length(unique(conds))) # Feature Gmeans
      condition.list <- vector("list", length(unique(conds)))    # list to store conditions
      
      for (i in 1:length(unique(conds)))
      {
        condition.list[[i]] <- which(conds == unique(conds)[i]) # Condition list
        feat.result[[i]] <- lapply( p[condition.list[[i]]], function(m) {
          apply(log2(m), 2, function(x){mean(x[feature.subset[[i]]])})
        })
      }
      set.rev <- unlist(feat.result, recursive=FALSE) # Unlist once to aggregate samples
      p.copy <- p
      for (i in 1:length(set.rev))
      {
        p.copy[[i]] <- as.data.frame(p.copy[[i]])
        p[[i]] <- apply(log2(p.copy[[i]]),1, function(x){ x - (set.rev[[i]])})
        p[[i]] <- t(p[[i]])
      }
      l2p <- p    # Save the set in order to generate the aldex.clr variable
    } else if (is.vector(feature.subset)){
      # Default ALDEx2, iqlr, user defined, lvha
      # denom[1] is put in explicitly for the user-defined denominator case
      if (has.BiocParallel){
        if (denom[1] != "median"){
          l2p <- bplapply( p, function(m) {
            apply( log2(m), 2, function(col) { col - mean(col[feature.subset]) } )
          })
        } else if (denom[1] == "median"){
          l2p <- bplapply( p, function(m) {
            apply( log2(m), 2, function(col) { col - median(col[feature.subset]) } )
          })
        }
        names(l2p) <- names(p)
      }
      else{
        if (denom[1] != "median"){
          l2p <- lapply( p, function(m) {
            apply( log2(m), 2, function(col) { col - mean(col[feature.subset]) } )
          })
        } else if (denom[1] == "median"){
          l2p <- lapply( p, function(m) {
            apply( log2(m), 2, function(col) { col - median(col[feature.subset]) } )
          })
        }
      }
    }  else {
      message("the denominator is not recognized, use a different denominator")
    }
    
    # sanity check on data
    for ( i in 1:length(l2p) ) {
      if ( any( ! is.finite( l2p[[i]] ) ) ) stop("non-finite log-frequencies were unexpectedly computed")
    }
    if (verbose == TRUE) message("transformation complete")
  }
  
  return(new("aldex.clr",reads=reads,mc.samples=mc.samples,conds=conds,denom=feature.subset,verbose=verbose,useMC=useMC,dirichletData=p,analysisData=l2p, scaleSamps = scale_samples))
}


# Function for running this scale model
CLR.mod <- function(d, n,gamma, seq.depth, pval = 0.05, mc.samples = 500){
  dd <- length(d)/2
  truth1 <- !(d[1:dd]==d[(dd+1):(2*dd)])##testing if the mean is different
  truth2 <- (1:dd)[truth1]##Locations of the differences
  
  fdr <- rep(NA,length(n))
  for(i in 1:length(n)){
    dat <- create_true_abundances(d, n=n[i])
    rdat <- resample_data(dat, seq.depth=seq.depth)
    sample.totals=rowSums(dat[,-1])
    
    countdata <- t(dat[,-1,drop=F])
    colnames(countdata) <- paste0("n", 1:ncol(countdata))

    clr <- new.clr(countdata, as.character(dat$Condition), gamma = gamma, mc.samples = mc.samples)
    ttest <- aldex.ttest(clr)
    pos <- rownames(ttest %>% filter(we.eBH < 0.05))
    tp <- sum(c("Taxa3", "Taxa4", "Taxa15", "Taxa20") %in% pos)
    fp <- length(pos) - tp
    fdr[i] <- fp/(fp+tp)
  }
  fdr <- ifelse(is.na(fdr), 0, fdr)
  return(list(fdr = fdr))
}# end of function

n <- c(25,50,100, 500, 750, 1000)
mod.small <- CLR.mod(d,n,.1,5000)
mod.med <- CLR.mod(d,n,.5,5000)

FDR.small <- mod.small$fdr
FDR.med <- mod.med$fdr

# Plotting the fdr

label <- rep(NA, length(n)*2)
label[c(1, 1 + length(n))] <- c("gamma == 0.1", "gamma == 0.5")
graph.df <- data.frame("N" = rep(n,2), "val" = c(FDR.small, FDR.med), "label" = label, "type" = c(rep("0.1", length(n)), rep("0.5", length(n))))

ggplot(graph.df, aes(x = N, y = val, group = type, color = type, label = label)) + 
  geom_line(aes(linetype = type), lwd = 1.25) +
  scale_color_npg()+
  scale_linetype_manual(values=c("twodash", "longdash", "dotted")) +
  theme_bw() +
  ylab("False Discovery Rate") +
  xlab("Sample Size") +
  theme(legend.position = "none") + 
  geom_label_repel(parse = TRUE,
                   na.rm = TRUE,
                   size = 5,
                   direction="y",
                   nudge_y = .05) +
  ylim(c(0,1)) +
  theme(text=element_text(size=20))

ggsave(file.path("results", "inc-fdr.pdf"), height = 4, units = "in", width = 7)


#-------------------------------------------------------------------------------
