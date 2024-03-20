library(ALDEx2)
library(tidyverse)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggsci)
library(DESeq2)
library(limma)
library(edgeR)
library(baySeq)
set.seed(12345)

# The directory change assumes you are using the R project.
setwd("selex")

# Selex data example

# Data can be loaded directly from ALDEx2
data(selex)
conds <- c(rep("NS",7), rep("S",7))

# Computing the values of \theta^\perp------------------------------------------

# Value implied by the data
normalized <- apply(selex+.5, 2, FUN = function(x){x/sum(x)})
gm_selex <- apply(normalized, 2, FUN = function(x){mean(log2(x))})
mean(-1*mean(gm_selex[1:7])) - mean(-1*mean(gm_selex[8:14]))

## Value implied by the dirichlet samples
## CLR estimate
clr.samps <- aldex.clr(selex, conds, mc.samples = 1000, denom = "all")
dir.sams <- clr.samps@dirichletData

gm_imp <- rep(NA,1000)

for(i in 1:1000){
  tmp <- matrix(NA,nrow = 1600, ncol = 14)
  for(k in 1:14){
    tmp[,k] <- dir.sams[[k]][,i]
  }
  tmp_gm <- apply(tmp, 2, FUN = function(x){mean(log2(x))})
  gm_imp[i] <- mean(-1*mean(tmp_gm[1:7])) - mean(-1*mean(tmp_gm[8:14]))
}

mean(gm_imp)

## Reading in the true positives------------------------------------------------

## True positives
df_true <- read.delim("selex-truth.txt")
true_pos <- c(df_true$X)

#-------------------------------------------------------------------------------


# Selex analysis ---------------------------------------------------------------

# Selex resampling different methods

# Helper functions to run each of these methods
run_aldex2 <- function(dat, conds, denom="all", gamma= NULL, mc.samples = 500){
  aldex.fit <- aldex(dat, conds, denom=denom, mc.samples = mc.samples, gamma = gamma, effect = FALSE) %>% 
    as.data.frame() %>% 
    rownames_to_column("category") %>% 
    dplyr::select(category, we.eBH) %>% 
    mutate(padj=we.eBH) %>% 
    mutate(low=NA, high=NA) %>%
    mutate(sig = ifelse(padj < 0.05, TRUE, FALSE))
  return(aldex.fit)
}

run_deseq2 <- function(dat,conds){
  coldata <- matrix(conds, ncol = 1)
  colnames(coldata) <- "Condition"
  dds <- DESeqDataSetFromMatrix(countData=dat,
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

run_limma <- function(dat,conds){
  coldata <- ifelse(conds == "NS", 0, 1)
  coldata <- matrix(coldata, ncol = 1)
  design <- model.matrix(~0+conds)
  
  y <- DGEList(counts=dat,group=conds)
  y <- calcNormFactors(y, method = "TMM")
  v <- voom(y,design,plot=FALSE)
  fit <- lmFit(v,design)

  contr <- makeContrasts(condsNS - condsS, levels = colnames(coef(fit)))
  contr.fit <- contrasts.fit(fit, contr)
  contr.fit <- eBayes(contr.fit)
  et <- topTable(contr.fit,number =1600, p.value = 1) %>%
    rownames_to_column("category") %>%
    mutate(sig = ifelse(adj.P.Val < 0.05, TRUE, FALSE))
  return(et)
}

run_edger <- function(dat, conds){
  coldata <- ifelse(conds == "NS", 0, 1)
  coldata <- matrix(coldata, ncol = 1)
  
  y <- DGEList(counts=dat,group=c(coldata[,1]))
  y <- calcNormFactors(y)
  y <- estimateDisp(y)
  et <- exactTest(y)
  et <- topTags(et, n=1600, p.value = 1) %>%
    as.data.frame() %>%
    rownames_to_column("category") %>%
    mutate(sig = ifelse(FDR < 0.05, TRUE, FALSE))
  return(et)
}

run_baySeq <- function(dat,conds){
  coldata <- ifelse(conds == "NS", 0, 1)
  groups <- list(NDE=rep(1, length(conds)), 
                 DE=as.numeric(coldata))
  dds <- new("countData", data = as.matrix(dat), groups = groups, replicates = conds)
  libsizes(dds) <- getLibsizes(dds)
  dds <- getPriors.NB(dds, cl = NULL)
  dds <- getLikelihoods(dds,cl=NULL,bootStraps=3,verbose=FALSE)
  res <- topCounts(dds,group="DE", number = 1600)
  
  ##Merge to the data to see which taxa because the names are uninformative
  taxa_names <- dat %>%
    as.data.frame() %>%
    rownames_to_column("category") %>%
    join(res, by = names(dat), type = "right") %>%
    mutate(sig = ifelse(FWER.DE < 0.05, TRUE, FALSE))
  
  return(taxa_names)
  
}

 # Helper function to run all of the methods

runBenchmark_DF <- function(rdat, conds, taxa_truth, pval = 0.05, mc.samples = 500){
  tp <- rep(NA, 7)
  fp <- rep(NA, 7)
  tn <- rep(NA, 7)
  fn <- rep(NA, 7)
  afit <- run_aldex2(rdat, conds, mc.samples = mc.samples)
  afit.sig <- afit %>%
    filter(padj < 0.05)
  afit.nosig <- afit %>%
    filter(padj > 0.05)
  tp[1] <- sum(afit.sig$category %in% taxa_truth)
  fp[1] <- sum(!(afit.sig$category %in% taxa_truth))
  tn[1] <- sum(!(afit.nosig$category %in% taxa_truth))
  fn[1] <- sum((afit.nosig$category %in% taxa_truth))
  
  
  dfit <- run_deseq2(rdat, conds)
  dfit.sig <- dfit %>%
    filter(padj < 0.05)
  dfit.nosig <- dfit %>%
    filter(padj > 0.05)
  tp[2] <- sum(dfit.sig$category %in% taxa_truth)
  fp[2] <- sum(!(dfit.sig$category %in% taxa_truth))
  tn[2] <- sum(!(dfit.nosig$category %in% taxa_truth))
  fn[2] <- sum((dfit.nosig$category %in% taxa_truth))
  
  lfit <- run_limma(rdat, conds) 
  lfit.sig <- lfit %>%
    filter(adj.P.Val < 0.05)
  lfit.nosig <- lfit %>%
    filter(adj.P.Val > 0.05)
  tp[3] <- sum(lfit.sig$category %in% taxa_truth)
  fp[3] <- sum(!(lfit.sig$category %in% taxa_truth))
  tn[3] <- sum(!(lfit.nosig$category %in% taxa_truth))
  fn[3] <- sum((lfit.nosig$category %in% taxa_truth))
  
  bsfit <- run_baySeq(rdat, conds)
  bsfit.sig <- bsfit %>%
    filter(FWER.DE < 0.05)
  bsfit.nosig <- bsfit %>%
    filter(FWER.DE > 0.05)
  tp[4] <- sum(bsfit.sig$category %in% taxa_truth)
  fp[4] <- sum(!(bsfit.sig$category %in% taxa_truth))
  tn[4] <- sum(!(bsfit.nosig$category %in% taxa_truth))
  fn[4] <- sum((bsfit.nosig$category %in% taxa_truth))
  
  efit <- run_edger(rdat, conds)
  efit.sig <- efit %>%
    filter(FDR < 0.05)
  efit.nosig <- efit %>%
    filter(FDR > 0.05)
  tp[5] <- sum(efit.sig$category %in% taxa_truth)
  fp[5] <- sum(!(efit.sig$category %in% taxa_truth))
  tn[5] <- sum(!(efit.nosig$category %in% taxa_truth))
  fn[5] <- sum((efit.nosig$category %in% taxa_truth))
  
  tfit.mid <- run_aldex2(rdat, conds, mc.samples = mc.samples, gamma =  .5)
  tfit.mid.sig <- tfit.mid %>%
    filter(padj < 0.05)
  tfit.mid.nosig <- tfit.mid %>%
    filter(padj > 0.05)
  tp[6] <- sum(tfit.mid.sig$category %in% taxa_truth)
  fp[6] <- sum(!(tfit.mid.sig$category %in% taxa_truth))
  tn[6] <- sum(!(tfit.mid.nosig$category %in% taxa_truth))
  fn[6] <- sum((tfit.mid.nosig$category %in% taxa_truth))
  
  tfit.large <- run_aldex2(rdat, conds, mc.samples = mc.samples, gamma = 5)
  tfit.large.sig <- tfit.large %>%
    filter(padj < 0.05)
  tfit.large.nosig <- tfit.large %>%
    filter(padj > 0.05)
  tp[7] <- sum(tfit.large.sig$category %in% taxa_truth)
  fp[7] <- sum(!(tfit.large.sig$category %in% taxa_truth))
  tn[7] <- sum(!(tfit.large.nosig$category %in% taxa_truth))
  fn[7] <- sum((tfit.large.nosig$category %in% taxa_truth))
  
  res <- data.frame("tp" = tp, "fp" = fp, "tn" = tn, "fn"= fn, "method" = c("ALDEx2 (gamma = 0)", "DESeq2", "limma", "baySeq", "edgeR", "ALDEx2 (gamma = 0.5)", "ALDEx2 (gamma = 5)"))
  return(res)
}# end of function

# Running over sample sizes

n_to_test <-  c(5,7,10,15,20,25,30,40,50,60,70,80,90,100)
benchmark_df <- data.frame()
k <- 3 # Number of replicates

for(j in 1:k){
  for(i in 1:length(n_to_test)){
    # Resampling observations
    if(n_to_test[i] <= 7){
      samp1 <- sample(1:7, n_to_test[i], replace = FALSE)
      samp2 <- sample(8:14, n_to_test[i], replace = FALSE)
    } else{
      samp1 <- sample(1:7, n_to_test[i], replace = TRUE)
      samp2 <- sample(8:14, n_to_test[i], replace = TRUE)
    }
    
    # Selecting the samples and running the benchmark function
    samp_all <- c(samp1,samp2)
    conds2 <- c(rep("NS", n_to_test[i]), rep("S", n_to_test[i]))
    selex_exp <- selex[, samp_all]
    tmp <- runBenchmark_DF(selex_exp, conds2, taxa_truth = true_pos, mc.samples = 1000)
    tmp$n <- rep(n_to_test[i], nrow(tmp))
    tmp$k <- rep(j, nrow(tmp))
    benchmark_df <- rbind(benchmark_df, tmp)
    print(i)
  }
}

# Saving the results 
write.csv(benchmark_df, file.path("results", "benchmarking_by_method_results.csv"))

## Plotting typeI error
benchmark_df$typeI <- benchmark_df$fp/(benchmark_df$tn + benchmark_df$fp)

ggplot(benchmark_df, aes(x=n, y = typeI, group = method, color = method)) +
  #geom_line(aes(linetype = method), lwd = 0.5) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "loess", alpha = 0.05, linewidth = 1, span = 1) +
  scale_color_npg() +
  theme_bw() +
  geom_hline(yintercept = 0.05) +
  theme(legend.title = element_blank(),
        text = element_text(size=16)) +
  xlab("Sample Size") +
  ylim(c(0,1.01)) +
  ylab("Type-I Error Rate")
#geom_label_repel(aes(label = label),
#                 na.rm = TRUE) 

ggsave(file.path("results", "selex-typeI-by-method-by-size.pdf"), height = 7, units = "in", width = 10)

##Plotting sensitivity
benchmark_df$sensitivity <- (benchmark_df$tp + benchmark_df$fn)/(benchmark_df$tp + benchmark_df$fn + benchmark_df$tn + benchmark_df$fp)

ggplot(benchmark_df, aes(x=n, y = sensitivity, group = method, color = method)) +
  #geom_line(aes(linetype = method), lwd = 0.5) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "loess", alpha = 0.05, linewidth = 1, span = 1) +
  scale_color_npg() +
  theme_bw() +
  theme(legend.title = element_blank(),
        text = element_text(size=16)) +
  xlab("Sample Size") +
  ylim(c(0,1.01)) +
  ylab("Sensitivity")
#geom_label_repel(aes(label = label),
#                 na.rm = TRUE) 

ggsave(file.path("results", "selex-sensitivity-by-method-by-size.pdf"), height = 7, units = "in", width = 10)

#-------------------------------------------------------------------------------

# T-test example----------------------------------------------------------------

# Helper function to compare the two different methods
calc_diff <- function(t, df = 10){
  print("Average p-value for two sided tests...")
  p.old <- pt(abs(t), df = df, lower.tail = FALSE)*2
  print(mean(p.old))
  print("Average p-value when controlling for sign...")
  p.new <- 2*pt(t, df = df, lower.tail = FALSE)
  print(min(mean(p.new), 2-mean(p.new)))
}

# based on the aldex2 t.fast function
# trimmed down to be only what we need
t.fast.subset <- function(data, group){
  
  grp1 <- group == unique(group)[1]
  grp2 <- group == unique(group)[2]
  n1 <- sum(grp1)
  n2 <- sum(grp2)
  
  
  t <- multtest::mt.teststat(data, as.numeric(grp1), test = "t", nonpara = "n")
  s1 <- apply(data[, grp1], 1, sd)
  s2 <- apply(data[, grp2], 1, sd)
  df <- ((s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  return(list(t = t, df = df))
}

##selex example
data(selex)
conds <- c(rep("NS", 7), rep("S", 7))
clr <- aldex.clr(selex, conds, mc.samples=2000, denom="all", verbose=FALSE, gamma = 1)

# Most is taken right from aldex.ttest
# get dimensions, names, etc from the input data
smpl.ids <- getSampleIDs(clr)
feature.names <- getFeatureNames(clr)
feature.number <- numFeatures(clr)
mc.instances <- numMCInstances(clr)

conditions <- as.factor( conds )
levels     <- levels( conditions )

if ( length( conditions ) !=  numConditions(clr) ){
  stop(paste("mismatch btw 'length(conditions)' and 'length(names(clr))'. len(condtitions):",
             length(conditions),"len(names(clr)):",numConditions(clr)))}
if ( length( levels ) != 2 ) stop("only two condition levels are currently supported")

# generate the comparison sets from the condition levels
levels <- vector( "list", length( levels ) )
names( levels ) <- levels( conditions )
sets <- names(levels)
setAsBinary <- as.numeric(conditions == sets[1])
setA <- which(conditions == sets[1])
setB <- which(conditions == sets[2])

t.matrix <- as.data.frame(matrix(1, nrow = feature.number, 
                                 ncol = mc.instances))
df.matrix <- t.matrix # duplicate result container

# mc.i is the i-th Monte-Carlo instance
mc.all <- getMonteCarloInstances(clr)
for(mc.i in 1:mc.instances){
  
  # generate a matrix of i-th Monte-Carlo instance, columns are samples, rows are features
  t.input <- sapply(mc.all, function(y){y[, mc.i]})
  
  tmp <- t.fast.subset(t.input, setAsBinary)
  t.matrix[, mc.i] <- tmp$t
  df.matrix[, mc.i] <- tmp$df
  
}


t.matrix <- as.matrix(t.matrix)
df.matrix <- as.matrix(df.matrix)

# Plotting
df <- data.frame("t" = c(t.matrix[4,], t.matrix[635,], t.matrix[20,]), "group" = c(rep("Grp1", ncol(t.matrix)), rep("Grp2", ncol(t.matrix)), rep("Grp3", ncol(t.matrix))))

p <- ggplot(df, aes(x = t, group = group, color = group, fill = group)) +
  geom_density(alpha = 0.1, aes(linetype = group)) +
  theme_bw() +
  scale_fill_npg() +
  scale_color_npg() +
  scale_linetype_manual(values = c("dotted", "twodash", "solid")) +
  theme(legend.position = "none") + 
  xlab("Test Statistic") +
  ylab("Density")

p
ggsave(file.path("results", "tstat_densities.pdf"), height = 7, units = "in", width = 7)

# Using a helper function to caculate the changes
calc_diff(t.matrix[4,], df.matrix[4,])
calc_diff(t.matrix[635,], df.matrix[635,])
calc_diff(t.matrix[20,], df.matrix[20,])

#-------------------------------------------------------------------------------
