source("sim_utils.R")
library(biomformat)
library(abind)
library(stringr)
library(lmerTest)
library(merTools)
library(tidyverse)
library(ggplot2)
library(ggpattern)
library(brms)
library(ExtDist)

## Load sequence count data
biom_data <- read_hdf5_biom("oral_trimmed_deblur.biom")
Y <- do.call(rbind, biom_data$data)
row.names(Y) <- sapply(biom_data$rows, function(item){item$id})

## Get metadata, add average flowcount if needed
metadata <- read.csv("oral_trimmed_metadata.csv", sep="\t")
subset_metadata <- metadata[,c("HostSubject", "Timepoint.", "X.SampleID", "brushing_event", "flow.cell.5min.1", "flow.cell.5min.2")]
subset_metadata[,"flowcount"] <- round(rowMeans(subset_metadata[,c("flow.cell.5min.1", "flow.cell.5min.2")]))

## Add taxonomy info
tax <- read.csv("taxonomy.tsv", sep="\t")
inds <- match(row.names(Y), tax[,"Feature.ID"])
row.names(Y) <- tax[inds,"Taxon"]

## Collapse to Genus Level
Y_genus <- c()
genus <- unname(sapply(row.names(Y), function(x) {strsplit(x, split=";")[[1]][6]}))
unique_genus <- unique(genus[!is.na(genus)])
all_inds <- c()
for(genus_n in unique_genus) {
    genus_inds <- which(genus%in%genus_n)
    all_inds <- c(all_inds, genus_inds)
    if(length(genus_inds)>1) {
        new_row <- colSums(Y[genus_inds,])
    } else {
        new_row <- Y[genus_inds,]
    }
    Y_genus <- rbind(Y_genus, new_row)
}
row.names(Y_genus) <- unique_genus

## Collapse rows with less than a third of samples having counts
other <- colSums(Y[which(is.na(genus)),])
other <- other + colSums(Y_genus[rowSums(Y_genus>0)<12,])
Y_genus <- Y_genus[rowSums(Y_genus>0)>=12,]
other <- other + Y_genus[" g__",]
Y_genus <- Y_genus[row.names(Y_genus)!=" g__",]
Y_genus <- rbind(Y_genus, other)

## Build Metadata Matrix for Mixed Effects Model
meta_mat <- cbind(c(log2(subset_metadata[,"flow.cell.5min.1"]),
                    log2(subset_metadata[,"flow.cell.5min.2"])),
                  c(subset_metadata[,"Timepoint."],
                    subset_metadata[,"Timepoint."]),
                  c(subset_metadata[,"HostSubject"],
                    subset_metadata[,"HostSubject"]),
                  c(subset_metadata[,"brushing_event"],
                    subset_metadata[,"brushing_event"]))
colnames(meta_mat) <- c("flowcounts", "time_point", "Host", "brushing_event")
meta_mat <- data.frame(meta_mat)
meta_mat$brushing_event <- factor(meta_mat$brushing_event, levels=c("before", "after"))
meta_mat$flowcounts <- as.numeric(meta_mat$flowcounts)
meta_mat$time_point_me <- apply(meta_mat, 1, function(row) {
    if(row[2]%in%c(1,2)) {
        return("morning")
    } else {
        return("evening")
    }
})

## Mixed Effects Model with 95% Confidence Interval
# set.seed(54623)
# 
# mixed_effect_model <- brm(flowcounts~brushing_event*time_point_me+(1|Host), meta_mat, iter = 10000, warmup = 6000 )
# 
# # Drawing from the posterior predictive distribution
# scale_gs <- 2^(t(posterior_linpred(mixed_effect_model, ndraws = 1000, newdata = meta_mat[1:32,])))
# 
# ## Run Analyses
# data <- list()
# data$Y <- Y_genus
# data$X <- cbind(1, as.numeric(subset_metadata[,"brushing_event"]=="after"))
# 
# 
# ## Gold Standard Using CI & Aldex2 Sens Interval Test
# aldex2_gold <- run_aldex2(data, 1000, denom="all", gamma = scale_gs)

## Mixed Effects Model with 95% Confidence Interval
mixed_effect_model <- lmer(flowcounts~brushing_event+(1|Host)+(1|time_point_me), meta_mat)
z <- lm(flowcounts~brushing_event, data=meta_mat)
1-(2^(coef(z)[2]))
scale_ci <- confint(mixed_effect_model, level=0.95)[5,]
summary(mixed_effect_model)

set.seed(54623)
## Run Analyses
data <- list()
data$Y <- Y_genus
data$X <- cbind(1, as.numeric(subset_metadata[,"brushing_event"]=="after"))

## Gold Standard Using CI & Aldex2 Sens Interval Test
aldex2_gold <- run_indexa(data, 1000, scale_ci[1], scale_ci[2], denom="none")
aldex2_gold[(aldex2_gold[,c("ctt.pval.BH.adj")]<0.05)|(aldex2_gold[,c("gtt.pval.BH.adj")]<0.05),]
aldex2_gold[c(" g__Rothia", " g__Selenomonas", " g__Schwartzia"),]

# Bio informed model
# This replicates his results
# Just make this a uniform distribution
# Do not touch this
outside.int <- c(-6.6, log2(1.25)) # Top interval comes from Funahara et al
scale_bio <- matrix(1, nrow = 32, ncol = 1000)

for(n in seq(2,32,by = 2)){
  scale_bio[n,] <- runif(1000, min = 2^outside.int[1], max = 2^outside.int[2])
}

aldex2_bio <- run_aldex2(data, 1000, denom="all", gamma = scale_bio)


## Existing Methods
## DESeq2
deseq2_res <- run_deseq2(data, "poscounts")
## Limma
limma_res <- run_limma(data)
## edgeR
edger_res <- run_edger(data)
## ALDEx2
aldex2_res <- run_aldex2(data, 1000, "BH")
## ALDEx2 default scale
aldex2_ds_res <- run_aldex2(data, 1000, "BH", gamma = 0.5)


## False/True Positive/Negative
## Helper Function
get_tf_np <- function(obs_lfc, obs_pvals, gold_standard_lfc, tp_rows, tn_rows) {
    res <- rep("", length(obs_lfc))
    tp <- (sign(obs_lfc[tp_rows])==sign(gold_standard_lfc[tp_rows]))
    tp <- tp&(obs_pvals[tp_rows]<0.05)
    res[tp_rows][tp] <- "TP"
    res[tp_rows][which(obs_pvals[tp_rows]>=0.05)] <- "FN"
    res[tn_rows][which(obs_pvals[tn_rows]<0.05)] <- "FP"
    res[tn_rows][which(obs_pvals[tn_rows]>=0.05)] <- "TN"
    return(res)
}
tp_tn_fp_fn <- c()

tp_rows <- which(row.names(aldex2_gold)%in%row.names(aldex2_gold[aldex2_gold[,"gtt.pval.BH.adj"] < 0.05,]))
tn_rows <- which(row.names(aldex2_gold)%in%row.names(aldex2_gold[aldex2_gold[,"gtt.pval.BH.adj"] >= 0.05,]))

tp_tn_fp_fn <- cbind(tp_tn_fp_fn, get_tf_np(aldex2_gold[,"median.effect"],
                                            aldex2_gold[,"gtt.pval.BH.adj"],
                                            aldex2_gold[,"median.effect"], tp_rows, tn_rows))
tp_tn_fp_fn <- cbind(tp_tn_fp_fn, get_tf_np(deseq2_res[,"log2FoldChange"],
                                            deseq2_res[,"padj"],
                                            aldex2_gold[,"median.effect"], tp_rows, tn_rows))
tp_tn_fp_fn <- cbind(tp_tn_fp_fn, get_tf_np(limma_res[,"logFC"],
                                            limma_res[,"adj.P.Val"],
                                            aldex2_gold[,"median.effect"], tp_rows, tn_rows))
tp_tn_fp_fn <- cbind(tp_tn_fp_fn, get_tf_np(edger_res[,"logFC"],
                                            edger_res[,"FDR"],
                                            aldex2_gold[,"median.effect"], tp_rows, tn_rows))
tp_tn_fp_fn <- cbind(tp_tn_fp_fn, get_tf_np(aldex2_res[,"effect"],
                                            aldex2_res[,"we.eBH"],
                                            aldex2_gold[,"median.effect"], tp_rows, tn_rows))
tp_tn_fp_fn <- cbind(tp_tn_fp_fn, get_tf_np(aldex2_ds_res[,"effect"],
                                            aldex2_ds_res[,"we.eBH"],
                                            aldex2_gold[,"median.effect"], tp_rows, tn_rows))
tp_tn_fp_fn <- cbind(tp_tn_fp_fn, get_tf_np(aldex2_bio[,"effect"],
                                            aldex2_bio[,"we.eBH"],
                                            aldex2_gold[,"median.effect"], tp_rows, tn_rows))


colnames(tp_tn_fp_fn) <- c("Gold Standard", "DESeq2", "limma", "edgeR", "ALDEx2",
                           "ALDEx2 (   = 0.5)", "ALDEx2 (Biological Model)")
row.names(tp_tn_fp_fn) <- row.names(aldex2_gold)


pl <- tp_tn_fp_fn %>%
  as.data.frame() %>%
  rownames_to_column("Taxa") %>%
  mutate(Taxa = trimws(str_replace(Taxa, "g__", ""))) %>%
  filter(Taxa %in% c("Streptococcus", "Haemophilus", "Rothia", "[Prevotella]", "Selenomonas", "Schwartzia")) %>%
  mutate(Taxa = str_replace(Taxa, "\\[", "")) %>%
  mutate(Taxa = str_replace(Taxa, "\\]", "")) %>%
  pivot_longer(-Taxa, names_to = "Model", values_to = "sigcode") %>%
  mutate(Taxa=factor(Taxa, levels = c("Haemophilus", "Streptococcus", "Prevotella", "Rothia", "Schwartzia", "Selenomonas")), sigcode=factor(sigcode, 
                                           levels=c("TP", "TN", 
                                                    "FP", "FN"))) %>% 
  mutate(Model=factor(Model, levels = c("ALDEx2 (   = 0.5)","ALDEx2 (Biological Model)",  "ALDEx2", "limma", "edgeR", "DESeq2", "Gold Standard"))) %>% 
  ggplot(aes(x=Taxa, y=Model)) +
  geom_tile_pattern(aes(fill=sigcode, pattern = sigcode), color="darkgrey",pattern_fill = 'grey',pattern_colour  = 'black', pattern_density = 0.015, pattern_key_scale_factor = 0.5) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        legend.title=element_blank(),
        text = element_text(size=16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Model") +
  scale_pattern_manual(values = c(TP = "none", TN = "none", FP = "stripe", FN = "stripe")) +
  scale_fill_manual(values= c("grey", "white", "grey", "white"))

pl

ggsave("oral-microbiome-plot-2.pdf", pl, height = 6, width = 8)

