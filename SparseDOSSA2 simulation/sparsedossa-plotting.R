source("sim_utils.R")
library(ggplot2)
library(ExtDist)
library(gridExtra)
library(ggpubr)

set.seed(3207485)
no <- 10000
seq_depths <- runif(no, 50000, 200000)
data_other <- sim_dossa2(no, "Stool", seq_depths)
species_names <- row.names(data_other$Y)
keep_taxa <- names(which((rowSums(data_other$Y!=0)/no)>=0.25))
other_taxa <- names(which((rowSums(data_other$Y!=0)/no)<0.25))
print(length(keep_taxa))

## LFCs
lfcs <- abs(rt(65, 4, 1))
subset_species <- sample(keep_taxa, length(lfcs))
prevalences <- lfcs/2
tpi <- match(subset_species, keep_taxa)
fpi <- (1:length(keep_taxa))[-tpi]
true_lfcs <- rep(0, length(keep_taxa)+1)
true_lfcs[tpi] <- lfcs

pos <- length(tpi)
neg <- length(fpi)+1 ## other category for +1

rename_vec <- c("ALDEx2", "limma", "DESeq2", "edgeR", "ALDEx2 - Default Scale (Gamma = 0.5)", "ALDEx2 - Default Scale (Gamma = 1)", "ALDEx2 - Right Skew")
names(rename_vec) <- c("aldex2_padj", "limma_padj", "deseq2_padj", "edger_padj", "aldex2_sm_0.5_padj", "aldex2_sm_1_padj", "aldex2_br_padj")

lfc_sorted <- sort(true_lfcs, decreasing=T)
df <- data.frame(lfcs=lfc_sorted)
write.csv(df, file.path("output", "fig_2_lfcs.csv"))

fpr_df <- c()
power_df <- c()
reps <- 100
for(n in c(14,30,50,70,100,120,150,170,200,250,300)) {
    res_l <- readRDS(paste0("output/sparse_dossa_out_", reps, "_", n, "_BH.RDS"))
    for(m in colnames(res_l$all_padj[,,1])) {
        if(m=="aldex2_0_padj") {next}
        power_avg <- 0
        fdr_avg <- 0
        fp_avg <- 0
        tp_avg <- 0
        for(r in 1:(dim(res_l$all_padj)[3])) {
            col <- res_l$all_padj[,m,r]
            col[is.na(col)] <- 1
            tp <- sum(col[tpi] < 0.05)
            fp <- sum(col[fpi] < 0.05)
            tn <- sum(col[fpi] >= 0.05)
            fn <- sum(col[tpi] >= 0.05)
            if(is.na(tp)&is.na(fp)) {
                ;
            } else if((tp==0)&(fp==0)) {
                ;
            } else {
                tp_avg <- tp_avg + tp
                fp_avg <- fp_avg + fp
                power_avg <- power_avg + (tp/pos)
                fdr_avg <- fdr_avg + (fp/(fp+tp))
            }
        }
        print(c(m, n, fp_avg/reps, tp_avg/reps))
        fpr_df <- rbind(fpr_df, c(n, rename_vec[m], fdr_avg/reps))
        power_df <- rbind(power_df, c(n, rename_vec[m], power_avg/reps))
    }
}


colnames(fpr_df) <- c("n", "method", "fpr")
fpr_df <- data.frame(fpr_df)
fpr_df[,1] <- as.numeric(fpr_df[,1])
#fpr_df[,2] <- as.factor(fpr_df[,2])
fpr_df[,3] <- as.numeric(fpr_df[,3])

g1 <- ggplot(fpr_df, aes(x=n, y=fpr, color=method, linetype=method))
g1 <- g1 + geom_smooth(alpha=0, size=1.5)
g1 <- g1 + theme_bw() + ylim(0, 0.5)
g1 <- g1 + geom_hline(yintercept=0.05, linetype="dotted", alpha=0.5)
g1 <- g1 + scale_color_manual(values=c("#000000",   "#d62728", "#d62728", "#1f77b4", "#666666", "#9E9E9E", "#CCCCCC"))
g1 <- g1 + scale_linetype_manual(values=c("dotted", "solid", "dotdash", "dashed", "dotted", "dotted",
                                          "dotted")) 
g1 <- g1 + theme(legend.position="bottom")
g1 <- g1 + ylab("False Positive Rate") + xlab("Sample Size")
g1 <- g1 + theme(legend.title=element_blank())
g1 <- g1 + theme(legend.key.size = unit(1, "cm"))

ggsave("output/geom_smooth_plot_2_fpr_data.pdf", g1, height = 6, width = 8)

colnames(power_df) <- c("n", "method", "power")
power_df <- data.frame(power_df)
power_df[,1] <- as.numeric(power_df[,1])
#power_df[,2] <- as.factor(power_df[,2])
power_df[,3] <- as.numeric(power_df[,3])

g1 <- ggplot(power_df, aes(x=n, y=power, color=method, linetype=method))
g1 <- g1 + geom_smooth(alpha=0, size=1.5)
g1 <- g1 + theme_bw() + ylim(0, 0.8)
g1 <- g1 + geom_hline(yintercept=0.05, linetype="dotted", alpha=0.5)
g1 <- g1 + scale_color_manual(values=c("#000000",   "#d62728", "#d62728", "#1f77b4", "#666666", "#9E9E9E", "#CCCCCC"))
g1 <- g1 + scale_linetype_manual(values=c("dotted", "solid", "dotdash", "dashed", "dotted", "dotted",
                                          "dotted")) 
g1 <- g1 + theme(legend.position="bottom")
g1 <- g1 + ylab("Power") + xlab("Sample Size")
g1 <- g1 + theme(legend.title=element_blank())
g1 <- g1 + theme(legend.key.size = unit(1, "cm"))

ggsave("output/geom_smooth_plot_2_power_data.pdf", g1, height = 6, width = 8)



