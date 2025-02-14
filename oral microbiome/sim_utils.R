library(ALDEx2)
library(DESeq2)
library(limma)
library(edgeR)
library(baySeq)
library(latex2exp)
library(lemon)
library(abind)
library(INDExA)

create_true_abundances <- function(d, n) {
    W <- replicate(n/2, rpois(ncol(d), lambda=d[1,]))
    W <- cbind(W, replicate(n/2, rpois(ncol(d), lambda=d[2,])))
    X <- cbind(1, c(rep(0, n/2), rep(1, n/2)))
    tp <- which(d[1,]!=d[2,])
    return(list(W=W, X=X, tp=tp))
}

resample_data <- function(W, seq_depth) {
    Y <- apply(W, 2, function(col) {
        rmultinom(1, size=seq_depth, prob=col/sum(col))
    })
    return(Y)
}

run_aldex2 <- function(data, mc.samples, fdr.method, alpha=0.05, denom="all",
                       gamma=NULL) {
    conds <- data$X
    colnames(conds) <- c("Intercept", "Cond")
    res <- aldex(data$Y, conds[,2], mc.samples=mc.samples,
                 denom=denom, verbose=F, gamma=gamma)
    row_names <- row.names(data$Y)
    new_res <- c()
    for(n in row_names) {
        new_res <- rbind(new_res, res[n,])
    }
    row.names(new_res) <- row_names
    return(new_res)
}

run_deseq2 <- function(data, sfType="ratio") {
    counts <- data$Y
    
    colnames(counts) <- paste0("n", 1:ncol(counts))
    coldata <- data.frame(Condition=data$X[,2])
    rownames(coldata) <- paste0("n", 1:ncol(counts))
    dds <- DESeqDataSetFromMatrix(countData=counts,
                                  colData=coldata,
                                  design=~Condition)
    dds <- results(DESeq(dds, sfType=sfType))
    return(dds)
}

run_limma <- function(data) {
    counts <- data$Y
    colnames(counts) <- paste0("n", 1:ncol(counts))
    rownames(counts) <- paste0("Taxa", 1:nrow(counts))
    coldata <- data.frame(Condition=as.character(data$X[,2]))
    rownames(coldata) <- paste0("n", 1:ncol(counts))
    design <- model.matrix(~0+Condition, data=coldata)
    y <- DGEList(counts=counts)
    y <- calcNormFactors(y)
    v <- voom(y, design, plot=F)
    fit <- lmFit(v, design)
    contr <- makeContrasts(Condition1 - Condition0, levels=colnames(coef(fit)))
    contr.fit <- contrasts.fit(fit, contr)
    contr.fit <- eBayes(contr.fit)
    res <- topTable(contr.fit, number=nrow(counts), p.value = 1)[row.names(counts),]
    return(res)
}

run_edger <- function(data) {
    counts <- data$Y
    colnames(counts) <- paste0("n", 1:ncol(counts))
    rownames(counts) <- paste0("Taxa", 1:nrow(counts))
    coldata <- data.frame(Condition=as.character(data$X[,2]))
    rownames(coldata) <- paste0("n", 1:ncol(counts))
    design <- model.matrix(~coldata$Condition)
    y <- DGEList(counts=counts, group=c(coldata$Condition))
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design)
    et <- exactTest(y)
    res <- as.data.frame(topTags(et, n=nrow(counts), p.value=1)[row.names(counts),])
    return(res)
}

run_baySeq <- function(data) {
    counts <- data$Y
    colnames(counts) <- paste0("n", 1:ncol(counts))
    rownames(counts) <- paste0("Taxa", 1:nrow(counts))
    coldata <- data.frame(Condition=as.character(data$X[,2]))
    rownames(coldata) <- paste0("n", 1:ncol(counts))
    groups <- list(NDE=rep(1, nrow(coldata)), 
                   DE=as.numeric(coldata$Condition))
    dds <- new("countData", data=counts, groups=groups,
               replicates=as.character(coldata$Condition))
    libsizes(dds) <- getLibsizes(dds)
    dds <- getPriors.NB(dds, cl=NULL)
    dds <- getLikelihoods(dds, cl=NULL, bootStraps=3, verbose=F)
    res <- topCounts(dds, group="DE", number=20)
    taxa_names <- apply(res, 1, function(row) {
        for(taxa in row.names(counts)) {
            colz <- c("n1", "n2", "n3", "n4", "n5")
            if(all(counts[taxa, colz]==as.numeric(row[colz]))) {
                return(taxa)
            }
        }
    })
    row.names(res) <- unname(taxa_names)
    res <- res[row.names(counts),]
    return(res)
}

run_indexa <- function(data, mc.samples, epsilon_l, epsilon_u,
                       denom="none", equal.var=T) {
  indexa.obj <- indexa.mc(data$Y, mc.samples, denom=denom)
  res <- indexa.test(indexa.obj, data$X[,2], c(epsilon_l, epsilon_u),
                     equal.var=equal.var)
  row_names <- row.names(data$Y)
  new_res <- c()
  for(n in row_names) {
    if(!(n %in% row.names(res))) {
      new_res <- rbind(new_res, NA)
    } else {    
      new_res <- rbind(new_res, res[n,])
    }
  }
  row.names(new_res) <- row_names
  return(new_res)
}