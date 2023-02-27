# DATA ----

rm(list=ls())
set.seed(3)

library(GEOquery)
library(readxl)
library(ggplot2)
library(edgeR)
library(DESeq2)
library(dearseq)
library(limma)
library(verification)
library(cowplot)


London <- read_excel("GSE107991_edgeR_normalized_Berry_London.xlsx") # edgeR preprocessed
London_raw <- read_excel("GSE107991_Raw_counts_Berry_London.xlsx") # raw counts
genes <- London$Genes
London <- as.matrix(London[,-c(1:3)])


# metadata

GSE107991_metadata <- GEOquery::getGEO("GSE107991", GSEMatrix = FALSE)

get_info <- function(i){
  name <- GSE107991_metadata@gsms[[i]]@header$source_name_ch1
  name <- gsub("Active_TB", "ActiveTB", name)
  name <- gsub("Test_set", "TestSet", name)
  unlist(strsplit(name, split="_"))
}

infos <- sapply(1:length(GSE107991_metadata@gsms), FUN=get_info)
infos_df <- cbind.data.frame("SampleID" = names(GSE107991_metadata@gsms), t(infos))
rownames(infos_df) <- names(GSE107991_metadata@gsms)
colnames(infos_df)[-1] <- c("Cohort", "Location", "Set", "Status")

# groups

TB373 <- as.data.frame(read_excel("41467_2018_4579_MOESM3_ESM.xlsx"))
genes_373 <- TB373[-c(1:2),1]
TB_control <- genes_373[which(as.numeric(TB373[-c(1:2),5])<0.05)]
TB_LTBI <- genes_373[which(as.numeric(TB373[-c(1:2),8])<0.05)]
LTBI_control <- genes_373[which(as.numeric(TB373[-c(1:2),11])<0.05)]

group_London_TB <- which(infos_df$Status=="ActiveTB") # 21
group_London_LTBI <- which(infos_df$Status=="LTBI") # 21
group_London_Control <- which(infos_df$Status=="Control") # 12


### DIFFERENTIAL EXPRESSION ANALYSIS ----

Y <- ifelse(infos_df$Status=="ActiveTB",1,0)

# ActiveTB vs Control #

group <- rep(NA,length(infos_df$Status[-group_London_LTBI]))
group[which(infos_df$Status[-group_London_LTBI]=="ActiveTB")] <- rep(1,length(which(infos_df$Status[-group_London_LTBI]=="ActiveTB")))
group[which(infos_df$Status[-group_London_LTBI]=="Control")] <- rep(0,length(which(infos_df$Status[-group_London_LTBI]=="Control")))


London_raw <- read_excel("GSE107991_Raw_counts_Berry_London.xlsx")
genes_raw <- London_raw$Genes
London_raw <- London_raw[,-c(1:3)]
London_raw <- as.matrix(London_raw)
dgList <- DGEList(counts=London_raw,group=Y)
cpm <- cpm(dgList)
countCheck <- cpm > 2
keep <- which(rowSums(countCheck) >= 5)
dgList <- dgList[keep,]
dgList <- calcNormFactors(dgList, method="TMM")
dgList_group <- dgList[,-group_London_LTBI]
design <- model.matrix(~group)
dgList_group <- estimateGLMTrendedDisp(dgList_group, design=design)
fit <- glmFit(dgList_group, design)
lrt <- glmLRT(fit) # logFC with edgeR

# dearseq

res_dearseq <- dear_seq(exprmat = London[,-group_London_LTBI],
                                covariates = matrix(rep(1,ncol(London[,-group_London_LTBI])),ncol=1),
                                variables2test = matrix(group, ncol=1),
                                gene_based_weights=FALSE, which_weights="loclin",
                                which_test='permutation',n_perm=1000,
                                preprocessed=TRUE)$pval$adjPval
df <- data.frame(logFC=lrt$table$logFC,res_dearseq)
res <- rbind(subset(df, res_dearseq<.05 & logFC>1),subset(df, res_dearseq<.05 & logFC<(-1)))
genes_dearseq <- genes[as.numeric(rownames(res))]
length(genes_dearseq)

# limma-voom

voom_fit <- voom(dgList_group$counts, cbind(1,group))
voom_fit <- limma::lmFit(voom_fit,cbind(1,group))
voom_fit <- limma::eBayes(voom_fit)
res_voom <- p.adjust(voom_fit$p.value[,2], method = "BH")
df_voom <- data.frame(logFC=lrt$table$logFC,res_voom)
res <- rbind(subset(df_voom, res_voom<.05 & logFC>1),subset(df_voom, res_voom<.05 & logFC<(-1)))
genes_voom <- genes_raw[as.numeric(rownames(res))]
length(genes_voom)

# DESeq2

dds <- DESeqDataSetFromMatrix(countData = dgList_group$counts,
                              colData = cbind.data.frame("int" = 1, "phi" =  matrix(as.factor(group), ncol=1)), design=~1+phi)
dds <- DESeq(dds)
dds2 <- results(dds)
res_deseq <- dds2$padj
df_deseq <- data.frame(logFC=lrt$table$logFC,res_deseq)
res <- rbind(subset(df_deseq, res_deseq<.05 & logFC>1),subset(df_deseq, res_deseq<.05 & logFC<(-1)))
genes_deseq <- genes[as.numeric(rownames(res))]
length(genes_deseq)

# edgeR

ydge <- edgeR::estimateDisp(dgList_group, cbind(1,group), robust=TRUE)
fit <- glmQLFit(dgList_group, cbind(1,group), robust=TRUE, dispersion = ydge$trended.dispersion)
qlf <- glmQLFTest(fit)
pval_edger <- qlf$table[,4]
res_edger <- p.adjust(pval_edger, method = "BH")
df_edger <- data.frame(logFC=lrt$table$logFC,res_edger)
res <- rbind(subset(df_edger, res_edger<.05 & logFC>1),subset(df_edger, res_edger<.05 & logFC<(-1)))
genes_edger <- genes[as.numeric(rownames(res))]
length(genes_edger) 


# ActiveTB vs LTBI #

group <- rep(NA,length(infos_df$Status[-group_London_Control]))
group[which(infos_df$Status[-group_London_Control]=="ActiveTB")] <- rep(1,length(which(infos_df$Status[-group_London_Control]=="ActiveTB")))
group[which(infos_df$Status[-group_London_Control]=="LTBI")] <- rep(0,length(which(infos_df$Status[-group_London_Control]=="LTBI")))
dgList_group <- dgList[,-group_London_Control]
design <- model.matrix(~group)
dgList_group <- estimateGLMTrendedDisp(dgList_group, design=design)
fit <- glmFit(dgList_group, design)
lrt <- glmLRT(fit)

# dearseq

res_dearseq <- dear_seq(exprmat = London[,-group_London_Control],
                                covariates = matrix(rep(1,ncol(London[,-group_London_Control])),ncol=1),
                                variables2test = matrix(group, ncol=1),
                                gene_based_weights=FALSE, which_weights="loclin",
                                which_test='permutation', n_perm=1000,
                                preprocessed=TRUE)$pval$adjPval
df <- data.frame(logFC=lrt$table$logFC,res_dearseq)
res <- rbind(subset(df, res_dearseq<.05 & logFC>1),subset(df, res_dearseq<.05 & logFC<(-1)))
genes_dearseq_2 <- genes[as.numeric(rownames(res))]
length(genes_dearseq_2)

# limma-voom

voom_fit <- voom(dgList_group$counts, cbind(1,group))
voom_fit <- limma::lmFit(voom_fit,cbind(1,group))
voom_fit <- limma::eBayes(voom_fit)
res_voom <- p.adjust(voom_fit$p.value[,2], method = "BH")
df_voom <- data.frame(logFC=lrt$table$logFC,res_voom)
res <- rbind(subset(df_voom, res_voom<.05 & logFC>1),subset(df_voom, res_voom<.05 & logFC<(-1)))
genes_voom2 <- genes_raw[as.numeric(rownames(res))]
length(genes_voom2)

# DESeq2

dds <- DESeqDataSetFromMatrix(countData = dgList_group$counts,
                              colData = cbind.data.frame("int" = 1, "phi" =  matrix(as.factor(group), ncol=1)), design=~1+phi)
dds <- DESeq(dds)
dds2 <- results(dds)
res_deseq <- dds2$padj
df_deseq <- data.frame(logFC=lrt$table$logFC,res_deseq)
res <- rbind(subset(df_deseq, res_deseq<.05 & logFC>1),subset(df_deseq, res_deseq<.05 & logFC<(-1)))
genes_deseq2 <- genes[as.numeric(rownames(res))]
length(genes_deseq2)

# edgeR

ydge <- edgeR::estimateDisp(dgList_group, cbind(1,group), robust=TRUE)
fit <- glmQLFit(dgList_group, cbind(1,group), robust=TRUE, dispersion = ydge$trended.dispersion)
qlf <- glmQLFTest(fit)
pval_edger <- qlf$table[,4]
res_edger <- p.adjust(pval_edger, method = "BH")
df_edger <- data.frame(logFC=lrt$table$logFC,res_edger)
res <- rbind(subset(df_edger, res_edger<.05 & logFC>1),subset(df_edger, res_edger<.05 & logFC<(-1)))
genes_edger2 <- genes[as.numeric(rownames(res))]
length(genes_edger2)

# LTBI vs Control #

group <- rep(NA,length(infos_df$Status[-group_London_TB]))
group[which(infos_df$Status[-group_London_TB]=="LTBI")] <- rep(1,length(which(infos_df$Status[-group_London_TB]=="LTBI")))
group[which(infos_df$Status[-group_London_TB]=="Control")] <- rep(0,length(which(infos_df$Status[-group_London_TB]=="Control")))
dgList_group <- dgList[,-group_London_TB]
design <- model.matrix(~group)
dgList_group <- estimateGLMTrendedDisp(dgList_group, design=design)
fit <- glmFit(dgList_group, design)
lrt <- glmLRT(fit)

# dearseq

res_dearseq <- dear_seq(exprmat = London[,-group_London_TB],
                                covariates = matrix(rep(1,ncol(London[,-group_London_TB])),ncol=1),
                                variables2test = matrix(group, ncol=1),
                                gene_based_weights=FALSE, which_weights="loclin",
                                which_test='permutation', n_perm=1000,
                                preprocessed=TRUE)$pval$adjPval
df <- data.frame(logFC=lrt$table$logFC,res_dearseq)
res <- rbind(subset(df, res_dearseq<.05 & logFC>1),subset(df, res_dearseq<.05 & logFC<(-1)))
genes_dearseq_3 <- genes[as.numeric(rownames(res))]
length(genes_dearseq_3)

# 

genes_dearseq_all <- c(genes_dearseq,genes_dearseq_2)
genes_dearseq_all <- unique(genes_dearseq_all)
length(genes_dearseq_all)

# voom

voom_fit <- voom(dgList_group$counts, cbind(1,group))
voom_fit <- limma::lmFit(voom_fit,cbind(1,group))
voom_fit <- limma::eBayes(voom_fit)
res_voom <- p.adjust(voom_fit$p.value[,2], method = "BH")
df_voom <- data.frame(logFC=lrt$table$logFC,res_voom)
res <- rbind(subset(df_voom, res_voom<.05 & logFC>1),subset(df_voom, res_voom<.05 & logFC<(-1)))
genes_voom3 <- genes_raw[as.numeric(rownames(res))]
length(genes_voom3)

#

genes_voom_all <- c(genes_voom,genes_voom2,genes_voom3)
genes_voom_all <- unique(genes_voom_all)
length(genes_voom_all) # 402

length(intersect(genes_voom_all,genes_dearseq_all))
length(intersect(genes_voom_all,genes_373))

# DESeq2

dds <- DESeqDataSetFromMatrix(countData = dgList_group$counts,
                              colData = cbind.data.frame("int" = 1, "phi" =  matrix(as.factor(group), ncol=1)), design=~1+phi)
dds <- DESeq(dds)
dds2 <- results(dds)
res_deseq <- dds2$padj
df_deseq <- data.frame(logFC=lrt$table$logFC,res_deseq)
res <- rbind(subset(df_deseq, res_deseq<.05 & logFC>1),subset(df_deseq, res_deseq<.05 & logFC<(-1)))
genes_deseq3 <- genes[as.numeric(rownames(res))]
length(genes_deseq3) 

#

genes_deseq_all <- c(genes_deseq,genes_deseq2,genes_deseq3)
genes_deseq_all <- unique(genes_deseq_all)
length(genes_deseq_all)

length(intersect(genes_deseq_all,genes_dearseq_all))
length(intersect(genes_deseq_all,genes_373))

# edgeR

ydge <- edgeR::estimateDisp(dgList_group, cbind(1,group), robust=TRUE)
fit <- glmQLFit(dgList_group, cbind(1,group), robust=TRUE, dispersion = ydge$trended.dispersion)
qlf <- glmQLFTest(fit)
pval_edger <- qlf$table[,4]
res_edger <- p.adjust(pval_edger, method = "BH")
df_edger <- data.frame(logFC=lrt$table$logFC,res_edger)
res <- rbind(subset(df_edger, res_edger<.05 & logFC>1),subset(df_edger, res_edger<.05 & logFC<(-1)))
genes_edger3 <- genes[as.numeric(rownames(res))]
length(genes_edger3)

#

genes_edger_all <- c(genes_edger,genes_edger2,genes_edger3)
genes_edger_all <- unique(genes_edger_all)
length(genes_edger_all)

# Signature ----

dearseq_signature <- London[which(is.element(London$Genes,genes_dearseq_all)),1:3]
deseq_signature <- London[which(is.element(London$Genes,genes_deseq_all)),1:3]
voom_signature <- London[which(is.element(London$Genes,genes_voom_all)),1:3]
Singhania_signature <- London[which(is.element(London$Genes,genes_373)),1:3]

write.csv(dearseq_signature,file="dearseq_signature.csv")
write.csv(deseq_signature,file="deseq_signature.csv")
write.csv(voom_signature,file="voom_signature.csv")
write.csv(Singhania_signature,file="Singhania_signature.csv")


# Venn Diagramm ----

library(VennDiagram)

venn.diagram(list(dearseq = genes_dearseq_all, edgeR=genes_373, DESeq2 = genes_deseq_all, limma_voom = genes_voom_all),
             filename = "Venn_4_new.png", alpha=0.2,
             col= c("#FF9900","#003399","#009900","#440154FF"),
             fill= c("#FF9900","#003399","#009900","#440154FF"))

venn.diagram(list(TB_Control = genes_dearseq, TB_LTBI=genes_dearseq_2),
             filename = "dearseq_venn_new.png", alpha=0.2,
             col= c("red","black"),
             fill= c("red","black"))


### dearseq vs edgeR (373 genes) Figure 4 ----

genes_private_edgeR <- setdiff(genes_373,genes_dearseq_all)
genes_private_dearseq <- setdiff(genes_dearseq_all,genes_373)

X_dearseq <- London[which(is.element(genes,genes_private_dearseq)),]
X_dearseq <- t(X_dearseq)
X_373 <- London[which(is.element(genes,genes_private_edgeR)),]
X_373 <- t(X_373)

Y <- rep(NA,length(infos_df$Status))
Y[which(infos_df$Status=="ActiveTB")] <- rep(1,length(which(infos_df$Status=="ActiveTB")))
Y[which(infos_df$Status=="LTBI")] <- rep(0,length(which(infos_df$Status=="LTBI")))
Y[which(infos_df$Status=="Control")] <- rep(0,length(which(infos_df$Status=="Control")))

# Brier score boxplot

temp_373 <- matrix(0,length(group),dim(X_373)[2])

for (i in 1:dim(X_373)[1]){
  for (j in 1:dim(X_373)[2]){
    fit_373 <- glm(Y~X,data = data.frame(Y=as.factor(Y[-i]),X = X_373[-i,j]),family=binomial())
    prob_373 <- predict(fit_373,newdata=data.frame(X=X_373[i,j]), type="response")
    temp_373[i,j] <- prob_373
  }
}

BS_rl_373 <- NULL
for (i in 1:dim(X_373)[2]){
  BS_rl_373[i] <- brier(Y, temp_373[,i])$bs
}

temp_dearseq <- matrix(0,length(group),dim(X_dearseq)[2])
for (i in 1:dim(X_dearseq)[1]){
  for (j in 1:dim(X_dearseq)[2]){
    fit_dearseq <- glm(Y~X,data = data.frame(Y=as.factor(Y[-i]),X = X_dearseq[-i,j]),family=binomial())
    prob_dearseq <- predict(fit_dearseq,newdata=data.frame(X=X_dearseq[i,j]), type="response")
    temp_dearseq[i,j] <- prob_dearseq
  }
}

BS_rl_dearseq <- NULL
for (i in 1:dim(X_dearseq)[2]){
  BS_rl_dearseq[i] <- brier(Y, temp_dearseq[,i])$bs
}

genes_common <- intersect(genes_373,genes_dearseq_all) 
X_common <- London[which(is.element(genes,genes_common)),]
X_common<- t(X_common)
temp_common <- matrix(0,length(group),dim(X_common)[2])

for (i in 1:dim(X_common)[1]){
  for (j in 1:dim(X_common)[2]){
    fit_common <- glm(Y~X,data = data.frame(Y=as.factor(Y[-i]),X = X_common[-i,j]),family=binomial())
    prob_common <- predict(fit_common,newdata=data.frame(X=X_common[i,j]), type="response") # On predit le ieme decoupage
    temp_common[i,j] <- prob_common
  }
}

BS_rl_common <- NULL
for (i in 1:dim(X_common)[2]){
  BS_rl_common[i] <- brier(Y, temp_common[,i])$bs
}

data <- data.frame(method=factor(c(rep("dearseq",length(BS_rl_dearseq)),rep("edgeR",length(BS_rl_373)))),Brier.Score=c(BS_rl_dearseq,BS_rl_373))
data$group <- factor(data$method, levels = c("dearseq", "edgeR"))

my_comparisons <- list(c("dearseq", "edgeR"))
g1 <- ggboxplot(data, x = "method", y = "Brier.Score",
                color = "method")+  theme(legend.position='none') +
  scale_color_manual("Method", values = c("#FF9900",viridis::viridis(n=4)[2], grDevices::heat.colors(n=2)[3])) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + ylab("Brier Score")

# Brier score of each gene called significant using `dearseq` and `edgeR`.

df <- data.frame(Brier.score = c(BS_rl_dearseq,BS_rl_373,BS_rl_common),
                 rank=rank( c(BS_rl_dearseq,BS_rl_373,BS_rl_common),ties.method = "random"),
                 method=as.factor(c(rep("dearseq",length(BS_rl_dearseq)),rep("edgeR",length(BS_rl_373)),
                                    rep("common",length(BS_rl_common)))))
#summary(df)
df$method <- factor(as.character(df$method), ordered = TRUE, levels=c("common", "edgeR", "dearseq"))
df <- df[order(df$method),]

g2 <- ggplot(df, aes(x=rank, y=Brier.score, color=method)) +
  theme_minimal() +
  scale_color_manual("Method", values = c("#CCCCCC",viridis::viridis(n=4)[2],"#FF9900", grDevices::heat.colors(n=2)[3])) +
  geom_point(size=1.5,aes(alpha=method)) + 
  theme(axis.title = element_text(size=16), axis.text = element_text(size=12), legend.title = element_text(size=16), legend.text = element_text(size=14)) +
  guides(color = guide_legend(override.aes = list(alpha=1))) +
  scale_alpha_manual("Method", values=rep(0.5, 3), breaks=c("common","edgeR", "dearseq")) + 
  geom_rug(alpha=0.5)  + 
  ylab("Brier Score")
g2


# p-values

pval_dearseq <- rep(0,dim(X_dearseq)[2])
for (j in 1:dim(X_dearseq)[2]){
  rl_dearseq <- glm(Y~X,data = data.frame(Y=as.factor(Y),X = X_dearseq[,j]),family=binomial())
  pval_dearseq[j] <- coef(summary(rl_dearseq))[2,4]
}

pval_edger <- rep(0,dim(X_373)[2])
for (j in 1:dim(X_373)[2]){
  rl_edger <- glm(Y~X,data = data.frame(Y=as.factor(Y),X = X_373[,j]),family=binomial())
  pval_edger[j] <- coef(summary(rl_edger))[2,4]
}

pval_common <- rep(0,dim(X_common)[2])
for (j in 1:dim(X_common)[2]){
  rl_common <- glm(Y~X,data = data.frame(Y=as.factor(Y),X = X_common[,j]),family=binomial())
  pval_common[j] <- coef(summary(rl_common))[2,4]
}

df <- data.frame(method=factor(c(rep("common",length(pval_common)),rep("dearseq",length(pval_dearseq)),rep("edgeR",length(pval_edger)))),
                 p.values=c(pval_common,pval_dearseq,pval_edger),rank=rank(c(pval_common,pval_dearseq,pval_edger),
                                                                           ties.method = "random"))
#summary(df)
df$method <- factor(as.character(df$method), ordered = TRUE, levels=c("common", "edgeR", "dearseq"))
df <- df[order(df$method),]

g3 <- ggplot(df, aes(x=rank, y=p.values, color=method)) +
  theme_minimal() +
  scale_color_manual("Method", values = c("#CCCCCC",viridis::viridis(n=4)[2],"#FF9900", grDevices::heat.colors(n=2)[3])) +
  geom_point(size=1.5,aes(alpha=method)) + 
  theme(axis.title = element_text(size=16), axis.text = element_text(size=12), legend.title = element_text(size=16), legend.text = element_text(size=14)) +
  guides(color = guide_legend(override.aes = list(alpha=1))) +
  scale_alpha_manual("Method", values=rep(0.5, 3), breaks=c("common","edgeR", "dearseq")) + 
  geom_rug(alpha=0.5) + geom_hline(yintercept=0.05,col="red") + 
  ylab("p-values") 
g3

plot_grid(g1, g2, g3, labels=c("A", "B", "C"), ncol = 1, nrow = 3)

ggsave(filename = "dearseq_edger_new.pdf",dpi = 400, units = "in", width = 8.5, height = 11.5, device="pdf")


### dearseq vs voom / dearseq vs DESeq2 ----

## dearseq vs voom

genes_private_voom <- setdiff(genes_voom_all,genes_dearseq_all)
genes_private_dearseq <- setdiff(genes_dearseq_all,genes_voom_all)

X_dearseq <- London[which(is.element(genes,genes_private_dearseq)),]
X_dearseq <- t(X_dearseq)
X_voom <- London[which(is.element(genes,genes_private_voom)),]
X_voom <- t(X_voom)

# Brier score boxplot

temp_voom <- matrix(0,length(group),dim(X_voom)[2])

for (i in 1:dim(X_voom)[1]){
  for (j in 1:dim(X_voom)[2]){
    fit_voom <- glm(Y~X,data = data.frame(Y=as.factor(Y[-i]),X = X_voom[-i,j]),family=binomial())
    prob_voom <- predict(fit_voom,newdata=data.frame(X=X_voom[i,j]), type="response") # On predit le ieme decoupage
    temp_voom[i,j] <- prob_voom
  }
}

BS_rl_voom <- NULL
for (i in 1:dim(X_voom)[2]){
  BS_rl_voom[i] <- brier(Y, temp_voom[,i])$bs
}

temp_dearseq <- matrix(0,length(group),dim(X_dearseq)[2])
for (i in 1:dim(X_dearseq)[1]){
  for (j in 1:dim(X_dearseq)[2]){
    fit_dearseq <- glm(Y~X,data = data.frame(Y=as.factor(Y[-i]),X = X_dearseq[-i,j]),family=binomial())
    prob_dearseq <- predict(fit_dearseq,newdata=data.frame(X=X_dearseq[i,j]), type="response") # On predit le ieme decoupage
    temp_dearseq[i,j] <- prob_dearseq
  }
}

BS_rl_dearseq <- NULL
for (i in 1:dim(X_dearseq)[2]){
  BS_rl_dearseq[i] <- brier(Y, temp_dearseq[,i])$bs
}

genes_common <- intersect(genes_voom_all,genes_dearseq_all)
X_common <- London[which(is.element(genes,genes_common)),]
X_common<- t(X_common)
temp_common <- matrix(0,length(group),dim(X_common)[2])

for (i in 1:dim(X_common)[1]){
  for (j in 1:dim(X_common)[2]){
    fit_common <- glm(Y~X,data = data.frame(Y=as.factor(Y[-i]),X = X_common[-i,j]),family=binomial())
    prob_common <- predict(fit_common,newdata=data.frame(X=X_common[i,j]), type="response") # On predit le ieme decoupage
    temp_common[i,j] <- prob_common
  }
}

BS_rl_common <- NULL
for (i in 1:dim(X_common)[2]){
  BS_rl_common[i] <- brier(Y, temp_common[,i])$bs
}

data <- data.frame(method=factor(c(rep("common",length(BS_rl_common)),rep("limma-voom",length(BS_rl_voom)))),Brier.Score=c(BS_rl_common,BS_rl_voom))
data$group <- factor(data$method, levels = c("common", "limma-voom"))

my_comparisons <- list(c("limma-voom", "common"))
p1 <- ggboxplot(data, x = "method", y = "Brier.Score",
                color = "method")+  theme(legend.position='none')+
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + ylab("Brier score") +
  scale_color_manual("method", values = c('#CCCCCC',viridis::viridis(n=4)[1], grDevices::heat.colors(n=2)[3]))
p1

# Brier score of each gene called significant using `dearseq` and `voom`.

df <- data.frame(Brier.score = c(BS_rl_voom,BS_rl_common,BS_rl_dearseq),
                 rank=rank(c(BS_rl_voom,BS_rl_common, BS_rl_dearseq),ties.method = "random"),
                 method=as.factor(c(rep("limma-voom",length(BS_rl_voom)),
                                    rep("common",length(BS_rl_common)),rep("dearseq",length(BS_rl_dearseq)))))

df$method <- factor(as.character(df$method), ordered = TRUE, levels=c("common", "limma-voom", "dearseq"))
df <- df[order(df$method),]

p2 <- ggplot(df, aes(x=rank, y=Brier.score, color=method)) +
  theme_minimal() +
  scale_color_manual("Method", values = c("#CCCCCC",viridis::viridis(n=4)[1],"#FF9900", grDevices::heat.colors(n=2)[3])) +
  geom_point(size=1.5,aes(alpha=method)) + 
  theme(axis.title = element_text(size=16), axis.text = element_text(size=12), legend.title = element_text(size=16), legend.text = element_text(size=14)) +
  guides(color = guide_legend(override.aes = list(alpha=1))) +
  scale_alpha_manual("Method", values=rep(0.5, 3), breaks=c("common","limma-voom", "dearseq")) + 
  geom_rug(alpha=0.5)  + 
  ylab("Brier Score")
p2


# p-values

pval_common <- rep(0,dim(X_common)[2])
for (j in 1:dim(X_common)[2]){
  rl_common <- glm(Y~X,data = data.frame(Y=as.factor(Y),X = X_common[,j]),family=binomial())
  pval_common[j] <- coef(summary(rl_common))[2,4]
}

pval_voom <- rep(0,dim(X_voom)[2])
for (j in 1:dim(X_voom)[2]){
  rl_voom <- glm(Y~X,data = data.frame(Y=as.factor(Y),X = X_voom[,j]),family=binomial())
  pval_voom[j] <- coef(summary(rl_voom))[2,4]
}

pval_dearseq <- rep(0,dim(X_dearseq)[2])
for (j in 1:dim(X_dearseq)[2]){
  rl_dearseq <- glm(Y~X,data = data.frame(Y=as.factor(Y),X = X_dearseq[,j]),family=binomial())
  pval_dearseq[j] <- coef(summary(rl_dearseq))[2,4]
}

df <- data.frame(method=factor(c(rep("common",length(pval_common)),rep("limma-voom",length(pval_voom)),rep("dearseq",length(pval_dearseq)))),
                 p.values=c(pval_common,pval_voom,pval_dearseq),rank=rank(c(pval_common,pval_voom,pval_dearseq),
                                                                          ties.method = "random"))

df$method <- factor(as.character(df$method), ordered = TRUE, levels=c("common", "limma-voom", "dearseq"))
df <- df[order(df$method),]

p3 <- ggplot(df, aes(x=rank, y=p.values, color=method)) +
  theme_minimal() +
  scale_color_manual("Method", values = c("#CCCCCC",viridis::viridis(n=4)[1],"#FF9900", grDevices::heat.colors(n=2)[1])) +
  geom_point(size=1.5,aes(alpha=method)) + 
  theme(axis.title = element_text(size=16), axis.text = element_text(size=12), legend.title = element_text(size=16), legend.text = element_text(size=14)) +
  guides(color = guide_legend(override.aes = list(alpha=1))) +
  scale_alpha_manual("Method", values=rep(0.5, 3), breaks=c("common","limma-voom", "dearseq")) + 
  geom_rug(alpha=0.5)  + geom_hline(yintercept=0.05,col="red") +
  ylab("p-values")
p3


## dearseq vs deseq2

genes_private_deseq <- setdiff(genes_deseq_all,genes_dearseq_all)
genes_private_dearseq <- setdiff(genes_dearseq_all,genes_deseq_all)

length(genes_private_dearseq)
length(genes_private_deseq)

X_dearseq <- London[which(is.element(genes,genes_private_dearseq)),]
X_dearseq <- t(X_dearseq)
X_deseq <- London[which(is.element(genes,genes_private_deseq)),]
X_deseq <- t(X_deseq)

# Brier score boxplot

temp_deseq <- matrix(0,length(group),dim(X_deseq)[2])

for (i in 1:dim(X_deseq)[1]){
  for (j in 1:dim(X_deseq)[2]){
    fit_deseq <- glm(Y~X,data = data.frame(Y=as.factor(Y[-i]),X = X_deseq[-i,j]),family=binomial())
    prob_deseq <- predict(fit_deseq,newdata=data.frame(X=X_deseq[i,j]), type="response") # On predit le ieme decoupage
    temp_deseq[i,j] <- prob_deseq
  }
}

BS_rl_deseq <- NULL
for (i in 1:dim(X_deseq)[2]){
  BS_rl_deseq[i] <- brier(Y, temp_deseq[,i])$bs
}

genes_common <- intersect(genes_deseq_all,genes_dearseq_all) # 234
X_common <- London[which(is.element(genes,genes_common)),]
X_common<- t(X_common)
temp_common <- matrix(0,length(group),dim(X_common)[2])

for (i in 1:dim(X_common)[1]){
  for (j in 1:dim(X_common)[2]){
    fit_common <- glm(Y~X,data = data.frame(Y=as.factor(Y[-i]),X = X_common[-i,j]),family=binomial())
    prob_common <- predict(fit_common,newdata=data.frame(X=X_common[i,j]), type="response") # On predit le ieme decoupage
    temp_common[i,j] <- prob_common
  }
}

BS_rl_common <- NULL
for (i in 1:dim(X_common)[2]){
  BS_rl_common[i] <- brier(Y, temp_common[,i])$bs
}

data <- data.frame(method=factor(c(rep("DESeq2",length(BS_rl_deseq)),rep("common",length(BS_rl_common)))),Brier.Score=c(BS_rl_deseq,BS_rl_common))
data$group <- factor(data$method, levels = c("DESeq2","common"))

my_comparisons <- list(c("DESeq2", "common"))
p4 <- ggboxplot(data, x = "method", y = "Brier.Score",
                color = "method")+  theme(legend.position='none') + scale_color_manual("method", values = c('#CCCCCC',viridis::viridis(n=4)[3], grDevices::heat.colors(n=2)[3])) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + ylab("Brier score")

# Brier score of each gene called significant using `dearseq` and `deseq`.

df <- data.frame(Brier.score = c(BS_rl_deseq,BS_rl_common),
                 rank=rank( c(BS_rl_deseq,BS_rl_common),ties.method = "random"),
                 method=as.factor(c(rep("DESeq2",length(BS_rl_deseq)),
                                    rep("common",length(BS_rl_common)))))

df$method <- factor(as.character(df$method), ordered = TRUE, levels=c("common", "DESeq2", "dearseq"))
df <- df[order(df$method),]

p5 <- ggplot(df, aes(x=rank, y=Brier.score, color=method)) +
  theme_minimal() +
  scale_color_manual("Method", values = c("#CCCCCC",viridis::viridis(n=4)[3],"#FF9900", grDevices::heat.colors(n=2)[3])) +
  geom_point(size=1.5,aes(alpha=method)) + 
  theme(axis.title = element_text(size=16), axis.text = element_text(size=12), legend.title = element_text(size=16), legend.text = element_text(size=14)) +
  guides(color = guide_legend(override.aes = list(alpha=1))) +
  scale_alpha_manual("Method", values=rep(0.5, 3), breaks=c("common","DESeq2", "dearseq")) + 
  geom_rug(alpha=0.5)  + 
  ylab("Brier Score")
p5


# Marginal p-values

pval_deseq <- rep(0,dim(X_deseq)[2])
for (j in 1:dim(X_deseq)[2]){
  rl_deseq <- glm(Y~X,data = data.frame(Y=as.factor(Y),X = X_deseq[,j]),family=binomial())
  pval_deseq[j] <- coef(summary(rl_deseq))[2,4]
}

pval_common <- rep(0,dim(X_common)[2])
for (j in 1:dim(X_common)[2]){
  rl_common <- glm(Y~X,data = data.frame(Y=as.factor(Y),X = X_common[,j]),family=binomial())
  pval_common[j] <- coef(summary(rl_common))[2,4]
}

df <- data.frame(method=factor(c(rep("common",length(pval_common)),rep("DESeq2",length(pval_deseq)))),
                 p.values=c(pval_common,pval_deseq),rank=rank(c(pval_common,pval_deseq),
                                                              ties.method = "random"))

# plot

df$method <- factor(as.character(df$method), ordered = TRUE, levels=c("common", "DESeq2", "dearseq"))
df <- df[order(df$method),]

p6 <- ggplot(df, aes(x=rank, y=p.values, color=method)) +
  theme_minimal() +
  scale_color_manual("Method", values = c("#CCCCCC",viridis::viridis(n=4)[3],"#FF9900", grDevices::heat.colors(n=2)[3])) +
  geom_point(size=1.5,aes(alpha=method)) + 
  theme(axis.title = element_text(size=16), axis.text = element_text(size=12), legend.title = element_text(size=16), legend.text = element_text(size=14)) +
  guides(color = guide_legend(override.aes = list(alpha=1))) +
  scale_alpha_manual("Method", values=rep(0.5, 3), breaks=c("common","DESeq2", "dearseq")) + 
  geom_rug(alpha=0.5) + geom_hline(yintercept=0.05,col="red") +
  ylab("p-values")
p6

plot_grid(p1, p4, p2, p5, p3, p6, labels=c("A", "D", "B", "E", "C", "F"), ncol = 2, nrow = 3)

ggsave(filename = "dearseq_deseq_voom_new.pdf",dpi = 400, units = "in", width = 14.5, height = 14.5, device="pdf")


### limma-voom / dearseq / edgeR (373 genes) / DESeq2 ----

# edgeR (373 genes)

X_373 <- London[which(is.element(genes,genes_373)),] 
X_373 <- t(X_373)
temp_373 <- matrix(0,length(group),dim(X_373)[2])

for (i in 1:dim(X_373)[1]){
  for (j in 1:dim(X_373)[2]){
    fit_373 <- glm(Y~X,data = data.frame(Y=as.factor(Y[-i]),X = X_373[-i,j]),family=binomial())
    prob_373 <- predict(fit_373,newdata=data.frame(X=X_373[i,j]), type="response") # On predit le ieme decoupage
    temp_373[i,j] <- prob_373
  }
}

BS_rl_373 <- NULL
for (i in 1:dim(X_373)[2]){
  BS_rl_373[i] <- brier(Y, temp_373[,i])$bs
}

# dearseq

X_dearseq <- London[which(is.element(genes,genes_dearseq_all)),]
X_dearseq <- t(X_dearseq)
temp_dearseq <- matrix(0,length(group),dim(X_dearseq)[2])

for (i in 1:dim(X_dearseq)[1]){
  for (j in 1:dim(X_dearseq)[2]){
    fit_dearseq <- glm(Y~X,data = data.frame(Y=as.factor(Y[-i]),X = X_dearseq[-i,j]),family=binomial())
    prob_dearseq <- predict(fit_dearseq,newdata=data.frame(X=X_dearseq[i,j]), type="response") # On predit le ieme decoupage
    temp_dearseq[i,j] <- prob_dearseq
  }
}

BS_rl_dearseq <- NULL
for (i in 1:dim(X_dearseq)[2]){
  BS_rl_dearseq[i] <- brier(Y, temp_dearseq[,i])$bs
}

# voom

X_voom <- London[which(is.element(genes,genes_voom_all)),]
X_voom<- t(X_voom)
temp_voom <- matrix(0,length(group),dim(X_voom)[2])

for (i in 1:dim(X_voom)[1]){
  for (j in 1:dim(X_voom)[2]){
    fit_voom <- glm(Y~X,data = data.frame(Y=as.factor(Y[-i]),X = X_voom[-i,j]),family=binomial())
    prob_voom <- predict(fit_voom,newdata=data.frame(X=X_voom[i,j]), type="response") # On predit le ieme decoupage
    temp_voom[i,j] <- prob_voom
  }
}

BS_rl_voom <- NULL
for (i in 1:dim(X_voom)[2]){
  BS_rl_voom[i] <- brier(Y, temp_voom[,i])$bs
}

# deseq2

X_deseq <- London[which(is.element(genes,genes_deseq_all)),]
X_deseq<- t(X_deseq)
temp_deseq <- matrix(0,length(group),dim(X_deseq)[2])

for (i in 1:dim(X_deseq)[1]){
  for (j in 1:dim(X_deseq)[2]){
    fit_deseq <- glm(Y~X,data = data.frame(Y=as.factor(Y[-i]),X = X_deseq[-i,j]),family=binomial())
    prob_deseq <- predict(fit_deseq,newdata=data.frame(X=X_deseq[i,j]), type="response") # On predit le ieme decoupage
    temp_deseq[i,j] <- prob_deseq
  }
}

BS_rl_deseq <- NULL
for (i in 1:dim(X_deseq)[2]){
  BS_rl_deseq[i] <- brier(Y, temp_deseq[,i])$bs
}

# plot 

data <- data.frame(method=factor(c(rep("DESeq2",length(BS_rl_deseq)),rep("edgeR",length(BS_rl_373)),rep("limma-voom",length(BS_rl_voom)),rep("dearseq",length(BS_rl_dearseq)))),
                   Brier.Score=c(BS_rl_deseq,BS_rl_373,BS_rl_voom,BS_rl_dearseq))
data$method <- factor(data$method, levels = c("DESeq2","edgeR","limma-voom", "dearseq"))


my_comparisons <- list( c("dearseq", "edgeR"), c("dearseq", "limma-voom"), c("edgeR", "limma-voom"), c("edgeR", "DESeq2"), c("dearseq", "DESeq2"), c("DESeq2", "limma-voom") )
ggboxplot(data, x = "method", y = "Brier.Score", color = "method") +
  theme(legend.position='none')+ scale_color_manual("Method", labels = c("DESeq2" ,"limma-voom", "edgeR", "dearseq"),
                                                    values = c(viridis::viridis(n=4)[c(3,2,1)],"#FF9900", grDevices::heat.colors(n=4)[3])) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + ylab("Brier score")


ggsave(filename = "comparisons_new.pdf",dpi = 400, units = "in", width = 8.5, height = 6.5)



