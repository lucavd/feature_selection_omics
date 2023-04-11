library(MASS)
library(GEOquery)
library(edgeR)
library(DESeq2)
library(limma)
library(readxl)
library(dearseq)
library(fitdistrplus)
library(reshape2)
library(ggplot2)

# # data set
London_raw<- read_excel("GSE107991_Raw_counts_Berry_London.xlsx")

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
#
group_London_TB <- which(infos_df$Status=="ActiveTB") # 21
group_London_LTBI <- which(infos_df$Status=="LTBI") # 21
group_London_Control <- which(infos_df$Status=="Control") # 12

Y <- rep(NA,length(infos_df$Status))
Y[which(infos_df$Status=="ActiveTB")] <- rep(1,length(which(infos_df$Status=="ActiveTB")))
Y[which(infos_df$Status=="LTBI")] <- rep(0,length(which(infos_df$Status=="LTBI")))
Y[which(infos_df$Status=="Control")] <- rep(0,length(which(infos_df$Status=="Control")))


genes_raw <- London_raw$Genes
London_raw <- London_raw[,-c(1:3)]
London_raw <- as.matrix(London_raw)

dgList <- DGEList(counts=London_raw,group=Y)
cpm <- cpm(dgList)
countCheck <- cpm > 2
keep <- which(rowSums(countCheck) >= 5)
London_raw <- London_raw[keep,] 

London_H0 <- London_raw[,-group_London_TB]

# MLE

# MLE_size_0 <- rep(NA,nrow(London_H0))
# MLE_mu_0 <- rep(NA,nrow(London_H0))
# 
# for (i in 1:nrow(London_H0)){
#   tryCatch({
#     print(i)
#     fit_0 <- fitdist(as.integer(London_H0[i,]), "nbinom")$estimate
#     MLE_size_0[i] <- fit_0[1]
#     MLE_mu_0[i] <- fit_0[2]
#   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
# }

# Genes simulation

load("MLE_size_0_save.RData")
load("MLE_mu_0_save.RData")

MLE_mu_0 <- MLE_mu_0_save
MLE_size_0 <- MLE_size_0_save

# Missing Values

ind_NA <- which(is.na(MLE_mu_0))
MLE_mu_0<- MLE_mu_0[-ind_NA]
MLE_size_0 <- MLE_size_0[-ind_NA]

# Generate 10 000 genes

sampling_matrix_gene <- function(size_sample,n_genes,n_H1){
  res_0 <- matrix(0,n_genes,size_sample)
  res_1 <- matrix(0,n_genes,size_sample)
  pert <- rep(0,length(n_genes))
  
  for (i in 1:n_genes){
    
    mu_i <- MLE_mu_0[i]
    size_i <- MLE_size_0[i]
    pert <- runif(size_sample,0.3,0.8)
    
    if (size_sample==4){
      while( (length(unique(res_0[i,1:(size_sample/2)]))==1) & (length(unique(res_0[i,((size_sample/2)+1):size_sample]))==1)){
        
        gene_i <- rnbinom(n = size_sample, size = size_i, mu=mu_i)
        res_0[i,] <- gene_i 
        
        if (i%%2==0){
          res_1[i,] <- round(gene_i + (gene_i*pert))
        }
        
        else{
          res_1[i,] <- round(gene_i - (gene_i*pert))
        }
      }
    }
    
    else{
      while((length(unique(res_0[i,1:(size_sample/2)]))<4) & (length(unique(res_0[i,((size_sample/2)+1):size_sample]))<4)){
        
        gene_i <- rnbinom(n = size_sample, size = size_i, mu=mu_i)
        res_0[i,] <- gene_i 
        
        if (i%%2==0){
          res_1[i,] <- round(gene_i + (gene_i*pert))
        }
        else{
          res_1[i,] <- round(gene_i - (gene_i*pert))
        }
      }
    }
  }
  
  res <- res_0
  res[(n_genes-(n_H1-1)):n_genes,((size_sample/2)+1):size_sample] <- res_1[(n_genes-(n_H1-1)):n_genes,
                                                                           ((size_sample/2)+1):size_sample]
  return(res)
}

# Truth : 9000 non DE and 1000 DE genes

truth_bin <- c(rep(0,9500),rep(1,500))

# Functions

FDR <- function(x){
  if (sum(x<0.05)==0){
    res <- 0
  }
  if (sum(x<0.05)!=0){
    res <- sum(x[which(truth_bin==0)]<0.05)/sum(x<0.05)
  }
  
  return(res)
}


type1_error <- function(x){
  if (sum(x<0.05)==0){
    res <- 0
  }
  if (sum(x<0.05)!=0){
    res <- sum(x[which(truth_bin==0)]<0.05)/length(which(truth_bin==0))
  }
  return(res)
}

TDR <- function(x){
  if (sum(x<0.05)==0){
    res <- 0
  }
  if (sum(x<0.05)!=0){
    res <- sum(x[which(truth_bin==1)]<0.05)/sum(x<0.05)
  }
  
  return(res)
}

stat_power <- function(x){
  if (sum(x<0.05)==0){
    res <- 0
  }
  if (sum(x<0.05)!=0){
    res <- sum(x[which(truth_bin==1)]<0.05)/length(which(truth_bin==1))
  }
  return(res)
}

# Simulation

size <- c(4,8,16,50,100,150,200,300)
fdrs <- matrix(NA,length(size),5)
type1 <- matrix(NA,length(size),5)
tdrs <- matrix(NA,length(size),5)
pwr <- matrix(NA,length(size),5)


for (n in 1:length(size)){
  
  # Genes Matrix: 10000 genes of which 1000 are DE
  genes_matrix <- sampling_matrix_gene(size[n],10000,500)
  
  # 50% of 0 and 50% of 1
  condition<- c(rep(0,ncol(genes_matrix)/2),rep(1,ncol(genes_matrix)/2))
  
  # dearseq
  
  res_dearseq <- dear_seq(exprmat = genes_matrix,
                          covariates = matrix(rep(1,ncol(genes_matrix)),ncol=1),
                          variables2test = matrix(condition, ncol=1),
                          gene_based_weights=FALSE, which_weights="loclin",
                          which_test='asymptotic', progressbar = TRUE,
                          preprocessed=FALSE)
  
  
  res_dearseq_perm <- dear_seq(exprmat = genes_matrix,
                               covariates = matrix(rep(1,ncol(genes_matrix)),ncol=1),
                               variables2test = matrix(condition, ncol=1),
                               gene_based_weights=FALSE, which_weights="loclin",
                               which_test='permutation', progressbar = TRUE,
                               preprocessed=FALSE)
  
  # voom
  
  voom_fit <- voom(genes_matrix, cbind(1,condition))
  voom_fit <- limma::lmFit(voom_fit,cbind(1,condition))
  voom_fit <- limma::eBayes(voom_fit)
  res_voom <- data.frame(raw=voom_fit$p.value[,2],adj=p.adjust(voom_fit$p.value[,2], method = "BH"))
  
  # DESeq2
  
  dds <- DESeqDataSetFromMatrix(countData = genes_matrix,
                                colData = cbind.data.frame("int" = 1, "phi" =  matrix(as.factor(condition), ncol=1)),
                                design=~1+phi)
  
  dds <- DESeq(dds )
  dds2 <- results(dds, cooksCutoff=FALSE )
  res_deseq <- dds2
  
  # edgeR
  
  ydge <- edgeR::estimateDisp(genes_matrix, cbind(1,condition), robust=TRUE)
  fit <- glmQLFit(genes_matrix, cbind(1,condition), robust=TRUE, dispersion = ydge$trended.dispersion)
  qlf <- glmQLFTest(fit)
  res_edger <- data.frame(raw=qlf$table[,4],padj=p.adjust(qlf$table[,4], method = "BH"))
  
  # p-values data frame
  
  pvs <- data.frame(deseq=res_deseq$pvalue,
                    edger=res_edger$raw,
                    voom=res_voom$raw,
                    dearseq=res_dearseq$pvals$rawPval,
                    dearseq_perm=res_dearseq_perm$pvals$rawPval)
  
  # FDR
  fdrs[n,] <- apply(apply(pvs, 2, p.adjust, "BH"), 2, function(x){FDR(x)})
  # Type-1 error
  type1[n,] <- apply(pvs, 2, function(x){type1_error(x)})
  # TDR
  tdrs[n,] <- apply(apply(pvs, 2, p.adjust, "BH"), 2, function(x){TDR(x)})
  # Power
  pwr[n,] <- apply(pvs, 2, function(x){stat_power(x)})
  
}

df <- data.frame(setting=rep("nb",40),n=rep(size,5),
                 method=factor(as.character(rep(c("voom","edger","deseq","vs","vsp"),each=8)),ordered=TRUE),
                 fdr=matrix(fdrs,ncol=1),
                 ti_err= matrix(type1,ncol=1),
                 tpr=matrix(tdrs,ncol=1),
                 pwr= matrix(pwr,ncol=1))


############## PLOT ########################

res_nb_singhania <- data.frame(setting=rep("nb",40),n=rep(size,5),
                  method=factor(as.character(rep(c("voom","edger","deseq","vs","vsp"),each=8)),ordered=TRUE),
                  fdr=matrix(fdrs,ncol=1),
                  ti_err= matrix(type1,ncol=1),
                  tpr=matrix(tdrs,ncol=1),
                  pwr= matrix(pwr,ncol=1))
