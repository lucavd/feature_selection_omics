library(tcgsaseq)
library(dplyr)
library(edgeR)
library(DESeq2)
library(tidyr)
library(stringr)
library(ggplot2)

pn <- 0
str_remove <- function (string, pattern) {
  str_replace(string, pattern, "")
}

#--------------------------------
## sim fns
#---------------------------------
sim_nb_data <- function(n = 100,
                        beta = 2,
                        nGenes = 1000,
                        re_sd = 1,
                        gene_sd = 1,
                        p0 = 0.9,
                        two_group = TRUE) {

  mu_x <- rexp(n)*10
  # mu_x <- rnorm(n)
  x <- cbind(1, rnorm(n, mean = mu_x))
  if (two_group) {
    z <- sample(0:1, size = n, replace = TRUE)
  } else z <- rnorm(n)


  b_0 <- replicate(nGenes, rnorm(n))
  b_1 <- replicate(nGenes, rnorm(n, sd = re_sd))
  gene_f <- rnorm(nGenes)*gene_sd
  mu <- y <- b_0
  int <- 1000

  n.h1 <- round(nGenes*(1-p0))
  for (ii in 1:ncol(b_0)) {
    if (ii > n.h1) {
      mu[,ii] <- int + b_0[,ii] + rowSums(x)
    } else mu[,ii] <- int + b_0[,ii] + rowSums(x) + (b_1[,ii] + beta + gene_f[ii])*(z*x[,2])

    mu[,ii] <- ifelse(mu[,ii] < 0, 0, mu[,ii])
    y[,ii] <- rnbinom(n, mu = mu[,ii], size = rexp(n)) + 1
  }
  indiv <- 1:n
  list(x = x, z = z, y = y, indiv = indiv)
}

sim_nonlin_data <- function(n = 250,
                            beta = 2,
                            nGenes = 1000,
                            rho = 0.5,
                            re_sd = 1,
                            gene_sd = 1,
                            p0 = 0.9) {
  # browser()
  x <- cbind(1, matrix(rnorm(n, mean = 100, sd = 50), n, 1))
  u <- rexp(nGenes, rate = 1/100)
  e <- t(replicate(nGenes, rnorm(n)))*u + u
  grp <- rnorm(n)

  b.1 <- matrix(rnorm(n, sd = re_sd), n, 1)
  bb.1 <- t(t(b.1) + beta)
  n.h1 <- round(nGenes*(1-p0))
  gene_f <- rnorm(n.h1, sd = gene_sd)

  esd <- 0.01*apply(e, 1, sd)
  b_0 <- matrix(rnorm(n*nGenes, sd = 0.01*esd), n, nGenes) %>% t
  eta <- t(e + b_0) + rowMeans(x)
  eta <- cbind(eta[,1:n.h1] + matrix(grp * bb.1, nrow(eta), n.h1) +
                 matrix(grp, nrow(eta), n.h1)*matrix(gene_f, nrow(eta), n.h1, byrow = TRUE),
               eta[,-(1:n.h1)]) +
    rnorm(n*nGenes, sd = rexp(n*nGenes))
  eta <- eta*rowMeans(eta)/1000

  y <- cbind(eta[,1:n.h1] + matrix(grp * bb.1, nrow(eta), n.h1) +
               matrix(grp, nrow(eta), n.h1)*matrix(gene_f, nrow(eta), n.h1, byrow = TRUE),
             eta[,-(1:n.h1)]) +
    exp(rnorm(n*nGenes, sd = rexp(n*nGenes)))
  y <- ifelse(y <= 0, 1e-7, y)
  indiv <- 1:n
  y <- pmin(pmax(ceiling(y), 1), 1e9)

  list(x = x, z = grp, y = y, indiv = indiv)
}

#----------------------
## setting up parms
#-----------------------
p1 <- 0.3
nGenes <- 10000

  b <- -75
  r <- 0
  g <- 0
  n <- 10
  m <- 'nb'

    sim_data <- sim_nb_data(n = n,
                            beta = b,
                            nGenes = nGenes,
                            re_sd = r,
                            gene_sd = g,
                            p0 = 1-p1,
                            two_group = FALSE)

  y <- sim_data$y
  x <- sim_data$x
  z <- sim_data$z %>% as.matrix
  design_r <- cbind(x, z)


  vsi_res <- try(varseq(exprmat = t(y), covariates = x, variables2test = z, which_test = 'asymptotic', doPlot = FALSE, gene_based_weights = FALSE), silent = TRUE)
  if (class(vsi_res) != 'try-error') {
    vsi_p <- data.frame(
      gene = 1:nGenes,
      vsi_pval = vsi_res$pvals$rawPval,
      vsi_pval_adj = vsi_res$pvals$adjPval
    )
  } else {
    vsi_p <- data.frame(
      gene = 1:nGenes,
      vsi_pval = NA,
      vsi_pval_adj = NA
    )
  }

  vsn_res <- varseq(exprmat = t(y), covariates = x, variables2test = z, which_test = 'asymptotic', doPlot = FALSE, which_weights = 'none')
  vsn_p <- data.frame(
    gene = 1:nGenes,
    vs_pval = vsn_res$pvals$rawPval,
    vs_pval_adj = vsn_res$pvals$adjPval
  )
  # if (n <= 100) {
  #system.time(
    vsp_res <- varseq(exprmat = t(y), covariates = x, variables2test = z, which_test = 'permutation', doPlot = FALSE, which_weights = 'none',
                                n_perm = 1000)
  #  )
  vsp_p <- data.frame(
    gene = 1:nGenes,
    vsp_pval = vsp_res$pvals$rawPval,
    vsp_pval_adj = vsp_res$pvals$FDR
  )
  qplot(vsp_res$pvals$rawPval, vsp_res$pvals$adjPval) + geom_abline(slope=1, intercept=0, col="red") +xlim(0,1) + ylim(0,1)
  sum(vsp_res$pvals$FDR<0.05)
  nreject <- min(sum(vsp_res$pvals$adjPval < 0.05), 1, na.rm=TRUE)
  eFDR <- sum(vsp_res$pvals$adjPval[-(1:(nGenes*p1))] < 0.05)/nreject
  eTDR <- sum(vsp_res$pvals$adjPval[1:(nGenes*p1)] < 0.05)/nreject
  cat("5% FDR:", eFDR, " TDR:", eTDR, "\n")
  plot(y = sapply(seq(0, 1, by=0.001), function(x){max(sum(vsp_res$pvals$adjPval[-(1:(nGenes*p1))] < x)/sum(vsp_res$pvals$adjPval < x), 0, na.rm=TRUE)}),
      x = seq(0, 1, by=0.001),
      type = "l", xlab = "Nominal FDR level", ylab = "Empirical FDR", col = "red", lwd = 2,
      ylim = c(0,1))
  abline(a = 0, b = 1, lty = 2)


  ydge <- edgeR::DGEList(counts=t(y))
  ydge <- edgeR::calcNormFactors(ydge)
  ydge <- edgeR::estimateDisp(ydge, design_r, robust=TRUE)
  edger_fit <- edgeR::glmQLFit(ydge, design = design_r)
  edger_res <- edgeR::glmQLFTest(edger_fit, coef = 3)
  edger_p <- data.frame(
    gene = 1:nGenes,
    edger_pval = edger_res$table$PValue
  ) %>%
    mutate(edger_pval_adj = p.adjust(edger_pval, method = 'fdr'))

  vv <- limma::voom(ydge,design_r)
  voom_fit <- limma::lmFit(vv,design_r)
  voom_fit <- limma::eBayes(voom_fit)
  voom_p <- data.frame(
    gene = 1:nGenes,
    voom_pval = voom_fit$p.value[,3]
  ) %>%
    mutate(voom_pval_adj = p.adjust(voom_pval, method = 'fdr'))

  # if (n <=150) {
  dsq_df <- cbind.data.frame("indiv"=as.factor(sim_data$indiv),
                             "time"=as.numeric(z),
                             "x"=as.numeric(x[,2]))
  y_dsq <- DESeq2::DESeqDataSetFromMatrix(countData = t(y),
                                          colData = dsq_df,
                                          design = ~  time + x)


  res_dsq <- try(DESeq2::DESeq(y_dsq, test="LRT", reduced = ~ x), silent = TRUE)
  if (class(res_dsq) == 'try-error') {
    dsq_p <-
      data.frame(
        gene = 1:nGenes,
        dsq_pval = NA,
        dsq_pval_adj = NA)
  } else {
    dsq_p <-
      data.frame(
        gene = 1:nGenes,
        dsq_pval = results(res_dsq)$pvalue
      ) %>%
      mutate(dsq_pval_adj = p.adjust(dsq_pval, method = 'fdr'))
  }
  # }


  all_p <-
    vsn_p %>%
    inner_join(vsi_p) %>%
    inner_join(voom_p) %>%
    inner_join(edger_p)
  # if (n <= 150) {
  all_p <- all_p %>%
    inner_join(dsq_p)
  # }
  # if (n <= 100) {
  all_p <- all_p %>%
    inner_join(vsp_p)
  # }

  n1 <- if_else(b == 0 & r == 0 & g == 0, 0, round(nGenes*p1))
  res_p <- all_p %>%
    select(
      -contains('adj')
    ) %>%
    gather(method, pval, -gene) %>%
    group_by(method) %>%
    summarise(
      ti_err = mean(pval[gene > n1] < 0.05),
      pwr = mean(pval[gene <= n1] < 0.05)
    ) %>%
    mutate(method = str_remove(method, '_pval'))
  adj_res_p <- all_p %>%
    select(gene, #vs = vs_pval_adj, voom = voom_pval_adj,
           #edger = edger_pval_adj#, dsq = dsq_pval_adj
           contains('adj')
    ) %>%
    gather(method, adj_pval, -gene) %>%
    group_by(method) %>%
    summarise(
      fdr = sum(adj_pval[gene > n1] < 0.05)/sum(adj_pval < 0.05),
      tpr = mean(adj_pval[gene <= n1] < 0.05)
    ) %>%
    mutate(n = n,
           p1 = p1,
           setting = m,
           beta = b,
           re_sd = r,
           gene_effect_sd = g)%>%
    mutate(method = str_remove(method, '_pval_adj'))
  res_p <- res_p %>%
    inner_join(adj_res_p)
  # })
  res_p %>% print

  qplot(vsp_p$vsp_pval, vsp_p$vsp_pval_adj)
