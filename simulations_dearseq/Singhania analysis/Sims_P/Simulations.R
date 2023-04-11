library(readxl)
library(tidyverse)
library(janitor)
###effects data
dt <- read_excel("Effects_tb_vs_notb.xlsx")
dt<-clean_names(dt)

###overall data

db_counts_ens <- readRDS("db_counts_ens.rds")


cond<-db_counts_ens$tb_status%in%"Control"

ds_contr<-db_counts_ens[cond,]
ds_cases<-db_counts_ens[!cond,]
ds<-db_counts_ens

means_contr<-apply(ds_contr[,3:ncol(ds_contr)],2,mean,na.rm=T)
dss_contr<-apply(ds_contr[,3:ncol(ds_contr)],2,sd,na.rm=T)
p_contr<-apply(ds_contr[,3:ncol(ds_contr)],2,
                 function(x)sum(x==0, na.rm = T)/length(x))

means_cases<-apply(ds_cases[,3:ncol(ds_cases)],2,mean,na.rm=T)
dss_cases<-apply(ds_cases[,3:ncol(ds_cases)],2,sd,na.rm=T)
p_cases<-apply(ds_cases[,3:ncol(ds_cases)],2,
               function(x)sum(x==0, na.rm = T)/length(x))

means<-apply(ds[,3:ncol(ds)],2,mean,na.rm=T)
dss<-apply(ds[,3:ncol(ds)],2,sd,na.rm=T)
p<-apply(ds[,3:ncol(ds)],2,
               function(x)sum(x==0, na.rm = T)/length(x))

n<-100


###negbin var=mu+mu^2/size
#size=-mu^2/(mu-var)

do_gen_dat<-function(x){
  if(names(means_cases[x])%in%dt$x1){
  ##cases
  Y_c <- rnbinom(n,mu=means_cases[x],
          size = -means_cases[x]^2/(means_cases[x]-dss_cases[x]^2) )
  
  if(all(is.na(Y_c))){
    Y_c[is.na(Y_c)]<-0} else {
     NULL}
  
  U_c <- sample(c(0, 1), size = n, prob = c(p_cases[x], 1-p_cases[x]), replace = TRUE)
  v_c<-U_c*Y_c
  
  #controls
  
  Y_co <- rnbinom(n,mu=means_contr[x],
                 size = -means_contr[x]^2/(means_contr[x]-dss_contr[x]^2) )
  
  if(all(is.na(Y_co))){
    Y_co[is.na(Y_co)]<-0} else {
      NULL}
  
  
  U_co <- sample(c(0, 1), size = n, prob = c(p_contr[x], 1-p_contr[x]), replace = TRUE)
  v_co<-U_co*Y_co
  
  c(v_c,v_co)
  }else{
    
    Y_c <- rnbinom(n*2,mu=means[x],
                   size = -means[x]^2/(means[x]-dss[x]^2) )
    
    if(all(is.na(Y_c))){
      Y_c[is.na(Y_c)]<-0} else {
        NULL}
    
    U_c <- sample(c(0, 1), size = n*2, prob = c(p[x], 1-p[x]), replace = TRUE)
    v_c<-U_c*Y_c
    
    
    
  }
}


d_sim<-data.frame(sapply(1:length(means_cases),do_gen_dat))
colnames(d_sim)<-names(means_cases)
d_sim$group<-c(rep(1,n),
               rep(0,n))

# Wilcoxon

do_pval<-function(x){
  p_val<-wilcox.test(d_sim[1:n,x],d_sim[n+1:n*2,x])
  p_val$p.value  
}

pvals<-sapply(1:(ncol(d_sim)-1),do_pval)
names(pvals)<-colnames(d_sim)[1:(ncol(d_sim)-1)]

pv_s<-names(pvals)[pvals<0.05][complete.cases(names(pvals)[pvals<0.05])]

nomi<-colnames(d_sim)[colnames(d_sim)%in%dt$x1]

length(nomi[nomi%in%pv_s])/length(nomi)


