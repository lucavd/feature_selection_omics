

library(readxl)
library(janitor)

###effects data
dt <- read_excel("Effects_tb_vs_notb.xlsx")
dt<-clean_names(dt)


set.seed(666)
p<-0.1 ###This is the baseline TB prob
pg<-0.1 ####probability to obserbe the gene
b0<-log(p/(1-p)) ### intercept baseline log odds
n<-1000


b1 = dt$log2_fold_change[x] ## gene effects along db

x1 = rbinom(n,1,pg)


z = b0 + b1*x1        # linear combination 
pr = 1/(1+exp(-z))         # pass through an inv-logit function
y = rbinom(1000,1,pr)      # bernoulli response variable TB yes no




####noise variables

#continous
z1<- runif(n, 0:100)

z2<- rnorm(n)

###binari
zk<- rbinom(n,1,0.5)