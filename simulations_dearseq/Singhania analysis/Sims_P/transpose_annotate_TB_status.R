library(tidyverse)

# get sample type

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

group_London_TB <- which(infos_df$Status=="ActiveTB") # 21
group_London_LTBI <- which(infos_df$Status=="LTBI") # 21
group_London_Control <- which(infos_df$Status=="Control") # 12

Y <- rep(NA,length(infos_df$Status))
Y[which(infos_df$Status=="ActiveTB")] <- rep(1,length(which(infos_df$Status=="ActiveTB")))
Y[which(infos_df$Status=="LTBI")] <- rep(2,length(which(infos_df$Status=="LTBI")))
Y[which(infos_df$Status=="Control")] <- rep(0,length(which(infos_df$Status=="Control")))

# import data

db_counts_raw <- readxl::read_excel(here::here(
  'GSE107991_Raw_counts_Berry_London.xlsx')) |> 
  select(-Gene_biotype)


# transpose with ENSEMBL ID
# I Gene_names non sono unici, teniamo ENSEMBL ID
db_counts_ens <- db_counts_raw |> 
  select(-Gene_name) |> 
  pivot_longer(-Genes) |> 
  pivot_wider(names_from = "Genes", values_from = "value") |> 
  
  mutate(tb_status = Y) |> 
  relocate(tb_status, .after = name) |> 
  filter(tb_status != "2") |> 
  mutate(tb_status = ifelse(tb_status == 0, "Control", "TB"))

# transpose with Gene Name
# db_counts_gene <- db_counts_raw |> 
#   select(-Genes) |> 
#   pivot_longer(-Gene_name) |> 
#   pivot_wider(names_from = "Gene_name", values_from = "value")
