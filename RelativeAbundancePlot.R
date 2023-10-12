## libraries
library(tidyverse)
library(tidyr)
library(phyloseq)
library(ggplot2)
library(stringr)

##### Read in input data ####
seqtab <- read.csv("Pardis_input/seqtab.csv", header = T, row.names = 1)
taxa <- read.csv("Pardis_input/taxatab.csv", header = T, row.names = 1)
meta <- read.csv(".csv", header = T, row.names = 1)

# combine seqtab and taxa in correct order
# fake meta
meta <- as.data.frame(unlist(str_split(rownames(seqtab),"_"))[seq(1,5*length(rownames(seqtab)),5)]); rownames(meta) <- rownames(seqtab); colnames(meta) <- "sampno"
meta <- meta %>% cbind(repl=c(rep("A",49),rep("B",100),rep("C","55")))

newseqtab <- t(seqtab) %>% cbind(taxa[colnames(seqtab),])
newseqtab <- newseqtab %>% gather( sample, count, `1_S1_L001_R1_filt.fastq.qz`:`PK-60_S60_L001_R1_filt.fastq.qz`,factor_key=T)

plot_df <- newseqtab %>% filter(!is.na(Order)) %>%
  group_by(sample, Order) %>%
  # Add abundances within each order and sample
  summarize_at("count", sum) %>%
  # Get order proportions within each sample
  mutate(Percent = count / sum(count)*100) %>%
  cbind(meta[plot_df$sample,])

ggplot(plot_df, aes(x = sample, y = Percent, fill = Order))+
  geom_bar(stat = "identity")+ theme_classic()+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

av_df <- newseqtab %>% filter(!is.na(Order)) %>%
  cbind(meta[newseqtab[!is.na(newseqtab$Order),]$sample,])  %>%
  group_by(repl, Order) %>%
  summarize_at("count", sum) %>%
  # Get order proportions within each sample
  mutate(Percent = count / sum(count)*100)
  
ggplot(av_df,
       aes(x = repl, y = Percent, fill = Order))+
  geom_bar(stat = "identity", colour="black")+ theme_classic()+
  theme(axis.ticks.x = element_blank())

