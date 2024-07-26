library(tidyverse)
library(plotly)
library(scales)
hla = read_tsv('MHC_vars_anno_full.vep.csv')
hla <- hla[!(is.na(hla$beta)),]
hla <- hla[!(is.na(hla$SYMBOL)),]
hla$p_value <- as.double(hla$p_value)
hla$log10pval <- -log10(hla$p_value)

p <- ggplot(hla, aes(x=POS, y = log10pval, color=SYMBOL)) +
  geom_point(size=0.5) +
  xlab("Chr6") +
  ggtitle('MHC GWAS') +
  geom_hline(yintercept=8, linetype="dashed", color = "red")+
  theme(legend.position="none",axis.title.x=element_text(hjust=0))+
  scale_x_continuous(labels = label_comma())
p


