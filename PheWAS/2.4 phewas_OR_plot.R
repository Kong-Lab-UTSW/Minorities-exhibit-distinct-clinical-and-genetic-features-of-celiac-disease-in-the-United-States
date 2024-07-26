library(tidyverse)

#load in phewas results
phecode_res_fin <- read_csv('~/Downloads/phecode_results_full_v2_2.csv')
#find OR and confidence interval from data
phecode_res_fin$OR_conf_int_1 <- exp(phecode_res_fin$conf_int_1)
phecode_res_fin$OR_conf_int_2 <- exp(phecode_res_fin$conf_int_2)
phecode_res_fin$OR <- exp(phecode_res_fin$beta_ind)
#remove celiac disease
phecode_res_fin <- phecode_res_fin %>% filter(phecode!=557.1) %>% arrange(OR)
phecode_res_fin$yAxis <- length(phecode_res_fin$description):1
#order by highest to lowest OR
phecode_res_fin$description <- fct_inorder(phecode_res_fin$description)

#select p<0.001
phecode_res_fin_2 <- phecode_res_fin %>% filter(p_value<=0.001)
#remove redundant conditions
phecode_res_fin_3 <- phecode_res_fin_2 %>% filter(!description %in% c('Chronic lymphocytic thyroiditis','Thyroiditis','Morbid obesity','Overweight, obesity and other hyperalimentation'
, 'Iron deficiency anemias, unspecified or not due to blood loss', 'Hypothyroidism NOS'))

#Odds ratio plot
p <- ggplot(phecode_res_fin_3, aes(x = OR, y = description,color = group)) 

p + geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = OR_conf_int_2, xmin = OR_conf_int_1), size = .5, height = 
                   .2, color = "gray50") +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),text=element_text(size=12, family="Arial"), axis.title.x = element_text(hjust=1)) +
  scale_x_continuous(breaks = seq(0.5,2,0.5) ) +
  scale_color_manual(values=c(unique(phecode_res_fin_2$color))) +
  ylab("Condition") +
  xlab("Odds ratio") +
  annotate(geom = "text", y =1.1, x = 1.4, label ="p < 0.001", size = 3.5, hjust = 0) +
  ggtitle("odds ratios for\nsignificant phewas conditions")

