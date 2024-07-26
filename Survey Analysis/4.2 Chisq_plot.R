setwd("E:/AllofUS_project/AOU survey data/Survey analysis")

library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)

#All survey data
df <- read_csv('20240202_CeD_survey.csv')

colnames(df2)

df <- as.data.table(df)


class(df$question_concept_id)

unique(df$question_concept_id)

df$Padj[df$Padj == ""] <- NA

#find unique questions
df[, .(total_unique_count = uniqueN(question_concept_id),
                           count_unique_Padj_less_than_00001 = sum(!duplicated(question_concept_id) & as.numeric(format(Padj, scientific = FALSE)) < 0.000001)), by = .(survey)]

df[, .(total_unique_count = uniqueN(question_concept_id),
  count_unique_Padj_less_than_00001 = sum(!duplicated(question_concept_id) & !is.na(as.numeric(format(Padj, scientific = FALSE))) & as.numeric(format(Padj, scientific = FALSE)) < 0.000001)),
  by = .(survey)]

#filter for p<0.05 and p<1e-25 
P05 <- df[, .(total_unique_count = uniqueN(question_concept_id),
       count_unique_Padj_less_than_00001 = sum(!duplicated(question_concept_id) & !is.na(as.numeric(format(Padj, scientific = FALSE))) & as.numeric(format(Padj, scientific = FALSE)) < 0.05)),
   by = .(survey)]

P_25_2 <- df[, .(total_unique_count = uniqueN(question_concept_id),
              count_unique_Padj_less_than_00001 = sum(!duplicated(question_concept_id) & 
              !is.na(as.numeric(format(Padj, scientific = FALSE))) & 
                -log(as.numeric(format(Padj, scientific = FALSE))) >25)),
          by = .(survey)]

write_csv(P05, "padj04_variable.csv")
write_csv(P_25_2, "padj_25_variable.csv")

table(is.na(df$question_concept_id))

table(is.na(df$Padj))

#filter out skipped responses and failed chi-square tests
df2 <- df %>% 
  filter(df$celiac >= 20 & !is.na(Padj) & 
           as.numeric(format(Padj, scientific = FALSE)) < 0.000001 &
           df$answer != "PMI: Skip" &
           df$answer != "PMI: Dont Know" &
           df$healthy >= 20) 


df3 <- df2[, .(unique_answer_count = uniqueN(answer_concept_id)), by = question_concept_id]

df4 <- df2[, .(min_OR = min(OR, na.rm = TRUE), max_OR = max(OR, na.rm = TRUE)), by = question_concept_id]

#Odds ratio; get the max enriched OR or the min protective OR for each question,as well as the most enriched/protected response 

df4$min_OR[df4$min_OR > 1] <- NA

df4$max_OR[df4$max_OR < 1] <- NA

df5 <- df4[, comparison_result := ifelse(-log2(min_OR) > log2(max_OR), 'min', 'max'), by = question_concept_id]


df6 <- df5[, comparison_result := ifelse(is.na(comparison_result), ifelse(is.na(min_OR), 'max', 'min'), comparison_result), by = question_concept_id]

df7<- df6[, final_OR := ifelse(comparison_result == 'max', max_OR, min_OR), by = question_concept_id]

Padj <- df2 %>%
  select(question_concept_id, question, Padj) %>%
  distinct()

df8 <- left_join(df7, Padj, by = "question_concept_id")

answer <- df2 %>% 
  select(question_concept_id, question, answer, OR, survey) %>%
  distinct()

df9 <- left_join(df8, answer, by = c("question_concept_id", "final_OR" ="OR","question"="question"))



## cap OR values to make plot look better

df$log_OR_capped <- pmin(pmax(log(df$OR), -2.5), 2.5)


#df9
df9 <- df9 %>% 
  filter(-log10(Padj)<200)

write_csv(df9, "survey_plotting.csv")


###Final plot for use in the figure
unique_survey_values <- unique(df9$survey)


ggplot(df9, aes(x = log_OR_capped, y = -log10(Padj),
                fill = survey)) +
  geom_point(shape = 21, size = 3) +
  geom_hline(yintercept = 25, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = log(2), linetype = "dashed", color = "green") +
  geom_vline(xintercept = -log(2), linetype = "dashed", color = "green") +
  labs(x = "Log Odds Ratio (log OR)", y = "-log10(Padj)", title = "Scatter Plot of log OR vs -log10(Padj)") +
  coord_cartesian(ylim = c(0, 200)) +
  scale_fill_manual(values = rainbow(length(unique_survey_values))) +
  guides(fill = guide_legend(title = "Survey"))+
  theme_minimal()+
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )


