library(dplyr)
library(tidyr)

# No swaps
perm_clean <- results_df %>%
  mutate(Test = "Permutation",
         Scenario = "No swaps")

# 12 swaps
perm_12 <- results_df_12 %>%
  mutate(Test = "Permutation",
         Scenario = "12 swaps")

# 16 swaps
perm_16 <- results_df_16 %>%
  mutate(Test = "Permutation",
         Scenario = "16 swaps")

perm_all <- bind_rows(perm_clean, perm_12, perm_16)


df_all   # contains K-type and H-type

df_all_tests <- bind_rows(df_all, perm_all)


df_all_tests <- df_all_tests %>%
  mutate(PValue_adj = ifelse(PValue == 0, 1e-16, PValue),
         logP = -log10(PValue_adj))

df_all_tests$Scenario <- factor(df_all_tests$Scenario,
                                levels = c("No swaps",
                                           "12 swaps",
                                           "16 swaps"))

df_all_tests$Test <- factor(df_all_tests$Test,
                            levels = c("Permutation",
                                       "K-type",
                                       "H-type"))



p_all_tests <- ggplot(df_all_tests,
                      aes(x = Model,
                          y = logP,
                          fill = Test)) +
  
  geom_col(position = position_dodge(width = 0.75),
           width = 0.65,
           color = "black",
           linewidth = 0.4) +
  
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             linewidth = 0.8,
             color = "blue") +
  
  scale_fill_manual(values = c("Permutation" = "grey20",
                               "K-type"     = "grey50",
                               "H-type"     = "grey80")) +
  
  facet_wrap(~ Scenario, nrow = 1) +
  
  labs(
    x = "",
    y = expression(-log[10](italic(p_value))),
    fill = ""
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid = element_blank()
  )

print(p_all_tests)
