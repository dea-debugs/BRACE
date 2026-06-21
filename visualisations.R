# EDA - visualisation


set.seed(123)

# Get data
in_path <- "C:/Users/parki/OneDrive/Desktop/Honours/Datasets/brca_metabric/"
df <- read.csv(paste0(in_path, "dataset_genes.csv"))

#path = "C:/Users/parki/OneDrive/Desktop/Honours/Datasets/GSE20685/"
#df <- read.csv(paste0(path, "dataset.csv"))

# Option 2
library(survival)
library(survminer)

km_fit <- survfit(Surv(time, event) ~ 1, data = df)

ggsurvplot(
  km_fit,
  data = df,
  conf.int = TRUE,
  risk.table = TRUE,     # shows how many are still at risk
  xlab = "Time",
  ylab = "Survival probability",
  ggtheme = theme_minimal(base_size=16),
  risk.table.font.size = 5
)

# Option 3
df$id <- seq_len(nrow(df))

ggplot(df, aes(y = id, x = time, colour = factor(event))) +
  geom_point() +
  scale_colour_manual(
    values = c("0" = "grey70", "1" = "black"),
    labels = c("Censored", "Event"),
    name = "Status"
  ) +
  labs(
    x = "Time",
    y = "Individual",
    title = "Individual Follow-up Times"
  ) +
  theme_minimal(base_size=16) +
  theme(axis.text.y = element_blank())

# Treatment histogram
df$treatment <- factor(
  df$treatment,
  levels = c(0, 1),
  labels = c("Untreated", "Treated")
)

ggplot(df, aes(x = treatment)) +
  geom_bar(fill = "grey50") +
  geom_text(
    stat = "count",
    aes(label = ..count..),
    vjust = -0.5,
    size = 10
  ) +
  labs(
    x = "Treatment group",
    y = "Count",
    title = "METABRIC Chemotherapy Treatment Assignment"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    plot.title = element_text(size = 24, face = "bold"),
    axis.title = element_text(size = 24),
    axis.text  = element_text(size = 20)
  ) +
  ylim(0, max(table(df$treatment)) * 1.15)

# Survival histogram
ggplot(df, aes(x = time, fill = factor(event))) +
  geom_histogram(
    bins = 30,
    alpha = 0.5,
    position = "identity"
  ) +
  geom_vline(
    xintercept = 150,
    linetype = "dashed",
    linewidth = 1
  ) +
  scale_fill_manual(
    values = c("0" = "steelblue", "1" = "firebrick"),
    labels = c("Alive / Censored", "Death"),
    name = "Status"
  ) +
  labs(
    x = "Time (months)",
    y = "Count",
    title = "Distribution of Follow-up Times by Survival Status"
  ) +
  theme_minimal(base_size = 16)
