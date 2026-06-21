# Check significance of the recommendation using difference in survival outcomes
metric_significance <- function(df_patient, cutoff = NULL) {
  
  # If cutoff is provided, truncate time and censor events
  if (!is.null(cutoff)) {
    df_use <- df_patient %>%
      mutate(
        time_use  = pmin(time, cutoff),
        event_use = ifelse(time <= cutoff, event, 0)
      )
    time_var  <- "time_use"
    event_var <- "event_use"
    
    
    # Else, use full data
  } else {
    df_use <- df_patient
    time_var  <- "time"
    event_var <- "event"
  }
  
  # Log-rank test
  logrank <- survdiff(
    as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ followed")),
    data = df_use
  )
  
  # P-value
  p_value <- 1 - pchisq(logrank$chisq, df = length(logrank$n) - 1)
  
  # Message
  if (is.null(cutoff)) {
    message("Does following our recommendation improve survival (full follow-up)?")
  } else {
    message(paste0(
      "Does following our recommendation improve survival up to ",
      cutoff,
      " months?"
    ))
  }
  
  message(paste0("P-value: ", signif(p_value, 4)))
  
  invisible(p_value)
}

# Check Concordance Index
# This answers the following questions.
# Is the model’s internal ordering meaningful?
# Does this hold in short-term vs long-term horizons?
# Do patients estimated to benefit more from treatment tend to live longer?
metric_hte_concordance <- function(df_patient, cutoff = NULL) {
  
  # Truncate follow-up
  if (!is.null(cutoff)) {
    df_use <- df_patient %>%
      mutate(
        time_use  = pmin(time, cutoff),
        event_use = ifelse(time <= cutoff, event, 0)
      )
    time_var  <- "time_use"
    event_var <- "event_use"
    
    # Else, use full data
  } else {
    df_use <- df_patient
    time_var  <- "time"
    event_var <- "event"
  }
  
  # Treated group C-index
  c_treated <- concordance(
    as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ est_effect")),
    data = df_use %>% filter(treatment == 1)
  )$concordance
  
  # Control group C-index
  c_control <- concordance(
    as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ est_effect")),
    data = df_use %>% filter(treatment == 0)
  )$concordance
  
  # Messages
  if (is.null(cutoff)) {
    message("HTE concordance index (full follow-up):")
  } else {
    message(paste0(
      "HTE concordance index up to ",
      cutoff,
      " months:"
    ))
  }
  
  message(paste0("  Treated: ", round(c_treated, 3)))
  message(paste0("  Control: ", round(c_control, 3)))
  
  invisible(
    list(
      treated = c_treated,
      control = c_control
    )
  )
}

# Plot survival by group
metric_survival_plot <- function(df_patient, lines = c(NULL)) {
  
  par(mar = c(5.5, 6, 4, 2) + 0.1)
  
  df_patient$followed <- factor(df_patient$followed, levels = c("N", "Y"))
  
  # calculate survival
  km_fit <- survfit(
    Surv(time, event) ~ followed
    , data = df_patient
  )
  
  # get order of labels
  strata_labels <- names(km_fit$strata)
  n_strata <- length(strata_labels)
  
  plot(
    km_fit,
    col = 1:n_strata,
    lwd = 2,
    xlab = "Time (months)",
    ylab = "Survival probability",
    #main = paste("Survival –", HTE_METHOD, "-", FS_METHOD),
    #main = paste("Survival: PAM50 + CausalPFN"),
    #sub  = "Trained on METABRIC Luminal A, evaluated on METABRIC Luminal A",
    cex.lab = 1.8,     # axis labels
    cex.axis = 1.8,    # tick labels
    cex.main = 2,    # title
  )
  
  #mtext(
  #  "Trained on METABRIC Luminal A, evaluated on METABRIC Luminal A",
  #  side = 3,
  #  line = 0.8,    # distance below main title
  #  cex = 1.5
  #)
  
  # optional vertical lines
  if (!is.null(lines)) {
    
    y = 0.7
    for (i in lines) {
      
      # plot vertical line
      abline(v = i, lty = 2, lwd = 2)
      
      # add label
      text(
        x = i
        , y = y
        , labels = paste0(i, " months")
        , pos = 4
        , cex = 1.8
      )
      
      y <- y+0.07
    }
  }
  
  # add legend
  #legend(
  #  "bottomleft",
  #  legend = strata_labels,
  #  col = 1:n_strata,
  #  lwd = 2,
  #  title = "Followed our recommendation?",
  #  cex = 1.8,        # legend text
  #  pt.cex = 2,
  #  title.cex = 2
  #)
}

metric_coverage <- function(y_test, y_pred) {
  y_test <- factor(y_test, levels = c(0, 1))
  y_pred  <- factor(y_pred,  levels = c(0, 1))
  
  #fourfoldplot(as.table(confusionMatrix(y_test, y_pred)),
  #             color = c("yellow", "pink"),
  #             main = "Confusion Matrix")
  
  library(ggplot2)
  
  cm <- as.data.frame(confusionMatrix(data = y_pred, reference = y_test)$table)
  
  cm$Followed <- ifelse(
    (cm$Prediction == 0 & cm$Reference == 0) | (cm$Prediction == 1 & cm$Reference == 1),
    "Followed",
    "Not Followed"
  )
  
  print(cm)
  
  ggplot(cm, aes(x = Prediction, y = Reference, fill = Followed)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Freq), size = 6) +
    scale_fill_manual(values = c("Followed" = "cadetblue2", "Not Followed" = "orange")) +
    labs(x = "Recommended", y = "Actual", title = "Confusion Matrix", fill = "") +
    theme_minimal(base_size = 16) +
    theme(
      axis.text = element_text(size = 14),
      plot.title = element_text(size = 18, face = "bold")
    )
  
}