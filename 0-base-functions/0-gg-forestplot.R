##############################################
##############################################
# Documentation: gg_fp
# Usage: gg_fp(result_df, age, interaction_flag, binary_outcome)
# Description: forest plot for visualize single-mediator analysis outcomes


gg_fp <- function(result_df, modifier_flag=0, age, interaction_flag, binary_outcome) {
  df = result_df %>% filter(interaction == interaction_flag, age_group == age)
  df$use_modifier = modifier_flag
  
  df$label <- ifelse(df$use_modifier == 0, paste0(df$age_group, ": ", df$mediator, " → ", df$outcome), 
                     paste0(df$age_group, ": ", df$mediator, " → ", df$outcome, " (", df$modifier, ")"))
  
  long_df <- data.frame(
    label = rep(df$label, times = 3),
    effect_size = c(df$ACME_average, df$ADE_average, df$total_effect),
    lower_CI = c(
      df$ACME_average_lower_CI,
      df$ADE_average_lower_CI,
      df$total_effect_lower_CI
    ),
    upper_CI = c(
      df$ACME_average_upper_CI,
      df$ADE_average_upper_CI,
      df$total_effect_upper_CI
    ),
    effect_type = rep(c("ACME", "ADE", "Total Effect"), each = nrow(df))
  )
  long_df <- long_df[order(long_df$label, long_df$effect_type),]
  
  adjustment_factor <- 1
  long_df$y_position <- 4*adjustment_factor * (as.numeric(as.factor(long_df$label)) - 1)
  long_df$y_position <- long_df$y_position + adjustment_factor +
    ifelse(long_df$effect_type == "Total Effect",-adjustment_factor,
           ifelse(long_df$effect_type == "ADE", 0, adjustment_factor))
  
  y_breaks <-
    seq(3, max(long_df$y_position), by = 4 * adjustment_factor)
  
  unique_breaks <-
    unique(long_df$y_position[long_df$effect_type == "ADE"])
  
  plot1 <- ggplot(long_df, aes(y = y_position, x = effect_size)) +
    geom_point(aes(color = effect_type), shape = 15, size = 2) +
    geom_errorbarh(aes(
      xmin = lower_CI,
      xmax = upper_CI,
      color = effect_type
    ),
    height = 0.1) +
    geom_vline(xintercept = ifelse(binary_outcome==1, 1, 0), linetype = "longdash") +
    geom_hline(yintercept = max(long_df$y_position) + 4,
               color = "black") +
    geom_hline(yintercept = y_breaks,
               color = "grey",
               linetype = "dotted") +
    scale_y_continuous(name = "",
                       labels = unique(long_df$label),
                       breaks = unique_breaks) +
    scale_color_manual(values = c(
      "ACME" = "#FF6347",
      "ADE" = "#3CB371",
      "Total Effect" = "#1E90FF"
    )) +
    xlab(ifelse(interaction_flag == 1, "interacted", "non-interacted")) +
    ylab(" ") +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.y = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.text.y = element_text(size = 12, colour = "black"),
      axis.ticks.y = element_blank(),
      axis.text.x.bottom = element_text(size = 12, colour = "black"),
      axis.title.x = element_text(size = 12, colour = "black"),
      legend.position = c(1, 0.06),
      legend.background = element_rect(fill='transparent')
    ) +
    ggtitle("")
  
  table_base <- ggplot(long_df, aes(y = y_position)) +
    ylab(NULL) + xlab(" ") +
    theme(
      plot.title = element_text(size = 12, 
                                hjust = 6,
                                vjust = -20),
      axis.text.x = element_text(color = "white", size = 10),
      ## This is used to help with alignment
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none",
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_blank()
    )
  
  ## Effect size table
  tab1 <- table_base +
    geom_text(aes(
      y = y_position,
      x = 1,
      label = sprintf("%0.4f", effect_size)
    ), size = 3.5) +
    ggtitle("Effect size") +
    geom_hline(yintercept = max(long_df$y_position) + 4,
               color = "black") +
    theme(plot.title = element_text(
      size = 12,
      hjust = 0.5,
      vjust = -10
    ))
  
  ## 95% CI table
  tab2 <- table_base +
    geom_text(aes(
      y = y_position,
      x = 1,
      label = sprintf("(%0.4f, %0.4f)", lower_CI, upper_CI)
    ), size = 3.5) +
    ggtitle("95% CI") +
    geom_hline(yintercept = max(long_df$y_position) + 4,
               color = "black") +
    theme(plot.title = element_text(
      size = 12,
      hjust = 0.5,
      vjust = -10
    ))
  
  grid.arrange(plot1, tab1, tab2, widths = c(4.5, 1, 1))
  
  # cowplot::plot_grid(plot1, tab1, tab2,
  #   nrow = 1,
  #   align = "h",
  #   rel_widths = c(4,1,1)
  # )
}
