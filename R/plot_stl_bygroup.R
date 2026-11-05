plot_stl_bygroup <- function(df, mod, name, years = c(1985, 2000, 2015), plot_legend = FALSE, interactions = FALSE){
  
  # predict for visualisation
  preddata <- ggpredict(
    mod, 
    terms = c("logQMD", paste0("year [", years[1], ", ", years[2], ", ", years[3], "]")), 
    full.data = TRUE) |>  # full.data = TRUE to include random effects. full.data = FALSE ignores group-specific random effects and gives predictions for an "average" group.
    as_tibble()
  
  # significance of year in subtitle
  out <- summary(mod)
  
  if (interactions){
    
    df_caption <- tibble(
      variable = c("year", "qmd_year"),
      estimate = out$coefficients[c("scale(year)", "scale(logQMD):scale(year)"), "Estimate"],
      pvalue = out$coefficients[c("scale(year)", "scale(logQMD):scale(year)"), "Pr(>|t|)"]
      ) |> 
      mutate(
        estimate_lab = round(estimate, 3),
        pval_lab = ifelse(pvalue < 0.001, "< 0.001 ***", paste0(round(pvalue, 3)))
      )
    
    subtitle <- bquote(
      italic(n) == .(as.character(nobs(mod))) ~ ~ ~
        italic(p)(yr) == .(df_caption$pval_lab[1]) ~ ~ ~
        italic(p)(QMD %*% yr) == .(df_caption$pval_lab[2])
    )
    
  } else {
    
    df_caption <- tibble(
      variable = c("year"),
      estimate = out$coefficients[c("scale(year)"), "Estimate"],
      pvalue = out$coefficients[c("scale(year)"), "Pr(>|t|)"]
    ) |> 
      mutate(
        estimate_lab = round(estimate, 3),
        pval_lab = ifelse(pvalue < 0.001, "< 0.001 ***", paste0(round(pvalue, 3)))
      )
    
    subtitle <- bquote(
      italic(n) == .(as.character(nobs(mod))) ~ ~ ~
        italic(p)[yr] == .(df_caption$pval_lab[1])
    )
  }
  
  # panel for final plot
  ggplot() + 
    geom_point(
      data = df, 
      aes(x = logQMD, y = logDensity),
      alpha = 0.5, 
      size = 1, 
      col = "darkgrey", 
      shape = 16, 
      inherit.aes = FALSE
    ) + 
    geom_ribbon(
      data = preddata, 
      aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group), 
      alpha = 0.2, 
      show.legend = TRUE
    ) + 
    geom_smooth(
      data = preddata, 
      aes(x = x, y = predicted, color = group), 
      method = "lm", 
      fullrange = FALSE, 
      size = .6, 
      se = FALSE
    ) +
    labs(
      x = expression(ln(QMD)),
      y = expression(ln(italic(N))), 
      title = name,
      subtitle = subtitle,
      color  = "Year"
    ) +
    scale_color_manual("Year",
                       breaks = c(as.character(years[1]), as.character(years[2]), as.character(years[3])), 
                       values = c(viridis(3)[2], viridis(3)[3], viridis(3)[1])) +
    scale_fill_manual("Year",
                      breaks = c(as.character(years[1]), as.character(years[2]), as.character(years[3])), 
                      values = c(viridis(3)[2], viridis(3)[3], viridis(3)[1])) +
    theme_classic() +  
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 12), 
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 10), 
      legend.title = element_text(size = 10),
      plot.title = element_text(size = 12),
      plot.subtitle = element_text(size = 10),
      legend.key = element_rect(fill = NA, color = NA),
      legend.position = ifelse(plot_legend, "bottom", "none"),
      plot.caption = element_text(vjust = -1),
      plot.title.position = "plot"
    ) +
  scale_x_continuous(limits = c(2,4.5), breaks = seq(2,4,1)) +
  scale_y_continuous(limits = c(2.9,9.3), breaks = seq(4,8,2))
  
}
