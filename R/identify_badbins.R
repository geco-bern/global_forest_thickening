identify_badbins <- function(df, binwidth = 10){

  df_binned <- df |> 
    mutate(yearbin = floor(year / binwidth) * binwidth) %>%
    group_by(yearbin) %>%
    nest() |> 
    arrange(yearbin) |> 
    mutate(nobs = purrr::map_int(data, ~nrow(.)))
    # filter(nobs >= 10)

  # compare pairwise between bins
  pairs <- expand_grid(b1 = df_binned$yearbin, b2 = df_binned$yearbin) %>%
    filter(b1 < b2) %>%
    mutate(
      test = map2_dbl(b1, b2, ~ {
        x <- df_binned$data[df_binned$yearbin == .x][[1]]$logQMD
        y <- df_binned$data[df_binned$yearbin == .y][[1]]$logQMD
        
        # test measures the distance of the two medians, normalised by the inter-quartile range
        abs(median(x) - median(y)) / diff(quantile(x, probs = c(0.25, 0.75)))

      })
    )

  # count how often each bin is “different”
  bin_scores <- pairs %>%
    mutate(diff = test > 1.0) %>%
    pivot_longer(c(b1, b2), values_to = "bin") %>%
    group_by(bin) %>%
    summarise(diff_rate = mean(diff), .groups = "drop")

  # Remove bins at head (old years) that differ from the majority
  tmp <- bin_scores %>%
    filter(diff_rate > 2/3) %>%   # differs from >75% of bins
    pull(bin)

  remove_smaller_equal_yearbin <- ifelse(
    length(tmp) > 0, 
    max(tmp),
    -9999)

  out <- df_binned |> 
    mutate(badqmdbin = yearbin <= remove_smaller_equal_yearbin) |> 
    ungroup() |> 
    select(-yearbin, -nobs) |> 
    unnest(data)

  # # visualise bad bin removal
  # df_binned |> 
  #   unnest(data) |> 
  #   mutate(yearbin = as.factor(yearbin)) |> 
  #   ggplot(aes(yearbin, logQMD, fill = badqmdbin)) +
  #   geom_boxplot()

  return(out)

}