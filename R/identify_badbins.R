identify_badbins <- function(df, dataset_name, binwidth = 10){

  message(unique(df$dataset))

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
        
        # # testing differences of logQMD distributions in pairs of year bins
        # ks.test(x, y)$p.value

        # test measures the distance of the two medians, normalised by the inter-quartile range
        abs(median(x) - median(y)) / diff(quantile(x, probs = c(0.45, 0.65)))

      })
    )

  if (nrow(pairs) == 1 && pairs$test > 1.0){
    # drop all data -- excessive QMD shift between the two available bins
    return(
      list(
        df = df_binned |> 
          unnest(data) |> 
          mutate(badqmdbin = TRUE) |> 
          ungroup(),
        gg = NA
      )
    )
    
  } else {

    # count how often each bin is “different”
    bin_scores <- pairs %>%
      mutate(diff = test > 1.0) %>%   # for normalised median comparisons
      # mutate(diff = test < 1e-6) %>%  # for ks test
      pivot_longer(c(b1, b2), values_to = "bin") %>%
      group_by(bin) %>%
      summarise(diff_rate = mean(diff), .groups = "drop") |> 

      # cut to bins that
      # keep only low == TRUE
      # split into consecutive blocks
      arrange(bin) %>%
      mutate(
        low = diff_rate <= 2/3,
        grp = cumsum(lag(!low, default = TRUE))
      )

    # keep the longest block
    longest_group <- bin_scores %>%
      filter(low) %>%
      group_by(grp) %>%
      summarise(len = n()) |> 
      arrange(desc(len)) |> 
      ungroup() |> 
      slice(1) |> 
      pull(grp)

    if (length(longest_group) > 0){

      keep_bins <- bin_scores |> 
        filter(grp == longest_group & low) |> 
        pull(bin)

      df <- df_binned |> 
        unnest(data) |> 
        mutate(badqmdbin = !(yearbin %in% keep_bins)) |> 
        ungroup()

      # visualise bad bin removal
      gg <- df |> 
        mutate(yearbin = as.factor(yearbin)) |> 
        ggplot(aes(yearbin, logQMD, fill = badqmdbin)) +
        geom_boxplot() +
        labs(
          x = "10-year bin",
          title = dataset_name
        ) +
        scale_fill_manual(values = c("grey40", "tomato")) +
        theme_minimal() +
        theme(
          legend.position = "none"
        )

    } else {
      # all are bad bins
      df <- df_binned |> 
        unnest(data) |> 
        mutate(badqmdbin = TRUE) |> 
        ungroup()

      gg <- df |> 
        mutate(yearbin = as.factor(yearbin)) |> 
        ggplot(aes(yearbin, logQMD, fill = "tomato")) +
        geom_boxplot() +
        labs(
          x = "10-year bin",
          title = dataset_name
        ) +
        theme_minimal() +
        theme(
          legend.position = "none"
        )

    }

  }

  return(
    list(
      df = df,
      gg = gg
    )
  )

}