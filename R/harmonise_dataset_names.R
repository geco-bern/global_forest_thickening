# Canonical dataset identifiers for processing, tables and plot labels.
# Keep the original raw-data filenames and distinct BCI/SCBI sites unchanged.
harmonise_dataset_names <- function(dataset) {
  aliases <- c(
    greece_stand = "fep_gre_stand",
    france = "inrae_lessem",
    forst = "fvabw",
    incds = "unitbv",
    lwf_tree = "lwf",
    npvbw = "nbw",
    nwfva_tree = "nwfva",
    tuzvo_tree = "tuzvo",
    ul_tree = "ul",
    palanam = "palanan",
    df_forestplots = "forestplots",
    df_rainfor = "rainfor",
    rainfor_plots = "rainfor"
  )
  dataset <- as.character(dataset)
  replacement <- unname(aliases[dataset])
  replace(dataset, !is.na(replacement), replacement[!is.na(replacement)])
}
