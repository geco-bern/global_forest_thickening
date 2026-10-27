# Dataset-specific SI figures, using the saved harmonised compilation before
# the common analysis filters. Run from the project root with:
# Rscript analysis/si_figures/01_plot_dataset_histograms.R
# Optional argument: path to an alternative harmonised df_all.rds snapshot.
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

root <- here::here()
source(file.path(root, "R/harmonise_dataset_names.R"))
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args)) args[[1]] else file.path(root, "data/inputs/df_all.rds")
figure_dir <- file.path(root, "manuscript/figures/dataset_histograms")
tex_dir <- file.path(root, "manuscript/si/dataset_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tex_dir, recursive = TRUE, showWarnings = FALSE)

# Order matches the constituent-dataset descriptions in the SI. These are
# canonical identifiers shared with Table S1 and the processing scripts.
datasets <- c(
  "nfi_spain", "nfi_sweeden", "fia_us", "nfi_switzerland", "nfi_norway",
  "efm_swi", "uholka", "fep_gre_stand", "inrae_lessem", "bnp", "czu",
  "fvabw", "iberbas", "unitbv", "lwf", "nbw", "nfr_swi",
  "nwfva", "tuzvo", "ul", "unito", "urk", "wuls",
  "luquillo", "scbi", "palanan", "serc", "wytham", "pasoh", "mudumalai",
  "forestplots", "aus_plots", "rainfor_plots"
)
d <- readRDS(input) |>
  mutate(dataset = harmonise_dataset_names(dataset))
required <- c("dataset", "plotID", "year", "management_since_census1_yrs")
stopifnot(all(required %in% names(d)), is.numeric(d$year),
          is.numeric(d$management_since_census1_yrs))
if (!setequal(unique(d$dataset), datasets)) {
  stop("Dataset identifiers differ from the SI mapping; review the mapping first.")
}
if (anyNA(d$plotID)) stop("Missing plot identifiers require resolution before counting plots.")

# An interval is a plot attribute, not one observation per census. Do not
# silently average inconsistent values across censuses. NA/NaN are missing.
one_interval <- function(x) {
  values <- unique(x[is.finite(x)])
  if (length(values) > 1L) {
    stop("Management interval varies within a plot; resolve before plotting.")
  }
  if (length(values)) values[[1]] else NA_real_
}
span <- function(year) {
  year <- year[is.finite(year)]
  if (length(year)) diff(range(year)) else NA_real_
}
plots <- d |>
  group_by(dataset, plotID) |>
  summarise(
    n_census = n_distinct(year[is.finite(year)]),
    n_undated_records = sum(!is.finite(year)),
    record_length = span(year),
    management_years = one_interval(management_since_census1_yrs),
    .groups = "drop"
  )

fmt <- function(x) format(x, big.mark = ",", trim = TRUE, scientific = FALSE)
integer_breaks <- function(limits) {
  b <- pretty(limits, n = 4)
  b[abs(b - round(b)) < 1e-8 & b >= 0]
}
# Use integer-year bins and a small number of readable ticks. A constant
# distribution remains a histogram with one occupied one-year bin.
histogram <- function(x, title, subtitle, xlabel, colour) {
  stopifnot(length(x) > 0L, all(is.finite(x)))
  width <- max(1, ceiling(diff(range(x)) / 20))
  lower <- floor(min(x) / width) * width - width / 2
  upper <- ceiling(max(x) / width) * width + width / 2
  limits <- c(lower, upper)
  if (length(unique(x)) == 1L) limits <- x[[1]] + c(-1.5, 1.5)
  p <- ggplot(data.frame(value = x), aes(value)) +
    geom_histogram(binwidth = width, boundary = -width / 2,
                   fill = colour, colour = "white", linewidth = 0.2) +
    scale_x_continuous(limits = limits, breaks = function(lim) {
      if (length(unique(x)) == 1L) return(x[[1]] + c(-1, 0, 1))
      pretty(lim, n = 4)
    }, expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(breaks = integer_breaks, labels = fmt,
                       expand = expansion(mult = c(0, 0.08))) +
    labs(title = title, subtitle = subtitle, x = xlabel, y = "Number of plots") +
    theme_minimal(base_size = 9, base_family = "sans") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "#e6e6e6", linewidth = 0.25),
          axis.text = element_text(colour = "#333333"),
          plot.title = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(size = 8, margin = margin(b = 7)),
          axis.title = element_text(size = 9),
          plot.margin = margin(5, 7, 5, 5))
  # Verify that axis limits/bin boundaries have not discarded any plots.
  stopifnot(sum(ggplot_build(p)$data[[1]]$count) == length(x))
  p
}

summaries <- list()
for (id in datasets) {
  p <- filter(plots, dataset == id)
  lengths <- p$record_length[is.finite(p$record_length)]
  # 999 is a categorical marker, not an observed elapsed time. Preserve it
  # in the audit counts and captions, but exclude it from the numeric axis.
  marker <- is.finite(p$management_years) & p$management_years == 999
  management <- p$management_years[is.finite(p$management_years) & !marker]
  missing_management <- sum(!is.finite(p$management_years))
  subtitle <- paste0(fmt(nrow(p)), " plots; ", fmt(sum(p$n_census)), " plots \u00d7 census")
  a <- histogram(lengths, "A  Record length", subtitle,
                 "Last minus first census (years)", "#376b7e")
  if (length(management)) {
    b <- histogram(management, "B  Management history",
      paste0(fmt(length(management)), " plots with numeric intervals"),
      "Years since last management\nat first census", "#a66b36")
    fig <- a + b + plot_layout(nrow = 1)
  } else {
    fig <- a
  }
  ggsave(file.path(figure_dir, paste0(id, ".pdf")), fig,
         width = if (length(management)) 6.3 else 3.15,
         height = 2.5, units = "in", device = grDevices::cairo_pdf, bg = "white")

  tex_id <- gsub("_", "\\_", id, fixed = TRUE)
  caption <- paste0("Record-length distribution for \\texttt{", tex_id,
    "}, before the common analysis filters. Each plot contributes once; ",
    "plots $\\times$ census counts distinct dated plot--year observations.")
  single <- sum(p$n_census == 1L)
  if (single) caption <- paste0(caption, " The ", fmt(single),
    " plots with one dated census have a record length of zero.")
  undated <- sum(p$n_undated_records)
  no_year <- sum(p$n_census == 0L)
  if (undated) caption <- paste0(caption, " There are ", fmt(undated),
    " undated records, excluded from the census count and record-length calculation; ",
    fmt(no_year), " plots have no dated census and are absent from panel A but included in the plot count.")
  if (length(management)) {
    caption <- paste0(caption, " Panel B shows the stored interval since last management ",
      "at the first census, once per plot; ", fmt(missing_management),
      " plots have missing intervals")
    if (any(marker)) caption <- paste0(caption, " and ", fmt(sum(marker)),
      " carry the non-numeric-history marker 999")
    caption <- paste0(caption, " and are omitted from panel B.")
    if (any(management < 0)) caption <- paste0(caption,
      " Negative values are retained and indicate a recorded intervention after the first census.")
    if (id == "nfi_switzerland") caption <- paste0(caption,
      " Swiss NFI values follow the special screening calculation described in the text.")
  } else if (any(marker)) {
    caption <- paste0(caption, " No numeric management histogram is shown: ",
      fmt(sum(marker)), " plots carry marker 999 (long absence of intervention/never managed), ",
      "which is not treated as a measured interval; ", fmt(missing_management),
      " plots have missing intervals.")
  } else {
    caption <- paste0(caption, " Management intervals are unavailable; panel B is omitted.")
  }
  tex <- c("% Generated by analysis/si_figures/01_plot_dataset_histograms.R; do not edit manually.",
    "\\begin{figure}[H]", "\\centering",
    paste0("\\includegraphics[width=", if (length(management)) "\\textwidth" else "0.5\\textwidth",
           "]{../figures/dataset_histograms/", id, ".pdf}"),
    paste0("\\caption{", caption, "}"),
    paste0("\\label{fig:dataset-", gsub("_", "-", id), "}"), "\\end{figure}")
  writeLines(tex, file.path(tex_dir, paste0(id, ".tex")))
  summaries[[id]] <- data.frame(dataset = id, n_plots = nrow(p),
    n_plot_census = sum(p$n_census), n_single_census = single,
    n_undated_records = undated, n_plots_without_year = no_year,
    n_management_numeric = length(management), n_management_missing = missing_management,
    n_management_marker999 = sum(marker), n_management_negative = sum(management < 0),
    management_panel = length(management) > 0)
}
audit <- bind_rows(summaries)
write.csv(audit, file.path(figure_dir, "summary.csv"), row.names = FALSE)
writeLines(c(
  "Generated by analysis/si_figures/01_plot_dataset_histograms.R",
  paste("Input:", normalizePath(input)),
  paste("Input MD5:", unname(tools::md5sum(input))),
  paste("Input modified (UTC):", format(file.info(input)$mtime, tz = "UTC", usetz = TRUE)),
  "Scope: saved harmonised compilation, before common analysis filters.",
  "Dataset labels harmonised with R/harmonise_dataset_names.R; BCI and SCBI remain distinct.",
  "Plot key: dataset + plotID. Census count: distinct finite calendar years per plot.",
  "Record length: max(year) - min(year); single-census plots have length zero.",
  "Management: one unique finite interval per plot; 999 excluded as a marker.",
  "Missing years and negative intervals are reported explicitly in captions.",
  capture.output(sessionInfo())
), file.path(figure_dir, "provenance.txt"))
print(audit, row.names = FALSE)
