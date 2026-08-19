#' Run plot_cov with the specified arguments

library(argparser, quietly = TRUE)
library(here)
library(rmarkdown)

here::i_am("scripts/run_plot.R")

arg_parser <- argparser::arg_parser("Plot coverage for a given amplicon") |>
  argparser::add_argument("out_base", help = "Output directory") |>
  argparser::add_argument("samplesheet", help = "Path to sample sheet for run") |>
  argparser::add_argument("threshold", help = "Minimum coverage to pass QC") |>
  argparser::add_argument("maxdepth", help = "Maximum depth to plot") |>
  argparser::add_argument("--region_file",
                          help = "Path to the region file to be used (optional, only for a single reference file)",
                          default = "")

args <- argparser::parse_args(arg_parser)

rmarkdown::render(here::here("scripts", "plot_cov.Rmd"),
                  output_format = "html_document",
                  output_file = paste0(args$out_base, "/plot_cov.html"),
                  params = list("threshold" = args$threshold,
                                "maxDepth" = args$maxdepth,
                                "path" = args$out_base,
                                "samplesheet" = args$samplesheet,
                                "region_file" = args$region_file))