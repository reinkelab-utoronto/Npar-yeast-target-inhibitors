#!/usr/bin/env Rscript
# Install the CRAN packages needed to reproduce the figures.
# Run once:  Rscript scripts/00_install_packages.R
pkgs <- c("readxl", "tidyverse", "ggplot2", "patchwork",
          "ggnewscale", "ggrepel", "multcomp")
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}
# Figs. 1 and 2 need ggplot2 >= 4.0 (position_jitterdodge() dodging by group)
# and patchwork >= 1.3 (free(type = "label", side = "l")).
stopifnot(packageVersion("ggplot2") >= "4.0.0",
          packageVersion("patchwork") >= "1.3.0")
message("All packages installed.")
