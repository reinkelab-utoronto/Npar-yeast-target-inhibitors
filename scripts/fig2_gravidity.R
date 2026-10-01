# =============================================================
# JIP Figure - Gravidity (2 panels)
# Full-width (190 mm) grouped bar graphs of % gravidity per drug
# across 3 concentrations (10, 40, 100 uM), 3 reps each.
#   Panel A: infected animals   (control = Infected+DMSO)
#   Panel B: uninfected animals (control = DMSO)
# Same visual style as the infection figure: shared colours,
# replicate points, Dunnett stats, alternating green column bands,
# and a rotated drug-name strip under each panel.
#
# FIX: Panel B's first column (alignment placeholder for Panel A's
# Uninfected+DMSO) previously rendered a stray "." mark. The
# placeholder is now a plain empty string "" with no label drawn,
# so the slot holds its position but renders nothing.
# =============================================================

# ---- Packages ----
# Load multcomp/MASS BEFORE tidyverse so dplyr's select()/filter() win masking.
library(multcomp)
library(readxl)
library(tidyverse)
library(patchwork)

# ---- File paths ----
infected_file   <- "gravidity infected.xlsx"
uninfected_file <- "gravidity uninfected.xlsx"

# ---- Alignment placeholder (empty slot in Panel B under Panel A's
#      Uninfected+DMSO). Plain empty string = nothing to render. ----
placeholder <- ""

# ---- Shared style constants (match the infection figure) ----
conc_levels <- c("10 \u00b5M", "40 \u00b5M", "100 \u00b5M")
conc_fill   <- c("10 \u00b5M" = "#bdd7e7",
                 "40 \u00b5M" = "#6baed6",
                 "100 \u00b5M" = "#2171b5")
band_fill   <- "#e6f4ea"   # light-green alternating column dividers

base_theme <- theme_classic(base_size = 8) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    axis.text.y  = element_text(size = 7, margin = margin(r = 1)),
    axis.title   = element_text(size = 8),
    axis.title.y = element_text(size = 8, vjust = -0.2),
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 7),
    legend.key.size = unit(3, "mm"),
    plot.margin  = margin(2, 2, 0, 2)
  )

strip_theme <- theme_void(base_size = 8) +
  theme(
    plot.margin  = margin(0, 2, 0, 2),
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 7),
    legend.key.size = unit(3, "mm")
  )

# =============================================================
# Helper: read a 9-value (3 conc x 3 rep) gravidity file and
# reshape to long + summary, with a chosen control ordered first.
# =============================================================
read_gravidity <- function(path, controls) {
  raw <- read_excel(path, col_names = FALSE)

  dat <- raw %>%
    rename(drug = 1) %>%
    filter(!is.na(drug)) %>%
    mutate(drug = trimws(as.character(drug))) %>%   # strip stray whitespace
    filter(!is.na(suppressWarnings(as.numeric(.[[2]]))))

  val_cols <- paste0("v", 1:9)
  colnames(dat)[2:10] <- val_cols
  dat <- dat %>% mutate(across(all_of(val_cols), as.numeric))

  col_to_conc <- tibble(
    col  = val_cols,
    conc = rep(conc_levels, each = 3),
    rep  = rep(1:3, times = 3)
  )

  long <- dat %>%
    pivot_longer(all_of(val_cols), names_to = "col", values_to = "pct") %>%
    left_join(col_to_conc, by = "col") %>%
    mutate(conc = factor(conc, levels = conc_levels))

  # ordering: control(s) first, then alphabetical
  all_drugs  <- unique(long$drug)
  rest       <- sort(setdiff(all_drugs, controls))
  drug_order <- c(intersect(controls, all_drugs), rest)

  long <- long %>% mutate(drug = factor(drug, levels = drug_order))
  summ <- long %>%
    group_by(drug, conc) %>%
    summarise(mean = mean(pct), sd = sd(pct), .groups = "drop")

  list(long = long, summ = summ, drug_order = drug_order)
}

# =============================================================
# Helper: Dunnett stats (one-way ANOVA per concentration, each
# drug vs the control), returning star positions per drug x conc.
# =============================================================
sig_symbol <- function(p) {
  ifelse(is.na(p), "",
  ifelse(p < 0.0001, "****",
  ifelse(p < 0.001,  "***",
  ifelse(p < 0.01,   "**",
  ifelse(p < 0.05,   "*", "ns")))))
}

dunnett_stats <- function(long, summ, drug_order, control_name,
                          exclude = NULL, dodge_w = 0.8) {
  n_conc <- length(conc_levels)
  conc_offset <- tibble(
    conc = factor(conc_levels, levels = conc_levels),
    dx   = ((seq_len(n_conc) - 0.5) / n_conc - 0.5) * dodge_w
  )

  stats_list <- lapply(conc_levels, function(cc) {
    d <- long %>%
      filter(conc == cc, !drug %in% exclude) %>%
      mutate(drug = factor(drug))
    if (!(control_name %in% levels(d$drug)) || nlevels(d$drug) < 2) return(NULL)
    d$drug <- relevel(d$drug, ref = control_name)
    fit <- aov(pct ~ drug, data = d)
    tryCatch({
      dn <- summary(glht(fit, linfct = mcp(drug = "Dunnett")))
      tibble(conc = cc,
             drug = sub(" - .*$", "", names(dn$test$coefficients)),
             p    = as.numeric(dn$test$pvalues))
    }, error = function(e) NULL)
  })
  stats_df <- bind_rows(stats_list)
  if (nrow(stats_df) == 0) return(stats_df)

  stats_df %>%
    mutate(sig  = sig_symbol(p),
           conc = factor(conc, levels = conc_levels),
           drug = factor(drug, levels = drug_order)) %>%
    left_join(summ, by = c("drug", "conc")) %>%
    left_join(conc_offset, by = "conc") %>%
    mutate(xpos  = as.numeric(drug) + dx,
           ystar = pmin(mean + sd, 100) + 1.5)
}

# =============================================================
# Helper: build a grouped-bar panel (bars + points + stars).
# =============================================================
make_bar_panel <- function(summ, long, stats_df, control_name, ytitle,
                            xlevels) {
  p <- ggplot(summ, aes(x = drug, y = mean, fill = conc)) +
    geom_col(position = position_dodge(width = 0.8),
             width = 0.7, colour = "black", linewidth = 0.2) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                  position = position_dodge(width = 0.8),
                  width = 0.3, linewidth = 0.25) +
    geom_point(data = long, aes(x = drug, y = pct, group = conc),
               position = position_jitterdodge(jitter.width = 0.12,
                                               dodge.width = 0.8, seed = 1),
               shape = 21, size = 0.7, stroke = 0.2,
               fill = "white", colour = "black", inherit.aes = FALSE) +
    scale_fill_manual(values = conc_fill, name = "Concentration") +
    scale_x_discrete(limits = xlevels, drop = FALSE) +  # shared A/B positions
    scale_y_continuous(limits = c(0, 112), breaks = seq(0, 100, 25),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(y = ytitle, x = NULL) +
    base_theme +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x  = element_blank(),
          plot.margin  = margin(2, 2, 1, 2))

  if (nrow(stats_df) > 0) {
    star_dat <- stats_df %>%
      filter(drug != control_name, sig != "ns", sig != "") %>%
      mutate(sig_v = vapply(sig,
              function(s) paste(strsplit(s, "")[[1]], collapse = "\n"),
              character(1)))
    if (nrow(star_dat) > 0) {
      p <- p + geom_text(data = star_dat,
                         aes(x = xpos, y = ystar, label = sig_v),
                         hjust = 0.5, vjust = 0, lineheight = 0.6,
                         size = 2.7, inherit.aes = FALSE)
    }
  }
  p
}

# =============================================================
# Helper: build the rotated drug-name strip (with green bands).
#   bottom = TRUE extends the band downward to cover hanging names
#   at the very bottom of the figure.
#   The alignment placeholder ("") is never given a label, so the
#   empty Panel B slot renders nothing.
# =============================================================
make_drug_strip <- function(drug_order, bottom = FALSE, bands = TRUE) {
  band_idx <- which(seq_along(drug_order) %% 2 == 0)
  band_dat <- tibble(xmin = band_idx - 0.5, xmax = band_idx + 0.5)

  if (bottom) {
    band_dat <- band_dat %>% mutate(ymin = -3, ymax = 1)
    ann <- tibble(drug = factor(drug_order, levels = drug_order)) %>%
      mutate(lab = ifelse(as.character(drug) == placeholder, "",
                          as.character(drug)))
    p <- ggplot(ann, aes(x = drug, y = 1))
    if (bands)
      p <- p + geom_rect(data = band_dat,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE, fill = band_fill)
    p +
      geom_text(aes(label = lab),
                size = 2.5, angle = 90, hjust = 1, vjust = 0.5) +
      scale_x_discrete(limits = drug_order, expand = c(0, 0.5)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      coord_cartesian(clip = "off") +
      labs(y = NULL) +
      strip_theme +
      theme(plot.margin = margin(0, 2, 22, 2))
  } else {
    ann <- tibble(drug = factor(drug_order, levels = drug_order)) %>%
      mutate(lab = ifelse(as.character(drug) == placeholder, "",
                          as.character(drug)))
    p <- ggplot(ann, aes(x = drug, y = 0.5))
    if (bands)
      p <- p + geom_rect(data = band_dat,
                aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                inherit.aes = FALSE, fill = band_fill)
    p +
      geom_text(aes(label = lab),
                size = 2.5, angle = 90, hjust = 0.5, vjust = 0.5) +
      scale_x_discrete(limits = drug_order, expand = c(0, 0.5)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      coord_cartesian(clip = "off") +
      labs(y = NULL) +
      strip_theme +
      theme(plot.margin = margin(0, 2, 0, 2))
  }
}

# =============================================================
# Panel A - infected animals (control = Infected+DMSO)
# =============================================================
A <- read_gravidity(infected_file,
                     controls = c("Uninfected+DMSO", "Infected+DMSO"))
A_stats <- dunnett_stats(A$long, A$summ, A$drug_order,
                         control_name = "Infected+DMSO",
                         exclude = "Uninfected+DMSO")
if (nrow(A_stats) > 0)
  write.csv(A_stats %>% dplyr::select(drug, conc, p, sig),
            "gravidity_panelA_stats.csv", row.names = FALSE)

p_barA   <- make_bar_panel(A$summ, A$long, A_stats,
                           control_name = "Infected+DMSO",
                           ytitle = "% Gravidity",
                           xlevels = A$drug_order)
p_drugA  <- make_drug_strip(A$drug_order, bottom = FALSE, bands = FALSE)

# =============================================================
# Panel B - uninfected animals (control = DMSO)
# =============================================================
B <- read_gravidity(uninfected_file, controls = c("DMSO"))

# ---- Align the two panels on a common x-axis (robust, name-based) ----
# Panel A's order is the master: [Uninfected+DMSO, Infected+DMSO, <alpha>].
# For each Panel A position we place: the placeholder where Panel A has
# Uninfected+DMSO (absent in B), Panel B's "DMSO" where Panel A has
# Infected+DMSO, and otherwise the same drug name (which Panel B shares).
# Building by name (not by blindly copying positions) guarantees every
# Panel B drug sits under its namesake in Panel A even if sorting differs.
B_data_drugs <- as.character(unique(B$long$drug))

B_levels <- vapply(A$drug_order, function(d) {
  if (d == "Uninfected+DMSO") return(placeholder)   # no Panel B counterpart
  if (d == "Infected+DMSO")  return("DMSO")          # B's control slot
  if (d %in% B_data_drugs)   return(d)               # shared drug
  return(d)                                          # keep slot even if absent
}, character(1), USE.NAMES = FALSE)

# Sanity check: every Panel B drug must appear in B_levels, else warn.
missing_in_levels <- setdiff(B_data_drugs, B_levels)
if (length(missing_in_levels) > 0)
  warning("Panel B drugs not placed on shared axis: ",
          paste(missing_in_levels, collapse = ", "))

B$long <- B$long %>% mutate(drug = factor(as.character(drug), levels = B_levels))
B$summ <- B$summ %>% mutate(drug = factor(as.character(drug), levels = B_levels))
B_control <- "DMSO"

B_stats <- dunnett_stats(B$long, B$summ, B_levels,
                         control_name = B_control,
                         exclude = placeholder)
if (nrow(B_stats) > 0)
  write.csv(B_stats %>% dplyr::select(drug, conc, p, sig),
            "gravidity_panelB_stats.csv", row.names = FALSE)

p_barB   <- make_bar_panel(B$summ, B$long, B_stats,
                           control_name = B_control,
                           ytitle = "% Gravidity",
                           xlevels = B_levels)
p_drugB  <- make_drug_strip(B_levels, bottom = TRUE)

# =============================================================
# Assemble: Panel A (bars + names) over Panel B (bars + names)
# =============================================================
figure <- p_barA / p_drugA / p_barB / p_drugB +
  plot_layout(heights = c(6, 1.7, 6, 1.7), guides = "collect") +
  plot_annotation(tag_levels = list(c("A", "", "B", ""))) &
  theme(legend.position = "right",
        panel.spacing = unit(0, "pt"),
        plot.tag = element_text(face = "bold", size = 10))

# =============================================================
# Save - full JIP width (190 mm)
# =============================================================
ggsave("figure_gravidity.pdf", figure,
       width = 190, height = 150, units = "mm")
ggsave("figure_gravidity.png", figure,
       width = 190, height = 150, units = "mm", dpi = 600)

message("Saved figure_gravidity.pdf and figure_gravidity.png")