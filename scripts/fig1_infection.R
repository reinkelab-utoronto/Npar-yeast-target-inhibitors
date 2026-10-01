# =============================================================
# JIP Figure - Panel A
# Full-width (190 mm) grouped bar graph of % infected per drug
# across 3 concentrations (10, 40, 100 uM), 3 reps each,
# with target annotations + N. parisii identity heatmap strip.
# =============================================================

# ---- Packages ----
# NOTE: load multcomp/MASS BEFORE tidyverse so that dplyr's select(),
# filter(), etc. take masking precedence over MASS's versions.
library(multcomp)    # Dunnett's post-hoc test (glht); loads MASS
library(readxl)
library(tidyverse)
library(patchwork)   # to stack bar panel + heatmap strip
library(ggnewscale)  # second fill scale for the "Not present" NA legend

# ---- File paths ----
infection_file    <- "all infection.xlsx"
severe_file       <- "severe infection.xlsx"
conservation_file <- "conservation data.xlsx"

# =============================================================
# 1. Read & reshape infection data
# =============================================================
# Layout: a header row with concentration labels spanning groups
# of 3 columns, then one row per drug with 9 numeric values.
# We read with no header logic and assign columns ourselves.

raw <- read_excel(infection_file, col_names = FALSE)

# Identify the data rows: first column = drug name, next 9 = values.
# Drop any fully-blank / header rows by requiring a name in col 1
# and at least one numeric value in cols 2:10.
dat <- raw %>%
  rename(drug = 1) %>%
  filter(!is.na(drug)) %>%
  # keep only rows whose 2nd column parses as a number (data rows)
  filter(!is.na(suppressWarnings(as.numeric(.[[2]]))))

# Assign the 9 measurement columns
val_cols <- paste0("v", 1:9)
colnames(dat)[2:10] <- val_cols
dat <- dat %>% mutate(across(all_of(val_cols), as.numeric))

# Map the 9 columns -> 3 concentrations x 3 reps
conc_levels <- c("10 \u00b5M", "40 \u00b5M", "100 \u00b5M")
col_to_conc <- tibble(
  col  = val_cols,
  conc = rep(conc_levels, each = 3),
  rep  = rep(1:3, times = 3)
)

long <- dat %>%
  pivot_longer(all_of(val_cols), names_to = "col", values_to = "pct") %>%
  left_join(col_to_conc, by = "col") %>%
  mutate(conc = factor(conc, levels = conc_levels))

# Summarise: mean +/- SD across the 3 reps per drug x concentration
summ <- long %>%
  group_by(drug, conc) %>%
  summarise(mean = mean(pct), sd = sd(pct), .groups = "drop")

# =============================================================
# 2. Read conservation / target data
# =============================================================
cons <- read_excel(conservation_file)

# Standardise the columns we need (names from the file)
cons_clean <- cons %>%
  transmute(
    drug      = .data[["drug"]],
    target    = .data[["Target"]],
    np_name   = .data[["N. parisii name"]],
    identity  = suppressWarnings(as.numeric(.data[["N. parisii identity"]]))
  )

# =============================================================
# 3. Drug ordering: Uninfected+DMSO, Infected+DMSO, then the rest
#    alphabetically. Applied to every layer.
# =============================================================
controls   <- c("Uninfected+DMSO", "Infected+DMSO")
all_drugs  <- unique(summ$drug)
rest       <- sort(setdiff(all_drugs, controls))
drug_order <- c(intersect(controls, all_drugs), rest)

summ       <- summ %>% mutate(drug = factor(drug, levels = drug_order))
long       <- long %>% mutate(drug = factor(drug, levels = drug_order))
cons_clean <- cons_clean %>%
  filter(drug %in% drug_order) %>%
  mutate(drug = factor(drug, levels = drug_order)) %>%
  arrange(drug)

# =============================================================
# 3b. Statistics: one-way ANOVA per concentration across drugs,
#     followed by Dunnett's post-hoc comparing each drug to the
#     Infected+DMSO control AT THE SAME CONCENTRATION.
#     Produces one significance star per drug x concentration.
# =============================================================
control_name   <- "Infected+DMSO"
exclude_from_stats <- "Uninfected+DMSO"   # not a treatment; don't test

stats_list <- lapply(conc_levels, function(cc) {
  d <- long %>%
    filter(conc == cc, drug != exclude_from_stats) %>%
    mutate(drug = factor(drug))

  if (!(control_name %in% levels(d$drug)) || nlevels(d$drug) < 2) return(NULL)

  # Drop drug groups with zero variance only if ALL groups are constant;
  # Dunnett needs residual variance > 0 overall.
  d$drug <- relevel(d$drug, ref = control_name)

  fit <- aov(pct ~ drug, data = d)

  tryCatch({
    dunnett <- summary(glht(fit, linfct = mcp(drug = "Dunnett")))
    comps   <- names(dunnett$test$coefficients)        # "drug - Infected+DMSO"
    tibble(
      conc = cc,
      drug = sub(" - .*$", "", comps),
      p    = as.numeric(dunnett$test$pvalues)
    )
  }, error = function(e) NULL)
})

stats_df <- bind_rows(stats_list)

sig_symbol <- function(p) {
  ifelse(is.na(p), "",
  ifelse(p < 0.0001, "****",
  ifelse(p < 0.001,  "***",
  ifelse(p < 0.01,   "**",
  ifelse(p < 0.05,   "*", "ns")))))
}

if (nrow(stats_df) > 0) {
  # numeric dodge offset per concentration so each star sits exactly
  # over the center of its bar. Bars use dodge width 0.8 across 3 groups.
  n_conc      <- length(conc_levels)
  dodge_w     <- 0.8
  # centers of the 3 dodged bars, relative to the drug's integer x.
  conc_offset <- tibble(
    conc = factor(conc_levels, levels = conc_levels),
    dx   = ((seq_len(n_conc) - 0.5) / n_conc - 0.5) * dodge_w
  )

  stats_df <- stats_df %>%
    mutate(
      sig  = sig_symbol(p),
      conc = factor(conc, levels = conc_levels),
      drug = factor(drug, levels = drug_order)
    ) %>%
    left_join(summ, by = c("drug", "conc")) %>%
    left_join(conc_offset, by = "conc") %>%
    mutate(
      xpos  = as.numeric(drug) + dx,            # explicit bar-center x
      ystar = pmin(mean + sd, 100) + 1.5        # star y-position
    )

  write.csv(stats_df %>% dplyr::select(drug, conc, p, sig),
            "panelA_stats_dunnett.csv", row.names = FALSE)
}


# Build the multi-line annotation label shown under each drug:
#   Target (yeast)
#   N. parisii name
ann <- cons_clean %>%
  mutate(
    np_name = ifelse(is.na(np_name) | np_name == "N/A", "\u2013", np_name),
    target  = ifelse(is.na(target), "\u2013", target),
    lab     = paste0(target, "\n", np_name)
  )

# Some drugs in the bar data may lack a conservation entry; keep them
# but with blank annotation so axes stay aligned.
ann_full <- tibble(drug = factor(drug_order, levels = drug_order)) %>%
  left_join(ann, by = "drug") %>%
  mutate(
    lab      = ifelse(is.na(lab), "\u2013\n\u2013", lab),
    identity = identity
  )

# =============================================================
# 4. Colours / theme
# =============================================================
conc_fill <- c("10 \u00b5M" = "#bdd7e7",
               "40 \u00b5M" = "#6baed6",
               "100 \u00b5M" = "#2171b5")

base_theme <- theme_classic(base_size = 8) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    axis.text.y  = element_text(size = 7, margin = margin(r = 0)),
    axis.title   = element_text(size = 8),
    axis.title.y = element_text(size = 8, margin = margin(r = -3)),
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 7),
    legend.key.size = unit(3, "mm"),
    plot.margin  = margin(2, 2, 0, 2)
  )

# Tidy the annotation fields for display
ann_full <- ann_full %>%
  mutate(
    target_disp = ifelse(is.na(target) | target == "\u2013", "", target),
    npname_disp = ifelse(is.na(np_name) | np_name == "\u2013", "", np_name)
  )


# =============================================================
# 4b. Multi-gene label handling
# =============================================================
# Some compounds have more than one relevant gene (e.g. the 5-FC
# activation pathway FCY2/FCY1/FUR1). Because the annotation text is
# rotated 90 degrees and the strips use clip = "off", a long single-line
# label overflows its strip and collides with the rows above and below.
#
# wrap_genes() splits on / , or | and puts each gene on its own line.
# For rotated text, extra lines stack SIDEWAYS, so the vertical extent
# is set by the longest single gene name rather than the whole string.
# fit_size() then shrinks any label whose longest line is still wide.
wrap_genes <- function(x) {
  vapply(x, function(v) {
    if (is.na(v) || v == "" || v == "\u2013") return(v)
    paste(trimws(strsplit(v, "\\s*[/,|]\\s*")[[1]]), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

# Longest line in a (possibly wrapped) label, in characters.
longest_line <- function(x) {
  vapply(x, function(v) {
    if (is.na(v) || v == "") return(0L)
    max(nchar(strsplit(v, "\n")[[1]]))
  }, integer(1), USE.NAMES = FALSE)
}

# Shrink only labels that are still too long; everything else keeps base size.
fit_size <- function(x, base, max_chars) {
  n <- longest_line(x)
  pmin(base, base * max_chars / pmax(n, max_chars))
}

ann_full <- ann_full %>%
  mutate(
    target_disp = wrap_genes(target_disp),
    npname_disp = wrap_genes(npname_disp),
    target_size = fit_size(target_disp, base = 2.4,  max_chars = 6),
    npname_size = fit_size(npname_disp, base = 2.16, max_chars = 10)
  )

# =============================================================
# 5. Panel A - bar graph (top)
# =============================================================
# Alternating vertical guide bands so the eye can follow each drug
# column down through the bars and the annotation rows. Every other
# drug gets a faint grey band. Defined here so all panels can reuse it.
band_idx  <- which(seq_along(drug_order) %% 2 == 0)   # shade even columns
band_dat  <- tibble(
  xmin = band_idx - 0.5,
  xmax = band_idx + 0.5
)
band_fill <- "#e6f4ea"   # light green column dividers (distinct from NA grey)

make_bands <- function(ymin = -Inf, ymax = Inf) {
  geom_rect(data = band_dat,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = band_fill)
}

p_bar <- ggplot(summ, aes(x = drug, y = mean, fill = conc)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7, colour = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(width = 0.8),
                width = 0.3, linewidth = 0.25) +
  # individual biological replicates (3 points per bar)
  geom_point(data = long,
             aes(x = drug, y = pct, group = conc),
             position = position_jitterdodge(jitter.width = 0.12,
                                             dodge.width = 0.8,
                                             seed = 1),
             shape = 21, size = 0.7, stroke = 0.2,
             fill = "white", colour = "black",
             inherit.aes = FALSE) +
  scale_fill_manual(values = conc_fill, name = "Concentration") +
  scale_y_continuous(limits = c(0, 112),
                     breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(y = "% Infected", x = NULL) +
  base_theme +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x  = element_blank(),
        plot.margin  = margin(2, 2, 1, 2))   # minimal gap below bars

# Significance stars from Dunnett's test (skip the control itself & "ns").
if (exists("stats_df") && nrow(stats_df) > 0) {
  star_dat <- stats_df %>%
    filter(drug != control_name, sig != "ns", sig != "")
  if (nrow(star_dat) > 0) {
    # Stack the significance characters vertically by inserting newlines
    # (e.g. "****" -> "*\n*\n*\n*"). Upright text centers cleanly on x
    # with hjust = 0.5, avoiding rotated-text anchor ambiguity.
    star_dat <- star_dat %>%
      mutate(sig_v = vapply(
        sig,
        function(s) paste(strsplit(s, "")[[1]], collapse = "\n"),
        character(1)
      ))
    p_bar <- p_bar +
      geom_text(data = star_dat,
                aes(x = xpos, y = ystar, label = sig_v),
                hjust = 0.5, vjust = 0, lineheight = 0.6,
                size = 2.7, inherit.aes = FALSE)
  }
}

# =============================================================
# 6. Annotation rows beneath the bars
#    Each piece of info gets its OWN horizontal row so nothing
#    overlaps:
#      row 3 (top)    : yeast target  (italic)
#      row 2 (middle) : N. parisii name
#      row 1 (bottom) : identity heatmap tile + number
#    Drug names are the rotated x-axis tick labels at the very
#    bottom of this strip.
# =============================================================

# Shared discrete x scale for every strip (identical ordering/expansion
# so the columns line up perfectly across panels).
x_scale <- scale_x_discrete(limits = drug_order, expand = c(0, 0.5))

# A theme for the text/heatmap strips: NO y axis at all, NO x text.
strip_theme <- theme_void(base_size = 8) +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 7.5,
                               margin = margin(r = -3)),
    plot.margin  = margin(0, 2, 0, 2),
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 7),
    legend.key.size = unit(3, "mm")
  )

# ---- Drug name row (top of the annotation block), rotated 90 deg ----
# Names centered within the strip so the green band fully covers them.
p_drug <- ggplot(ann_full, aes(x = drug, y = 0.5)) +
  make_bands() +
  geom_text(aes(label = as.character(drug)),
            size = 2.5, angle = 90, hjust = 0.5, vjust = 0.5) +
  x_scale +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(y = NULL) +
  strip_theme +
  theme(plot.margin = margin(0, 2, 0, 2))

# ---- Yeast target row (italic gene names), rotated 90 deg ----
p_target <- ggplot(ann_full, aes(x = drug, y = 0.5)) +
  make_bands() +
  geom_text(aes(label = target_disp, size = target_size),
            fontface = "italic", angle = 90, hjust = 0.5, vjust = 0.5,
            lineheight = 0.85) +
  scale_size_identity() +
  x_scale +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(y = "Yeast target") +
  strip_theme

# ---- N. parisii ortholog row, rotated 90 deg ----
p_npname <- ggplot(ann_full, aes(x = drug, y = 0.5)) +
  make_bands() +
  geom_text(aes(label = npname_disp, size = npname_size),
            angle = 90, hjust = 0.5, vjust = 0.5, lineheight = 0.85) +
  scale_size_identity() +
  x_scale +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(y = expression(italic("N. parisii") * " ortholog")) +
  strip_theme

# ---- Identity heatmap row (bottom of the annotation block). No drug
#      names here anymore -- they live on p_drug at the top. ----
# Separate the rows that have an identity value from those that don't,
# so the NA ("Not present") tiles get their own legend key.
heat_val <- ann_full %>% filter(!is.na(identity))
heat_na  <- ann_full %>%
  filter(is.na(identity), !drug %in% controls) %>%   # skip DMSO controls
  mutate(na_lab = "Not present")

p_heat <- ggplot(mapping = aes(x = drug, y = 1)) +
  # NA tiles first, with a manual discrete fill -> own legend entry
  geom_tile(data = heat_na, aes(fill = na_lab),
            colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = c("Not present" = "grey90"),
                    name = NULL) +
  ggnewscale::new_scale_fill() +
  # value tiles with the continuous orange gradient
  geom_tile(data = heat_val, aes(fill = identity),
            colour = "white", linewidth = 0.4) +
  geom_text(data = heat_val,
            aes(label = sprintf("%.0f", identity)),
            size = 2.1,
            colour = ifelse(heat_val$identity > 55, "white", "black")) +
  scale_fill_gradient(low = "#fff5eb", high = "#d94801",
                      name = expression(italic("N. parisii") * " identity (%)"),
                      limits = c(0, 100)) +
  x_scale +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "% Identity") +
  strip_theme

# =============================================================
# 6b. PANEL B - severe infection, single concentration (100 uM)
# =============================================================
# Layout of severe infection.xlsx: a header row ("100 uM" over the
# group), then one row per drug with 3 replicate values.
sev_raw <- read_excel(severe_file, col_names = FALSE)

sev_dat <- sev_raw %>%
  rename(drug = 1) %>%
  filter(!is.na(drug)) %>%
  filter(!is.na(suppressWarnings(as.numeric(.[[2]]))))

sev_val_cols <- paste0("v", 1:3)
colnames(sev_dat)[2:4] <- sev_val_cols
sev_dat <- sev_dat %>% mutate(across(all_of(sev_val_cols), as.numeric))

sev_conc <- "100 \u00b5M"   # single concentration

sev_long <- sev_dat %>%
  pivot_longer(all_of(sev_val_cols), names_to = "col", values_to = "pct") %>%
  mutate(conc = factor(sev_conc, levels = conc_levels),
         drug = factor(drug, levels = drug_order))

sev_summ <- sev_long %>%
  group_by(drug, conc) %>%
  summarise(mean = mean(pct), sd = sd(pct), .groups = "drop")

# ---- Stats: one-way ANOVA + Dunnett vs Infected+DMSO (same design) ----
sev_stats <- {
  d <- sev_long %>%
    filter(drug != exclude_from_stats) %>%
    mutate(drug = factor(drug))
  if (control_name %in% levels(d$drug) && nlevels(d$drug) >= 2) {
    d$drug <- relevel(d$drug, ref = control_name)
    fit <- aov(pct ~ drug, data = d)
    tryCatch({
      dn <- summary(glht(fit, linfct = mcp(drug = "Dunnett")))
      tibble(
        drug = sub(" - .*$", "", names(dn$test$coefficients)),
        p    = as.numeric(dn$test$pvalues)
      )
    }, error = function(e) NULL)
  } else NULL
}

if (!is.null(sev_stats) && nrow(sev_stats) > 0) {
  sev_stats <- sev_stats %>%
    mutate(
      sig  = sig_symbol(p),
      conc = factor(sev_conc, levels = conc_levels),
      drug = factor(drug, levels = drug_order)
    ) %>%
    left_join(sev_summ, by = c("drug", "conc")) %>%
    mutate(
      xpos  = as.numeric(drug),                 # single bar -> centered on tick
      ystar = pmin(mean + sd, 100) + 1.5
    )
  write.csv(sev_stats %>% dplyr::select(drug, conc, p, sig),
            "panelB_stats_dunnett.csv", row.names = FALSE)
}

# ---- Panel B bar plot (single 100 uM bar per drug) ----
p_barB <- ggplot(sev_summ, aes(x = drug, y = mean)) +
  geom_col(width = 0.7, fill = conc_fill[["100 \u00b5M"]],
           colour = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.3, linewidth = 0.25) +
  geom_point(data = sev_long, aes(x = drug, y = pct),
             position = position_jitter(width = 0.12, height = 0, seed = 1),
             shape = 21, size = 0.7, stroke = 0.2,
             fill = "white", colour = "black", inherit.aes = FALSE) +
  scale_y_continuous(limits = c(0, 112),
                     breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(y = "% Severe infection", x = NULL) +
  base_theme +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x  = element_blank(),
        plot.margin  = margin(2, 2, 1, 2))

# stars (vertically stacked, centered) on Panel B
if (!is.null(sev_stats) && nrow(sev_stats) > 0) {
  sev_star <- sev_stats %>%
    filter(drug != control_name, sig != "ns", sig != "") %>%
    mutate(sig_v = vapply(sig,
            function(s) paste(strsplit(s, "")[[1]], collapse = "\n"),
            character(1)))
  if (nrow(sev_star) > 0) {
    p_barB <- p_barB +
      geom_text(data = sev_star,
                aes(x = xpos, y = ystar, label = sig_v),
                hjust = 0.5, vjust = 0, lineheight = 0.6,
                size = 2.7, inherit.aes = FALSE)
  }
}

# ---- Panel B drug-name strip (repeated under Panel B) ----
# Bands extended downward (ymin < 0) with clip off so the green fully
# covers the longest rotated drug name, which overflows below the strip.
band_dat_B <- band_dat %>% mutate(ymin = -3, ymax = 1)
p_drugB <- ggplot(ann_full, aes(x = drug, y = 1)) +
  geom_rect(data = band_dat_B,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = band_fill) +
  geom_text(aes(label = as.character(drug)),
            size = 2.5, angle = 90, hjust = 1, vjust = 0.5) +
  x_scale +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(y = NULL) +
  strip_theme +
  theme(plot.margin = margin(0, 2, 22, 2))   # bottom room for hanging names


# =============================================================
# 7. Keep each y-axis title close to its own axis
# =============================================================
# patchwork normally aligns the widths of left-side plot decorations
# across vertically stacked plots. Because the annotation-strip titles
# are horizontal and relatively long, that alignment pushes the y-axis
# titles for Panels A and B far away from their axes.
#
# free(type = "label", side = "l") frees ONLY the left-side axis-title
# alignment. The plotting panels remain aligned, so all drug columns
# still line up vertically.
p_bar_close    <- patchwork::free(p_bar,    type = "label", side = "l")
p_target_close <- patchwork::free(p_target, type = "label", side = "l")
p_npname_close <- patchwork::free(p_npname, type = "label", side = "l")
p_heat_close   <- patchwork::free(p_heat,   type = "label", side = "l")
p_barB_close   <- patchwork::free(p_barB,   type = "label", side = "l")

panelA <- p_bar_close / p_drug / p_target_close / p_npname_close / p_heat_close +
  plot_layout(heights = c(6, 1.8, 1.5, 1.3, 0.7), guides = "collect") &
  theme(legend.position = "right",
        legend.justification = "top",   # pin legend to top of its region
        panel.spacing = unit(0, "pt"))

# =============================================================
# 7b. Assemble the FULL figure (Panel A block + Panel B block).
#     Both bar graphs share identical x ordering, so the drug
#     columns line up vertically across the whole figure.
# =============================================================
figure <- p_bar_close / p_drug / p_target_close / p_npname_close / p_heat_close /
          p_barB_close / p_drugB +
  plot_layout(heights = c(6, 1.7, 1.3, 1.1, 0.5,   # Panel A block
                          6, 1.7),                    # Panel B block
              guides = "collect") +
  plot_annotation(tag_levels = list(c("A", "", "", "", "", "B", ""))) &
  theme(legend.position = "right",
        legend.justification = "top",   # pin collected legend to the top
        panel.spacing = unit(0, "pt"),
        plot.tag = element_text(face = "bold", size = 10))

# =============================================================
# 8. Save - full JIP width (190 mm)
# =============================================================
ggsave("figure_panelAB.pdf", figure,
       width = 190, height = 188, units = "mm")
ggsave("figure_panelAB.png", figure,
       width = 190, height = 188, units = "mm", dpi = 600)

message("Saved figure_panelAB.pdf and figure_panelAB.png")