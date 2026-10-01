#!/usr/bin/env Rscript
# =============================================================
# JIP Figure - Target conservation and inhibitory activity
#
# Panel A: N. parisii ortholog identity vs % inhibition of
#          infection, for the 12 compounds whose S. cerevisiae
#          target has a detectable N. parisii ortholog.
# Panel B: % inhibition by ortholog status, with a third group
#          in which the previously reported anti-microsporidia
#          compounds are removed.
#
# Style matched to fig1_AB_3.R.
# =============================================================

# ---- Locale ----
# Column headers and axis labels contain a micro sign; R must be in a UTF-8
# locale to render it. Harmless no-op where the locale is already UTF-8.
suppressWarnings(invisible(
  Sys.setlocale("LC_CTYPE", if (.Platform$OS.type == "windows")
    "English_United States.utf8" else "C.UTF-8")))

# ---- Packages ----
library(readxl)
library(tidyverse)
library(patchwork)
library(ggrepel)      # NEW dependency vs fig1_AB_3.R

# ---- File paths ----
infection_file    <- "Data_S1.xlsx"
conservation_file <- "conservation data.xlsx"
infection_sheet   <- "Figure 1A"

# ---- Analysis settings ----
CONC_PRIMARY <- 100      # uM
# Compounds already reported to inhibit microsporidia, removed in the third
# group of Panel B to ask whether the ortholog-present group still shows
# elevated activity once the known actives are set aside.
KNOWN_ACTIVE <- c("Fumagillin", "Dexrazoxane", "Albendazole")

controls     <- c("Uninfected+DMSO", "Infected+DMSO")
control_name <- "Infected+DMSO"

# =============================================================
# 1. Read & reshape
# =============================================================
rep_cols <- function(df, conc) {
  nm  <- names(df)
  sel <- startsWith(nm, paste0(conc, " ")) & grepl("rep", nm, fixed = TRUE)
  if (!any(sel)) stop("No replicate columns found for ", conc, " uM")
  nm[sel]
}

infection <- read_excel(infection_file, sheet = infection_sheet)
cons      <- read_excel(conservation_file)

summ <- infection %>%
  pivot_longer(-drug, names_to = "col", values_to = "pct") %>%
  filter(col %in% rep_cols(infection, CONC_PRIMARY)) %>%
  group_by(drug) %>%
  summarise(mean = mean(pct), sd = sd(pct), .groups = "drop")

ctrl <- summ$mean[summ$drug == control_name]
summ <- summ %>% mutate(inhibition = 100 * (ctrl - mean) / ctrl)

dat <- cons %>%
  transmute(drug     = .data[["drug"]],
            target   = .data[["Target"]],
            identity = suppressWarnings(as.numeric(.data[["N. parisii identity"]]))) %>%
  left_join(summ, by = "drug") %>%
  filter(!drug %in% controls) %>%
  mutate(ortholog_present = !is.na(identity),
         known_active     = drug %in% KNOWN_ACTIVE,
         lab              = sprintf('"%s"~(italic("%s"))', drug, target))

stopifnot(all(!is.na(dat$inhibition)))
write.csv(dat, "conservation_figure_data.csv", row.names = FALSE)

con <- dat %>% filter(ortholog_present)      # Panel A analysis set
n_obs <- nrow(con)

# =============================================================
# 2. Panel A statistics
# =============================================================
pear <- cor.test(con$identity, con$inhibition, method = "pearson")
spear <- suppressWarnings(cor.test(con$identity, con$inhibition,
                                   method = "spearman", exact = FALSE))
fit  <- lm(inhibition ~ identity, data = con)
fs   <- summary(fit)
slope   <- coef(fit)[2]
slope_p <- fs$coefficients[2, 4]

# Sensitivity: collapse compounds that share a target (TOP2) to one point
con_uniq <- con %>% group_by(target) %>%
  summarise(identity = first(identity), inhibition = mean(inhibition),
            .groups = "drop")
pear_uniq <- cor.test(con_uniq$identity, con_uniq$inhibition)

statA <- tibble(
  analysis     = c("All compounds with ortholog", "One point per unique target"),
  n            = c(n_obs, nrow(con_uniq)),
  pearson_r    = c(unname(pear$estimate), unname(pear_uniq$estimate)),
  ci_lo        = c(pear$conf.int[1], pear_uniq$conf.int[1]),
  ci_hi        = c(pear$conf.int[2], pear_uniq$conf.int[2]),
  pearson_p    = c(pear$p.value, pear_uniq$p.value),
  spearman_rho = c(unname(spear$estimate), NA),
  spearman_p   = c(spear$p.value, NA))
write.csv(statA, "panelA_correlation_stats.csv", row.names = FALSE)

# =============================================================
# 3. Panel B groups and statistics
# =============================================================
grp_levels <- c("No\northolog", "Ortholog\npresent",
                "Ortholog present,\nexcluding known inhibitors")

panelB <- bind_rows(
  dat %>% filter(!ortholog_present)              %>% mutate(group = grp_levels[1]),
  dat %>% filter(ortholog_present)               %>% mutate(group = grp_levels[2]),
  dat %>% filter(ortholog_present, !known_active) %>% mutate(group = grp_levels[3])
) %>% mutate(group = factor(group, levels = grp_levels))

# NOTE: group 3 is a subset of group 2, so the two are not independent and a
# global test across all three groups would be invalid. Each group is instead
# compared to the "No ortholog" group with a Wilcoxon rank-sum test.
wt <- function(a, b) suppressWarnings(wilcox.test(
  panelB$inhibition[panelB$group == a],
  panelB$inhibition[panelB$group == b], exact = FALSE))

test12 <- wt(grp_levels[1], grp_levels[2])
test13 <- wt(grp_levels[1], grp_levels[3])

grp_summary <- panelB %>% group_by(group) %>%
  summarise(n = n(), mean = mean(inhibition), sd = sd(inhibition),
            median = median(inhibition), .groups = "drop")
write.csv(grp_summary, "panelB_group_summary.csv", row.names = FALSE)
write.csv(tibble(comparison = c(paste(grp_levels[1], "vs", grp_levels[2]),
                                paste(grp_levels[1], "vs", grp_levels[3])),
                 test = "Wilcoxon rank-sum",
                 p = c(test12$p.value, test13$p.value)),
          "panelB_stats.csv", row.names = FALSE)

# =============================================================
# 4. Colours / theme (matched to fig1_AB_3.R)
# =============================================================
conc_fill <- c("10 \u00b5M" = "#bdd7e7", "40 \u00b5M" = "#6baed6",
               "100 \u00b5M" = "#2171b5")
pt_fill   <- conc_fill[["100 \u00b5M"]]

base_theme <- theme_classic(base_size = 8) +
  theme(
    axis.text.x  = element_text(size = 7),
    axis.text.y  = element_text(size = 7, margin = margin(r = 0)),
    axis.title   = element_text(size = 8),
    axis.title.y = element_text(size = 8, margin = margin(r = 1)),
    legend.position = "none",
    plot.margin  = margin(2, 2, 2, 2)
  )

# =============================================================
# 5. Panel A - identity vs inhibition
# =============================================================
statA_lab <- sprintf(
  paste0("atop(italic(r)=='%.2f'*','~italic(R)^2=='%.2f'*','~italic(p)=='%.2f',",
         "rho=='%.2f'*','~rho^2=='%.2f'*','~italic(p)=='%.2f')"),
  unname(pear$estimate), unname(pear$estimate)^2, pear$p.value,
  unname(spear$estimate), unname(spear$estimate)^2, spear$p.value)

p_cor <- ggplot(con, aes(x = identity, y = inhibition)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = "grey35", linewidth = 0.4) +
  geom_point(shape = 21, size = 2.0, fill = pt_fill,
             colour = "black", stroke = 0.25) +
  geom_text_repel(aes(label = lab), parse = TRUE, size = 1.95,
                  segment.size = 0.2, segment.colour = "grey50",
                  min.segment.length = 0, box.padding = 0.4,
                  point.padding = 0.15, max.overlaps = Inf, seed = 1) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.04, vjust = 1.25,
           label = statA_lab, parse = TRUE, size = 2.3, lineheight = 1.1) +
  scale_x_continuous(limits = c(25, 85), breaks = seq(30, 80, 10),
                     expand = expansion(mult = 0.02)) +
  scale_y_continuous(limits = c(-12, 148), breaks = seq(0, 100, 25),
                     expand = expansion(mult = 0.01)) +
  labs(x = expression(italic("N. parisii") * " ortholog identity to " *
                        italic("S. cerevisiae") * " target (%)"),
       y = sprintf("%% Inhibition of infection (%d \u00b5M)", CONC_PRIMARY)) +
  base_theme

# =============================================================
# 6. Panel B - inhibition by ortholog status
# =============================================================
# Mean +/- 1 SD, matching the summary statistic used in Figs. 1 and 2.
mean_sd <- function(x) {
  data.frame(y = mean(x), ymin = mean(x) - sd(x), ymax = mean(x) + sd(x))
}
fmt_p <- function(p) if (p < 0.001) "italic(p)<'0.001'" else
  sprintf("italic(p)=='%.3f'", p)

y_top  <- 118      # first bracket
y_top2 <- 136      # second bracket
brk    <- 4        # bracket tick height

bracket <- function(x1, x2, y) {
  data.frame(x = c(x1, x1, x2, x2), y = c(y - brk, y, y, y - brk))
}
b12 <- bracket(1, 2, y_top)
b13 <- bracket(1, 3, y_top2)

p_grp <- ggplot(panelB, aes(x = group, y = inhibition)) +
  # SD whiskers, then the group mean as a horizontal bar
  stat_summary(fun.data = mean_sd, geom = "errorbar",
               width = 0.25, linewidth = 0.25, colour = "grey20") +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = "errorbar", width = 0.42,
               linewidth = 0.45, colour = "grey20") +
  geom_point(position = position_jitter(width = 0.13, height = 0, seed = 1),
             shape = 21, size = 1.6, fill = pt_fill,
             colour = "black", stroke = 0.25) +
  geom_path(data = b12, aes(x, y), inherit.aes = FALSE,
            linewidth = 0.25, colour = "grey20") +
  geom_path(data = b13, aes(x, y), inherit.aes = FALSE,
            linewidth = 0.25, colour = "grey20") +
  annotate("text", x = 1.5, y = y_top  + 2, label = fmt_p(test12$p.value),
           parse = TRUE, size = 2.3, vjust = 0) +
  annotate("text", x = 2.0, y = y_top2 + 2, label = fmt_p(test13$p.value),
           parse = TRUE, size = 2.3, vjust = 0) +
  annotate("text", x = seq_along(grp_levels), y = -13,
           label = sprintf("italic(n)=='%d'", grp_summary$n),
           parse = TRUE, size = 2.2, colour = "grey30") +
  scale_y_continuous(limits = c(-16, 152), breaks = seq(0, 100, 25),
                     expand = expansion(mult = 0.01)) +
  labs(x = NULL,
       y = sprintf("%% Inhibition of infection (%d \u00b5M)", CONC_PRIMARY)) +
  base_theme +
  theme(axis.text.x = element_text(size = 6.5, lineheight = 0.9),
        plot.margin = margin(2, 12, 2, 2))

# =============================================================
# 7. Assemble & save
# =============================================================
figure <- p_cor + p_grp +
  plot_layout(widths = c(1.25, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 10))

ggsave("figure_conservation_two_panel.pdf", figure,
       width = 190, height = 95, units = "mm")
ggsave("figure_conservation_two_panel.png", figure,
       width = 190, height = 95, units = "mm", dpi = 600)

# =============================================================
# 8. Console report
# =============================================================
cat("\n--- PANEL A: identity vs inhibition ------------------------\n")
print(as.data.frame(statA), row.names = FALSE, digits = 3)
cat(sprintf("\nLinear model: inhibition = %.2f + %.2f * identity\n",
            coef(fit)[1], slope))
cat(sprintf("  R^2 = %.3f, slope p = %.3f\n", fs$r.squared, slope_p))

cat("\n--- PANEL B: inhibition by ortholog status -----------------\n")
print(as.data.frame(grp_summary), row.names = FALSE, digits = 3)
cat(sprintf("\nWilcoxon rank-sum vs 'No ortholog':\n"))
cat(sprintf("  vs Ortholog present:                    p = %.4f\n", test12$p.value))
cat(sprintf("  vs Ortholog present, excluding known:   p = %.4f\n", test13$p.value))
cat(sprintf("\nCompounds removed in group 3: %s\n",
            paste(KNOWN_ACTIVE, collapse = ", ")))

message("\nSaved figure_conservation_two_panel.pdf and .png")
