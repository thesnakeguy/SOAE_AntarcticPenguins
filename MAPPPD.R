
######################################################################################################
## This script draws penguin data from https://ipt.biodiversity.aq/resource?r=mapppd_count_data and
## performs data visualization as well as statitical tests for (recent) population trends
######################################################################################################

library(MASS)
library(ggplot2)
library(tidyverse)
library(mgcv)
library(emmeans)
library(lme4)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)   
library(scales)
library(ggrepel)
library(dplyr)


###################################
### 1. Preparing data for analysis ####
###################################

# Draw data in
if (!file.exists("MAPPPD_AllCounts.csv")) {
  download.file("https://www.penguinmap.com/mapppd/DownloadAll",
                "MAPPPD_AllCounts.csv",
                mode = "wb")
}

pengs_raw <- readr::read_csv("MAPPPD_AllCounts.csv", show_col_types = FALSE)

# Filter and clean
pengs_clean <- pengs_raw %>%
  mutate(species = case_when(
    str_detect(common_name, "adelie") ~ "Adelie",
    str_detect(common_name, "chinstrap") ~ "Chinstrap",
    str_detect(common_name, "gentoo") ~ "Gentoo",
    str_detect(common_name, "macaroni penguin") ~ "Macaroni penguin",
    str_detect(common_name, "king penguin") ~ "King penguin",
    str_detect(common_name, "emperor penguin") ~ "Emperor penguin"
  ))

# Check what count types you're working with
table(pengs_clean$count_type)

# See which sites have consistent monitoring
monitoring_summary <- pengs_clean %>%
  group_by(cammlr_region, site_name, species) %>%
  summarize(
    n_years = n_distinct(year),
    first_year = min(year),
    last_year = max(year),
    year_span = last_year - first_year + 1,
    coverage = n_years / year_span,
    # proportion of years monitored
    .groups = "drop"
  ) %>%
  arrange(desc(n_years))

# Focus on well-monitored sites since 2000
coverage_threshold <- 0.5
n_yrs_req <- (as.numeric((format(Sys.Date(), "%Y")))-2000)*coverage_threshold

well_monitored_summary <- pengs_clean %>%
  filter(year >= 2000) %>%
  group_by(cammlr_region, site_name, species) %>%
  summarize(
    n_years = n_distinct(year),
    first_year = min(year),
    last_year = max(year),
    year_span = last_year - first_year + 1,
    coverage = n_years / year_span
  ) %>%
  filter(n_years > n_yrs_req) %>%
  arrange(desc(n_years))

# Join back to get analysis dataset
pengs_analysis <- pengs_clean %>%
  right_join(well_monitored_summary, by = c("cammlr_region", "site_name", "species"))

# Aggregate to site-year-species level (just to be sure, if a site would have been visited twice in one year)
pengs_yearly <- pengs_analysis %>%
  group_by(cammlr_region, site_name, site_id, species, year, count_type, longitude_epsg_4326, latitude_epsg_4326) %>%
  summarize(
    mean_count = mean(penguin_count, na.rm = TRUE),
    max_count = max(penguin_count, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

###################################
### 2. Bubble map grid  ####
###################################

bubble_data <- filter(pengs_yearly, year > 2000)

# --- Trends + latest population size 
calculate_trends <- function(data) {
  data %>%
    group_by(
      site_id, site_name,
      longitude_epsg_4326, latitude_epsg_4326,
      species,       
      count_type
    ) %>%
    summarise(
      slope        = coef(lm(mean_count ~ year))[2],
      intercept    = coef(lm(mean_count ~ year))[1],
      latest_count = mean_count[which.max(year)],
      .groups = "drop"
    )
}

trend_data <- calculate_trends(bubble_data) %>%
  mutate(
    longitude_epsg_4326 = as.numeric(longitude_epsg_4326),
    latitude_epsg_4326  = as.numeric(latitude_epsg_4326)
  )

last_year <- max(bubble_data$year)

trend_data_sf <- st_as_sf(
  trend_data,
  coords = c("longitude_epsg_4326", "latitude_epsg_4326"),
  crs = 4326
) %>%
  st_transform(crs = 3031)

antarctica <- ne_countries(
  continent = "Antarctica",
  returnclass = "sf",
  scale = "medium"
) %>%
  st_transform(crs = 3031)

size_cap    <- quantile(trend_data_sf$latest_count, 0.95, na.rm = TRUE)
size_limits <- c(0, size_cap)

slope_cap    <- quantile(abs(trend_data_sf$slope), 0.95, na.rm = TRUE)
slope_limits <- c(-slope_cap, slope_cap)  # symmetric around 0

species_list    <- sort(unique(trend_data_sf$species))
count_type_list <- sort(unique(trend_data_sf$count_type))

get_species_extent <- function(sp, pad_frac = 0.35, min_half_extent = 400000) {
  sub <- trend_data_sf %>% filter(species == sp)
  bb  <- st_bbox(sub)
  
  half_x <- max((bb["xmax"] - bb["xmin"]) / 2, min_half_extent)
  half_y <- max((bb["ymax"] - bb["ymin"]) / 2, min_half_extent)
  cx <- mean(c(bb["xmin"], bb["xmax"]))
  cy <- mean(c(bb["ymin"], bb["ymax"]))
  
  # Make it square (so coord_sf doesn't distort) using the larger half-extent,
  # then pad by pad_frac for breathing room around the outermost points.
  half <- max(half_x, half_y) * (1 + pad_frac)
  
  list(
    xlim = c(cx - half, cx + half),
    ylim = c(cy - half, cy + half)
  )
}

# --- Build one column (a stack of count_type panels) per species 
make_species_column <- function(sp) {
  ext <- get_species_extent(sp)
  sub <- trend_data_sf %>% filter(species == sp)
  
  ggplot() +
    geom_sf(data = antarctica, fill = "grey88", color = "grey50", linewidth = 0.2) +
    geom_sf(
      data = sub,
      aes(size = latest_count, fill = slope),
      shape = 21, color = "black", stroke = 0.3, alpha = 0.85
    ) +
    scale_size_area(
      max_size = 14,
      limits = size_limits,
      oob = scales::squish,
      labels = function(x) ifelse(x >= size_cap, paste0(scales::comma(x), "+"), scales::comma(x)),
      name = "Population size\n(latest count, capped at 95th pct)"
    ) +
    scale_fill_gradient2(
      low = "red", mid = "grey90", high = "blue",
      midpoint = 0, limits = slope_limits,
      oob = scales::squish,
      name = "Trend (slope,\ncapped at 95th pct)"
    ) +
    coord_sf(crs = st_crs(3031), xlim = ext$xlim, ylim = ext$ylim) +
    facet_wrap(~count_type, ncol = 1, strip.position = "right") +
    geom_text_repel(
      data = sub,
      aes(label = site_name, geometry = geometry),
      stat = "sf_coordinates",
      segment.color = "grey60",
      min.segment.length = 1,
      max.overlaps = 15,
      show.legend = FALSE
    ) +
    labs(title = sp) +
    theme_antarctic() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_line(color = "grey93"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      strip.text = element_text(face = "bold")
    )
}

species_plots <- lapply(species_list, make_species_column)

final_plot <- wrap_plots(species_plots, nrow = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste0("Penguin count bubble map - trend 2000 until ", last_year),
    theme = theme(plot.title = element_text(size = 20, face = "bold"), 
                  legend.position = "right")
  )

print(final_plot)

ggsave(
  "Penguin_bubble_map.svg",
  final_plot,
  width = 20,     
  height = 12,    
  dpi = 300
)


###################################
### 3. Raw site-level trajectories (GAM modeling)  ####
###################################

line_data <- pengs_yearly %>%
  filter(year > 2000, !is.na(mean_count), mean_count > 0, count_type == "nests") %>%
  arrange(species, site_id, year)

# ------------------------------------------------------------------
# One GAM per site -> smooth trajectory + average % change per year
#
# Model:  log(mean_count) ~ s(year, k), fitted with REML.
#   - the smooth lets the trajectory bend (crash/recovery, plateau, ...)
#   - the log scale keeps it multiplicative, same as the log-linear version
#
# The headline number is the *average compound annual rate* implied by the
# smooth between the first and last survey year:
#     (exp((f(end) - f(start)) / n_years) - 1) * 100
# Its 95% interval comes from posterior simulation of the GAM coefficients
# (Bayesian covariance matrix, incl. smoothing-parameter uncertainty), so it
# is a proper uncertainty band for the *smooth*, not just for a straight line.
# ------------------------------------------------------------------
min_years <- 6      # fewer survey years than this -> site is not modelled
k_max     <- 6      # cap on basis dimension; keeps curves sensible for short series
n_sim     <- 2000   # posterior draws for the CI

fit_site_gam <- function(d) {
  d <- arrange(d, year)
  k <- min(k_max, n_distinct(d$year) - 1)
  m <- gam(log(mean_count) ~ s(year, k = k), data = d, method = "REML")

  # predict on a complete yearly grid (interpolates across survey gaps)
  yrs <- seq(min(d$year), max(d$year))
  Xp  <- predict(m, data.frame(year = yrs), type = "lpmatrix")
  eta_hat <- as.numeric(Xp %*% coef(m))              # smooth on log scale

  set.seed(1)
  beta <- MASS::mvrnorm(n_sim, coef(m), vcov(m, unconditional = TRUE))
  eta  <- Xp %*% t(beta)                             # years x draws (log scale)

  span <- max(yrs) - min(yrs)
  rate_draws <- (exp((eta[length(yrs), ] - eta[1, ]) / span) - 1) * 100
  ci <- quantile(rate_draws, c(0.025, 0.975), names = FALSE)

  list(
    trend = tibble(
      pct_per_yr = (exp((eta_hat[length(yrs)] - eta_hat[1]) / span) - 1) * 100,
      lo  = ci[1],
      hi  = ci[2],
      edf = summary(m)$edf          # ~1 = essentially linear; larger = more wiggly
    ),
    curve = tibble(
      year = yrs,
      fit  = exp(eta_hat),
      lo   = exp(apply(eta, 1, quantile, 0.025)),
      hi   = exp(apply(eta, 1, quantile, 0.975))
    )
  )
}

fits <- line_data %>%
  group_by(species, site_id, site_name) %>%
  filter(n_distinct(year) >= min_years) %>%
  nest() %>%
  ungroup() %>%
  mutate(fit = map(data, fit_site_gam))

trend_tbl <- fits %>%
  transmute(species, site_id, site_name, trend = map(fit, "trend")) %>%
  unnest(trend) %>%
  mutate(
    direction = case_when(lo > 0 ~ "Increasing",
                          hi < 0 ~ "Decreasing",
                          TRUE   ~ "No clear trend"),
    direction = factor(direction, c("Increasing", "No clear trend", "Decreasing")),
    panel     = paste0(site_name, "\n", sprintf("%+.1f%%/yr", pct_per_yr))
  )

curve_tbl <- fits %>%
  transmute(species, site_id, site_name, curve = map(fit, "curve")) %>%
  unnest(curve) %>%
  left_join(select(trend_tbl, species, site_id, direction, panel),
            by = c("species", "site_id"))

dir_cols <- c("Increasing"     = "#1b9e77",
              "No clear trend" = "grey55",
              "Decreasing"     = "#d95f02")

plot_data <- line_data %>%
  inner_join(select(trend_tbl, species, site_id, direction, panel, pct_per_yr),
             by = c("species", "site_id"))

# ------------------------------------------------------------------
# One panel per site, own log y-axis, GAM smooth + 95% band
# ------------------------------------------------------------------
plot_species <- function(pts, crv, ncol = 4) {
  ggplot() +
    geom_line(data = pts, aes(year, mean_count), colour = "grey75", linewidth = 0.4) +
    geom_ribbon(data = crv, aes(year, ymin = lo, ymax = hi, fill = direction),
                alpha = 0.18) +
    geom_line(data = crv, aes(year, fit, colour = direction), linewidth = 0.9) +
    geom_point(data = pts, aes(year, mean_count), colour = "grey25", size = 0.9) +
    facet_wrap(vars(panel), ncol = ncol, scales = "free_y") +
    scale_y_log10(labels = label_number(scale_cut = cut_short_scale()),
                  breaks = log_breaks(n = 3)) +
    scale_x_continuous(breaks = seq(2000, 2030, by = 10)) +
    scale_colour_manual(values = dir_cols, drop = FALSE) +
    scale_fill_manual(values = dir_cols, drop = FALSE) +
    labs(title = unique(pts$species), x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title       = element_text(face = "bold", size = 12),
      strip.text       = element_text(hjust = 0, size = 8.5, lineheight = 1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92"),
      legend.position  = "none"
    )
}

species_names <- unique(plot_data$species)
ncol_panels   <- 4

plots <- map(species_names, ~ plot_species(filter(plot_data, species == .x),
                                           filter(curve_tbl, species == .x),
                                           ncol = ncol_panels))
heights <- map_dbl(species_names,
                   ~ ceiling(n_distinct(filter(plot_data, species == .x)$site_id) / ncol_panels))

site_plot <- wrap_plots(plots, ncol = 1, heights = heights) +
  plot_annotation(
    title    = "Nest-count trajectories by site",
    subtitle = paste0("Log-scaled y-axis (free per site); line = site-level GAM smooth with 95% band; ",
                      "panel label = average % change per year"),
    caption  = "Green = significant increase, orange = significant decrease, grey = no clear trend",
    theme    = theme(plot.title = element_text(face = "bold", size = 14))
  )

# ggsave("site_level_trends.png", site_plot, width = 12, height = 14, dpi = 300)

###################################
### 4. Recent trends: mixed model  ####
###################################

# ------------------------------------------------------------------
# a. parameterization
# ------------------------------------------------------------------

# Statistical test (5 most recent yrs vs 5 prior to that)
latest           <- max(pengs_yearly$year[!is.na(pengs_yearly$max_count)])
recent_years     <- (latest - 4):latest
historical_years <- (latest - 9):(latest - 5)   # window midpoints are exactly 5 years apart

# Options
min_sites         <- 5      # random-effect variance is unreliable with very few sites
paired_sites_only <- TRUE   # only sites with data in BOTH windows -> a true "change" estimate
use_year_re       <- TRUE   # (1 | year): shared good/bad years; stops CIs being too narrow
                            # (set FALSE if models struggle to converge)

mm_data <- pengs_yearly %>%
  filter(!is.na(max_count), year %in% c(historical_years, recent_years)) %>%
  mutate(
    period = factor(if_else(year %in% recent_years, "recent", "historical"),
                    levels = c("historical", "recent")),   # explicit reference level
    year_f = factor(year)
  )

if (paired_sites_only) {
  mm_data <- mm_data %>%
    group_by(cammlr_region, species, count_type, site_id) %>%
    filter(n_distinct(period) == 2) %>%
    ungroup()
}

# ------------------------------------------------------------------
# b. Negative binomial mixed model per region x species x count_type
# ------------------------------------------------------------------
ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

fit_nb <- possibly(function(d) {
  if (use_year_re) {
    glmer.nb(max_count ~ period + (1 | site_id) + (1 | year_f), data = d, control = ctrl)
  } else {
    glmer.nb(max_count ~ period + (1 | site_id), data = d, control = ctrl)
  }
}, otherwise = NULL)   # a failed fit returns NULL instead of killing the whole pipeline

mm_fits <- mm_data %>%
  group_by(cammlr_region, species, count_type) %>%
  nest() %>%
  ungroup() %>%
  mutate(
    n_sites          = map_int(data, ~ n_distinct(.x$site_id)),
    n_obs_historical = map_int(data, ~ sum(.x$period == "historical")),
    n_obs_recent     = map_int(data, ~ sum(.x$period == "recent"))
  ) %>%
  filter(n_sites >= min_sites) %>%
  mutate(model = map(data, fit_nb))

failed <- mm_fits %>% filter(map_lgl(model, is.null)) %>%
  select(cammlr_region, species, count_type, n_sites)
if (nrow(failed) > 0) { message("Models that failed to fit:"); print(failed) }

# ------------------------------------------------------------------
# c. Summaries: % change + 95% CI, and the same effect as an annual rate
# ------------------------------------------------------------------
dir_cols <- c("Increasing"     = "#1b9e77",
              "No clear trend" = "grey55",
              "Decreasing"     = "#d95f02")

mm_tbl <- mm_fits %>%
  filter(!map_lgl(model, is.null)) %>%
  mutate(
    singular = map_lgl(model, isSingular),   # TRUE = a random-effect variance collapsed to ~0
    est = map(model, ~ broom.mixed::tidy(.x, effects = "fixed", conf.int = TRUE,
                                         conf.method = "Wald") %>%
                filter(term == "periodrecent") %>%
                select(estimate, conf.low, conf.high, p.value))
  ) %>%
  select(-data) %>%
  unnest(est) %>%
  mutate(
    pct_change = (exp(estimate)  - 1) * 100,    # recent vs historical, typical site
    pct_lo     = (exp(conf.low)  - 1) * 100,
    pct_hi     = (exp(conf.high) - 1) * 100,
    # same scale as the GAM figure (compound % per year); midpoints are 5 years apart
    pct_per_yr_equiv = (exp(estimate / 5) - 1) * 100,
    significant = p.value < 0.05,
    direction = case_when(pct_lo > 0 ~ "Increasing",
                          pct_hi < 0 ~ "Decreasing",
                          TRUE       ~ "No clear trend"),
    direction = factor(direction, c("Increasing", "No clear trend", "Decreasing"))
  )

mm_tbl %>% select(-model) %>% print(n = Inf)

# ------------------------------------------------------------------
# d. Visualisation 
# ------------------------------------------------------------------
plot_trends <- mm_tbl %>%
  mutate(em_df = map(model, ~ as.data.frame(emmeans(.x, "period", type = "response")))) %>%
  select(-model) %>%
  unnest(em_df) %>%
  mutate(
    period    = factor(period, levels = c("historical", "recent")),
    sig_label = paste0(sprintf("%+.1f%%", pct_change), if_else(significant, " **", ""))
  )

x_labels <- c(historical = paste(range(historical_years), collapse = "–"),
              recent     = paste(range(recent_years),     collapse = "–"))

plot_region <- function(region) {
  d <- filter(plot_trends, cammlr_region == region)

  ggplot(d, aes(x = period, y = response,
                group = interaction(species, count_type), colour = direction)) +
    geom_line(linewidth = 2, alpha = 0.8) +
    geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1, colour = "grey30") +
    geom_point(size = 5) +
    geom_text(data = filter(d, period == "recent"), aes(label = sig_label),
              hjust = -0.2, vjust = 0.5, size = 5, fontface = "bold", show.legend = FALSE) +
    scale_x_discrete(labels = x_labels) +
    scale_y_continuous(labels = label_comma()) +
    scale_colour_manual(values = dir_cols, drop = FALSE) +
    facet_grid(count_type ~ species, scales = "free_y") +
    theme_antarctic() +
    theme(
      strip.text    = element_text(face = "bold", size = 18),
      axis.title    = element_text(face = "bold"),
      panel.spacing = unit(2, "lines"),
      legend.position = "none",
      plot.title    = element_text(size = 22, face = "bold"),
      plot.subtitle = element_text(size = 14, color = "gray40")
    ) +
    labs(
      title    = paste0("Penguin population trends: region ", region),
      subtitle = paste0("Model-adjusted counts for a typical site, ",
                        x_labels[["recent"]], " vs ", x_labels[["historical"]],
                        " (5 vs 5 years). '**' denotes p < 0.05"),
      y = "Estimated count per site (NB mixed model)",
      x = ""
    ) +
    expand_limits(x = 2.6)
}

all_plots <- map(set_names(unique(plot_trends$cammlr_region)), plot_region)

# ------------------------------------------------------------------
# e. Save: iterate through the list and save
# ------------------------------------------------------------------

for (region_name in names(all_plots)) {
  file_name <- paste0("penguins_recent_change_", region_name, ".svg")
  
  ggsave(
    filename = file_name,
    plot = all_plots[[region_name]],
    width = 14,
    height = 10,
    dpi = 300
  )
  
  message("Saved: ", file_name)
}
