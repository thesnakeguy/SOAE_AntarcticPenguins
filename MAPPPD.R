######################################################################################################
## This script draws penguin data from https://ipt.biodiversity.aq/resource?r=mapppd_count_data and 
## performs data visualization as well as statitical tests for (recent) population trends
######################################################################################################

library(MASS)
library(ggplot2)
library(tidyverse)
library(mgcv)
library(emmeans)

###################################
### Preparing data for analysis
###################################

# Draw data in
if (!file.exists("MAPPPD_AllCounts.csv")) {
  download.file(
    "https://www.penguinmap.com/mapppd/DownloadAll",
    "MAPPPD_AllCounts.csv",
    mode = "wb"
  )
}

pengs_raw <- readr::read_csv("MAPPPD_AllCounts.csv", show_col_types = FALSE)

# Filter and clean
pengs_clean <- pengs_raw %>%
  filter(
    cammlr_region %in% c("48.1", "88.1", "58.4.2"),
    common_name %in% c("adelie penguin", "chinstrap penguin", "gentoo penguin"),
    !is.na(penguin_count),
    penguin_count > 0  
  ) %>%
  mutate(
    species = case_when(
      str_detect(common_name, "adelie") ~ "Adelie",
      str_detect(common_name, "chinstrap") ~ "Chinstrap",
      str_detect(common_name, "gentoo") ~ "Gentoo"
    )
  )

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
    coverage = n_years / year_span,  # proportion of years monitored
    .groups = "drop"
  ) %>%
  arrange(desc(n_years))

# Focus on well-monitored sites (e.g., at least 5 years of data)
well_monitored <- monitoring_summary %>%
  filter(n_years >= 5)

# Join back to get analysis dataset
pengs_analysis <- pengs_clean %>%
  semi_join(well_monitored, by = c("cammlr_region", "site_name", "species"))

# Aggregate to site-year-species level (just to be sure, if a site would have been visited twice in one year)
pengs_yearly <- pengs_clean %>%
  group_by(cammlr_region, site_name, site_id, species, year, count_type) %>%
  summarize(
    mean_count = mean(penguin_count, na.rm = TRUE),
    max_count = max(penguin_count, na.rm = TRUE),  
    n_obs = n(),
    .groups = "drop"
  )

###################################
### Recent trends
###################################

### Fitting a negative binomial statistical model

# Statistical test (last 6-3 yrs vs last 3 yrs)
latest <- max(pengs_yearly$year)
sixyrsago <- as.numeric(latest - 6)
threeyrsago <- as.numeric(latest - 3)

recent_test_mixed <- pengs_yearly %>%
  filter(year >= sixyrsago) %>%
  mutate(period = ifelse(year >= threeyrsago, "recent", "historical")) %>%
  group_by(cammlr_region, species, count_type) %>%
  filter(n_distinct(period) == 2) %>%
  nest() %>%
  mutate(
    # Negative binomial mixed model with site as random effect
    model = map(data, ~glmer.nb(max_count ~ period + (1|site_id), data = .x)),
    
    # Extract fixed effect for period
    fixed_effects = map(model, ~summary(.x)$coefficients),
    p_value = map_dbl(fixed_effects, ~.["periodrecent", "Pr(>|z|)"]),
    coefficient = map_dbl(fixed_effects, ~.["periodrecent", "Estimate"]),
    
    n_sites = map_int(data, ~n_distinct(.x$site_id)),
    n_obs_recent = map_dbl(data, ~sum(.x$period == "recent")),
    n_obs_historical = map_dbl(data, ~sum(.x$period == "historical"))
  ) %>%
  dplyr::select(-data, -fixed_effects) %>%
  mutate(
    significant = p_value < 0.05,
    direction = ifelse(coefficient > 0, "increase", "decrease"),
    pct_change = (exp(coefficient) - 1) * 100
  )

print(recent_test_mixed, n = Inf)

### Visualization
# 1. Calculate adjusted means for all groups using your nested model results
# (Assuming your nested list 'recent_test_mixed' still has the 'model' column)
plot_trends <- recent_test_mixed %>%
  mutate(
    em_df = map(model, ~as.data.frame(emmeans(.x, "period", type = "response")))
  ) %>%
  unnest(em_df)
# Prepare the labels and significance markers
selected_region <- 48.1
plot_data_final <- plot_trends %>%
  filter(cammlr_region == selected_region) %>%
  mutate(
    # Create significance label
    sig_label = ifelse(significant, paste0(round(pct_change, 1), "% **"), 
                       paste0(round(pct_change, 1), "%")),
    # Order periods correctly on X-axis
    period = factor(period, levels = c("historical", "recent"))
  )

ggplot(plot_data_final, aes(x = period, y = response, group = interaction(species, count_type))) +
  geom_line(aes(color = species), linewidth = 2, alpha = 0.8) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1, color = "gray30") +
  geom_point(aes(color = species), size = 5) +
  geom_text(data = filter(plot_data_final, period == "recent"),
            aes(label = sig_label), 
            hjust = -0.2, vjust = 0.5, size = 5, fontface = "bold") +
  scale_x_discrete(labels = c(paste(sixyrsago, "-", threeyrsago), paste(threeyrsago, "-", latest))) +
  facet_grid(count_type ~ species, scales = "free_y") +
  theme_minimal(base_size = 16) + 
  theme(
    strip.text = element_text(face = "bold", size = 18),
    axis.title = element_text(face = "bold"),
    panel.spacing = unit(2, "lines"),
    legend.position = "none",
    plot.title = element_text(size = 22, face = "bold"),
    plot.subtitle = element_text(size = 14, color = "gray40")
  ) +
  labs(
    title = paste0("Penguin Population Trends: Region ", selected_region),
    subtitle = "Model-adjusted counts comparing last 10yrs vs last 3yrs. '**' denotes p < 0.05",
    y = "Estimated Mean Count (Mixed Model)",
    x = ""
  ) +
  expand_limits(x = 2.5) # Make room for the text labels
