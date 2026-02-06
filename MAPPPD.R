library(dplyr)
######################################################################################################
## This script draws penguin data from https://ipt.biodiversity.aq/resource?r=mapppd_count_data and 
## performs data visualization as well as statitical tests for (recent) population trends
######################################################################################################


library(ggplot2)
library(tidyverse)
library(mgcv)

###################################
### Preparing data for analysis
###################################


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
### Fitting the data to a GAM model
###################################

### 1. Complete overview of the trends in every site

site_trends_gam <- pengs_yearly %>%
  group_by(cammlr_region, site_name, species, count_type) %>%
  filter(n_distinct(year) >= 5) %>%
  nest() %>% # nest allows you to create a database with a list column
  mutate(
    model_gam = map(data, ~gam(max_count ~ s(year, k = 4), # k equals the number of smoothing splines
                               data = .x, 
                               family = nb(),
                               method = "REML")),
    gam_summary = map(model_gam, summary),
    p_value_gam = map_dbl(gam_summary, ~.x$s.table[1, "p-value"]),
    edf = map_dbl(gam_summary, ~.x$s.table[1, "edf"]),
    
    # Determine trend direction from predictions
    predictions = map2(model_gam, data, ~{
      years <- range(.y$year)
      pred_start <- predict(.x, newdata = data.frame(year = years[1]), type = "response")
      pred_end <- predict(.x, newdata = data.frame(year = years[2]), type = "response")
      tibble(start = pred_start, end = pred_end, change = pred_end - pred_start)
    })
  ) %>%
  unnest(predictions) %>%  # Now unnest to get start, end, change
  mutate(
    pct_change = (change / start) * 100,
    trend_direction = case_when(
      p_value_gam >= 0.05 ~ "stable",
      change > 0 ~ "increasing",
      change < 0 ~ "decreasing"
    ),
    is_nonlinear = edf > 1.5
  ) %>%
  select(-data, -model_gam, -gam_summary)  # Clean up list columns

# View results
head(site_trends_gam)

### 2. Summary by region and species

trend_summary <- site_trends_gam %>%
  group_by(cammlr_region, species, count_type, trend_direction) %>%
  summarize(
    n_sites = n(),
    n_nonlinear = sum(is_nonlinear),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = trend_direction, values_from = c(n_sites, n_nonlinear), values_fill = 0)

print(trend_summary)

# Visualization
ggplot(pengs_yearly, aes(x = year, y = max_count, color = species)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "gam", 
              method.args = list(family = nb()),  # negative binomial
              formula = y ~ s(x, k = 5),  # smoothing spline
              se = TRUE, 
              linewidth = 1.2) +
  facet_wrap(~cammlr_region, scales = "free_y") +
  scale_y_log10() +
  theme_minimal() +
  labs(
    title = "Penguin Population Trends by CCAMLR Region",
    subtitle = "GAM smoothed trends with negative binomial family",
    y = "Penguin Count (log scale)",
    x = "Year",
    color = "Species"
  ) +
  theme(legend.position = "bottom")

### 3. Individual site trajectories per region

ggplot(pengs_yearly %>% filter(cammlr_region == "48.1"), 
       aes(x = year, y = max_count, group = site_name)) +
  geom_line(alpha = 0.4) +
  geom_smooth(aes(group = species, color = species), 
              method = "gam",
              method.args = list(family = nb()),
              formula = y ~ s(x, k = 4),
              se = TRUE, 
              linewidth = 1.5) +
  facet_wrap(~species, scales = "free_y") +
  theme_bw() +
  labs(title = "Site-level GAM trends in CCAMLR 48.1",
       subtitle = "Allows for non-linear trends")

# # Alternative: Show which sites have non-linear vs linear trends
# ggplot(pengs_yearly %>% 
#          semi_join(site_trends_gam, by = c("cammlr_region", "site_name", "species", "count_type")),
#        aes(x = year, y = max_count)) +
#   geom_point(alpha = 0.3) +
#   geom_smooth(data = . %>% 
#                 inner_join(site_trends_gam %>% filter(is_nonlinear), 
#                            by = c("cammlr_region", "site_name", "species", "count_type")),
#               method = "gam",
#               method.args = list(family = nb()),
#               formula = y ~ s(x, k = 4),
#               color = "blue",
#               se = TRUE) +
#   geom_smooth(data = . %>% 
#                 inner_join(site_trends_gam %>% filter(!is_nonlinear), 
#                            by = c("cammlr_region", "site_name", "species", "count_type")),
#               method = "lm",
#               color = "red",
#               linetype = "dashed",
#               se = TRUE) +
#   facet_wrap(~cammlr_region + species, scales = "free_y") +
#   theme_minimal() +
#   labs(title = "Blue = Non-linear trends, Red dashed = Linear trends")

### 4. Recent period 

pengs_yearly %>%
  mutate(period = ifelse(year >= 2015, "2015-2025", "Before 2015")) %>%
  group_by(cammlr_region, species, period) %>%
  summarize(median_count = median(max_count, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = species, y = median_count, fill = period)) +
  geom_col(position = "dodge") +
  facet_wrap(~cammlr_region, scales = "free_y") +
  theme_minimal() +
  labs(title = "Recent vs Historical Median Counts",
       y = "Median Count", fill = "Period")

# Statistical test (last 10 yrs vs last 3 yrs)
recent_test <- pengs_yearly %>%
  filter(year >= 2014) %>%
  mutate(period = ifelse(year >= 2022, "recent", "historical")) %>%
  group_by(cammlr_region, species) %>%
  filter(n_distinct(period) == 2) %>%
  nest() %>%
  mutate(
    test = map(data, ~wilcox.test(max_count ~ period, data = .x)),
    p_value = map_dbl(test, "p.value"),
    statistic = map_dbl(test, "statistic"),
    n_recent = map_dbl(data, ~sum(.x$period == "recent")),
    n_historical = map_dbl(data, ~sum(.x$period == "historical"))
  ) %>%
  select(-data, -test) %>%
  mutate(significant = p_value < 0.05)

print(recent_test)
