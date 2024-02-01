library(fs)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(hubUtils)

library(here)
setwd(here())

# read in "starting point" data sets
target_data_raw <- read.csv(
  "internal-data-raw/target-data-orig/truth-Incident Hospitalizations.csv")

locations_raw <- read.csv(
  "internal-data-raw/auxiliary-data-orig/locations.csv")


# save target data in time series format
target_data_ts <- target_data_raw[, c("date", "location", "value")] |>
  filter(location %in% locations_raw$location)
write_csv(target_data_ts,
          file = "target-data/flu-hospitalization-time-series.csv")


# create auxiliary-data/locations.csv

# Exploration used to identify thresholds in terms of rate per 100k population
locations <- locations_raw[, c("location", "abbreviation", "location_name",
                               "population")]

max_obs <- target_data_ts |>
  filter(location == "US") |>
  pull(value) |>
  max()

us_pop <- locations |>
  filter(location == "US") |>
  pull(population)
max_rate <- max_obs / (us_pop / 100000)

# max_rate is 7.92685
# for illustrative purposes, we will set thresholds at 2.5, 5, and 7.5
# hospital admissions per 100,000 population.
thresholds <- c(2.5, 5, 7.5)
for (i in seq_along(thresholds)) {
  locations[[paste0("threshold_", i)]] <-
    floor(locations[["population"]] / 100000 * thresholds[i])
}

# check with a plot
ggplot() +
  geom_line(
    data = target_data_ts |>
      mutate(date = as.Date(date)) |>
      filter(date >= "2022-10-01", date <= "2023-07-01") |>
      left_join(locations, by = "location") |>
      filter(!is.na(abbreviation)),
    mapping = aes(x = date, y = value)
  ) +
  geom_hline(
    data = locations,
    mapping = aes(yintercept = threshold_1)) +
  geom_hline(
    data = locations,
    mapping = aes(yintercept = threshold_2)) +
  geom_hline(
    data = locations,
    mapping = aes(yintercept = threshold_3)) +
  facet_wrap( ~ abbreviation, scales = "free_y") +
  theme_bw()

# save
write_csv(locations,
          file = "auxiliary-data/locations.csv")


# observed target values, suitable for use in scoring forecasts
tasks_config <- hubUtils::read_config(".", "tasks")
task_ids <- tasks_config$rounds[[1]]$model_tasks[[1]]$task_ids

target_data_target_values <- expand.grid(
    location = task_ids$location$optional,
    reference_date = as.Date(task_ids$reference_date$optional),
    horizon = as.integer(task_ids$horizon$optional)) |>
  mutate(
    target_end_date = as.character(reference_date + 7 * horizon)
  ) |>
  left_join(
    target_data_ts,
    by = join_by(location, target_end_date == date)
  ) |>
  arrange(location, reference_date, horizon)

write_csv(target_data_target_values,
          file = "target-data/flu-hospitalization-value-target-values.csv")


target_data_cat_target_values <- target_data_target_values |>
  left_join(locations, by = "location") |>
  mutate(
    output_type_id = case_when(
      value <= threshold_1 ~ "low",
      value <= threshold_2 ~ "moderate",
      value <= threshold_3 ~ "high",
      TRUE ~ "very high"
    ),
    value = 1
  ) |>
  select(location, reference_date, horizon, target_end_date, output_type_id,
         value)

# plot to double check
target_data_combined <- target_data_target_values |>
  left_join(target_data_cat_target_values |> select(-value))

ggplot(
  data = target_data_combined |>
    dplyr::filter(location %in% c("US", "06", "25"))) +
  geom_point(mapping = aes(x = target_end_date, y = value,
                           color = output_type_id)) +
  facet_grid(rows = vars(horizon), cols = vars(location), scales = "free") +
  theme_bw()

write_csv(target_data_cat_target_values,
          file = "target-data/flu-hospitalization-category-target-values.csv")
