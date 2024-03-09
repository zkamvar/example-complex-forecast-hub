library(dplyr)
library(readr)
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
          file = "target-data/time-series.csv")


# create auxiliary-data/locations.csv
locations <- locations_raw[, c("location", "abbreviation", "location_name",
                               "population")]

# add threholds on hospital admissions per 100,000 population defining intensity
# category boundaries
thresholds <- c(2.5, 5, 7.5)
for (i in seq_along(thresholds)) {
  locations[[paste0("threshold_", i)]] <-
    floor(locations[["population"]] / 100000 * thresholds[i])
}

write_csv(locations, file = "auxiliary-data/locations.csv")

# check with a plot
if (interactive()) {
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
}


# observed target values, suitable for use in scoring forecasts
tasks_config <- hubUtils::read_config(".", "tasks")

get_expanded_tasks_outputs <- function(target_block) {
  task_ids <- target_block$task_ids
  output_meta <- target_block$output_type

  tasks_expanded <- expand.grid(
      location = task_ids$location$optional,
      reference_date = as.Date(task_ids$reference_date$optional),
      horizon = as.integer(task_ids$horizon$optional),
      target = task_ids$target$optional) |>
    mutate(
      target_end_date = as.character(reference_date + 7 * horizon),
      .after = horizon
    )

  outputs_expanded <- purrr::map(
    names(output_meta),
    function(output_type) {
      output_type_meta <- output_meta[[output_type]]
      data.frame(
        output_type = output_type,
        output_type_id = output_type_meta$output_type_id |> unlist() |> unname()
      )
    }) |>
    purrr::list_rbind()

  cross_join(tasks_expanded, outputs_expanded)
}

t1 <- get_expanded_tasks_outputs(tasks_config$rounds[[1]]$model_tasks[[1]])
t2 <- get_expanded_tasks_outputs(tasks_config$rounds[[1]]$model_tasks[[2]])
t3 <- get_expanded_tasks_outputs(tasks_config$rounds[[1]]$model_tasks[[3]])

# bind then subset to reduce sensitivity to order of the tasks in the tasks.json
# file in case that changes later
target_rows <- rbind(t1, t2, t3)

obs_target_values_inc <- target_rows |>
  filter(target == "wk inc flu hosp") |>
  left_join(
    target_data_ts,
    by = join_by(location, target_end_date == date)
  )

obs_rate_cat <- target_data_ts |>
  left_join(locations, by = "location") |>
  mutate(
    target = "wk flu hosp rate category",
    output_type_id = case_when(
      value <= threshold_1 ~ "low",
      value <= threshold_2 ~ "moderate",
      value <= threshold_3 ~ "high",
      TRUE ~ "very high"
    ),
    value = 1
  )

obs_target_values_rate_cat <- target_rows |>
  filter(target == "wk flu hosp rate category") |>
  left_join(
    obs_rate_cat |> select(location, date, output_type_id, value),
    by = join_by(location, target_end_date == date, output_type_id)) |>
  mutate(
    value = ifelse(is.na(value), 0, 1)
  )

obs_rates <- target_data_ts |>
  left_join(locations, by = "location") |>
  mutate(rate = value / (population / 100000))

obs_target_values_rate <- target_rows |>
  filter(target == "wk flu hosp rate") |>
  left_join(
    obs_rates |> select(location, date, rate),
    by = join_by(location, target_end_date == date)) |>
  mutate(
    value = as.integer(as.numeric(output_type_id) >= rate)
  ) |>
  select(-rate)


target_data_complete <- bind_rows(
  obs_target_values_inc,
  obs_target_values_rate_cat,
  obs_target_values_rate)

write_csv(target_data_complete,
          file = "target-data/target-values-complete.csv")


target_data_distinct <- bind_rows(
  target_data_complete |>
    filter(target == "wk inc flu hosp") |>
    mutate(reference_date = "NA", horizon = "NA",
           output_type = "NA", output_type_id = "NA") |>
    distinct(),
  target_data_complete |>
    filter(target == "wk flu hosp rate category") |>
    mutate(reference_date = "NA", horizon = "NA") |>
    distinct(),
  target_data_complete |>
    filter(target == "wk flu hosp rate") |>
    mutate(reference_date = "NA", horizon = "NA") |>
    distinct()
)

write_csv(target_data_distinct,
          file = "target-data/target-values-distinct.csv")


target_data_obs_bin_only <- bind_rows(
  target_data_complete |>
    filter(target == "wk inc flu hosp") |>
    mutate(reference_date = "NA", horizon = "NA",
           output_type = "NA", output_type_id = "NA") |>
    distinct(),
  target_data_complete |>
    filter(target == "wk flu hosp rate category", value > 0) |>
    mutate(reference_date = "NA", horizon = "NA") |>
    distinct(),
  target_data_complete |>
    filter(target == "wk flu hosp rate", value > 0) |>
    mutate(reference_date = "NA", horizon = "NA",
           numeric_oti = as.numeric(output_type_id)) |>
    distinct() |>
    group_by(location, reference_date, horizon, target_end_date, target,
             output_type) |>
    slice_min(numeric_oti) |>
    ungroup() |>
    select(-numeric_oti)
)

write_csv(target_data_obs_bin_only,
          file = "target-data/target-values-obs-bin-only.csv")



# plot to double check
if (interactive()) {
  target_data_combined <- target_data_target_values |>
    left_join(target_data_cat_target_values_complete |> filter(value == 1) |> select(-target, -value))

  ggplot(
    data = target_data_combined |>
      dplyr::filter(location %in% c("US", "06", "25"))) +
    geom_point(mapping = aes(x = target_end_date, y = value,
                            color = output_type_id)) +
    facet_grid(rows = vars(horizon), cols = vars(location), scales = "free") +
    theme_bw()
}

# target_data <- dplyr::bind_rows(
#   target_data_target_values,
#   target_data_cat_target_values_complete
# ) |>
#   relocate("output_type_id", .before = "value")

# write_csv(target_data,
#           file = "target-data/target-values.csv")
