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
  filter(location %in% locations_raw$location) |>
  rename(observation = value)
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
      observation <= threshold_1 ~ "low",
      observation <= threshold_2 ~ "moderate",
      observation <= threshold_3 ~ "high",
      TRUE ~ "very high"
    ),
    observation = 1
  )

obs_target_values_rate_cat <- target_rows |>
  filter(target == "wk flu hosp rate category") |>
  left_join(
    obs_rate_cat |> select(location, date, output_type_id, observation),
    by = join_by(location, target_end_date == date, output_type_id)) |>
  mutate(
    observation = ifelse(is.na(observation), 0, 1)
  )

obs_rates <- target_data_ts |>
  left_join(locations, by = "location") |>
  mutate(rate = observation / (population / 100000))

obs_target_values_rate <- target_rows |>
  filter(target == "wk flu hosp rate") |>
  left_join(
    obs_rates |> select(location, date, rate),
    by = join_by(location, target_end_date == date)) |>
  mutate(
    observation = as.integer(as.numeric(output_type_id) >= rate)
  ) |>
  select(-rate)


target_data_complete <- bind_rows(
  obs_target_values_inc,
  obs_target_values_rate_cat,
  obs_target_values_rate)

# set output_type_id to "NA" for the quantile output_type,
# then remove duplicate rows created by that action
target_data_distinct <- bind_rows(
  target_data_complete |>
    filter(output_type == "quantile") |>
    mutate(output_type_id = "NA") |>
    distinct(),
  target_data_complete |>
    filter(output_type != "quantile")
)

# drop reference_date and horizon columns
target_data_distinct <- target_data_distinct |>
  select(-reference_date, -horizon) |>
  distinct()

write_csv(target_data_distinct,
          file = "target-data/target-observations.csv")
