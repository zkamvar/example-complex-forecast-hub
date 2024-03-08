library(fs)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(distfromq)
library(lubridate)

library(here)
setwd(here())

set.seed(42)

get_median_forecasts_from_q <- function(df) {
  result <- df |>
    dplyr::filter(abs(output_type_id - 0.5) < 0.0001) |>
    dplyr::mutate(
      output_type = "median",
      output_type_id = "NA"
    )

  return(result)
}

get_mean_forecasts_from_q <- function(df) {
  result <- df |>
    dplyr::filter(output_type == "quantile") |>
    dplyr::group_by(location, reference_date, horizon, target_end_date,
                    target) |>
    dplyr::summarize(
      value = distfromq::make_r_fn(ps = output_type_id, qs = value)(1e5) |>
                mean(),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      output_type = "mean",
      output_type_id = "NA"
    )

  return(result)
}

get_cat_probs <- function(output_type_id, value, threshold_1, threshold_2,
                          threshold_3) {
  p_fn <- distfromq::make_p_fn(ps = as.numeric(output_type_id),
                               qs = as.numeric(value))
  boundary_probs <- c(0, p_fn(c(threshold_1, threshold_2, threshold_3)), 1)
  cat_probs <- diff(boundary_probs)
  cat_probs <- pmin(cat_probs, 1)
  cat_probs <- pmax(cat_probs, 0)
  cat_probs <- cat_probs / sum(cat_probs)
  return(list(cat_probs))
}

get_cat_forecasts_from_q <- function(df, locations) {
  result <- df |>
    left_join(locations, by = "location") |>
    group_by(location, reference_date, horizon, target_end_date) |>
    summarize(
      target = "wk flu hosp rate category",
      value = get_cat_probs(output_type_id, value,
                            threshold_1[1], threshold_2[1], threshold_3[1]),
      output_type = "pmf",
      output_type_id = list(c("low", "moderate", "high", "very high")),
      .groups = "drop"
    ) |>
    unnest(cols = c(output_type_id, value)) |>
    select(location, reference_date, horizon, target_end_date, target,
           output_type, output_type_id, value)

  return(result)
}

get_cdf_values <- function(output_type_id, value, x) {
  p_fn <- distfromq::make_p_fn(ps = as.numeric(output_type_id),
                               qs = as.numeric(value))
  return(list(p_fn(x)))
}

get_rate_cdf_forecasts_from_q <- function(df, locations) {
  cdf_x_values <- seq(from = 0.25, to = 25, by = 0.25)
  result <- df |>
    left_join(locations, by = "location") |>
    mutate(rate = value / (population / 100000)) |>
    group_by(location, reference_date, horizon, target_end_date) |>
    summarize(
      target = "wk flu hosp rate",
      value = get_cdf_values(output_type_id, rate, cdf_x_values),
      output_type = "cdf",
      output_type_id = list(as.character(cdf_x_values)),
      .groups = "drop"
    ) |>
    unnest(cols = c(output_type_id, value)) |>
    select(location, reference_date, horizon, target_end_date, target,
           output_type, output_type_id, value)

  return(result)
}


locations <- read_csv("auxiliary-data/locations.csv")

files <- Sys.glob("internal-data-raw/model-output-orig/*/*")

for (f in files) {
  df <- read.csv(f) |>
    transmute(
      location = location,
      reference_date = as.Date(forecast_date) + 5,
      target_end_date = as.Date(target_end_date),
      target = "wk inc flu hosp",
      output_type = type,
      output_type_id = quantile,
      value = round(value)
    ) |>
    mutate(
      horizon = as.integer((target_end_date - reference_date) / 7),
      .before = target_end_date
    ) |>
    filter(
      location %in% locations$location,
      output_type == "quantile"
    )

  df <- dplyr::bind_rows(
    df %>% dplyr::mutate(output_type_id = as.character(output_type_id)),
    get_median_forecasts_from_q(df),
    get_mean_forecasts_from_q(df),
    get_cat_forecasts_from_q(df, locations),
    get_rate_cdf_forecasts_from_q(df, locations)
  )

  path_parts <- fs::path_split(f)[[1]][2:4]
  path_parts[[1]] <- "model-output"
  path_parts[[3]] <- paste0(df$reference_date[1], "-",
                            path_parts[[2]], ".csv")
  dir_out <- fs::path_join(path_parts[1:2])
  if (!dir.exists(dir_out)) {
    dir.create(dir_out, recursive = TRUE)
  }
  f_out <- fs::path_join(path_parts)
  write_csv(df, file = f_out)
}
