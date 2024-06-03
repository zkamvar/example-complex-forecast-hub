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

#' apply the schaake shuffle to matrices of simulated and observed values
schaake_shuffle <- function(sim, obs) {
  stopifnot(all(dim(sim) == dim(obs)))

  # apply shuffle, iterating over columns. See Eq. (9) - (12) in
  # Clark, Martyn, et al. "The Schaake shuffle: A method for reconstructing
  # space–time variability in forecasted precipitation and temperature fields."
  # Journal of Hydrometeorology 5.1 (2004): 243-262.
  result <- matrix(NA, nrow = nrow(sim), ncol = ncol(sim))
  for (j in seq_len(ncol(sim))) {
    result[order(obs[, j]), j] <- sort(sim[, j])
  }

  return(result)
}


#' Apply the Schaake shuffle to vector of simulated values x, with the
#' ground truth being draws from an AR(rho) process.
#'
#' This is more like using a Gaussian copula than the Schaake shuffle.
schaake_shuffle_ar <- function(value, horizon, sample_index, rho) {
  H <- max(horizon)
  n <- max(sample_index)

  # construct x, a matrix of provided values
  x <- matrix(NA, nrow = n, ncol = H)
  x[cbind(sample_index, horizon)] <- value

  # simulate n replicates from an AR(rho) process,
  # H observations for each replicate
  innovations <- matrix(rnorm(n * H, sd = 1), nrow = n, ncol = H)
  y_0 <- rnorm(n, mean = 0, sd = sqrt(1 / (1 - rho^2)))
  y <- matrix(NA, nrow = n, ncol = H)
  y[, 1] <- rho * y_0 + innovations[, 1]
  for (i in seq(from = 2, to = H)) {
    y[, i] <- rho * y[, i - 1] + innovations[, i]
  }

  # apply Schaake shuffle
  x_ordered <- schaake_shuffle(sim = x, obs = y)

  # return values from x_ordered, matching input horion and sample_index
  return(x_ordered[cbind(sample_index, horizon)])
}

get_sample_forecasts_from_q <- function(df, n = 100, rho = 0.9) {
  # set up a data frame with marginal samples for each date/location,
  # joint across horizon
  df <- df |>
    dplyr::filter(output_type == "quantile") |>
    dplyr::group_by(location, reference_date, horizon, target_end_date,
                    target) |>
    # draw samples from marginal distributions by location/date/horizon
    dplyr::summarize(
      value = list(
        distfromq::make_r_fn(ps = output_type_id, qs = value)(n) |>
          round()
      ),
      sample_index = list(seq_len(n)),
      .groups = "drop"
    ) |>
    tidyr::unnest(cols = c(value, sample_index)) |>
    dplyr::group_by(location, reference_date) |>
    # induce dependence across horizons for each location/date combo
    dplyr::mutate(
      value = schaake_shuffle_ar(value, horizon + 1, sample_index, rho)
    ) |>
    dplyr::ungroup() |>
    # get to desired columns: output_type, output_type_id, no sample_index
    dplyr::mutate(
      output_type = "sample",
      output_type_id = as.character(
        sample_index +
          n * (as.integer(factor(paste0(location, reference_date))) - 1)
      )
    ) |>
    dplyr::select(-sample_index) |>
    dplyr::ungroup()

  return(df)
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
           output_type, output_type_id, value) |>
    mutate(value = pmax(pmin(value, 1), 0))

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
    get_sample_forecasts_from_q(df),
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
