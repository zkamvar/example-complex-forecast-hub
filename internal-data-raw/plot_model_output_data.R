library(hubUtils)
library(hubData)
library(dplyr)
library(ggplot2)
library(here)
setwd(here())

df <- connect_hub(".") |>
  collect()

df_sub <- df |> filter(
  model_id == "Flusight-baseline",
  location == "US")

ggplot() +
  geom_line(
    data = df_sub |>
      filter(output_type == "cdf") |>
      mutate(output_type_id = as.numeric(output_type_id)),
    mapping = aes(x = output_type_id, y = value, group = model_id)
  ) +
  facet_grid(vars(reference_date), vars(horizon))
