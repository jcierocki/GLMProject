require(tidyverse)
require(knitr)
require(kableExtra)
require(magrittr)
require(stargazer)

rm(list = ls())

## data and summary

save_tex <- function(df, caption, label, filename, first_italic = T, dir = "output") {
  filepath <- file.path(dir, sprintf("%s.tex", filename))
  
  df |> 
    kable(digits = 2, format = "latex", caption = caption, label = label) |>
    kable_styling(position = "center", latex_options = "HOLD_position") |>
    row_spec(0, bold = TRUE, background = "#ffd966") |>
    column_spec(1, italic = first_italic, border_left = T) |>
    column_spec(ncol(df), border_right = T) |>
    write_file(filepath)
}

save_eps <- function(plotting_expr, filename, dir = "output") {
  setEPS()
  file.path(dir, sprintf("%s.eps", filename)) |> postscript()
  rlang::enquo(plotting_expr) |> rlang::eval_tidy()
  dev.off()
}

df <- read_table("data/16-victim.txt", skip = 21, col_names = c("index", "resp", "race")) |>
  select(-index) |>
  mutate(race = as.factor(str_remove_all(race, '"')))

df |> 
  group_by(race) |>
  summarise(mu = mean(resp), sigma = sd(resp), count = n()) %T>%
  save_tex(
    caption = "Dataset summary statistics",
    label = "dataset_summary",
    filename = "dataset_summary"
  ) |> 
  kable(format = "pipe")

df_crosscount <- df |> 
  group_by(race, resp) |>
  summarise(count = n())

df_crosscount %T>%
  save_tex(
    caption = "Dataset cross frequencies",
    label = "dataset_crosscount",
    filename = "dataset_scrosscount"
  ) |> 
  kable(format = "pipe")

## poisson reg

poisson_model <- glm(resp ~ race, family = "poisson", data = df)
summary(poisson_model)

##

exp(coef(poisson_model)[2])
mean(df$resp[df$race == "black"])/mean(df$resp[df$race == "white"])

##

fitted_means <- poisson_model |> 
  predict(newdata = tibble(race = c("white", "black"))) |>
  exp() |>
  set_names(c("white", "black"))

fitted_means

df_crosscount |> 
  pivot_wider(names_from = "race", values_from = "count", values_fill = 0) |>
  mutate(
    black_pred = dpois(resp, lambda = fitted_means[1]) * sum(df$race =="black"),
    white_pred = dpois(resp, lambda = fitted_means[1]) * sum(df$race =="white"),
  ) %T>%
  save_tex(
    caption = "Observed vs predicted counts",
    label = "obs_vs_pred",
    filename = "obs_vs_pred"
  ) |>
  kable(format = "pipe", digits = 2)

## rootogram

countreg::rootogram(poisson_model)
countreg::rootogram(poisson_model) |>
  save_eps("rootgram_poisson")
