require(tidyverse)
require(knitr)
require(kableExtra)
require(magrittr)
require(stargazer)
require(MASS, exclude = "select")
require(DHARMa)

rm(list = ls())

set.seed(1234)

## data and summary

save_table_tex <- function(df, caption, label, filename, first_italic = T, digits = 4, dir = "output", add_layers = identity, col_names = NA) {
  filepath <- file.path(dir, sprintf("%s.tex", filename))
  add_layers <- rlang::as_function(add_layers)
  
  if (!all(is.na(col_names))) {
    assertthat::assert_that(length(col_names) == ncol(df))
  }
  
  df |> 
    kable(
      digits = digits, 
      format = "latex", 
      caption = caption, 
      label = label,
      col.names = col_names
    ) |>
    kable_styling(position = "center", latex_options = "HOLD_position") |>
    row_spec(0, hline_after = T) |>
    column_spec(1, italic = first_italic) |> 
    add_layers() |>
    write_file(filepath)
}

save_last_plot_eps <- function(filename, dir = "output") {
  captured_plot <- recordPlot()

  setEPS()
  file.path(dir, sprintf("%s.eps", filename)) |> postscript()
  replayPlot(captured_plot)
  dev.off()
}

save_output <- function(expr, filename, dir = "output") {
  filepath <- file.path(dir, sprintf("%s.txt", filename))
  captured_output <- capture.output(expr) %T>% cat(sep = "\n")
  
  # captured_output[5:length(captured_output)] |> 
  captured_output |>
    write_lines(filepath)
}

####

df <- read_table("data/16-victim.txt", skip = 21, col_names = c("index", "resp", "race")) |>
  select(-index) |>
  mutate(race = as.factor(str_remove_all(race, '"')))

df |> 
  group_by(race) |>
  summarise(mu = mean(resp), sigma = sd(resp), count = n()) %T>%
  save_table_tex(
    caption = "Dataset summary statistics",
    label = "dataset_summary",
    filename = "dataset_summary"
  ) |>
  kable(format = "pipe")

df_crosscount <- df |> 
  group_by(race, resp) |>
  summarise(count = n())

df_crosscount %T>%
  save_table_tex(
    caption = "Dataset cross frequencies",
    label = "dataset_crosscount",
    filename = "dataset_scrosscount"
  ) |>
  kable(format = "pipe")

## poisson reg

poisson_model <- glm(resp ~ race, family = "poisson", data = df)
summary(poisson_model) |>
  save_output("poisson_summary")

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
  ) |> 
  select(
    resp, black, black_pred, white, white_pred
  ) %T>%
  save_table_tex(
    caption = "Observed vs predicted counts",
    label = "obs_vs_pred",
    filename = "obs_vs_pred",
    col_names = c("victims", "obs", "pred", "obs", "red"),
    add_layers = ~ add_header_above(.x, c(" " = 1, "Black" = 2, "White" = 2))
  ) |>
  kable(format = "pipe", digits = 2)

## tests

# poisson_model$df.residual
# summary(poisson_model)$df.residual
# dim(df)[1] - length(coef(poisson_model))

# poisson_model$deviance
# summary(poisson_model)$deviance

sum(residuals(poisson_model, type = "pearson")^2) %>% 
  c(., pchisq(., df = poisson_model$df.residual, lower.tail = F)) |>
  set_names(c("test statistic", "p-value"))

poisson_model$deviance %>% 
  c(., pchisq(., df = poisson_model$df.residual, lower.tail = F)) |>
  set_names(c("test statistic", "p-value"))

anova(poisson_model, test = "Chisq")

poisson_model |> car::Anova(test = "LR", type = 3)
poisson_model |> car::Anova(test = "Wald", type = 3)

## DHARMa

sim_resid_poisson <- poisson_model |> 
  simulateResiduals(n = 1000, plot = T, seed = 1234)
  
save_last_plot_eps("poisson_resid_diag_dharma")

sim_resid_poisson |> 
  testOutliers(type = "bootstrap", nBoot = 100) |> 
  save_output("poisson_resid_outlier_test")

save_last_plot_eps("poisson_resid_outlier_test_plot")

sim_resid_poisson |> 
  testDispersion() |> 
  save_output("poisson_resid_dispersion_test")

save_last_plot_eps("poisson_resid_dispersion_test_plot")

sim_resid_poisson |> 
  testUniformity(plot = F) |>
  save_output("poisson_resid_uniformity_test")

## tests agregated

bind_rows(
  sum(residuals(poisson_model, type = "pearson")^2) %>% tibble(
    statistic = .,
    pval = pchisq(., df = poisson_model$df.residual, lower.tail = F)
  ),
  anova(poisson_model, test = "Chisq")["race", c("Resid. Dev", "Pr(>Chi)")] |> 
    set_names(c("statistic", "pval")),
  car::Anova(poisson_model, test = "LR", type = 3)["race", c("LR Chisq", "Pr(>Chisq)")] |> 
    set_names(c("statistic", "pval")),
  car::Anova(poisson_model, test = "Wald", type = 3)["race", c("Chisq", "Pr(>Chisq)")] |> 
    set_names(c("statistic", "pval")),
  testOutliers(sim_resid_poisson, type = "bootstrap", nBoot = 100, plot = F)[c("estimate", "p.value")] |> 
    set_names(c("statistic", "pval")),
  testDispersion(sim_resid_poisson, plot = F)[c("statistic", "p.value")] |> 
    set_names(c("statistic", "pval")),
  testUniformity(sim_resid_poisson, plot = F)[c("statistic", "p.value")] |> 
    set_names(c("statistic", "pval"))
) |>
  mutate(
    test = c("Pearson", "Deviance", "LR", "Wald", "Bootstrap Outliers", "Dispersion", "K-S Uniformity"),
    .before = 1
  ) %T>% 
  save_table_tex(
    "Poisson regression test",
    "poisson_reg_tests",
    "poisson_reg_tests"
  ) |>
  kable(format = "pipe", digits = 4)

## rootogram

countreg::rootogram(poisson_model, main = "Poisson regression", ylab = "Square root of frequency")
save_last_plot_eps("rootgram_poisson")

## negative binomial model

neg_bin_model <- glm.nb(resp ~ race, data = df)
summary(neg_bin_model)

## fitted counts for Negative Binomial GLM

fitted_means_nb <- neg_bin_model |> 
  predict(newdata = tibble(race = c("white", "black"))) |>
  exp() |>
  set_names(c("white", "black"))

fitted_means_nb 

## dispersion parameter

neg_bin_model$theta

## estimated variance for the count

fitted_means_nb + fitted_means_nb^2 * (1/neg_bin_model$theta)

## quasi-likelihood model

quasilik_model <- glm(resp ~ race, data = df, family = quasipoisson)
summary(quasilik_model)

## model comparison

anova(poisson_model, neg_bin_model, quasilik_model)

stargazer(poisson_model, neg_bin_model, quasilik_model, type = "text")

stargazer(poisson_model, neg_bin_model, quasilik_model, type = "latex") |>
  capture.output(file = "output/stargazer_comparison.tex")
