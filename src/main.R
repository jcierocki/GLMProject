require(tidyverse)
require(knitr)
require(kableExtra)
require(magrittr)
require(stargazer)
require(MASS, exclude = "select")
require(DHARMa)

rm(list = ls())

set.seed(1234)

options(knitr.table.format = "pipe", digits = 4)

## functions

save_table_tex <- function(df, caption, label, filename = NULL, first_italic = T, digits = 4, dir = "output", add_layers = identity, col_names = NA) {
  if (is.null(filename)) {
    filename <- label
  }
  
  filepath <- file.path(dir, sprintf("%s.tex", filename))
  add_layers <- rlang::as_function(add_layers)
  
  if (!all(is.na(col_names))) {
    assertthat::assert_that(length(col_names) == ncol(df))
  }
  
  df |> 
    kbl(
      digits = digits, 
      format = "latex", 
      caption = caption, 
      label = label,
      col.names = col_names,
      booktabs = T
    ) |>
    kable_styling(position = "center", latex_options = c("striped")) |>
    # row_spec(0, hline_after = T) |>
    column_spec(1, italic = first_italic) |> 
    add_layers() |>
    save_kable(filepath)
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

glm_tests_combined <- function(glm_model, sim_resid = NULL, save_plots_model_name = NULL) {
  tests_df <- bind_rows(
    sum(residuals(glm_model, type = "pearson")^2) %>% tibble(
      statistic = .,
      pval = pchisq(., df = glm_model$df.residual, lower.tail = F)
    ),
    anova(glm_model, test = "Chisq")["race", c("Resid. Dev", "Pr(>Chi)")] |> 
      set_names(c("statistic", "pval")),
    car::Anova(glm_model, test = "LR", type = 3)["race", c("LR Chisq", "Pr(>Chisq)")] |> 
      set_names(c("statistic", "pval")),
    car::Anova(glm_model, test = "Wald", type = 3)["race", c("Chisq", "Pr(>Chisq)")] |> 
      set_names(c("statistic", "pval")),
  ) |>
    mutate(
      test = c("Pearson", "Deviance", "LR", "Wald"),
      .before = 1
    )
  
  save_dharma_plots <- !is.null(save_plots_model_name)
  
  if (is.null(sim_resid)) {
    sim_resid <- tryCatch({
      glm_model |> simulateResiduals(n = 1000, plot = save_dharma_plots, seed = 1234)
    }, error = function(e) NULL)
  }
  
  if (!is.null(sim_resid)) {
    save_plots_model_name <- save_plots_model_name |>
      str_to_lower() |> 
      str_remove_all("[:punct:]") |> 
      str_replace_all("\\s", "_")
    
    sprintf("%s_reg_simulated_residuals_plot", save_plots_model_name) |> 
      save_last_plot_eps()
    
    dharma_test_funs <- list(
      "Bootstrap Outliers" = ~ testOutliers(.x, type = "bootstrap", nBoot = 100, plot = save_dharma_plots)[c("estimate", "p.value")], 
      "Dispersion" = ~ testDispersion(.x, plot = save_dharma_plots)[c("statistic", "p.value")], 
      "K-S Uniformity" = ~ testUniformity(.x, plot = save_dharma_plots)[c("statistic", "p.value")], 
      "Zero Inflation" = ~ testZeroInflation(.x, plot = save_dharma_plots)[c("statistic", "p.value")]
    )
    
    dharma_plots_filenames = names(dharma_test_funs) |>
      str_to_lower() |> 
      str_remove_all("[:punct:]") |> 
      str_replace_all("\\s", "_") %>%
      sprintf("%s_reg_%s_plot", save_plots_model_name, .)
    
    dharma_tests_results_df <- dharma_test_funs |> 
      set_names(dharma_plots_filenames) |>
      imap_dfr(~ {
        out <- rlang::as_function(.x)(sim_resid) |>
          set_names(c("statistic", "pval"))
        
        if (save_dharma_plots) {
          save_last_plot_eps(.y)
        }
        
        out
      }) |>
        mutate(
          test = names(dharma_test_funs),
          .before = 1
        )
    
    tests_df <- bind_rows(
      tests_df,
      dharma_tests_results_df
    )
  }
  
  tests_df
}

glm_RR_table <- function(glm_model) { ## only poisson supported now
  estimates <- stats::coef(glm_model)
  
  bind_cols(
    parameter = names(estimates),
    RR = estimates,
    stats::confint(glm_model)
  ) |> mutate_if(is.numeric, exp)
}

## data and summary

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

glm_RR_table(poisson_model) %T>%
  save_table_tex(
    "Poisson Regression Risk Ratios",
    "poisson_reg_RR"
  ) |> 
  kable(format = "pipe", digits = 2)

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

# sum(residuals(poisson_model, type = "pearson")^2) %>% 
#   c(., pchisq(., df = poisson_model$df.residual, lower.tail = F)) |>
#   set_names(c("test statistic", "p-value"))
# 
# poisson_model$deviance %>% 
#   c(., pchisq(., df = poisson_model$df.residual, lower.tail = F)) |>
#   set_names(c("test statistic", "p-value"))
# 
# anova(poisson_model, test = "Chisq")
# 
# poisson_model |> car::Anova(test = "LR", type = 3)
# poisson_model |> car::Anova(test = "Wald", type = 3)

## tests inc. DHARMa

glm_tests_combined(poisson_model, save_plots_model_name = "poisson") %T>%
  save_table_tex(
    "Poisson regression tests results",
    "poisson_reg_tests"
  ) |>
  kable(format = "pipe", digits = 4)

## rootogram

countreg::rootogram(poisson_model, main = "Poisson regression", ylab = "Square root of frequency")
save_last_plot_eps("rootgram_poisson")

## negative binomial model

neg_bin_model <- glm.nb(resp ~ race, data = df)
summary(neg_bin_model)

glm_RR_table(neg_bin_model) %T>%
  save_table_tex(
    "Negative Binomial Regression Risk Ratios",
    "neg_bin_reg_RR"
  ) |> 
  kable(format = "pipe", digits = 2)

glm_tests_combined(neg_bin_model, save_plots_model_name = "negative binomial") %T>% 
  save_table_tex(
    "Negative Binomial regression test",
    "negbin_reg_tests",
    "negbin_reg_tests"
  ) |>
  kable(format = "pipe", digits = 4)

countreg::rootogram(neg_bin_model, main = "Negative Binomial regression", ylab = "Square root of frequency")
save_last_plot_eps("rootgram_negbin")

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

glm_RR_table(quasilik_model) %T>%
  save_table_tex(
    "Quasi Poisson Regression Risk Ratios",
    "quasipoisson_reg_RR"
  ) |> 
  kable(digits = 2)

glm_tests_combined(quasilik_model, save_plots_model_name = "quasi poisson") %T>% 
  save_table_tex(
    "Quasi Poisson regression test",
    "quasipoisson_reg_tests"
  ) |>
  kable()

## model comparison

anova(poisson_model, neg_bin_model, quasilik_model) |>
  save_output("anova_model_comparison")

stargazer(poisson_model, neg_bin_model, quasilik_model, type = "text")

stargazer(poisson_model, neg_bin_model, quasilik_model, type = "latex") |>
  capture.output(file = "output/stargazer_comparison.tex")
