library(ggplot2)
library(patchwork)
library(purrr)
library(furrr)
library(densityratio)
library(dplyr)
library(synthpop)



plan(multisession, workers = future::availableCores() - 1)

set.seed(123)


dlaplace <- function(x, mu = 0, sd = 1) exp(-abs(x-mu)/(sd / sqrt(2))) / (2*(sd / sqrt(2))) 
rlaplace <- function(n, mu = 0, sd = 1) {
  p <- runif(n)
  b <- sd / sqrt(2)
  mu - b * sign(p - 0.5) * log(1 - 2*abs(p - 0.5))
}
dratio_lap_norm <- function(x, mu = 0, sd = 1) {
  dnorm(x, mu, sd) / dlaplace(x, mu, sd)
}

dratio_lnorm_norm <- function(x, mu = 0, sd = 1) {
  mean_rlnorm <- log(mu^2 / sqrt(mu^2 + sd^2))
  sd_rlnorm <- sqrt(log(1 + sd^2 / mu^2))
  dnorm(x, mu, sd) / dlnorm(x, mean_rlnorm, sd_rlnorm)
}

dratio_t_norm <- function(x, mu = 0, sd = 1) {
  df <- 2 / (1 - 1/sd^2)
  dnorm(x - mu, 0, sd) / dt(x - mu, df)
}

dratio_norm_norm <- function(x, mu = 0, sd = 1) {
  dnorm(x, mu, sd) / dnorm(x, mu, sd)
}

mu <- 1
sd <- sqrt(2)


ggplot() +
  geom_function(aes(col = "A"), 
                fun = dlaplace, args = list(mu = mu, sd = sd),
                size = 0.5) +
  geom_function(aes(col = "B"), 
                fun = dnorm, args = list(mean = mu, sd = sd),
                size = 0.5) +
  xlim(-3, 5) +
  ylim(0, 1.2) +
  scale_color_brewer(labels = c("A" = "Laplace", "B" = "Normal"),
                     palette = "Set2",
                     name = "") +
  ylab("Density") +
  theme_minimal() +
  theme(legend.position = "bottom",
        text = element_text(family = "LM Roman 10"),
        legend.background = element_blank()) +
  ggplot() +
  geom_function(aes(col = "A"), 
                fun = dlnorm, args = list(meanlog = log(mu^2 / sqrt(mu^2 + sd^2)), 
                                          sdlog = sqrt(log(1 + sd^2/mu^2))),
                size = 0.5) +
  geom_function(aes(col = "B"), 
                fun = dnorm, args = list(mean = mu, sd = sd),
                size = 0.5) +
  xlim(-3, 5) +
  ylim(0, 1.2) +
  scale_color_brewer(labels = c("A" = "Log-normal", "B" = "Normal"),
                     palette = "Set2",
                     name = "") +
  ylab("Density") +
  theme_minimal() +
  theme(legend.position = "bottom",
        text = element_text(family = "LM Roman 10"),
        legend.background = element_blank()) +
  ggplot() +
  geom_function(aes(col = "A"), 
                fun = ~dt(.x - mu, df = 2 * sd^2 / (sd^2 - 1)),
                size = 0.5) +
  geom_function(aes(col = "B"), 
                fun = dnorm, args = list(mean = mu, sd = sd),
                size = 0.5) +
  xlim(-3, 5) +
  ylim(0, 1.2) +
  scale_color_brewer(labels = c("A" = expression(italic(lst)), "B" = "Normal"),
                     palette = "Set2",
                     name = "") +
  ylab("Density") +
  theme_minimal() +
  theme(legend.position = "bottom",
        text = element_text(family = "LM Roman 10"),
        legend.background = element_blank()) +
  ggplot() +
  geom_function(aes(col = "A"), 
                fun = dnorm, args = list(mean = mu, sd = sd),
                size = 0.5) +
  geom_function(aes(col = "B"), 
                fun = dnorm, args = list(mean = mu, sd = sd),
                size = 0.5) +
  xlim(-3, 5) +
  ylim(0, 1.2) +
  scale_color_brewer(labels = c("A" = "Normal", "B" = "Normal"),
                     palette = "Set2",
                     name = "") +
  ylab("Density") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.background = element_blank(),
        text = element_text(family = "LM Roman 10"))

ggsave(filename = "files/densities.png", device = "png", dpi = 500,
       width = 7, height = 4)


nsim <- 200
n <- 250
mu <- 1
sd <- sqrt(2)

x_eval <- seq(-3, 5, length.out = 10000) |> matrix()

plot_dr <- function(x_eval, predict_list, dr_fun, dr_fun_args) {
  nsim <- length(predict_list)
  ggplot() +
    geom_line(mapping = aes(x = rep(x_eval, nsim),
                            y = unlist(predict_list),
                            group = rep(1:nrow(x_eval), each = nsim)),
              alpha = 0.1, col = "#014c61") +
    stat_function(fun = dr_fun, args = dr_fun_args) +
    theme_minimal() +
    ylab(expression(hat(italic(r)))) +
    xlab(NULL)
}

dat_lap_norm <- purrr::map(1:nsim, ~list(laplace = rlaplace(n, mu, sd),
                                         norm = rnorm(n, mu, sd)))

dat_lnorm_norm <- purrr::map(1:nsim, ~list(lnorm = rlnorm(n,
                                                          log(mu^2 / sqrt(mu^2 + sd^2)), 
                                                          sqrt(log(1 + sd^2/mu^2))),
                                           norm = rnorm(n, mu, sd)))

dat_t_norm <- purrr::map(1:nsim, ~list(t = rt(n, 2/(1 - 1/sd^2)) + mu,
                                       norm = rnorm(n, mu, sd)))

dat_norm_norm <- purrr::map(1:nsim, ~list(norm1 = rnorm(n, mu, sd),
                                          norm2 = rnorm(n, mu, sd)))

lap_norm_dr <- future_map(dat_lap_norm, ~ulsif(.x$norm, 
                                               .x$laplace, 
                                               ncenters = 100, 
                                               nsigma = 10,
                                               nlambda = 10,
                                               progressbar = FALSE),
                          .options = furrr_options(seed = TRUE))

lnorm_norm_dr <- future_map(dat_lnorm_norm, ~ulsif(.x$norm,
                                                   .x$lnorm, 
                                                   ncenters = 100, 
                                                   nsigma = 10,
                                                   nlambda = 10,
                                                   progressbar = FALSE),
                            .options = furrr_options(seed = TRUE))

t_norm_dr <- future_map(dat_t_norm, ~ulsif(.x$norm,
                                           .x$t, 
                                           ncenters = 100, 
                                           nsigma = 10,
                                           nlambda = 10,
                                           progressbar = FALSE),
                        .options = furrr_options(seed = TRUE))

norm_norm_dr <- future_map(dat_norm_norm, ~ulsif(.x$norm1,
                                                 .x$norm2, 
                                                 ncenters = 100, 
                                                 nsigma = 10,
                                                 nlambda = 10,
                                                 progressbar = FALSE),
                           .options = furrr_options(seed = TRUE))

lap_norm_predict   <- map(lap_norm_dr, ~predict(.x, newdata = x_eval))
lnorm_norm_predict <- map(lnorm_norm_dr, ~ predict(.x, newdata = x_eval))
t_norm_predict     <- map(t_norm_dr, ~predict(.x, newdata = x_eval))
norm_norm_predict  <- map(norm_norm_dr, ~ predict(.x, newdata = x_eval))

list(
  plot_dr(x_eval, lap_norm_predict, dratio_lap_norm, list(mu = mu, sd = sd)) +
    ggtitle(expression(italic(P[Normal])/italic(P[Laplace]))) +
    ylim(-1, 4) +
    theme(text = element_text(family = "LM Roman 10", size = 9)),
  plot_dr(x_eval, lnorm_norm_predict, dratio_lnorm_norm, list(mu = mu, sd = sd)) +
    ggtitle(expression(italic(P[Normal])/italic(P[`Log-normal`]))) +
    ylim(-1, 10) +
    theme(text = element_text(family = "LM Roman 10", size = 9)),
  plot_dr(x_eval, t_norm_predict, dratio_t_norm, list(mu = mu, sd = sd)) +
    ggtitle(expression(italic(P[Normal])/italic(P[lst]))) +
    ylim(-1, 4) +
    theme(text = element_text(family = "LM Roman 10", size = 9)),
  plot_dr(x_eval, norm_norm_predict, dratio_norm_norm, list(mu = mu, sd = sd)) +
    ggtitle(expression(italic(P[Normal])/italic(P[Normal]))) +
    ylim(-1, 4) +
    theme(text = element_text(family = "LM Roman 10", size = 9))
) |>
  patchwork::wrap_plots(ncol = 2, byrow = TRUE)


ggsave(filename = "files/density-ratios.png", device = "png", dpi = 500,
       width = 7, height = 4)


load("files/cps5000.RData")

df <- cps |>
  select(-csp) |>
  mutate(educ = as.numeric(as.character(educ)), 
         educ = case_when(educ < 39 ~ 1,
                          educ < 40 ~ 2,
                          educ < 44 ~ 3,
                          educ < 47 ~ 4) |>
           factor(labels = c("NoHS", "HS", "AoBD", "MoH")),
         race = factor(race, labels = c("White", "Non-white", "Non-white", "Non-white")),
         marital = factor(marital, labels = c("Married", "Married", "Separated", 
                                              "Separated", "Widowed", "Single", 
                                              "WidowedORDivorced")),
         sex = factor(sex, labels = c("Male", "Female")))

method.ini <- c("norm", "norm", "norm", "polyreg", "polyreg", "logreg", "logreg", "norm")
visit <- c(7, 1, 2, 3, 4, 5, 6, 8)
m <- 5
synlist <- list(
  unadj = syn(
    df, 
    m = m, 
    method = method.ini, 
    visit.sequence = visit, 
    seed = 1234, 
    print.flag = FALSE
  ), 
  trans = syn(
    df |> mutate(across(c(income, tax, age, ss), ~ .x^{1/3})),
    m = m, 
    method = method.ini,
    visit.sequence = visit, 
    seed = 1234,
    print.flag = FALSE
  ),
  semi = syn(
    df |> mutate(across(c(income, tax, age, ss), ~ .x^{1/3})),
    m = m,
    method = method.ini,
    visit.sequence = visit,
    semicont = list(tax = "0", ss = "0"),
    seed = 1234,
    print.flag = FALSE
  )
)

synlist$trans$syn <- synlist$trans$syn |>
  map(~ .x |> mutate(across(c(income, tax, age, ss), ~.x^3)))
synlist$semi$syn <- synlist$semi$syn |>
  map(~ .x |> mutate(across(c(income, tax, age, ss), ~.x^3)))


comb_df <- bind_rows(
  Real = df,
  Naive = bind_rows(synlist$unadj$syn), 
  Transformed = bind_rows(synlist$trans$syn),
  `Semi-continuous` = bind_rows(synlist$semi$syn),
  .id = "Data"
) |>
  select(Data,
         Age = age, 
         Income = income, 
         Tax = tax, 
         `Social security` = ss) |>
  mutate(Data = factor(Data, levels = c("Real", "Naive", "Transformed", "Semi-continuous")),
         RealSyn = ifelse(Data == "Real", 1, 2) |> factor(labels = c("Real", "Synthetic")),
         across(c(Age, Income, Tax, `Social security`), ~abs(.x)^{1/3} * sign(.x))) |>
  tidyr::pivot_longer(cols = c(Age, Income, Tax, `Social security`), names_to = "Variable")

purrr::map(c("Age", "Income", "Social security", "Tax"), ~
             ggplot(comb_df |> filter(Variable == .x), aes(x = value, fill = RealSyn, after_stat(density))) +
             geom_histogram(col = "black", bins = 20) +
             scale_fill_brewer(palette = "Set2") + 
             facet_wrap(~Data, ncol = 5) + 
             theme_minimal() +
             ylab(.x) +
             theme(legend.position = "none", 
                   axis.title.x = element_blank(), 
                   strip.text.x = element_text(size = 8),
                   text = element_text(family = "LM Roman 10"))) |>
  patchwork::wrap_plots(nrow = 4)

ggsave(filename = "files/syn-vars.png", device = "png", dpi = 500,
       width = 7, height = 4)

vars <- c("age", "income", "ss", "tax")

PE <- synlist |>
  map(function(x) {
    map_dbl(vars, function(var) {
      map_dbl(x$syn, ~ulsif(.x[[var]], df[[var]], nsigma = 10, nlambda = 10, 
                            ncenters = 100, progressbar = FALSE) |>
                summary(test = FALSE) |>
                (\(x) x$PE)()
      ) |> mean()
    })
  }, .progress = TRUE)

PE_allvars <- synlist |>
  map(function(x) {
    map_dbl(x$syn, ~ulsif(
      .x |> mutate(across(everything(), as.numeric)),
      df |> mutate(across(everything(), as.numeric)),
      nsigma = 10, nlambda = 10, ncenters = 100) |>
        summary(test = FALSE) |>
        (\(x) x$PE)()
    ) |> mean()
  })

bind_rows(PE, PE_allvars) |>
  mutate(variable = factor(1:5, labels = c("Age", "Income", "Social security", "Tax", "All"))) |>
  tidyr::pivot_longer(cols = c(unadj, trans, semi)) |>
  mutate(`Synthesis method` = factor(name, levels = c("unadj", "trans", "semi"),
                                     labels = c("Naïve", "Transformed", "Semi-continuous"))) |>
  ggplot(aes(x = variable, y = value, col = `Synthesis method`, shape = `Synthesis method`, 
             group = `Synthesis method`)) +
  geom_point(size = 3) +
  # geom_line(alpha = 0.1, linetype = "dashed") +
  scale_y_log10() +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  ylab("Pearson divergence") +
  xlab("Variable") +
  theme(legend.position = "bottom",
        text = element_text(family = "LM Roman 10", size = 9))

ggsave(filename = "files/syn-PEs.png", device = "png", dpi = 500,
       width = 7, height = 4, bg = "white")
