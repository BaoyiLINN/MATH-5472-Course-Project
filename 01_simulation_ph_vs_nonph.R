## 01_simulation_ph_vs_nonph.R
## Simulation study: PH vs Non-PH scenarios
## Metrics: beta_hat distribution + time-averaged Brier + RMST bias

rm(list = ls())
set.seed(2025)

library(survival)
library(ggplot2)
library(dplyr)
library(tidyr)

############################################
## 1. Data-generating functions: PH and Non-PH scenarios
############################################

simulate_PH <- function(n, lambda0 = 0.1, beta = log(2), cens_rate = 0.03) {
  Z <- rnorm(n, 0, 1)
  rate_event <- lambda0 * exp(beta * Z)
  T_event <- rexp(n, rate = rate_event)
  C <- rexp(n, rate = cens_rate)
  time <- pmin(T_event, C)
  status <- as.integer(T_event <= C)
  data.frame(time = time, status = status, Z = Z, scenario = "PH")
}

simulate_nonPH <- function(n, lambda0 = 0.1, 
                           beta1 = log(2), beta2 = -log(2),
                           t_cut = 10, cens_rate = 0.03) {
  Z <- rnorm(n, 0, 1)
  rate1 <- lambda0 * exp(beta1 * Z)
  rate2 <- lambda0 * exp(beta2 * Z)
  T1 <- rexp(n, rate = rate1)
  T_event <- ifelse(T1 <= t_cut,
                    T1,
                    t_cut + rexp(n, rate = rate2))
  C <- rexp(n, rate = cens_rate)
  time <- pmin(T_event, C)
  status <- as.integer(T_event <= C)
  data.frame(time = time, status = status, Z = Z, scenario = "NonPH")
}

############################################
## 2. Survival prediction + Brier + RMST utilities
############################################

## Cox PH survival prediction: using basehaz + linear predictor
predict_surv_cox <- function(fit, data, times) {
  base <- basehaz(fit, centered = FALSE)
  H0_t <- sapply(times, function(t0) {
    idx <- which(base$time <= t0)
    if (length(idx) == 0) return(0)
    base$hazard[max(idx)]
  })
  lp <- predict(fit, newdata = data, type = "lp")
  sapply(H0_t, function(Ht) exp(-Ht * exp(lp)))
}

## Log-normal AFT survival prediction
predict_surv_aft_lognormal <- function(fit, data, times) {
  lp <- predict(fit, newdata = data, type = "lp")  # μ(x)
  sigma <- fit$scale
  outer(lp, times, FUN = function(mu, t0) {
    1 - pnorm((log(t0) - mu) / sigma)
  })
}

## Brier(t): without censoring correction (censoring rate is low in simulation)
brier_curve <- function(time, status, surv_mat, times) {
  n <- length(time)
  sapply(seq_along(times), function(j) {
    t0 <- times[j]
    y <- as.integer(time <= t0 & status == 1)
    p <- 1 - surv_mat[, j]  # event probability
    mean((y - p)^2)
  })
}

## Mean Brier score over [0, tau]
compute_brier_mean <- function(time, status, surv_mat, times, tau = max(times)) {
  idx <- which(times <= tau)
  if (length(idx) < 1) return(NA_real_)
  tt <- times[idx]
  bs <- brier_curve(time, status, surv_mat, tt)
  mean(bs, na.rm = TRUE)
}

## KM-based RMST: ∫_0^τ S_KM(t) dt
compute_rmst_km <- function(time, status, tau) {
  fit <- survfit(Surv(time, status) ~ 1)
  t <- fit$time
  s <- fit$surv
  
  ## keep only points <= tau
  idx <- which(t <= tau)
  t_trim <- t[idx]
  s_trim <- s[idx]
  
  ## If max event time < tau, append a point at tau
  if (length(t_trim) == 0 || max(t_trim) < tau) {
    t_all <- c(0, t_trim, tau)
    s_all <- c(1, s_trim, ifelse(length(s_trim) == 0, 1, tail(s_trim, 1)))
  } else {
    t_all <- c(0, t_trim)
    s_all <- c(1, s_trim)
  }
  
  ## Trapezoidal integration for RMST
  sum(s_all[-length(s_all)] * diff(t_all))
}

## Model-predicted RMST: integrate S_hat(t|x) over time, then average across individuals
compute_rmst_pred <- function(surv_mat, times, tau) {
  if (is.null(dim(surv_mat))) return(NA_real_)
  if (nrow(surv_mat) == 0) return(NA_real_)
  
  idx <- which(times <= tau)
  if (length(idx) < 1) return(NA_real_)
  tt <- times[idx]
  
  ## Integration grid: 0, t1, ..., tk, tau
  if (max(tt) < tau) {
    tt_all <- c(0, tt, tau)
    S_tau <- surv_mat[, length(tt), drop = FALSE]
    S_all <- cbind(1, surv_mat[, idx, drop = FALSE], S_tau)
  } else {
    tt_all <- c(0, tt)
    S_all <- cbind(1, surv_mat[, idx, drop = FALSE])
  }
  
  ## RMST_i = ∑ S_i(t_j) Δt_j
  dt <- diff(tt_all)
  rmst_i <- S_all[, -ncol(S_all), drop = FALSE] %*% dt
  mean(rmst_i)
}

############################################
## 3. Single simulation: fit three models and compute metrics
############################################

fit_and_evaluate <- function(dat, 
                             beta_true_PH = log(2),
                             t_cut = 10,
                             eval_times = seq(2, 30, by = 1),
                             tau_rmst = max(eval_times)) {
  
  time   <- dat$time
  status <- dat$status
  Z      <- dat$Z
  n      <- nrow(dat)
  
  out_list <- list()
  
  ## 1) Cox PH model
  fit_cox <- try(coxph(Surv(time, status) ~ Z, data = dat), silent = TRUE)
  if (!inherits(fit_cox, "try-error")) {
    beta_hat <- coef(fit_cox)["Z"]
    se_hat   <- sqrt(vcov(fit_cox)["Z", "Z"])
    ci_low   <- beta_hat - 1.96 * se_hat
    ci_up    <- beta_hat + 1.96 * se_hat
    cover    <- as.integer(beta_true_PH >= ci_low & beta_true_PH <= ci_up)
    
    ## Survival prediction
    surv_mat <- try(predict_surv_cox(fit_cox, dat, eval_times), silent = TRUE)
    if (inherits(surv_mat, "try-error")) {
      surv_mat <- matrix(NA_real_, nrow = n, ncol = length(eval_times))
    }
    
    ## Brier & RMST
    brier_mean <- compute_brier_mean(time, status, surv_mat, eval_times, tau = tau_rmst)
    rmst_km    <- compute_rmst_km(time, status, tau_rmst)
    rmst_pred  <- compute_rmst_pred(surv_mat, eval_times, tau_rmst)
    rmst_bias  <- rmst_pred - rmst_km
    
    out_list[["CoxPH"]] <- data.frame(
      model      = "CoxPH",
      beta_hat   = beta_hat,
      se         = se_hat,
      cover      = cover,
      brier_mean = brier_mean,
      rmst_km    = rmst_km,
      rmst_pred  = rmst_pred,
      rmst_bias  = rmst_bias
    )
  }
  
  ## 2) Log-normal AFT model
  fit_aft <- try(survreg(Surv(time, status) ~ Z,
                         data = dat, dist = "lognormal"),
                 silent = TRUE)
  if (!inherits(fit_aft, "try-error")) {
    beta_hat <- coef(fit_aft)["Z"]
    se_hat   <- sqrt(vcov(fit_aft)["Z", "Z"])
    
    surv_mat <- try(predict_surv_aft_lognormal(fit_aft, dat, eval_times),
                    silent = TRUE)
    if (inherits(surv_mat, "try-error")) {
      surv_mat <- matrix(NA_real_, nrow = n, ncol = length(eval_times))
    }
    
    brier_mean <- compute_brier_mean(time, status, surv_mat, eval_times, tau = tau_rmst)
    rmst_km    <- compute_rmst_km(time, status, tau_rmst)
    rmst_pred  <- compute_rmst_pred(surv_mat, eval_times, tau_rmst)
    rmst_bias  <- rmst_pred - rmst_km
    
    out_list[["AFT"]] <- data.frame(
      model      = "AFT_lognormal",
      beta_hat   = beta_hat,
      se         = se_hat,
      cover      = NA_integer_,
      brier_mean = brier_mean,
      rmst_km    = rmst_km,
      rmst_pred  = rmst_pred,
      rmst_bias  = rmst_bias
    )
  }
  
  ## 3) TVC Cox: only for beta_hat comparison, skip prediction metrics
  fit_tvc <- try(
    coxph(
      Surv(time, status) ~ Z + tt(Z),
      data = dat,
      tt = function(z, t, ...) ifelse(t > t_cut, z, 0)
    ),
    silent = TRUE
  )
  
  if (!inherits(fit_tvc, "try-error")) {
    coefs <- coef(fit_tvc)
    beta1_hat <- coefs["Z"]  ## main effect for t <= t_cut
    out_list[["TVC"]] <- data.frame(
      model      = "TVC_step",
      beta_hat   = beta1_hat,
      se         = NA_real_,
      cover      = NA_integer_,
      brier_mean = NA_real_,
      rmst_km    = NA_real_,
      rmst_pred  = NA_real_,
      rmst_bias  = NA_real_
    )
  }
  
  bind_rows(out_list)
}

############################################
## 4. Main simulation loop
############################################

n_sim     <- 500
n_sample  <- 500
beta_true <- log(2)
t_cut     <- 10
eval_times <- seq(2, 30, by = 1)
tau_rmst  <- max(eval_times)

results_list <- list()

for (s in c("PH", "NonPH")) {
  cat("Running scenario:", s, "\n")
  for (b in 1:n_sim) {
    if (s == "PH") {
      dat <- simulate_PH(n_sample, lambda0 = 0.1,
                         beta = beta_true, cens_rate = 0.03)
    } else {
      dat <- simulate_nonPH(n_sample, lambda0 = 0.1,
                            beta1 = beta_true, beta2 = -beta_true,
                            t_cut = t_cut, cens_rate = 0.03)
    }
    tmp <- fit_and_evaluate(
      dat,
      beta_true_PH = beta_true,
      t_cut        = t_cut,
      eval_times   = eval_times,
      tau_rmst     = tau_rmst
    )
    tmp$scenario <- s
    tmp$sim      <- b
    results_list[[length(results_list) + 1]] <- tmp
  }
}

res <- bind_rows(results_list)

############################################
## 5. Metric summary
############################################

summary_effect <- res %>%
  group_by(scenario, model) %>%
  summarise(
    mean_beta   = mean(beta_hat, na.rm = TRUE),
    bias        = mean(beta_hat - beta_true, na.rm = TRUE),
    rmse        = sqrt(mean((beta_hat - beta_true)^2, na.rm = TRUE)),
    cover       = mean(cover, na.rm = TRUE),
    brier_mean  = mean(brier_mean, na.rm = TRUE),
    rmst_bias   = mean(rmst_bias, na.rm = TRUE),
    .groups     = "drop"
  )

print(summary_effect)

############################################
## 6. Figure 1: beta_hat distribution (PH vs NonPH)
############################################

dir.create("figures/simulation", recursive = TRUE, showWarnings = FALSE)

res_beta_plot <- res %>%
  filter(model %in% c("CoxPH", "AFT_lognormal", "TVC_step"))

gg_fig1 <- ggplot(res_beta_plot,
                  aes(x = model, y = beta_hat, fill = model)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.12, size = 0.5) +
  facet_wrap(~ scenario, nrow = 1) +
  geom_hline(yintercept = beta_true, linetype = "dashed", color = "red") +
  labs(
    title = expression(hat(beta)~"distributions under PH vs. Non-PH scenarios"),
    x = NULL,
    y = expression(hat(beta))
  ) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 11))

print(gg_fig1)
ggsave("figures/simulation/figure1_beta_distribution.pdf", gg_fig1, width = 8, height = 4)
ggsave("figures/simulation/figure1_beta_distribution.png", gg_fig1, width = 8, height = 4)

############################################
## 7. Figure 2: RMST bias (PH vs NonPH)
############################################

res_perf_plot <- summary_effect %>%
  filter(model %in% c("CoxPH", "AFT_lognormal")) %>%
  select(scenario, model, rmst_bias) %>%
  pivot_longer(
    cols      = c(rmst_bias),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = dplyr::recode(metric,
                                rmst_bias  = "RMST bias (model - KM)"))

## Compute p-values for RMST bias difference between CoxPH and AFT
p_values_rmst <- res %>%
  filter(model %in% c("CoxPH", "AFT_lognormal")) %>%
  select(scenario, model, rmst_bias) %>%
  drop_na() %>%
  group_by(scenario) %>%
  summarise(
    p_value = t.test(
      rmst_bias[model == "CoxPH"],
      rmst_bias[model == "AFT_lognormal"]
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_label = paste0("p = ", format(p_value, digits = 3, scientific = TRUE))
  )

## Place p-values above bars
p_values_rmst$y_pos <- res_perf_plot %>%
  group_by(scenario) %>%
  summarise(y = max(value, na.rm = TRUE)) %>%
  pull(y) * 1.05

gg_fig2 <- ggplot(res_perf_plot,
                  aes(x = model, y = value, fill = model)) +
  geom_col(position = position_dodge()) +
  facet_grid(metric ~ scenario, scales = "free_y") +
  geom_text(
    data = p_values_rmst,
    aes(x = 1.5, y = y_pos, label = p_label),
    inherit.aes = FALSE,
    size = 4
  ) +
  labs(
    title = "RMST bias under PH vs Non-PH scenarios",
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 11))


print(gg_fig2)
ggsave("figures/simulation/figure2_rmst.png", gg_fig2, width = 8, height = 5)
ggsave("figures/simulation/figure2_rmst.pdf", gg_fig2, width = 8, height = 5)

############################################
##  End of simulation
############################################

library(patchwork)

#-------------------------
# Extract PH panels
#-------------------------
fig1_PH  <- gg_fig1 + facet_wrap(~scenario, nrow = 1, scales = "free") + 
  scale_x_discrete() + 
  coord_cartesian() +
  theme(strip.text = element_blank())

fig1_PH  <- fig1_PH %+% subset(res_beta_plot, scenario == "PH")
ggsave("figures/simulation/figure_beta_PH.png", fig1_PH, width = 8, height = 7)

fig2_PH  <- gg_fig2 + facet_grid(metric ~ scenario, scales = "free_y") +
  theme(strip.text = element_blank())

fig2_PH  <- fig2_PH %+% subset(res_perf_plot, scenario == "PH")

fig2_PH <- res_perf_plot %>%
  filter(scenario == "PH") %>%
  ggplot(aes(x = model, y = value, fill = model)) +
  geom_col(position = position_dodge()) +
  facet_grid(metric ~ ., scales = "free_y") +
  geom_text(
    data = subset(p_values_rmst, scenario == "PH"),
    aes(x = 1.5, y = y_pos, label = p_label),
    inherit.aes = FALSE, size = 4
  ) +
  labs(
    x = NULL, y = NULL,
    title = "RMST bias under PH scenario"
  ) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 11))

print(fig2_PH)

ggsave("figures/simulation/figure_rmstbias_PH.png", 
       fig2_PH, width = 8, height = 7)

#-------------------------
# Extract NonPH panels
#-------------------------
fig1_nonPH <- gg_fig1 %+% subset(res_beta_plot, scenario == "NonPH") +
  facet_wrap(~scenario, nrow = 1, scales = "free") +
  theme(strip.text = element_blank())

fig2_nonPH <- gg_fig2 %+% subset(res_perf_plot, scenario == "NonPH") +
  facet_grid(metric ~ scenario, scales = "free_y") +
  theme(strip.text = element_blank())

fig2_nonPH <- res_perf_plot %>%
  filter(scenario == "NonPH") %>%
  ggplot(aes(x = model, y = value, fill = model)) +
  geom_col(position = position_dodge()) +
  facet_grid(metric ~ ., scales = "free_y") +
  geom_text(
    data = subset(p_values_rmst, scenario == "NonPH"),
    aes(x = 1.5, y = y_pos, label = p_label),
    inherit.aes = FALSE, size = 4
  ) +
  labs(
    x = NULL, y = NULL,
    title = "NonPH scenario"
  ) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 11))


combined_PH <- fig1_PH / fig2_PH + plot_layout(heights = c(1, 1.2))

combined_nonPH <- fig1_nonPH / fig2_nonPH + plot_layout(heights = c(1, 1.2))

ggsave("figures/simulation/figure_NonPH_combined.png", 
       combined_nonPH, width = 8, height = 7)
