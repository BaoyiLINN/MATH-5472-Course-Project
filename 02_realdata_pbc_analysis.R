## 02_realdata_pbc_analysis.R
## Real-data experiment using the PBC dataset
## - Cox PH
## - PH diagnostics (Schoenfeld residuals)
## - KM curves (bili high vs low)
## - AFT model
## - Time-varying Cox model
## - PEC prediction error (IBS) for Cox / AFT / TVC

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(splines)
library(pec)

## 0. Directory setup ---------------------------------------------------------

dir.create("figures/pbc", recursive = TRUE, showWarnings = FALSE)
dir.create("results/pbc", recursive = TRUE, showWarnings = FALSE)

## 1. Data loading and preprocessing ------------------------------------------

data(pbc, package = "survival")

df <- pbc %>%
  filter(!is.na(bili), !is.na(albumin), !is.na(protime)) %>%
  transmute(
    time    = time,
    status  = as.numeric(status == 2),  # 2 = death
    bili    = scale(bili),
    albumin = scale(albumin),
    protime = scale(protime)
  )

## 2. Cox PH model -------------------------------------------------------------

fit.cox <- coxph(Surv(time, status) ~ bili + albumin + protime, data = df)
print(summary(fit.cox))

## 3. PH test + Schoenfeld residuals ------------------------------------------

zph <- cox.zph(fit.cox)
print(zph)

# survminer::ggcoxzph returns a list on your machine; elements 1/2/3 correspond to each covariate
p_zph_list <- ggcoxzph(zph)

# Extract variable names from zph$table (exclude GLOBAL row)
vars <- rownames(zph$table)
vars <- vars[vars != "GLOBAL"]  # typically c("bili", "albumin", "protime")

# Safety check: ensure the list length matches
n_plot <- min(length(p_zph_list), length(vars))

for (i in seq_len(n_plot)) {
  p_i <- p_zph_list[[i]]
  if (inherits(p_i, "gg")) {
    fname_png <- paste0("figures/pbc/schoenfeld_", vars[i], ".png")
    fname_pdf <- paste0("figures/pbc/schoenfeld_", vars[i], ".pdf")
    ggsave(fname_png, p_i, width = 6, height = 4)
    ggsave(fname_pdf, p_i, width = 6, height = 4)
  }
}

## 4. KM curves (high vs low bilirubin) ---------------------------------------

df <- df %>%
  mutate(
    bili_group = ifelse(bili > median(bili), "High bilirubin", "Low bilirubin")
  )

fit.km <- survfit(Surv(time, status) ~ bili_group, data = df)

p_km <- ggsurvplot(
  fit.km, data = df,
  risk.table = FALSE,
  pval = TRUE,
  conf.int = FALSE,
  ggtheme = theme_minimal(base_size = 14),
  palette = c("#d95f02", "#1b9e77"),
  xlab = "Time (days)",
  ylab = "Overall survival probability",
  legend.title = "",
  legend.labs = c("High bilirubin", "Low bilirubin")
)

ggsave("figures/pbc/km_bilirubin_pbc.png", p_km$plot, width = 7, height = 5)
ggsave("figures/pbc/km_bilirubin_pbc.pdf", p_km$plot, width = 7, height = 5)

## 5. AFT model ---------------------------------------------------------------

fit.aft <- survreg(
  Surv(time, status) ~ bili + albumin + protime,
  data = df,
  dist = "lognormal"
)
print(summary(fit.aft))

## 6. Time-varying coefficient Cox model --------------------------------------

fit.tv <- coxph(
  Surv(time, status) ~ ns(time, df = 3) * bili + albumin + protime,
  data = df
)
print(summary(fit.tv))

## 7. Predictive performance evaluation: Cox / AFT / TVC -----------------------

set.seed(123)
idx   <- sample(seq_len(nrow(df)), size = floor(0.7 * nrow(df)))
train <- df[idx, ]
test  <- df[-idx, ]

fit.cox.t <- coxph(Surv(time, status) ~ bili + albumin + protime,
                   data = train, x = TRUE, y = TRUE)

fit.aft.t <- survreg(Surv(time, status) ~ bili + albumin + protime,
                     data = train, dist = "lognormal")

fit.tv.t <- coxph(
  Surv(time, status) ~ ns(time, df = 3) * bili + albumin + protime,
  data = train, x = TRUE, y = TRUE
)

## 7. PEC prediction error for Cox & TVC --------------------------------------

models <- list(
  Cox = fit.cox.t,
  TVC = fit.tv.t
)

times_grid <- seq(100, quantile(test$time, 0.9), length.out = 20)

pec_res <- pec(
  object     = models,
  formula    = Surv(time, status) ~ 1,
  data       = test,
  cens.model = "cox",
  times      = times_grid
)

saveRDS(pec_res, "results/pbc/pec_results.rds")

## --- Manual IBS computation (always works) ---
## ---------------------------------------------------
## Robust IBS calculation compatible with all pec outputs
## ---------------------------------------------------

compute_IBS <- function(pec_obj, model_name) {
  pe <- pec_obj$AppErr[[model_name]]
  times <- pec_obj$time
  
  # Case 1: pe is a matrix (n_test × n_times)
  if (is.matrix(pe)) {
    bs_mean <- colMeans(pe, na.rm = TRUE)
    
    # Case 2: pe is a numeric vector
  } else if (is.numeric(pe)) {
    bs_mean <- pe
    
    # Case 3: pe is a list of numeric vectors (older pec structure)
  } else if (is.list(pe)) {
    bs_mat <- do.call(rbind, pe)
    bs_mean <- colMeans(bs_mat, na.rm = TRUE)
    
  } else {
    stop("Unsupported AppErr format for model: ", model_name)
  }
  
  # Trapezoidal integration
  ibs <- sum(diff(times) * (head(bs_mean, -1) + tail(bs_mean, -1)) / 2)
  return(ibs)
}

## Compute IBS for all models
ibs_df <- data.frame(
  model = names(models),
  ibs   = sapply(names(models), function(m) compute_IBS(pec_res, m))
)

print(ibs_df)
write.csv(ibs_df, "results/pbc/ibs_values.csv", row.names = FALSE)

## --- PEC plot (standard visual output) ---
png("figures/pbc/pec_error_pbc.png", width = 7, height = 5, units = "in", res = 300)
plot(pec_res, xlab = "Time (days)", legend.x = "topright")
dev.off()

pdf("figures/pbc/pec_error_pbc.pdf", width = 7, height = 5)
plot(pec_res, xlab = "Time (days)", legend.x = "topright")
dev.off()
