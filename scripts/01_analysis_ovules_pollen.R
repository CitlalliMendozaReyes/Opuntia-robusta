############################################################
# 01_analysis_ovules_pollen.R
# Analysis of ovules, pollen, and P/O ratio
############################################################

message("\n[01] Loading packages...")

required_packages <- c("tidyverse", "emmeans", "broom")

missing <- required_packages[!required_packages %in% installed.packages()[, "Package"]]

if (length(missing) > 0) {
  stop(paste(
    "Missing packages:", paste(missing, collapse = ", "),
    "\nPlease install them before running the pipeline."
  ))
}

lapply(required_packages, library, character.only = TRUE)


############################################################
# 1. Load data
############################################################

message("[01] Loading ovules_pollen.csv...")

data_path <- file.path("data", "ovules_pollen.csv")

if (!file.exists(data_path)) {
  stop("File 'data/ovules_pollen.csv' not found. 
       Please download it from Zenodo and place it in the data/ folder.")
}

datos <- read_csv(data_path) %>%
  rename(
    Population = population,
    Plant      = plant,
    Fruit      = fruit,
    Ovules     = ovules,
    Pollen     = pollen_total
  )


############################################################
# 2. ANALYSIS A — Ovules per flower
############################################################

message("[01] Fitting ovule model...")

mod_ov <- lm(log(Ovules) ~ Population, data = datos)
emm_ov <- emmeans(mod_ov, ~ Population)
ovules_back <- summary(emm_ov, type = "response")

tabla_ovules <- tidy(mod_ov) %>%
  mutate(
    estimate_exp = exp(estimate),
    CI_low       = exp(estimate - 1.96 * std.error),
    CI_high      = exp(estimate + 1.96 * std.error)
  )

fig_ovules <- ggplot(datos, aes(x = Population, y = Ovules, color = Population)) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Number of ovules per flower",
    y = "Ovules per flower"
  )


############################################################
# 3. ANALYSIS B — Pollen per flower
############################################################

message("[01] Fitting pollen model...")

mod_po <- lm(log(Pollen) ~ Population, data = datos)
emm_po <- emmeans(mod_po, ~ Population)
pollen_back <- summary(emm_po, type = "response")

tabla_pollen <- tidy(mod_po) %>%
  mutate(
    estimate_exp = exp(estimate),
    CI_low       = exp(estimate - 1.96 * std.error),
    CI_high      = exp(estimate + 1.96 * std.error)
  )

fig_pollen <- ggplot(datos, aes(x = Population, y = Pollen, color = Population)) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Number of pollen grains per flower",
    y = "Pollen grains per flower"
  )


############################################################
# 4. ANALYSIS C — Cruden P/O ratio
############################################################

message("[01] Computing Cruden P/O ratio...")

datos <- datos %>%
  mutate(
    PO    = Pollen / Ovules,
    logPO = log(PO)
  )

mod_cruden <- lm(logPO ~ Population, data = datos)
emm_cruden <- emmeans(mod_cruden, ~ Population)
cruden_back <- summary(emm_cruden, type = "response")

tabla_cruden <- tidy(mod_cruden) %>%
  mutate(
    estimate_exp = exp(estimate),
    CI_low       = exp(estimate - 1.96 * std.error),
    CI_high      = exp(estimate + 1.96 * std.error)
  )

fig_logPO <- ggplot(datos, aes(x = Population, y = logPO, color = Population)) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "P/O ratio (Cruden model)",
    y = "log(P/O)"
  )

fig_PO <- ggplot(datos, aes(x = Population, y = PO, color = Population)) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "P/O ratio (original scale)",
    y = "P/O"
  )


############################################################
# 5. Save results
############################################################

message("[01] Saving results...")

dir.create("results", showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)

save(
  tabla_ovules, tabla_pollen, tabla_cruden,
  ovules_back, pollen_back, cruden_back,
  file = "results/01_analysis_ovules_pollen.RData"
)

ggsave("figs/ovules_per_flower.png", fig_ovules, width = 6, height = 5)
ggsave("figs/pollen_per_flower.png", fig_pollen, width = 6, height = 5)
ggsave("figs/logPO_ratio.png", fig_logPO, width = 6, height = 5)
ggsave("figs/PO_ratio.png", fig_PO, width = 6, height = 5)

message("[01] Ovule–pollen analysis completed.\n")
