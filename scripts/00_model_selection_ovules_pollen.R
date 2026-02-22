############################################################
# 00_model_selection_ovules_pollen.R
# Model selection for ovules and pollen counts
############################################################

message("\n[00] Loading packages...")

required_packages <- c("tidyverse", "MASS", "lme4", "MuMIn", "broom")

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

message("[00] Loading ovules_pollen.csv...")

data_path <- file.path("data", "ovules_pollen.csv")

if (!file.exists(data_path)) {
  stop("File 'data/ovules_pollen.csv' not found. 
       Please download it from Zenodo and place it in the data/ folder.")
}

datos <- read_csv(data_path) %>%
  rename(
    Poblacion = population,
    Planta    = plant,
    Fruto     = fruit,
    Ovulos    = ovules,
    Polen     = pollen_total
  )


############################################################
# 2. Overdispersion function
############################################################

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp  <- residuals(model, type = "pearson")
  sum(rp^2) / rdf
}


############################################################
# 3. Models for Ovules
############################################################

message("[00] Fitting ovule models...")

m1_ov <- glm(Ovulos ~ Poblacion, data = datos, family = poisson)
m2_ov <- glm(Ovulos ~ Poblacion, data = datos, family = quasipoisson)
m3_ov <- glm.nb(Ovulos ~ Poblacion, data = datos)
m4_ov <- glm(Ovulos ~ Poblacion, data = datos, family = Gamma(link = "log"))
m5_ov <- lm(log(Ovulos) ~ Poblacion, data = datos)
m6_ov <- glmer(Ovulos ~ Poblacion + (1 | Planta), data = datos, family = poisson)


############################################################
# 4. Models for Pollen
############################################################

message("[00] Fitting pollen models...")

m1_po <- glm(Polen ~ Poblacion, data = datos, family = poisson)
m2_po <- glm(Polen ~ Poblacion, data = datos, family = quasipoisson)
m3_po <- glm.nb(Polen ~ Poblacion, data = datos)
m4_po <- glm(Polen ~ Poblacion, data = datos, family = Gamma(link = "log"))
m5_po <- lm(log(Polen) ~ Poblacion, data = datos)
m6_po <- glmer(Polen ~ Poblacion + (1 | Planta), data = datos, family = poisson)


############################################################
# 5. AICc comparison
############################################################

message("[00] Computing AICc tables...")

AICc_ovulos <- AICc(m1_ov, m3_ov, m4_ov, m5_ov)
AICc_polen  <- AICc(m1_po, m3_po, m4_po, m5_po)


############################################################
# 6. Overdispersion
############################################################

message("[00] Checking overdispersion...")

sobredisp_ovulos <- c(
  Poisson = overdisp_fun(m1_ov),
  NegBin  = overdisp_fun(m3_ov)
)

sobredisp_polen <- c(
  Poisson = overdisp_fun(m1_po),
  NegBin  = overdisp_fun(m3_po)
)


############################################################
# 7. Clean tables for manuscript
############################################################

message("[00] Creating tidy tables...")

tabla_ovulos <- tidy(m3_ov) %>%
  mutate(
    OR     = exp(estimate),
    CI_low = exp(estimate - 1.96 * std.error),
    CI_high= exp(estimate + 1.96 * std.error)
  )

tabla_polen <- tidy(m4_po) %>%
  mutate(
    OR     = exp(estimate),
    CI_low = exp(estimate - 1.96 * std.error),
    CI_high= exp(estimate + 1.96 * std.error)
  )


############################################################
# 8. Save results
############################################################

dir.create("results", showWarnings = FALSE)

save(
  AICc_ovulos, AICc_polen,
  sobredisp_ovulos, sobredisp_polen,
  tabla_ovulos, tabla_polen,
  file = "results/00_model_selection_results.RData"
)

message("[00] Model selection completed.\n")
