############################################################
# 02_analysis_germination.R
# Germination analysis (quasi-binomial + GLMM)
############################################################

message("\n[02] Loading packages...")

required_packages <- c("tidyverse", "lme4", "emmeans", "broom.mixed")

missing <- required_packages[!required_packages %in% installed.packages()[, "Package"]]

if (length(missing) > 0) {
  stop(paste(
    "Missing packages:", paste(missing, collapse = ", "),
    "\nPlease install them before running the pipeline."
  ))
}

lapply(required_packages, library, character.only = TRUE)


############################################################
# 1. Load germination CSV
############################################################

message("[02] Loading germination.csv...")

data_path <- file.path("data", "germination.csv")

if (!file.exists(data_path)) {
  stop("File 'data/germination.csv' not found. 
       Please download it from Zenodo and place it in the data/ folder.")
}

germ <- read_csv(data_path)


############################################################
# 2. Build germ_master (one row per fruit)
############################################################

message("[02] Building germination summary per fruit...")

germ_master <- germ %>%
  group_by(plant, fruit, population, treatment) %>%
  summarise(
    sampled_seeds   = n(),
    germinated_total = sum(germinated),
    .groups = "drop"
  ) %>%
  mutate(
    germination_prop = germinated_total / sampled_seeds,
    failures = sampled_seeds - germinated_total
  )


############################################################
# 3. Quasi-binomial model
############################################################

message("[02] Fitting quasi-binomial model...")

modelo_quasi <- glm(
  cbind(germinated_total, failures) ~ population + treatment,
  data = germ_master,
  family = quasibinomial
)


############################################################
# 4. GLMM (random effect: plant)
############################################################

message("[02] Fitting GLMM (binomial, random plant)...")

modelo_glmm <- glmer(
  cbind(germinated_total, failures) ~ population + treatment + (1 | plant),
  data = germ_master,
  family = binomial
)


############################################################
# 5. AJB-style table (odds ratios + CI)
############################################################

message("[02] Creating tidy GLMM table...")

tabla_glmm <- summary(modelo_glmm)$coefficients %>%
  as.data.frame() %>%
  mutate(
    term   = rownames(.),
    OR     = exp(Estimate),
    CI_low = exp(Estimate - 1.96 * `Std. Error`),
    CI_high= exp(Estimate + 1.96 * `Std. Error`)
  ) %>%
  select(
    term, Estimate, SE = `Std. Error`, z = `z value`,
    p = `Pr(>|z|)`, OR, CI_low, CI_high
  )


############################################################
# 6. Save results
############################################################

message("[02] Saving results...")

dir.create("results", showWarnings = FALSE)

save(
  germ_master,
  modelo_quasi,
  modelo_glmm,
  tabla_glmm,
  file = "results/02_analysis_germination.RData"
)

write_csv(tabla_glmm, "results/germination_glmm_table.csv")

message("[02] Germination analysis completed.\n")
