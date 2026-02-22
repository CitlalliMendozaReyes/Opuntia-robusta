############################################################
# 03_PO_vs_germination.R
# Linking P/O ratio with germination patterns
############################################################

message("\n[03] Loading packages...")

required_packages <- c("tidyverse", "lme4", "emmeans")

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

message("[03] Loading datasets...")

data_ovp <- file.path("data", "ovules_pollen.csv")
data_germ <- file.path("data", "germination.csv")

if (!file.exists(data_ovp)) {
  stop("File 'data/ovules_pollen.csv' not found. 
       Please download it from Zenodo and place it in the data/ folder.")
}

if (!file.exists(data_germ)) {
  stop("File 'data/germination.csv' not found. 
       Please download it from Zenodo and place it in the data/ folder.")
}

datos <- read_csv(data_ovp)
germ  <- read_csv(data_germ)


############################################################
# 2. Build germ_master (one row per fruit)
############################################################

message("[03] Building germination summary per fruit...")

germ_master <- germ %>%
  group_by(plant, fruit, population, treatment) %>%
  summarise(
    sampled_seeds    = n(),
    germinated_total = sum(germinated),
    .groups = "drop"
  ) %>%
  mutate(
    germination_prop = germinated_total / sampled_seeds,
    failures         = sampled_seeds - germinated_total
  )


############################################################
# 3. Compute P/O ratio per flower
############################################################

message("[03] Computing P/O ratio...")

datos <- datos %>%
  mutate(
    PO    = pollen_total / ovules,
    logPO = log(PO)
  )


############################################################
# 4. Mean P/O per population
############################################################

PO_pob <- datos %>%
  group_by(population) %>%
  summarise(
    PO_mean    = mean(PO, na.rm = TRUE),
    logPO_mean = mean(logPO, na.rm = TRUE),
    n_PO       = n()
  )


############################################################
# 5. Mean germination per population
############################################################

germ_pob <- germ_master %>%
  group_by(population) %>%
  summarise(
    germ_mean = mean(germination_prop, na.rm = TRUE),
    n_germ    = n()
  )


############################################################
# 6. Merge both summaries
############################################################

message("[03] Merging P/O and germination summaries...")

pob_merged <- left_join(PO_pob, germ_pob, by = "population")


############################################################
# 7. Statistical comparisons
############################################################

message("[03] Running statistical models...")

# 7A. P/O between populations
modelo_PO <- lm(PO ~ population, data = datos)

# 7B. Germination between populations (GLMM)
modelo_germ <- glmer(
  cbind(germinated_total, failures) ~ population + (1 | plant),
  data = germ_master,
  family = binomial
)

emm_germ <- emmeans(modelo_germ, ~ population, type = "response")


############################################################
# 8. Save results
############################################################

message("[03] Saving results...")

dir.create("results", showWarnings = FALSE)

save(
  PO_pob, germ_pob, pob_merged,
  modelo_PO, modelo_germ, emm_germ,
  file = "results/03_PO_vs_germination.RData"
)

write_csv(pob_merged, "results/PO_vs_germination_summary.csv")

message("[03] P/O vs germination analysis completed.\n")
