############################################################
# 04_figures.R
# Final figure generation for the Opuntia robusta project
############################################################

message("\n[04] Loading packages...")

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
# 1. Load datasets
############################################################

message("[04] Loading datasets...")

data_ovp  <- file.path("data", "ovules_pollen.csv")
data_germ <- file.path("data", "germination.csv")

if (!file.exists(data_ovp)) {
  stop("File 'data/ovules_pollen.csv' not found.")
}

if (!file.exists(data_germ)) {
  stop("File 'data/germination.csv' not found.")
}

datos <- read_csv(data_ovp)
germ  <- read_csv(data_germ)


############################################################
# 2. Build derived datasets
############################################################

message("[04] Building germination master dataset...")

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

message("[04] Computing P/O ratios...")

datos <- datos %>%
  mutate(
    PO    = pollen_total / ovules,
    logPO = log(PO)
  )


############################################################
# 2.5 Recoding treatments and populations
############################################################

message("[04] Recoding treatments and populations...")

germ_master <- germ_master %>%
  mutate(
    treatment = str_trim(treatment),
    treatment = str_to_lower(treatment),
    treatment = recode(
      treatment,
      "man_self"  = "MS",
      "auto_self" = "AS",
      "ctrl"      = "CT",
      "cross"     = "CR",
      "sup"       = "SP"
    ),
    treatment = factor(treatment, levels = c("MS","AS","CT","CR","SP")),
    
    population = str_trim(population),
    population = toupper(population),
    population = recode(
      population,
      "SNT" = "Trioecious",
      "MR"  = "Hermaphrodite"
    ),
    population = factor(population, levels = c("Trioecious", "Hermaphrodite"))
  )


############################################################
# 3. GLMM models
############################################################

message("[04] Fitting GLMM models...")

modelo_glmm <- glmer(
  cbind(germinated_total, failures) ~ population + treatment + (1 | plant),
  data = germ_master,
  family = binomial
)

modelo_glmm_int <- glmer(
  cbind(germinated_total, failures) ~ population * treatment + (1 | plant),
  data = germ_master,
  family = binomial
)


############################################################
# 4. Create output directories
############################################################

dir.create("figs", showWarnings = FALSE)


############################################################
# FIGURE 1 — Ovules and pollen
############################################################

message("[04] Generating Figure 1...")

datos_long <- datos %>%
  pivot_longer(
    cols = c(ovules, pollen_total),
    names_to = "trait",
    values_to = "value"
  ) %>%
  mutate(trait = recode(trait,
                        "ovules" = "Ovules",
                        "pollen_total" = "Pollen"))

fig1 <- ggplot(datos_long, aes(x = population, y = value, fill = population)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  facet_wrap(~ trait, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Ovule number and pollen production in two populations",
    x = "Population",
    y = "Value"
  )

ggsave("figs/Figure1_ovules_pollen.pdf", fig1, width = 7, height = 5)


############################################################
# FIGURE 2 — P/O ratio
############################################################

message("[04] Generating Figure 2...")

fig2 <- ggplot(datos, aes(x = population, y = PO, fill = population)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Pollen–ovule ratio (P/O) in two populations",
    x = "Population",
    y = "P/O ratio"
  )

ggsave("figs/Figure2_PO_ratio.pdf", fig2, width = 7, height = 5)


############################################################
# FIGURE 3 — Predicted germination probability
############################################################

message("[04] Generating Figure 3...")

emm_pop <- emmeans(modelo_glmm, ~ population, type = "response") %>%
  as.data.frame()

fig3 <- ggplot(emm_pop, aes(x = population, y = prob, fill = population)) +
  geom_col(alpha = 0.7) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Predicted germination probability (GLMM without interaction)",
    x = "Population",
    y = "Germination probability"
  )

ggsave("figs/Figure3_germination_GLMM_no_interaction.pdf", fig3, width = 7, height = 5)


############################################################
# FIGURE 4 — Germination by treatment
############################################################

message("[04] Generating Figure 4...")

fig4 <- ggplot(germ_master,
               aes(x = treatment, y = germination_prop, fill = treatment)) +
  geom_violin(alpha = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "white") +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.15, color = "white") +
  facet_wrap(~ population) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  labs(
    title = "Germination proportion across pollination treatments",
    x = "Treatment",
    y = "Germination proportion"
  )

ggsave("figs/Figure4_germination_treatments.pdf", fig4, width = 8, height = 6)


############################################################
# FIGURE S1 — Heatmap of intra-population contrasts
############################################################

message("[04] Generating Supplementary Figure S1...")

emm_pop_int <- emmeans(modelo_glmm_int, ~ treatment | population)

contr_intra_int <- pairs(emm_pop_int) %>%
  as.data.frame() %>%
  mutate(
    population = as.character(population),
    contrast   = as.character(contrast),
    type       = "Intra"
  )

split_contrast <- function(x) {
  parts <- str_split(x, " - ", simplify = TRUE)
  tibble(group1 = parts[,1], group2 = parts[,2])
}

intra_df <- contr_intra_int %>%
  bind_cols(split_contrast(.$contrast)) %>%
  select(population, group1, group2, p.value, type) %>%
  mutate(
    population = factor(population,
                        levels = c("Trioecious", "Hermaphrodite")),
    group1 = factor(group1, levels = c("MS","AS","CT","CR","SP")),
    group2 = factor(group2, levels = c("MS","AS","CT","CR","SP"))
  ) %>%
  arrange(population)

figS1 <- ggplot(intra_df, aes(group1, group2, fill = p.value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", p.value)), size = 3) +
  scale_fill_distiller(
    palette = "Blues",
    direction = 1,
    name = "p-value",
    na.value = "white"
  ) +
  facet_wrap(~ population, nrow = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Pairwise contrasts of germination (p-values, interaction model)",
    x = "Treatment",
    y = "Treatment"
  ) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("figs/FigureS1_heatmap_intra_interaction.pdf", figS1, width = 10, height = 5)


############################################################
# DONE
############################################################

message("\n[04] All figures generated successfully.\n")
