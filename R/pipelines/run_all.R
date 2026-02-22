#!/usr/bin/env Rscript

cat("\n========================================\n")
cat("   Running full analysis pipeline\n")
cat("========================================\n\n")

if (!dir.exists("scripts")) {
  stop("The 'scripts/' directory was not found. Please run this script from the project root.")
}

run_step <- function(step_name, file) {
  cat(paste0("\n--- ", step_name, " ---\n"))
  tryCatch(
    {
      source(file)
      cat("✓ Completed successfully.\n")
    },
    error = function(e) {
      cat("✗ Error in ", file, ":\n", e$message, "\n", sep = "")
      stop("Pipeline aborted due to error.")
    }
  )
}

run_step("Model selection for ovules and pollen", "scripts/00_model_selection_ovules_pollen.R")
run_step("Analysis of ovules and pollen", "scripts/01_analysis_ovules_pollen.R")
run_step("Germination analysis", "scripts/02_analysis_germination.R")
run_step("P/O vs germination analysis", "scripts/03_PO_vs_germination.R")
run_step("Figure generation", "scripts/04_figures.R")

cat("\n========================================\n")
cat(" All analyses and figures completed.\n")
cat("========================================\n\n")
