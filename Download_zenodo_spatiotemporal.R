# Title: Fetch spatio-temporal hurdle model output from Zenodo

# Summary: Downloads the hurdle model raster stack from Zenodo using the zenodor package.
# - Installs zenodor from GitHub if missing.
# - Ensures Outputs/Spatio_temporal_prediction exists.
# - Skips download when the target .Rdata file is already present.

# ---- Package setup ----
required_package <- "zenodor"
if (!requireNamespace(required_package, quietly = TRUE)) {
  # Install from GitHub when zenodor is not already available locally
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("FRBCesab/zenodor", force = TRUE)
}
library(zenodor)

# ---- Configuration ----
record_id <- "18632166"
target_file_name <- "r_stack_Hurdle_0.05_0_1000_19992021.Rdata"
target_dir <- file.path("Outputs", "Spatio_temporal_prediction")
if (!dir.exists(target_dir)) {
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
}
target_path <- file.path(target_dir, target_file_name)

# ---- Download ----
if (file.exists(target_path)) {
  message("File already present: ", target_path)
} else {
  # Use the zenodor helper to grab the target file; fail fast if the helper is missing
  if (!"zen_download_files" %in% getNamespaceExports("zenodor")) {
    stop("zenodor::zen_download_files is not available. Please update zenodor.")
  }

  zenodor::zen_download_files(
    record_id = record_id,
    files = target_file_name,
    path = target_dir,
    progress = TRUE
  )

  if (file.exists(target_path)) {
    message("Download completed: ", target_path)
  } else {
    warning("Download attempted but the file is still missing: ", target_path)
  }
}
