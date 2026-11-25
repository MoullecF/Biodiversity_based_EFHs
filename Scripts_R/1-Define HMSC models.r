# ------------------------------------------------------------------------------
# Define HMSC models for fish species distribution and abundance
# ------------------------------------------------------------------------------

# Load libraries and set paths
source("Scripts_R/0-Load libraries.R")

localDir <- "."
ModelDir <- file.path(localDir, "Models")
DataDir <- file.path(localDir, "Inputs_HMSC")

# Load input data (S, X, Y, P) and traits
load(file = file.path(DataDir, "allData.RData"))
load(file = file.path(DataDir, "Traits.RData"))

# ------------------------------------------------------------------------------
# Filter species and align tables
# ------------------------------------------------------------------------------
source("Scripts_R/Keep_species_more_than_1percent_occ.R")

# Excluding rare taxa based on previous script output
load(file = file.path(DataDir, "selected_sp.RData"))
Y <- Y[, colnames(Y) %in% selected_sp]

# Remove hauls with no selected species across S, X, Y
S_rows <- rowSums(Y)
hauls_to_be_removed <- names(S_rows[S_rows == 0])

Y <- Y[!(rownames(Y) %in% hauls_to_be_removed), ]
S <- S[!(S$id %in% hauls_to_be_removed), ]
X <- X[!(rownames(X) %in% hauls_to_be_removed), ]

# Cast variables to factors
S$Year <- as.factor(S$Year)
S$id <- as.factor(S$id)

# ------------------------------------------------------------------------------
# Define fixed effects (environmental and anthropogenic covariates)
# ------------------------------------------------------------------------------

XFormula <- ~ poly(Depth, degree = 2, raw = TRUE) +
  poly(SBT, degree = 2, raw = TRUE) +
  poly(SST, degree = 2, raw = TRUE) +
  poly(SSS, degree = 2, raw = TRUE) +
  poly(Chla, degree = 2, raw = TRUE) +
  Fishing + Gravity
print(XFormula)

# ------------------------------------------------------------------------------
# Correlation matrix of covariates
# ------------------------------------------------------------------------------

cor_mat <- cor(X, method = "pearson", use = "pairwise.complete.obs")
cor_covariates <- ggcorrplot::ggcorrplot(cor_mat,
  hc.order = TRUE,
  type = "lower", lab = TRUE, lab_size = 3,
  colors = c("#2166AC", "white", "#B2182B")) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1), name = "Pearson\ncorrelation")

#ggplot2::ggsave(cor_covariates, filename = "Figures/Correlation_Matrix.png", width = 30, height = 15, units = "cm", dpi = 400)

# Variance inflation factor (VIF) to check for multicollinearity
usdm::vifstep(as.data.frame(X), th=5, size = 12648)
 
# ------------------------------------------------------------------------------
# Define trait effects (remove Caudal shape and subset to retained species)
# ------------------------------------------------------------------------------

Tr <- Tr[, -which(names(Tr) %in% c("Caudal_shape"))]
Tr_tmp <- Tr
Tr_tmp <- Tr_tmp[rownames(Tr_tmp) %in% colnames(Y), ]

TrFormula <- as.formula(paste(
  "~",
  paste(
    c("Trophic_level", "Min_size", "Max_size", "Growth_coefficient",
      "Trophic_guild", "Body_shape", "Water_column_position"),
    collapse = "+"
  )
))
print(TrFormula)

# ------------------------------------------------------------------------------
# Study design and spatial random effect setup
# ------------------------------------------------------------------------------

sum(duplicated(S[, c(2, 3)])) # Lat/Long duplicates

# Ensure unique coordinates by perturbing duplicates slightly
row.names(S) <- seq(1, nrow(S), 1)

dupl <- which(duplicated(S[, c("Long", "Lat")]))
for (i in 1:length(dupl)) {
  ds <- S[which(S$Long == S$Long[dupl[i]] & S$Lat == S$Lat[dupl[i]]), ]
  ds <- ds[-1, ]
  lon <- ds$Long[1]
  lat <- ds$Lat[1]
  for (n in 1:nrow(ds)) {
    lon <- lon + 0.000001
    lat <- lat + 0.000001
    S[dupl[i], "Long"] <- lon
    S[dupl[i], "Lat"] <- lat
  }
}

# Build study design table
studyDesign <- data.frame(
  sample = as.factor(S$id),
  year = as.factor(S$Year),
  grid.cell = as.factor(S$grid.id),
  grid.lon = as.factor(round(S$grid.lon, 5)),
  grid.lat = as.factor(round(S$grid.lat, 5))
)

# ------------------------------------------------------------------------------
# Random effects (temporal and spatial)
# ------------------------------------------------------------------------------

rL.year <- HmscRandomLevel(units = levels(studyDesign$year))
rL.year <- setPriors(rL.year, nfMin = 2, nfMax = 5)

xy <- data.frame(S[match(unique(S$grid.id), S$grid.id), c("grid.lon", "grid.lat")])
rownames(xy) <- unique(S$grid.id)

rL.location <- HmscRandomLevel(sData = xy, longlat = TRUE, sMethod = "NNGP", nNeighbours = 10)
rL.location <- setPriors(rL.location, nfMin = 2, nfMax = 5)

# ------------------------------------------------------------------------------
# Build models: presence/absence and abundance
# ------------------------------------------------------------------------------

# Presence/absence (probit)
Ypa <- Y
Ypa[is.na(Ypa)] <- 0
Ypa[Ypa > 0] <- 1

mpa <- Hmsc(
  Y = Ypa,
  XData = X,
  XFormula = XFormula,
  TrFormula = TrFormula,
  TrData = Tr_tmp,
  phyloTree = P,
  studyDesign = studyDesign,
  YScale = FALSE,
  XScale = TRUE,
  TrScale = TRUE,
  ranLevels = list(year = rL.year, grid.cell = rL.location),
  distr = "probit"
)

# Abundance (log-normal)
Yabu <- Y
Yabu[Yabu == 0] <- NA
Yabu <- log(Yabu)

mabu <- Hmsc(
  Y = Yabu,
  XData = X,
  XFormula = XFormula,
  TrFormula = TrFormula,
  TrData = Tr_tmp,
  phyloTree = P,
  studyDesign = studyDesign,
  YScale = TRUE,
  XScale = TRUE,
  TrScale = TRUE,
  ranLevels = list(year = rL.year, grid.cell = rL.location),
  distr = "normal"
)

# ------------------------------------------------------------------------------
# Save model objects
# ------------------------------------------------------------------------------

models <- list(mpa, mabu)
modelnames <- c("probit", "normal")

# save(models, modelnames, file = file.path(ModelDir, "unfitted_models"))