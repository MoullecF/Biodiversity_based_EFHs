nocc <- apply(Y > 0, 2, sum)
totabu <- apply(Y, 2, sum)
threshold <- floor(nrow(Y) * 0.01)

# Number of species kept
N_kept_sp <- sum(nocc > threshold)
print(paste0(N_kept_sp, " species will be retained"))

# Percentage of abundance conserved
percent_abund_kept <- sum(Y[, nocc > threshold]) / sum(Y) * 100
print(paste0("The ", N_kept_sp, " species retained account for ", round(percent_abund_kept,2)," % of all abundances")) # nolint

# Select only species with occurrence in at least 1% of grid cells
abu <- Y[, nocc > threshold]

selected_sp <- colnames(abu)
save(selected_sp, file = "Inputs_HMSC/selected_sp.RData")
