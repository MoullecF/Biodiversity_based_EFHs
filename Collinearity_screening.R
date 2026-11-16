###############################################################################
#
# Collinearity screening of predictors
#
###############################################################################

rm(list = ls())

library(dplyr)
library(ggplot2)
library(ggcorrplot)
library(car)
library(patchwork)

# -----------------------------------------------------------------------------
# 0. Settings and data
# -----------------------------------------------------------------------------
data_file <- "./Inputs_HMSC/Old_env_dataset.RData" # Dataset including sea bottom salinity and seabed substrate. 
out_dir   <- "./Outputs/Collinearity/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

env <- get(load(data_file))
# Remove any potential duplicate based on id (to be sure, keep only unique ids)
env <- env[!duplicated(env$id), ]

retained   <- c("Depth", "SST", "SBT", "SSS", "Chla", "Fishing", "Gravity")
candidates <- c(retained, "botSo") # continuous
substrate  <- "MSFD_BBHT" # categorical

# Convert substrate to factor
env[[substrate]] <- factor(env[[substrate]])

r_threshold   <- 0.7
vif_threshold <- 3

# -----------------------------------------------------------------------------
# 1. Pairwise Pearson correlations, before and after selection
# -----------------------------------------------------------------------------

cor_before <- cor(env[, candidates], method = "pearson", use = "pairwise.complete.obs")
cor_after  <- cor(env[, retained],   method = "pearson", use = "pairwise.complete.obs")

# Pairs above the threshold
flag_pairs <- function(cm, th) {
  idx <- which(abs(cm) > th & upper.tri(cm), arr.ind = TRUE)
  data.frame(var1 = rownames(cm)[idx[, 1]], var2 = colnames(cm)[idx[, 2]],
             r = round(cm[idx], 3), row.names = NULL)
}
cat("\nPairs with |r| >", r_threshold, "(all continuous candidates)\n")
print(flag_pairs(cor_before, r_threshold))
cat("\nPairs with |r| >", r_threshold, "(retained predictors)\n")
print(flag_pairs(cor_after, r_threshold))

# hc.order = FALSE so that both panels keep the same variable order
plot_cor <- function(cm, tag) {
  ggcorrplot::ggcorrplot(cm, hc.order = FALSE, type = "lower", lab = TRUE, lab_size = 3, colors = c("#2166AC", "white", "#B2182B")) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1), name = "Pearson\ncorrelation") +
    labs(tag = tag)
}
p_before <- plot_cor(cor_before, "A")
p_after  <- plot_cor(cor_after,  "B")

 ggsave(file.path("./Figures/Correlation_matrix_before_after.png"), 
        p_before + p_after + plot_layout(guides = "collect"),
        width = 34, height = 15, units = "cm", dpi = 400)

# -----------------------------------------------------------------------------
# 2. Variance Inflation Factors before and after selection
# -----------------------------------------------------------------------------

# VIF depends only on the predictors, not on the response: a dummy response is used to build the linear model that car::vif() needs.
vif_table <- function(vars, data) {
  d <- data[, vars, drop = FALSE]
  d$.y <- seq_len(nrow(d))
  v <- car::vif(lm(.y ~ ., data = d))
  if (is.matrix(v)) {
    # factor is present: GVIF
    data.frame(covariate = rownames(v), VIF = round(v[, "GVIF"], 2), Df = v[, "Df"],
               GVIF_scaled = round(v[, "GVIF^(1/(2*Df))"], 2), row.names = NULL)
  } else {
    data.frame(covariate = names(v), VIF = round(as.numeric(v), 2), Df = 1,
               GVIF_scaled = round(sqrt(as.numeric(v)), 2), row.names = NULL)
  }
}

vif_all        <- vif_table(c(candidates, substrate), env)
vif_continuous <- vif_table(candidates, env)
vif_retained   <- vif_table(retained, env)

print(vif_all);print(vif_continuous);print(vif_retained)

tab_vif <- full_join(
  vif_all      %>% dplyr::select(covariate, VIF_before = VIF, Df),
  vif_retained %>% dplyr::select(covariate, VIF_after = VIF),
  by = "covariate")

print(tab_vif)
write.csv(tab_vif, file.path(out_dir, "VIF_before_after.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Seabed substrate vs continuous predictors
# -----------------------------------------------------------------------------

assoc_substrate <- function(v, data, f) {
  x <- data[[v]]; g <- data[[f]]
  ok <- complete.cases(x, g); x <- x[ok]; g <- droplevels(g[ok])
  kw <- kruskal.test(x, g)
  n <- length(x); k <- nlevels(g)
  ss_tot <- sum((x - mean(x))^2)
  ss_btw <- sum(tapply(x, g, function(z) length(z) * (mean(z) - mean(x))^2))
  eta2 <- ss_btw / ss_tot
  data.frame(covariate = v, n = n,
             KW_H = round(unname(kw$statistic), 1), KW_df = unname(kw$parameter),
             KW_p = format.pval(kw$p.value, digits = 3),
             eta2 = round(eta2, 3), eta = round(sqrt(eta2), 3),
             epsilon2 = round((unname(kw$statistic) - k + 1) / (n - k), 3))
}

tab_substrate <- do.call(rbind, lapply(candidates, assoc_substrate, data = env, f = substrate))
tab_substrate <- tab_substrate[order(-tab_substrate$eta2), ]
print(tab_substrate, row.names = FALSE)
write.csv(tab_substrate, file.path(out_dir, "substrate_association.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Substrate classes: number of hauls and depth distribution
# -----------------------------------------------------------------------------

tab_classes <- env %>%
  group_by(.data[[substrate]]) %>%
  summarise(n_hauls = n(),
            depth_median = round(median(Depth), 1),
            depth_q25 = round(quantile(Depth, 0.25), 1),
            depth_q75 = round(quantile(Depth, 0.75), 1),
            pct_shallower_200m = round(100 * mean(Depth <= 200), 1),
            .groups = "drop") %>%
  mutate(pct_hauls = round(100 * n_hauls / sum(n_hauls), 2)) %>%
  arrange(desc(n_hauls))
print(as.data.frame(tab_classes))
write.csv(tab_classes, file.path(out_dir, "substrate_classes.csv"), row.names = FALSE)

# Plot depth by substrate class
p_depth <- ggplot(env, aes(x = .data[[substrate]], y = Depth)) +
  geom_boxplot(outlier.size = 0.3, fill = "grey85") +
  scale_y_reverse() +
  labs(x = "Seabed substrate (MSFD broad habitat type)", y = "Depth (m)") +
  theme_bw() +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 12),
        panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "Depth_by_substrate.png"), p_depth, width = 14, height = 11, units = "cm", dpi = 400)
