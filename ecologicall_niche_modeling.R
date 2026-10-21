#ECOLOGICAL NICHE MODELING: Aedes albopictus in Turkey
# Beste Yaprak - Graduation Thesis
# Method: MaxEnt | Climate + Population + Urban + OSM Roads
# ============================================================


# ============================================================
# SECTION 0: PACKAGES
# ============================================================

library(sf)
library(ggplot2)
library(terra)
library(dplyr)
library(maxnet)
library(tidyterra)
library(geodata)
library(corrplot)
library(osmdata)

# ============================================================
# SECTION 1: TURKEY BOUNDARIES
# ============================================================

turkey_provinces <- st_as_sf(gadm(country = "TUR", level = 1, path = tempdir()))
turkey_boundary  <- st_as_sf(gadm(country = "TUR", level = 0, path = tempdir()))

# ============================================================
# SECTION 2: GBIF OCCURRENCE DATA
# ============================================================
library(sf)
library(dplyr)

gbif_raw <- read.csv("aedes_albopictus_gbif.csv", sep = "\t", header = TRUE)
aedes_sub <- gbif_raw[!is.na(gbif_raw$decimalLongitude) & 
                        !is.na(gbif_raw$decimalLatitude), ]

aedes_sub <- aedes_sub[aedes_sub$decimalLongitude >= 4 & aedes_sub$decimalLongitude <= 47 &
                         aedes_sub$decimalLatitude >= 34 & aedes_sub$decimalLatitude <= 47, ]

# 3. Spatial thinning
aedes_sub <- aedes_sub[!duplicated(
  round(aedes_sub[, c("decimalLongitude", "decimalLatitude")], 2)), ]
aedes_pts <- st_as_sf(aedes_sub,
                      coords = c("decimalLongitude", "decimalLatitude"),
                      crs = 4326)

cat("Filtered # presence :", nrow(aedes_pts), "\n")

##------------------------OLD VERSION--------------------------

gbif_raw <- read.csv("aedes_albopictus_gbif.csv", sep = "\t", header = TRUE)

aedes_sub <- gbif_raw[!is.na(gbif_raw$decimalLongitude) &
                        !is.na(gbif_raw$decimalLatitude), ]
aedes_sub <- aedes_sub[aedes_sub$countryCode %in% c("TR", "IT", "GR", "BG", "GE"), ]
aedes_sub <- aedes_sub[!duplicated(
  round(aedes_sub[, c("decimalLongitude", "decimalLatitude")], 2)), ]

aedes_pts <- st_as_sf(aedes_sub,
                      coords = c("decimalLongitude", "decimalLatitude"),
                      crs = 4326)
cat("Occurrence points:", nrow(aedes_pts), "\n")


# ============================================================
# SECTION 3: DISTRIBUTION MAP
# ============================================================

turkey_only_pts <- aedes_pts[
  st_within(aedes_pts, turkey_boundary, sparse = FALSE)[,1], ]

turkey_coastal_buffer <- st_buffer(turkey_boundary, dist = 0.02)
turkey_only_pts_clean <- aedes_pts[st_intersects(aedes_pts, turkey_coastal_buffer, sparse = FALSE)[, 1], ]

occupied <- st_join(turkey_provinces, turkey_only_pts_clean, left = FALSE)
p_dist <- ggplot() +
  geom_sf(data = turkey_provinces, fill = "white", color = "gray85") +
  geom_sf(data = occupied, fill = "#8B0000", alpha = 0.3, color = "#8B0000") +
  geom_sf(data = turkey_only_pts_clean, color = "cyan", size = 1.2, alpha = 0.7) +
  coord_sf(xlim = c(26, 45), ylim = c(36, 42)) +
  labs(
    title    = NULL,
    subtitle = NULL,
    x        = NULL,
    y        = NULL,
    caption  = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.background = element_rect(fill = "aliceblue", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )
print(p_dist)
ggsave("aedes_distribution_map_NOCAPTIONS.png", plot = p_dist, width = 10, height = 6, dpi = 600)


# ============================================================
# SECTION 4: CLIMATE VARIABLES
# ============================================================
# BIO1 and BIO6 dropped after correlation analysis (r > 0.8 with BIO11)
# BIO11 retained: mean temp coldest quarter = key overwintering constraint
# Final: BIO11, BIO12, BIO14, BIO15

unlink(paste0(tempdir(), "/climate"), recursive = TRUE, force = TRUE)

tile_turkey <- worldclim_tile(
  var  = "bio",
  lon  = 35,  
  lat  = 39,   
  path = tempdir(),
  res  = 2.5    
)

# 2. REGION: Italy, Central Europe, Western Balkans
tile_italy <- worldclim_tile(
  var  = "bio",
  lon  = 12,   
  lat  = 42,   
  path = tempdir(),
  res  = 2.5    
)

# 3. merge
bioclim_regional <- terra::mosaic(tile_turkey, tile_italy)
# Crop to your regional extent
study_extent <- ext(4, 47, 34, 47)   
bioclim_regional <- crop(bioclim_regional, study_extent)

bio_selected <- bioclim_regional[[c(11, 12, 14, 15)]]
names(bio_selected) <- c("BIO11_MeanTempColdQ", "BIO12_AnnPrecip",
                         "BIO14_PrecipDriestM", "BIO15_PrecipSeasonality")
# ============================================================
# SECTION 9.1: MAXENT - CLIMATE-ONLY TRAINING (FUNDAMENTAL NICHE)
# ============================================================
cat("Extracting bioclimatic values for all cross-border occurrences...\n")

aedes_pts <- st_transform(aedes_pts, crs = st_crs(bio_selected))

# Climate only data set 
pres_vals_clim <- terra::extract(bio_selected, vect(aedes_pts), ID = FALSE)
set.seed(42)
# 10.000 random background points 
bg_vals_clim <- spatSample(bio_selected, 10000, "random", na.rm = TRUE, as.df = TRUE)

# data control
pres_vals_clim_clean <- pres_vals_clim[complete.cases(pres_vals_clim), ]
bg_vals_clim_clean   <- bg_vals_clim[complete.cases(bg_vals_clim), ]


p_vec_clim <- c(rep(1, nrow(pres_vals_clim_clean)), rep(0, nrow(bg_vals_clim_clean)))
mdata_clim <- rbind(pres_vals_clim_clean, bg_vals_clim_clean)
cat("Kontrol:\n")
library(pROC)

# ==============================================================================
#  5-Fold Cross-Validation 
# ==============================================================================
library(pROC)
cat("\n--- 5-Fold Cross-Validation ---\n")

set.seed(42)
k_folds <- 5

pres_folds <- sample(rep(1:k_folds, length.out = nrow(pres_vals_clim_clean)))
bg_folds   <- sample(rep(1:k_folds, length.out = nrow(bg_vals_clim_clean)))

cv_auc_scores <- numeric(k_folds)

for (i in 1:k_folds) {
  # 1. 80% train, 20% test for each fold
  train_pres <- pres_vals_clim_clean[pres_folds != i, ]
  test_pres  <- pres_vals_clim_clean[pres_folds == i, ]
  
  train_bg   <- bg_vals_clim_clean[bg_folds != i, ]
  test_bg    <- bg_vals_clim_clean[bg_folds == i, ]
  
  # combine data
  train_data <- rbind(train_pres, train_bg)
  train_p    <- c(rep(1, nrow(train_pres)), rep(0, nrow(train_bg)))
  
  test_data  <- rbind(test_pres, test_bg)
  test_p     <- c(rep(1, nrow(test_pres)), rep(0, nrow(test_bg)))
  
  fold_model <- maxnet(p = train_p, data = train_data)
  test_preds <- predict(fold_model, test_data, type = "cloglog")
  
  # Test AUC 
  roc_obj <- pROC::roc(response = test_p, predictor = as.numeric(test_preds), quiet = TRUE)
  cv_auc_scores[i] <- as.numeric(pROC::auc(roc_obj))
  
  cat(sprintf("Fold %d Test AUC: %.4f\n", i, cv_auc_scores[i]))
}

# mean and standart deviation for 5-fold 
cat(sprintf("\n>>> Mean Test AUC: %.4f ± %.4f <<<\n\n", mean(cv_auc_scores), sd(cv_auc_scores)))

# Fundamental Niche Model
aedes_climate_model <- maxnet(p = p_vec_clim, data = mdata_clim)

cat("Success! Pure climate model trained.\n")
cat("Total international presences used:", nrow(pres_vals_clim_clean), "\n")

library(pROC)

# 1. Training (Apparent) AUC 
pred_train <- predict(aedes_climate_model, mdata_clim, type = "cloglog")
roc_obj    <- roc(response = p_vec_clim, predictor = as.numeric(pred_train), quiet = TRUE)
train_auc  <- as.numeric(auc(roc_obj))

cat("\n================ MODEL PERFORMANCE ================\n")
cat(sprintf("Training (Apparent) AUC : %.4f\n", train_auc))
cat(sprintf("5-Fold Cross-Validation Test AUC: %.4f ± %.4f\n", mean(cv_auc_scores), sd(cv_auc_scores)))
cat("===================================================\n\n")

# 2. Permutation Importance
cat("Calculating Permutation Importance...\n")
calc_importance <- function(model, data, y_true, base_auc) {
  var_names <- colnames(data)
  imp_list <- numeric(length(var_names))
  names(imp_list) <- var_names
  
  set.seed(42)
  for (v in var_names) {
    shuffled_data <- data
    shuffled_data[[v]] <- sample(shuffled_data[[v]])
    pred_shuff <- predict(model, shuffled_data, type = "cloglog")
    shuff_auc <- as.numeric(auc(roc(response = y_true, predictor = as.numeric(pred_shuff), quiet = TRUE)))
    imp_list[v] <- max(0, base_auc - shuff_auc)
  }
  
  imp_percent <- (imp_list / sum(imp_list)) * 100
  return(sort(imp_percent, decreasing = TRUE))
}

var_importance <- calc_importance(aedes_climate_model, mdata_clim, p_vec_clim, train_auc)

cat("--- Permutation Importance (%) ---\n")
print(round(var_importance, 2))

# save it in CSV
df_importance <- data.frame(
  Variable = names(var_importance),
  Importance_Percent = round(as.numeric(var_importance), 2)
)
write.csv(df_importance, "Table_Variable_ImportanceNEW.csv", row.names = FALSE)

# 3. ROC Graph ( AUC of Test and Train as academic version )
png("Figure_Model_ROC_Curve_600dpiNEW.png", width = 1800, height = 1800, res = 300)
par(mar = c(4.5, 4.5, 2, 2))
plot(roc_obj, 
     col = "#2166AC", 
     lwd = 2.5, 
     main = "",
     xlab = "1 - Specificity (False Positive Rate)", 
     ylab = "Sensitivity (True Positive Rate)",
     cex.lab = 1.1, font.lab = 2)
abline(a = 0, b = 1, lty = 2, col = "gray50", lwd = 1.2)
dev.off()

# 4. Response Curves
png("Figure_Response_Curves_600dpiNEW.png", width = 2400, height = 2000, res = 300)
plot(aedes_climate_model, type = "cloglog", col = "darkred", lwd = 2)
dev.off()

cat("\nModel değerlendirme, ROC ve tepki eğrileri başarıyla kaydedildi!\n")


#==============================================================
#==============================================================
# MAP CLİMATE ONLY
#==============================================================
#==============================================================

bio_turkey <- terra::crop(bio_selected, turkey_boundary)
bio_turkey <- terra::mask(bio_turkey, turkey_boundary)

suit_climate <- predict(bio_turkey, aedes_climate_model, type = "cloglog", na.rm = TRUE)
names(suit_climate) <- "Suitability"
p_clim <- ggplot() +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "aliceblue", color = NA)) +
  geom_spatraster(data = suit_climate) +
  scale_fill_viridis_c(option = "inferno", na.value = "transparent",
                       name = "Climate\nSuitability", limits = c(0, 1)) +
  geom_sf(data = turkey_provinces, fill = "transparent",
          color = "white", linewidth = 0.15) +
  geom_sf(data = aedes_pts, color = "cyan", size = 0.9, alpha = 0.7) +
  coord_sf(xlim = c(26, 45), ylim = c(36, 42), expand = FALSE) +
  scale_x_continuous(labels = function(x) paste0(x, "°E")) +
  scale_y_continuous(labels = function(y) paste0(y, "°N")) +
  labs(
    title    = NULL,
    subtitle = NULL,
    x = NULL, y = NULL,
    caption  = NULL
  ) +
  theme(plot.title = element_text(face = "bold", size = 14),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_line(color = "white", linetype = "dashed",
                                        linewidth = 0.2))
print(p_clim)
ggsave("aedes_pure_climate_suitability_NOCAPTİONS.png", plot = p_clim,
       width = 12, height = 7, dpi = 600)
#NO POINTS
p_clim_no_pts <- ggplot() +
  geom_spatraster(data = suit_climate) +
  scale_fill_viridis_c(
    option   = "inferno", 
    na.value = "transparent",
    name     = "Climate\nSuitability", 
    limits   = c(0, 1)
  ) +
  geom_sf(data = turkey_provinces, fill = "transparent", color = "white", linewidth = 0.15) +
  coord_sf(xlim = c(26, 45), ylim = c(36, 42), expand = FALSE) +
  scale_x_continuous(labels = function(x) paste0(x, "°E")) +
  scale_y_continuous(labels = function(y) paste0(y, "°N")) +
  labs(title = NULL, subtitle = NULL, caption = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "white", linetype = "dashed", linewidth = 0.2),
    legend.position  = "right",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8),
    axis.text        = element_text(size = 8, color = "black")
  )

print(p_clim_no_pts)
ggsave("Figure3_Climate_Suitability_NoPoints_600dpiNEW.png", plot = p_clim_no_pts,
       width = 11, height = 6.5, dpi = 600)
# countries where the presence points obtained from
table(aedes_pts$countryCode[complete.cases(pres_vals_clim)])

nrow(pres_vals_clim_clean)
# climate summary of Italy points
summary(pres_vals_clim_clean[aedes_pts$countryCode[complete.cases(pres_vals_clim)] == "IT", ])
# 1. Appendix table 
country_counts <- as.data.frame(table(Country = aedes_pts$countryCode[complete.cases(pres_vals_clim)]))
colnames(country_counts) <- c("Country", "Occurrence_Count")
write.csv(country_counts, "Table_Appendix_Country_Occurrences.csv", row.names = FALSE)

# 2. save the summary 
it_clim_summary <- do.call(cbind, lapply(pres_vals_clim_clean[aedes_pts$countryCode[complete.cases(pres_vals_clim)] == "IT", ], summary))
write.csv(it_clim_summary, "Table_Appendix_Italy_Climate_Summary.csv")
# ============================================================
# after 9.1 run:SECTION 5: POPULATION DENSITY AND URBAN LAYERS (TURKEY ONLY)
# ============================================================

pop_global     <- geodata::population(year = 2020, res = 2.5, path = tempdir())


pop_turkey     <- terra::crop(pop_global, turkey_boundary)
pop_turkey     <- terra::mask(pop_turkey, turkey_boundary)
pop_turkey     <- terra::resample(pop_turkey, bio_turkey, method = "bilinear")
pop_turkey_log <- log1p(pop_turkey)
names(pop_turkey_log) <- "PopDensity_log"

# Threshold (> 150 person per km2)
# urban clasification via threshold (OECD/Eurostat standardı: ≥ 150 person/km2)
urban_threshold <- log1p(150)

urban_layer <- terra::classify(
  pop_turkey_log, 
  rcl = matrix(c(-Inf, urban_threshold, 0,  urban_threshold, Inf, 1), ncol = 3, byrow = TRUE)
)
names(urban_layer) <- "Urban_Binary"
cat("Population and urbanization layers created succesfully.\n")
# ============================================================
# SECTION 6: OSM ROAD NETWORK - DISTANCE TO MAJOR ROADS
# ============================================================
# SCIENTIFIC JUSTIFICATION:
# Ae. albopictus disperses passively through human vehicle transport.
# Eggs survive desiccation in used tyres and traded goods traveling
# along road networks (Eritja et al. 2017 Sci Rep; Kraemer et al. 2019).
# Motorways carry long-distance trade; trunk and primary roads carry
# regional movement. Together they define the dispersal corridor network.
# We convert road geometry to an inverse-distance raster so that pixels
# closer to roads receive higher values = higher introduction pressure.

# ============================================================
# FIXED SECTION 6: LOCAL HIGH RESOLUTION TRANSPORTATION NETWORK
# ============================================================
cat("reading Turkey transportation network doc (.shp)...\n")

roads_lines <- st_read("gis_osm_roads_free_1.shp")
names(roads_lines)
unique(roads_lines$fclass)
# filter for main roads (Kraemer,2015)
roads_lines <- roads_lines[roads_lines$fclass %in% c("motorway", "trunk", "primary"), ]

roads_lines <- st_transform(roads_lines, crs = st_crs(turkey_boundary))
roads_turkey <- st_intersection(roads_lines, turkey_boundary)

cat("Total road segments merged successfully:", nrow(roads_turkey), "\n")

# Convert to raster distance layer
roads_turkey_clean <- sf::st_make_valid(roads_turkey)
roads_vect         <- terra::vect(roads_turkey_clean)

# field = 1 and background = 0, bypass
roads_rast <- terra::rasterize(roads_vect, bio_turkey[[1]], field = 1, background = 0)
roads_rast[roads_rast == 0] <- NA

# BIOLOGICAL 15 KM BUFFER 
dist_to_roads    <- terra::distance(roads_rast)
dist_to_roads    <- terra::resample(dist_to_roads, bio_turkey[[1]], method = "bilinear")
dist_to_roads_km <- dist_to_roads / 1000

esik_mesafe   <- 15 
inv_road_norm <- 1 / (1 + (dist_to_roads_km / esik_mesafe)^2)
inv_road_norm        <- terra::mask(inv_road_norm, turkey_boundary)
names(inv_road_norm) <- "RoadProximity"
cat("Road proximity layer successfully created.\n")
# --- ROAD NETWORK VISUALIZATION MAP ---
# Produces a standalone figure showing the road layer for your methods section

p_roads <- ggplot() +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "aliceblue", color = NA)) +
  geom_spatraster(data = inv_road_norm) +
  scale_fill_viridis_c(option = "magma", na.value = "transparent",
                       name = "Road\nProximity\nIndex") +
  geom_sf(data = roads_turkey[roads_turkey$fclass == "motorway", ],
          color = "red", linewidth = 0.5, alpha = 0.9) +
  geom_sf(data = roads_turkey[roads_turkey$fclass == "trunk", ],
          color = "orange", linewidth = 0.3, alpha = 0.7) +
  geom_sf(data = roads_turkey[roads_turkey$fclass == "primary", ],
          color = "yellow", linewidth = 0.2, alpha = 0.5) +
  geom_sf(data = turkey_provinces, fill = "transparent",
          color = "white", linewidth = 0.2) +
  coord_sf(xlim = c(26, 45), ylim = c(36, 42), expand = FALSE) +
  scale_x_continuous(labels = function(x) paste0(x, "°E")) +
  scale_y_continuous(labels = function(y) paste0(y, "°N")) +
  labs(
    title    = NULL,
    subtitle = NULL,
    x = NULL, y = NULL,
    caption  = NULL
  ) +
  theme(
    plot.title   = element_text(face = "bold", size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(color = "gray20", linetype = "dashed",
                                    linewidth = 0.2)
  )
print(p_roads)
ggsave("aedes_road_proximity_layer2_NOCAPTIONS.png", plot = p_roads,
       width = 12, height = 7, dpi = 600)

# Save and export some data for reuse 
writeRaster(suit_climate,   "suit_climate_current.tif",   overwrite = TRUE)
writeRaster(inv_road_norm,  "inv_road_norm.tif",          overwrite = TRUE)
writeRaster(pop_turkey_log, "pop_turkey_log.tif",         overwrite = TRUE)
writeRaster(bio_turkey,     "bio_turkey_stack.tif",       overwrite = TRUE)
saveRDS(aedes_climate_model, "aedes_climate_model.rds")
writeRaster(suit_hybrid,   "suit_hybrid_current.tif",   overwrite = TRUE)

# sf objects
saveRDS(aedes_pts,        "aedes_pts_clean.rds")
saveRDS(turkey_only_pts_clean,  "turkey_only_pts_clean.rds")
saveRDS(turkey_provinces, "turkey_provinces.rds")

#BRINGING BACK
library(terra)
library(sf)
library(maxnet)

suit_climate        <- rast("suit_climate_current.tif")
suit_hybrid         <- rast("suit_hybrid_current.tif")
bio_turkey          <- rast("bio_turkey_stack.tif")

aedes_climate_model <- readRDS("aedes_climate_model.rds")
aedes_pts           <- readRDS("aedes_pts_clean.rds")
turkey_provinces    <- readRDS("turkey_provinces.rds")
# ============================================================
# SECTION 7: COMBINE FULL ENVIRONMENTAL STACK
# ============================================================
# Layer composition:
#   BIO11, BIO12, BIO14, BIO15  - climate (fundamental niche)
#   PopDensity_log              - human settlement density
#   Urban_Binary                - urban habitat classification
#   RoadProximity               - mechanistic dispersal corridor proxy
pop_turkey_log <- terra::resample(pop_turkey_log, bio_turkey, method = "bilinear")
urban_layer    <- terra::resample(urban_layer,    bio_turkey, method = "near")
inv_road_norm  <- terra::resample(inv_road_norm,  bio_turkey, method = "bilinear")

env_stack <- c(bio_turkey, pop_turkey_log, urban_layer, inv_road_norm)

cat("Final stack layers (Turkey Domain):\n", paste("-", names(env_stack), collapse = "\n"), "\n")

# ============================================================
# SECTION 8: VARIABLE CORRELATION CHECK
# ============================================================

bg_check   <- spatSample(bio_turkey, 5000, "random", na.rm = TRUE, as.df = TRUE)
cor_matrix <- cor(bg_check, use = "complete.obs")

png("variable_correlation_matrix.png", width = 2400, height = 2000, res = 300)
corrplot(cor_matrix, method = "color", type = "upper",
         tl.cex = 0.8, tl.col = "darkred",
         addCoef.col = "black", number.cex = 0.6,
         title = "Variable Correlation Matrix", mar = c(0, 0, 2, 0))
dev.off()

cat("\nRoad proximity correlations:\n")
print(round(cor_matrix["RoadProximity", ], 3))
# Road proximity is typically uncorrelated with climate variables
# because roads traverse all climate zones in Turkey.
# If any r > 0.7 appears, note it as a limitation in your report.
# ============================================================
# VARIABLE CORRELATION CHECK -performed in the beginings while determination 
# ============================================================
cat("Generating diagnostic correlation matrix including BIO1 and BIO6...\n")


bio_clim_diag <- bioclim_regional[[c(1, 6, 11, 12, 14, 15)]]
names(bio_clim_diag) <- c("BIO1_AnnualMeanTemp", "BIO6_MinTempColdestM", 
                          "BIO11_MeanTempColdQ", "BIO12_AnnPrecip", 
                          "BIO14_PrecipDriestM", "BIO15_PrecipSeasonality")

set.seed(42)
bg_clim_sample <- terra::spatSample(bio_clim_diag, 5000, "random", na.rm = TRUE, as.df = TRUE)
cor_clim_mat   <- cor(bg_clim_sample, use = "complete.obs")

png("Figure_Bioclim_Correlation_Matrix_600dpi.png", width = 2200, height = 2000, res = 300)
corrplot::corrplot(
  cor_clim_mat, 
  method = "color", 
  type = "upper",
  tl.cex = 0.8, 
  tl.col = "black", 
  addCoef.col = "black", 
  number.cex = 0.75,
  col = colorRampPalette(c("#B2182B", "#F7F7F7", "#2166AC"))(200),
  mar = c(0, 0, 1.5, 0)
)
dev.off()

# 2. full variable
pop_turkey_log_res <- resample(pop_turkey_log, bio_full_diagnostic, method = "bilinear")
urban_layer_res    <- resample(urban_layer,    bio_full_diagnostic, method = "near")
inv_road_norm_res  <- resample(inv_road_norm,  bio_full_diagnostic, method = "bilinear")

# Diagnostic Stack
diagnostic_stack <- c(bio_full_diagnostic, pop_turkey_log_res, urban_layer_res, inv_road_norm_res)

# 3. random 5000 sampling for pixel
bg_check_full <- spatSample(diagnostic_stack, 5000, "random", na.rm = TRUE, as.df = TRUE)

cor_matrix_full <- cor(bg_check_full, use = "complete.obs")

# 4. PNG
png("variable_correlation_matrix_full.png", width = 2800, height = 2400, res = 300)
corrplot(cor_matrix_full, 
         method = "color", 
         type = "upper",
         tl.cex = 0.7,             
         tl.col = "darkred",       
         addCoef.col = "black",    
         number.cex = 0.55,        
         title = "Full Variable Correlation Matrix (Diagnostic)", 
         mar = c(0, 0, 2, 0))
dev.off()
cat("\n--- DIAGNOSTIC CORRELATION CHECK ---\n")
cat("Correlation between BIO1 and BIO11:", round(cor_matrix_full["BIO1_AnnualMeanTemp", "BIO11_MeanTempColdQ"], 3), "\n")
cat("Correlation between BIO6 and BIO11:", round(cor_matrix_full["BIO6_MinTempColdestM", "BIO11_MeanTempColdQ"], 3), "\n")
cat("\nRoad proximity correlations across all variables:\n")
print(round(cor_matrix_full["RoadProximity", ], 3))
# ============================================================
# SECTION 9.2 (REVISED): REALIZED NICHE / HYBRID INVASION MAP
# Multi-Scale Deterministic Overlay (Kraemer et al. Framework)
# ============================================================
library(terra)
library(sf)

suit_climate   <- rast("suit_climate_current.tif")
pop_turkey_log <- rast("pop_turkey_log.tif")
inv_road_norm  <- rast("inv_road_norm.tif")

aedes_pts        <- readRDS("aedes_pts_clean.rds")
turkey_only_pts_clean  <- readRDS("turkey_only_pts_clean.rds")
turkey_provinces <- readRDS("turkey_provinces.rds")
# 1. population density normalization
pop_mm   <- minmax(pop_turkey_log)
pop_norm <- (pop_turkey_log - pop_mm[1]) / (pop_mm[2] - pop_mm[1])

# 2. Hybride formla : 50% climate + 25% pop + 25% road proximity
human_footprint <- (0.5 * pop_norm) + (0.5 * inv_road_norm)

suit_hybrid <- suit_climate * (0.5 + 0.5 * human_footprint)
suit_hybrid <- terra::mask(suit_hybrid, turkey_boundary)
names(suit_hybrid) <- "Suitability"

# 3. Figure 4: Current Realized Niche Map
p_hybrid <- ggplot() +
  geom_spatraster(data = suit_hybrid) +
  scale_fill_viridis_c(
    option   = "inferno", 
    na.value = "transparent",
    name     = "Invasion\nRisk", 
    limits   = c(0, 1)
  ) +
  geom_sf(data = turkey_provinces, fill = "transparent", color = "white", linewidth = 0.15) +
  geom_sf(data = turkey_only_pts_clean, color = "cyan", size = 0.9, alpha = 0.8) +
  coord_sf(xlim = c(26, 45), ylim = c(36, 42), expand = FALSE) +
  scale_x_continuous(labels = function(x) paste0(x, "°E")) +
  scale_y_continuous(labels = function(y) paste0(y, "°N")) +
  labs(title = NULL, subtitle = NULL, caption = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "white", linetype = "dashed", linewidth = 0.2),
    legend.position  = "right",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8),
    axis.text        = element_text(size = 8, color = "black")
  )

print(p_hybrid)

ggsave("Figure4_Current_Realized_Niche_Clean_NOCAPTIONS2_NEW.png", plot = p_hybrid,
       width = 11, height = 6.5, dpi = 600)
#NO POINTS
p_hybrid_no_pts <- ggplot() +
  geom_spatraster(data = suit_hybrid) +
  scale_fill_viridis_c(
    option   = "inferno", 
    na.value = "transparent",
    name     = "Invasion\nRisk", 
    limits   = c(0, 1)
  ) +
  geom_sf(data = turkey_provinces, fill = "transparent", color = "white", linewidth = 0.15) +
  coord_sf(xlim = c(26, 45), ylim = c(36, 42), expand = FALSE) +
  scale_x_continuous(labels = function(x) paste0(x, "°E")) +
  scale_y_continuous(labels = function(y) paste0(y, "°N")) +
  labs(title = NULL, subtitle = NULL, caption = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "white", linetype = "dashed", linewidth = 0.2),
    legend.position  = "right",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8),
    axis.text        = element_text(size = 8, color = "black")
  )

print(p_hybrid_no_pts)
ggsave("Figure4_Realized_Niche_NoPoints_600dpiNEW.png", plot = p_hybrid_no_pts,
       width = 11, height = 6.5, dpi = 600)
# ============================================================
# SECTION 10 (FORMER 12): FUTURE PROJECTIONS (SSP2-4.5 & SSP5-8.5)
# GCM: MPI-ESM1-2-HR (CMIP6) | Time: 2041-2060 (or 2061-2080)
# ============================================================

cat("DOWNLOADING...\n")

# 1. download CMIP6 selected future climatic layers (MPI-ESM1-2-HR)
future_ssp245_raw <- cmip6_world(model = 'MPI-ESM1-2-HR', ssp = '245',
                                 time = '2041-2060', var = 'bioc',
                                 res = 2.5, path = tempdir())

future_ssp585_raw <- cmip6_world(model = 'MPI-ESM1-2-HR', ssp = '585',
                                 time = '2041-2060', var = 'bioc',
                                 res = 2.5, path = tempdir())

# 2. Crop
prep_future_clim <- function(fut_raw) {
  fut_crop <- terra::crop(fut_raw, turkey_boundary)
  fut_mask <- terra::mask(fut_crop, turkey_boundary)
  
  # BIO11, BIO12, BIO14, BIO15
  fut_selected <- fut_mask[[c(11, 12, 14, 15)]]
  fut_resampled <- terra::resample(fut_selected, bio_turkey, method = "bilinear")
  names(fut_resampled) <- names(bio_turkey)
  return(fut_resampled)
}

fut_clim_245 <- prep_future_clim(future_ssp245_raw)
fut_clim_585 <- prep_future_clim(future_ssp585_raw)

# 3. MaxEnt Fundamental Niche
suit_clim_245 <- predict(fut_clim_245, aedes_climate_model, type = "cloglog", na.rm = TRUE)
suit_clim_585 <- predict(fut_clim_585, aedes_climate_model, type = "cloglog", na.rm = TRUE)

# 4. Future Hibrit Realized Niche (Multiplicative Formula)
# Kraemer socio-ecological: constant Transportation and Pop network
suit_hybrid_245 <- suit_clim_245 * (0.5 + 0.5 * human_footprint)
suit_hybrid_585 <- suit_clim_585 * (0.5 + 0.5 * human_footprint)

suit_hybrid_245 <- terra::mask(suit_hybrid_245, turkey_boundary)
suit_hybrid_585 <- terra::mask(suit_hybrid_585, turkey_boundary)

# 5.Niche shift (future - current)
change_245 <- suit_hybrid_245 - suit_hybrid
change_585 <- suit_hybrid_585 - suit_hybrid
names(change_245) <- "NicheChange"
names(change_585) <- "NicheChange"
# ============================================================
# FUTURE ABSOLUTE SUITABILITY
# ============================================================

plot_future_absolute <- function(suit_layer, filename) {
  p <- ggplot() +
    geom_spatraster(data = suit_layer) +
    scale_fill_viridis_c(
      option   = "inferno", 
      na.value = "transparent",
      name     = "Invasion\nRisk", 
      limits   = c(0, 1)
    ) +
    geom_sf(data = turkey_provinces, fill = "transparent", color = "white", linewidth = 0.15) +
    coord_sf(xlim = c(26, 45), ylim = c(36, 42), expand = FALSE) +
    scale_x_continuous(labels = function(x) paste0(x, "°E")) +
    scale_y_continuous(labels = function(y) paste0(y, "°N")) +
    labs(title = NULL, subtitle = NULL, caption = NULL, x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "aliceblue", color = NA),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.grid.major = element_line(color = "white", linetype = "dashed", linewidth = 0.2),
      legend.position  = "right",
      legend.title     = element_text(size = 9, face = "bold"),
      legend.text      = element_text(size = 8),
      axis.text        = element_text(size = 8, color = "black")
    )
  
  print(p)
  ggsave(filename, plot = p, width = 11, height = 6.5, dpi = 600)
}

# 1. Inferno Scale - comparative with current map
plot_future_absolute(suit_hybrid_245, "Figure6A_Future_Realized_SSP245_600dpiNEW.png")
plot_future_absolute(suit_hybrid_585, "Figure7A_Future_Realized_SSP585_600dpiNEW.png")

# ============================================================
# VISUALIZATION (FIGURE 6A, 6B, 7A, 7B)
# ============================================================
plot_shift <- function(change_layer, filename) {
  p <- ggplot() +
    geom_spatraster(data = change_layer) +
    scale_fill_gradient2(
      low      = "#B2182B",  
      mid      = "#F7F7F7",  
      high     = "#2166AC",  
      midpoint = 0,
      limits   = c(-0.6, 0.6),
      oob      = scales::squish,
      na.value = "transparent",
      name     = "Δ Suitability\n(Niche Shift)"
    ) +
    geom_sf(data = turkey_provinces, fill = "transparent", color = "gray40", linewidth = 0.2) +
    coord_sf(xlim = c(26, 45), ylim = c(36, 42), expand = FALSE) +
    scale_x_continuous(labels = function(x) paste0(x, "°E")) +
    scale_y_continuous(labels = function(y) paste0(y, "°N")) +
    labs(title = NULL, subtitle = NULL, caption = NULL, x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "aliceblue", color = NA),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.grid.major = element_line(color = "gray80", linetype = "dashed", linewidth = 0.2),
      legend.position  = "right",
      legend.title     = element_text(size = 9, face = "bold"),
      legend.text      = element_text(size = 8),
      axis.text        = element_text(size = 8, color = "black")
    )
  
  print(p)
  ggsave(filename, plot = p, width = 11, height = 6.5, dpi = 600)
}


# 2. Shift Map (RdBu Scale)
plot_shift(change_245, "Figure6B_Niche_Shift_SSP245_600dpiNEW.png")
plot_shift(change_585, "Figure7B_Niche_Shift_SSP585_600dpiNEW.png")

cat("Future projection maps completed.\n")
# ============================================================
# SECTION 11: LATE-CENTURY EXTREME RISK PROJECTIONS (2061-2080)
# GCM: MPI-ESM1-2-HR | Scenarios: SSP2-4.5 & SSP5-8.5
# ============================================================

cat("Getting climatic parameters for the nd of the century (2061-2080)...\n")

# 1. 2061-2080 download
late_ssp245_raw <- cmip6_world(model = 'MPI-ESM1-2-HR', ssp = '245',
                               time = '2061-2080', var = 'bioc',
                               res = 2.5, path = tempdir())

late_ssp585_raw <- cmip6_world(model = 'MPI-ESM1-2-HR', ssp = '585',
                               time = '2061-2080', var = 'bioc',
                               res = 2.5, path = tempdir())

# 2. prep_future_clim -masking
late_clim_245 <- prep_future_clim(late_ssp245_raw)
late_clim_585 <- prep_future_clim(late_ssp585_raw)

# 3. MaxEnt Fundamental Niche
late_suit_clim_245 <- predict(late_clim_245, aedes_climate_model, type = "cloglog", na.rm = TRUE)
late_suit_clim_585 <- predict(late_clim_585, aedes_climate_model, type = "cloglog", na.rm = TRUE)

# 4. Multiplicative Realized Niche 
late_hybrid_245 <- late_suit_clim_245 * (0.5 + 0.5 * human_footprint)
late_hybrid_585 <- late_suit_clim_585 * (0.5 + 0.5 * human_footprint)

late_hybrid_245 <- terra::mask(late_hybrid_245, turkey_boundary)
late_hybrid_585 <- terra::mask(late_hybrid_585, turkey_boundary)

# 5. Niche shift 2070 - today)
late_change_245 <- late_hybrid_245 - suit_hybrid
late_change_585 <- late_hybrid_585 - suit_hybrid
names(late_change_245) <- "NicheChange"
names(late_change_585) <- "NicheChange"

# ============================================================
# VISUALIZATION (FIGURE 8A, 8B, 9A, 9B)
# ============================================================

# SSP2-4.5 (2061-2080)
plot_future_absolute(late_hybrid_245, "Figure8A_LateCentury_Realized_SSP245_600dpiNEW.png")
plot_shift(late_change_245, "Figure8B_LateCentury_Shift_SSP245_600dpiNEW.png")

# SSP5-8.5 (2061-2080: EXTREME SCENARIO)
plot_future_absolute(late_hybrid_585, "Figure9A_LateCentury_Realized_SSP585_600dpiNEW.png")
plot_shift(late_change_585, "Figure9B_LateCentury_Shift_SSP585_600dpiNEW.png")

cat("Yüzyıl sonu ekstrem risk haritaları başarıyla oluşturuldu.\n")
# ==============================================================================
# GEOTIFF SAVING PROJECTIONS
# ==============================================================================
# 1. 2050 (2041-2060) 
writeRaster(suit_hybrid_245, "suit_hybrid_2050_ssp245.tif", overwrite = TRUE)
writeRaster(suit_hybrid_585, "suit_hybrid_2050_ssp585.tif", overwrite = TRUE)
writeRaster(change_245,      "niche_change_2050_ssp245.tif", overwrite = TRUE)
writeRaster(change_585,      "niche_change_2050_ssp585.tif", overwrite = TRUE)

# 2. 2070 (2061-2080) 
writeRaster(late_hybrid_245, "suit_hybrid_2070_ssp245.tif", overwrite = TRUE)
writeRaster(late_hybrid_585, "suit_hybrid_2070_ssp585.tif", overwrite = TRUE)
writeRaster(late_change_245, "niche_change_2070_ssp245.tif", overwrite = TRUE)
writeRaster(late_change_585, "niche_change_2070_ssp585.tif", overwrite = TRUE)

# 3. cmip6_world 
writeRaster(fut_clim_245,  "clim_processed_2050_ssp245.tif", overwrite = TRUE)
writeRaster(fut_clim_585,  "clim_processed_2050_ssp585.tif", overwrite = TRUE)
writeRaster(late_clim_245, "clim_processed_2070_ssp245.tif", overwrite = TRUE)
writeRaster(late_clim_585, "clim_processed_2070_ssp585.tif", overwrite = TRUE)

cat("ALL THE GeoTIFF LAYERS SAVED !\n")

# ============================================================
# SECTION 13 (OPTİMİZE): CONTINUOUS-TRANSITION GAP ANALYSIS
# ============================================================

# 1. Occupied 
turkey_pts_sf   <- st_transform(turkey_only_pts_clean, crs = 3857)
occupied_buffer <- st_buffer(turkey_pts_sf, dist = 20000)
occupied_poly   <- st_union(occupied_buffer)
occupied_poly   <- st_transform(occupied_poly, crs = 4326)

occupied_vect <- vect(occupied_poly)
occupied_mask <- rasterize(occupied_vect, suit_hybrid, field = 1, background = 0)
occupied_mask <- occupied_mask == 1

# 2. 
# threshold for lower risk: 0.18
# threshold for high risk : 0.35
risk_potential <- suit_hybrid >= 0.18  
risk_high      <- suit_hybrid >= 0.35  

# species is traveling : 15 km buffer (0.50 eşiği)
road_corridor  <- inv_road_norm >= 0.50

# 3. 4 step dynamic classification 
gap_raster <- rast(suit_hybrid)
gap_raster[] <- 0

# class:
# 1 = blue: occupied
# 2 = orange: secondary potential risk 
# 3 = red: surveillance (high Risk)
gap_raster[risk_potential & !occupied_mask] <- 2
gap_raster[risk_high & road_corridor & !occupied_mask] <- 3
gap_raster[risk_potential & occupied_mask] <- 1

gap_raster <- terra::mask(gap_raster, vect(turkey_boundary))
levels(gap_raster) <- data.frame(ID = 0:3, label = c("0", "1", "2", "3"))
# 4. lets paint it
p_gap_v2 <- ggplot() +
  geom_spatraster(data = gap_raster) +
  scale_fill_manual(
    values = c(
      "0" = "#F5F5F5",  
      "1" = "#2166AC",  
      "2" = "#FDAE61",  
      "3" = "#D7191C"   
    ),
    labels = c(
      "0" = "Climatically Unsuitable",
      "1" = "Documented Established Zone",
      "2" = "Secondary Risk / Dispersal Corridor",
      "3" = "High Surveillance Priority"
    ),
    na.translate = FALSE,
    na.value = "transparent",
    name = "Surveillance\nPriority"
  ) +
  geom_sf(data = turkey_provinces, fill = "transparent", color = "gray45", linewidth = 0.2) +
  coord_sf(xlim = c(26, 45), ylim = c(36, 42), expand = FALSE) +
  scale_x_continuous(labels = function(x) paste0(x, "°E")) +
  scale_y_continuous(labels = function(y) paste0(y, "°N")) +
  labs(title = NULL, subtitle = NULL, caption = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "gray85", linetype = "dashed", linewidth = 0.2),
    legend.position  = "right",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8),
    axis.text        = element_text(size = 8, color = "black")
  )

print(p_gap_v2)
ggsave("Figure10_Dispersal_Gap_Optimized_600dpiNEW.png", plot = p_gap_v2, width = 11, height = 6.5, dpi = 600)
# NOTE: new_scale_fill() requires the ggnewscale package
# install.packages("ggnewscale") and library(ggnewscale)
# if you want the two-color gap map. Otherwise simplify by
# plotting gap only without the high_risk overlay layer.


# ============================================================
# SECTION 14: QUANTITATIVE INVASION EXTENT & POPULATION EXPOSURE
# ============================================================

cat("calculating invasion regions and risk percentages...\n")

# 1. THRESHOLDS
# ecological threshold (0.28) and high risk threshold (0.50)
th_moderate <- 0.28
th_high     <- 0.50

# 2. list of materials gonna be investigated
layers_list <- list(
  "Current (Present)"       = suit_hybrid,
  "2050 SSP2-4.5 (MPI-ESM)" = suit_hybrid_245,
  "2050 SSP5-8.5 (MPI-ESM)" = suit_hybrid_585,
  "2070 SSP2-4.5 (MPI-ESM)" = late_hybrid_245,
  "2070 SSP5-8.5 (MPI-ESM)" = late_hybrid_585
)

# 3. terra::cellSize for real km2)
cell_area_km2 <- cellSize(suit_hybrid, unit = "km")
cell_area_km2 <- terra::mask(cell_area_km2, turkey_boundary)

total_area_tr <- as.numeric(global(cell_area_km2, "sum", na.rm = TRUE)[1, 1])

# 4. table calculation loop 
results_df <- data.frame(
  Scenario         = character(),
  Unsuitable_km2   = numeric(),
  ModerateRisk_km2 = numeric(),
  HighRisk_km2     = numeric(),
  TotalSuitable_km2= numeric(),
  TotalSuitable_Pct= numeric(),
  stringsAsFactors = FALSE
)

for (sc_name in names(layers_list)) {
  lyr <- layers_list[[sc_name]]
  
  # masking (?)
  m_unsuitable <- lyr < th_moderate
  m_moderate   <- lyr >= th_moderate & lyr < th_high
  m_high       <- lyr >= th_high
  m_total_suit <- lyr >= th_moderate
  
  # km2 
  area_unsuit <- as.numeric(global(cell_area_km2 * m_unsuitable, "sum", na.rm = TRUE)[1, 1])
  area_mod    <- as.numeric(global(cell_area_km2 * m_moderate, "sum", na.rm = TRUE)[1, 1])
  area_high   <- as.numeric(global(cell_area_km2 * m_high, "sum", na.rm = TRUE)[1, 1])
  area_suit   <- as.numeric(global(cell_area_km2 * m_total_suit, "sum", na.rm = TRUE)[1, 1])
  
  pct_suit    <- (area_suit / total_area_tr) * 100
  
  results_df <- rbind(results_df, data.frame(
    Scenario          = sc_name,
    Unsuitable_km2    = round(area_unsuit, 1),
    ModerateRisk_km2  = round(area_mod, 1),
    HighRisk_km2      = round(area_high, 1),
    TotalSuitable_km2 = round(area_suit, 1),
    TotalSuitable_Pct = round(pct_suit, 2)
  ))
}

cat("\n=================== INVASION RISK COVERAGE IN TURKEY ===================\n")
print(results_df)
write.csv(results_df, "Table_Invasion_Area_Summary.csv", row.names = FALSE)
cat("Özet tablo 'Table_Invasion_Area_Summary.csv' olarak kaydedildi.\n\n")

# ============================================================
# 5. FIGURE: AREA PERCENTAGE CHANGE
# ============================================================
library(tidyr)
library(dplyr)

# Long format
plot_data <- results_df %>%
  select(Scenario, ModerateRisk_km2, HighRisk_km2) %>%
  pivot_longer(cols = c("ModerateRisk_km2", "HighRisk_km2"),
               names_to = "RiskClass", values_to = "Area_km2") %>%
  mutate(
    Area_Pct = (Area_km2 / total_area_tr) * 100,
    RiskClass = factor(RiskClass, 
                       levels = c("HighRisk_km2", "ModerateRisk_km2"),
                       labels = c("High Risk (>0.50)", "Moderate / Secondary Risk (0.28 - 0.50)")),
    Scenario = factor(Scenario, levels = names(layers_list))
  )

p_area <- ggplot(plot_data, aes(x = Scenario, y = Area_Pct, fill = RiskClass)) +
  geom_bar(stat = "identity", width = 0.65, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("High Risk (>0.50)" = "#D7191C", 
                               "Moderate / Secondary Risk (0.28 - 0.50)" = "#FDAE61")) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, max(results_df$TotalSuitable_Pct) * 1.25),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = NULL,
    y = "Proportion of Turkey Land Area (%)",
    fill = "Risk Category"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.x = element_text(angle = 25, hjust = 1, face = "bold", size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.y = element_text(face = "bold", size = 10, margin = margin(r = 10)),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8.5)
  )

print(p_area)
ggsave("Figure11_Invasion_Area_Proportions_600dpiNEW.png", plot = p_area, width = 8.5, height = 5.5, dpi = 600)

cat("Grafik 'Figure11_Invasion_Area_Proportions_600dpi.png' olarak başarıyla kaydedildi.\n")
# ============================================================
# FINAL MODEL SUMMARY & RUN REPORT -THANK YOU FOR CHECKING :)
# ============================================================

cat("\n=================== FINAL MODEL SUMMARY ===================\n")
cat("Species            : Aedes albopictus (Asian tiger mosquito)\n")
cat("Study Area         : Turkey (Extensive Mediterranean/Regional calibration)\n")
cat("Occurrence Points  : Total Regional =", sum(p_vec_clim == 1), 
    "| Turkey Documented =", nrow(turkey_only_pts_clean), "\n")
cat("Background Points  : 10,000 pseudo-absences\n")
cat("Bioclimatic Covars : BIO11 (Coldest Quarter), BIO12 (Ann. Precip),\n")
cat("                     BIO14 (Driest Month), BIO15 (Precip Seasonality)\n")
cat("Human Footprint    : Population Density (log-transformed) +\n")
cat("                     OSM Major Road Network (15-km decay function)\n")
cat("Integration Method : Multiplicative Gatekeeper Model [Suit_Clim * (0.5 + 0.5 * HFP)]\n")
cat("Model Performance  : Training AUC =", sprintf("%.4f", train_auc), "(MaxEnt / maxnet)\n")
cat("GCM & Future Times : MPI-ESM1-2-HR | Mid-Century (2041-2060) & Late-Century (2061-2080)\n")
cat("Climate Scenarios  : SSP2-4.5 (Moderate) & SSP5-8.5 (High-emission/Extreme)\n")
cat("Generated Outputs  : \n")
cat("  - Figure 1: Raw Occurrence & Environmental Background\n")
cat("  - Figure 2: Human Footprint Infrastructure (Roads & Density)\n")
cat("  - Figure 3: Pure Climate Suitability (Fundamental Niche)\n")
cat("  - Figure 4: Current Realized Invasion Risk (Hybrid Model)\n")
cat("  - Figure ROC & Response: Model AUC (0.917) & Covariate Response Curves\n")
cat("  - Figure 6A/B: Mid-Century (2050) Projections & Niche Shift (SSP2-4.5)\n")
cat("  - Figure 7A/B: Mid-Century (2050) Projections & Niche Shift (SSP5-8.5)\n")
cat("  - Figure 8A/B & 9A/B: Late-Century (2070) Extreme Risk Scenarios\n")
cat("  - Figure 10: Realized Dispersal Gap & Surveillance Prioritization\n")
cat("===========================================================\n")