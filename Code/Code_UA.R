#Extraction des données

############## Préparation des données #######

library(readr)
library(sf)
library(deldir)
library(tidyverse)
library(dplyr)

generer_voronoi_climatique <- function(path_resultats_csv, path_zone_gpkg, coords = c("X", "Y"), epsg = 2154) {

  # Lire les données CSV
  resultats <- read_csv(path_resultats_csv) %>% drop_na()

  # Créer les points en objet sf
  points_sf <- st_as_sf(resultats, coords = coords, crs = epsg)

  # Charger la zone de découpage et la reprojeter
  zone_sf <- st_read(path_zone_gpkg, quiet = TRUE) %>%
    st_transform(crs = st_crs(points_sf))
  zone_union <- st_union(zone_sf)

  # Garder uniquement les points dans la zone
  points_sf <- points_sf[st_within(points_sf, zone_union, sparse = FALSE)[,1], ]

  # Fonction interne pour créer les polygones Voronoi
  voronoi <- function(points) {
    coords <- st_coordinates(points)
    del <- deldir(coords[,1], coords[,2])

    vor_polygons <- tile.list(del)
    polys <- lapply(vor_polygons, function(p) {
      coords <- cbind(p$x, p$y)
      coords <- rbind(coords, coords[1, ])  # Fermer le polygone
      st_polygon(list(coords))
    })

    st_sfc(polys, crs = st_crs(points))
  }

  # Générer les polygones de Voronoi
  vor_sf <- voronoi(points_sf)

  # Découper les polygones à la zone
  vor_sf <- st_intersection(vor_sf, zone_union)

  # Retourner les résultats sous forme sf
  resultats_voronoi <- st_sf(st_drop_geometry(points_sf), geometry = vor_sf)

  return(resultats_voronoi)
}


resultats_2010_2020_voronoi <- generer_voronoi_climatique(
  path_resultats_csv = "data/resultats_2010-2020.csv",
  path_zone_gpkg = "data/DEPT34_T15_sansmer.gpkg"
)

resultats_1990_2000_voronoi <- generer_voronoi_climatique(
  path_resultats_csv = "data/resultats_1990-2000.csv",
  path_zone_gpkg = "data/DEPT34_T15_sansmer.gpkg"
)

resultats_2000_2010_voronoi <- generer_voronoi_climatique(
  path_resultats_csv = "data/resultats_2000-2010.csv",
  path_zone_gpkg = "data/DEPT34_T15_sansmer.gpkg"
)

######################Classification clustgeo ##########
############## Choix des param ######

library(sf)
library(ClustGeo)
library(dplyr)

explorer_alpha_clustgeo <- function(voronoi_sf, K = 5,
                                    exclude_cols = c("CODE_STATION", "X", "Y"),
                                    alpha_range = seq(0, 1, 0.05)) {

  # Extraction et scaling des données climatiques
  clim_vars <- voronoi_sf %>%
    st_drop_geometry() %>%
    { .[, setdiff(names(.), exclude_cols)] } %>%
    select_if(is.numeric)

  D0 <- dist(scale(clim_vars))  # Distance climatique


  # Distance géographique (centroïdes)
  coords <- st_centroid(voronoi_sf$geometry)
  D1 <- st_distance(coords)
  D1 <- as.dist(as.matrix(D1))

  library(ggplot2)
  library(reshape2)
  library(viridis)

  station_ids <- as.character(1:nrow(voronoi_sf))

  # Conversion des matrices en format long
  D0_mat <- as.matrix(D0)
  rownames(D0_mat) <- colnames(D0_mat) <- station_ids
  D0_df <- melt(D0_mat)

  D1_mat <- as.matrix(D1)
  rownames(D1_mat) <- colnames(D1_mat) <- station_ids
  D1_df <- melt(D1_mat)

  # Facteurs ordonnés pour garantir un affichage correct
  D0_df$Var1 <- factor(D0_df$Var1, levels = station_ids)
  D0_df$Var2 <- factor(D0_df$Var2, levels = station_ids)
  D1_df$Var1 <- factor(D1_df$Var1, levels = station_ids)
  D1_df$Var2 <- factor(D1_df$Var2, levels = station_ids)

  # Plot D0 (dissimilarité climatique)
  p1 <- ggplot(D0_df, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma", direction = -1) +
    coord_equal() +
    labs(x = "Station", y = "Station", fill = "Dissimilarité") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, size = 6),
          axis.text.y = element_text(size = 6))

  # Plot D1 (dissimilarité géographique)
  p2 <- ggplot(D1_df, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma", direction = -1) +
    coord_equal() +
    labs(x = "Station", y = "Station", fill = "Distance (en m)") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, size = 6),
          axis.text.y = element_text(size = 6))

  print(p1)
  print(p2)

  # Analyse de l’effet d’alpha
  choix <- choicealpha(D0, D1, alpha_range, K, graph = TRUE)
  return(list(D0 = D0, D1 = D1, choix = choix))
}

resultats_alpha <- explorer_alpha_clustgeo(
  voronoi_sf = resultats_2010_2020_voronoi,
  K = 5,
  exclude_cols = c("CODE_STATION", "X", "Y"),
  alpha_range = seq(0, 1, 0.05)
)



resultats_alpha90 <- explorer_alpha_clustgeo(
  voronoi_sf = resultats_1990_2000_voronoi,
  K = 5,
  exclude_cols = c("CODE_STATION", "X", "Y"),
  alpha_range = seq(0, 1, 0.05)
)

resultats_alpha0010 <- explorer_alpha_clustgeo(
  voronoi_sf = resultats_2000_2010_voronoi,
  K = 5,
  exclude_cols = c("CODE_STATION", "X", "Y"),
  alpha_range = seq(0, 1, 0.05)
)


######## Classif param ajustés ######
library(ClustGeo)
library(ggplot2)
library(dplyr)
library(patchwork)
library(MetBrewer)

clustgeo_classification <- function(voronoi_sf, D0, D1, K = 5, alpha = 0.35) {


  # Clustering
  tree <- hclustgeo(D0, D1, alpha)
  clusters <- cutree(tree, K)
  voronoi_sf$clus <- as.factor(clusters)

  # Couleurs
  palette <- met.brewer("Archambault", K)

  # Carte des clusters
  p0 <- ggplot() +
    geom_sf(data = voronoi_sf, aes(fill = clus), color = "black") +
    scale_fill_manual(values = palette_prune, name = "Clusters") +
    theme_minimal()

  print(p0)

  return(voronoi_sf)
}


resultats_clustgeo <- clustgeo_classification(
  voronoi_sf = resultats_2010_2020_voronoi,
  D0 = resultats_alpha$D0,
  D1 = resultats_alpha$D1,
  K = 5,
  alpha = 0.15
)



resultats_clustgeo90 <- clustgeo_classification(
  voronoi_sf = resultats_1990_2000_voronoi,
  D0 = resultats_alpha90$D0,
  D1 = resultats_alpha90$D1,
  K = 5,
  alpha = 0.5
)

resultats_clustgeo0010 <- clustgeo_classification(
  voronoi_sf = resultats_2000_2010_voronoi,
  D0 = resultats_alpha0010$D0,
  D1 = resultats_alpha0010$D1,
  K = 5,
  alpha = 0.4
)

############ Transfo en shapefile + smoothing ##########

library(sf)
library(terra)
library(smoothr)
library(dplyr)
library(lwgeom)

smooth_zoning <- function(data,
                          x_col = "X",
                          y_col = "Y",
                          zone_col = "max_clus",
                          mask_path=NULL,
                          majority_window = 5,
                          simplify_tolerance = 100,
                          smooth_method = "ksmooth",
                          smoothness = 2,
                          max_hole_area = 1000) {

  # 1. Conversion en raster
  r <- rast(data[, c(x_col, y_col, zone_col)], type = "xyz")
  names(r) <- "zone"

  # 2. Filtre majoritaire
  if (majority_window > 1) {
    r_smooth <- focal(
      r,
      w = matrix(1, majority_window, majority_window),
      fun = function(x) {
        x <- x[!is.na(x)]
        if(length(x) == 0) return(NA)
        as.numeric(names(which.max(table(x))))
      }
    )
  } else {
    r_smooth <- r
  }

  # 3. Vectorisation initiale
  polys <- as.polygons(r_smooth, dissolve = TRUE) |>
    st_as_sf() |>
    rename(zone = 1) |>
    mutate(zone = as.factor(zone))

  # 4. Correction topologique initiale
  polys_repaired <- polys |>
    st_make_valid() |>
    group_by(zone) |>
    summarise() |>
    st_buffer(0)

  # 5. Simplification avant lissage
  polys_simplified <- st_simplify(
    polys_repaired,
    preserveTopology = TRUE,
    dTolerance = simplify_tolerance
  )

  # 6. Lissage
  if (!is.null(smooth_method)) {
    smoothed <- smooth(
      polys_simplified,
      method = smooth_method,
      smoothness = smoothness
    ) |>
      st_make_valid() |>
      st_buffer(0)

    # 7. Détection et comblement des trous aux intersections
    # Création d'un polygone englobant toutes les zones
    bbox <- st_as_sfc(st_bbox(smoothed)) |> st_buffer(10)

    # Trouver les espaces entre les zones (les trous)
    unioned <- st_union(smoothed)
    holes <- st_difference(bbox, unioned)

    # Filtrer seulement les petits trous (éviter les vrais espaces vides)
    if (length(holes) > 0) {
      holes <- holes[st_area(holes) < max_hole_area]

      # Si des trous sont trouvés, les combler avec la zone adjacente
      if (length(holes) > 0) {
        # Trouver la zone la plus proche pour chaque trou
        nearest <- st_nearest_feature(holes, smoothed)

        # Étendre les zones pour combler les trous
        smoothed_filled <- smoothed
        for (i in seq_along(holes)) {
          zone_to_expand <- nearest[i]
          smoothed_filled[zone_to_expand] <- st_union(
            smoothed_filled[zone_to_expand],
            holes[i]
          )
        }

        smoothed <- smoothed_filled |>
          st_make_valid() |>
          st_buffer(0)
      }
    }

    result <- smoothed
  } else {
    result <- polys_simplified
  }

  if(is.na(st_crs(result))) {
    result <- st_set_crs(result, 2154)  # RGF93/Lambert-93 pour la France
  }
  if (!is.null(mask_path)) {
    mask <- st_read(mask_path, quiet = TRUE) |>
      st_union() |>
      st_transform(2154)  # Assurez-vous que c'est le même CRS

    result <- st_intersection(result, mask) |>
      st_make_valid() |>
      st_collection_extract("POLYGON")
  }

  result <- result %>%
    mutate(zone = as.integer(as.character(zone)))

  return(result)
}


############## Random Forest ##############################

# foret sur chaque clus

resultats_rf_spatial <- function(df_clusters, grille, covariables, clusters, crs_code = "+init=epsg:2154") {
  library(randomForest)
  library(sp)
  library(rpart.plot) # pour le plot
  library(rpart)       # pour convertir l'arbre si nécessaire

  df_clusters$coord_X <- coordinates(df_clusters)[, 1]
  df_clusters$coord_Y <- coordinates(df_clusters)[, 2]
  grille$coord_X <- coordinates(grille)[, 1]
  grille$coord_Y <- coordinates(grille)[, 2]

  covariables_rf <- c(covariables, "coord_X", "coord_Y")
  resultats_final <- list()

  for (clus in clusters) {
    formule <- as.formula(paste(clus, "~", paste(covariables_rf, collapse = "+")))

    rf_model <- randomForest(formule, data = df_clusters, ntree = 1000)
    print(summary(rf_model))
    # Prédiction sur la grille
    pred_rf <- predict(rf_model, newdata = grille)

    # Stockage des résultats
    grille[[paste0("prediction_", clus)]] <- pred_rf
    resultats_final[[clus]] <- data.frame(grille@data, prediction = pred_rf)

    # ---- Nouveau : extraction d’un arbre et plot ----
    # Extraire un arbre (le premier ici)
    arbre <- getTree(rf_model, k = 1, labelVar = TRUE)

    # Convertir en objet rpart-like pour le plot (non direct depuis RF, donc on reconstruit un arbre pour illustrer)
    arbre_rpart <- rpart(formule, data = df_clusters, method = "class",
                         control = rpart.control(cp = 0.000001, minsplit = 2, maxdepth = 30))
    n_classes <- length(unique(df_clusters[[as.character(formule[[2]])]]))
    palette_magma <- list("#FEC287FF","#B63679FF")
    rpart.plot(arbre_rpart,
               # main = "Arbre de décision du cluster",
               box.palette = palette_magma,
               type = 5,
               branch = 1,
               uniform = TRUE,
               gap = 4,               # Augmente l'espace entre les nœuds
               leaf.round = 2,        # Boîtes avec coins plus arrondis (plus visibles)
               shadow.col = "gray",
               branch.lty = 5,
               cex = 1.2,             # Agrandit la taille des textes dans les boîtes
               fallen.leaves = TRUE,  # Positionne mieux les feuilles
               tweak = 1.2            # Ajuste la taille globale des nœuds (boîtes)
    )
  }

  # Ajouter aussi les coordonnées
  for (var in names(resultats_final)) {
    resultats_final[[var]]$X <- coordinates(grille)[, 1]
    resultats_final[[var]]$Y <- coordinates(grille)[, 2]
  }

  return(resultats_final)
}

resultats_rf_data <- resultats_rf_spatial(
  df_clusters = krige_data$df_clusters,  # df_clusters préparé
  grille = krige_data$grille,            # grille pour la prédiction
  covariables = c("mnt", "slope", "aspect", "distmer_lowres", "mrvbf"),  # Covariables utilisées pour la régression
  clusters = c("clus1", "clus2", "clus3", "clus4", "clus5")  # Liste des clusters
)

dfclusdens <- krige_data$df_clusters
dfgrilledens <- krige_data$grille

dfclusdens$source <- "stations"
dfgrilledens$source <- "grille"

# Combiner les deux jeux de données
df_combineddens <- rbind(
  dfclusdens[, c("mnt", "source")],
  dfgrilledens[, c("mnt", "source")]
)

# Palette actualisée
palette_prune1 <- c(
  "stations" = "#FEC287FF",  # Bleu-gris doux
  "grille"   = "#B63679FF"   # Vert sauge
)

# S'assurer que les niveaux sont bien définis
df_combineddens$source <- factor(df_combineddens$source,
                                 levels = c("stations", "grille"))

df_combineddens <- as.data.frame(df_combineddens)
# Graphe
ggplot(df_combineddens, aes(x = source, y = mnt, fill = source)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.05, outlier.shape = NA, alpha = 0.8) +
  scale_fill_manual(values = palette_prune1) +
  labs(title = " ",
       x = "Distribution",
       y = "MNT (Altitude)") +
  theme_minimal() +
  theme(legend.position = "none")



scale_fill_manual(values = palette_prune)



head( krige_data$df_clusters)


plot_rf <- function(resultats_rf) {
  library(ggplot2)
  library(dplyr)
  library(patchwork)

  plots <- list()

  for (var in names(resultats_rf)) {
    p <- ggplot(resultats_rf[[var]], aes(x = X, y = Y, fill = prediction)) +
      geom_tile() +
      scale_fill_viridis_c(option = "magma", guide = guide_colorbar(barwidth = 0.5, barheight = 3)) +
      coord_fixed() +
      theme_minimal() +
      labs(fill = var) +
      theme(
        legend.text = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8, angle = 45),
        plot.title = element_text(size = 10, hjust = 0.5)
      )

    plots[[var]] <- p
    print(p)
  }

  a <- wrap_plots(plots, ncol = 3)

  # Création d'un tableau combiné avec les prédictions
  res_pred <- resultats_rf$clus1 %>%
    select(X, Y, prediction) %>%
    rename(clus1 = prediction)

  res_pred$clus2 <- resultats_rf$clus2$prediction
  res_pred$clus3 <- resultats_rf$clus3$prediction
  res_pred$clus4 <- resultats_rf$clus4$prediction
  res_pred$clus5 <- resultats_rf$clus5$prediction

  # Calcul du cluster avec la prédiction maximale
  res_pred <- res_pred %>%
    rowwise() %>%
    mutate(
      max_pred = max(c(clus1, clus2, clus3, clus4, clus5)),
      max_clus = which.max(c(clus1, clus2, clus3, clus4, clus5))
    )

  # Carte des clusters prédits
  b <- ggplot(res_pred, aes(x = X, y = Y, fill = as.factor(max_clus))) +
    geom_tile() +
    scale_fill_manual(values = palette_prune, name = "Cluster") +
    coord_fixed() +
    theme_minimal() +
    theme(legend.position = "right")

  return(list(patch = a, rf_clusters = b))
}

plot_rf_data <- plot_rf(resultats_rf_data)

plot_rf_data$patch[1]

plot_rf_data$rf_clusters
plot_rf_data$patch + plot_rf_data$rf_clusters +
  plot_annotation(title = "Clusters prédits avec Random Forest",
                  theme = theme(plot.title = element_text(size = 12, hjust = 0.5)))

polys_rf_2010_2020 <- smooth_zoning(
  data = plot_rf_data$rf_clusters$data,
  majority_window = 5,
  mask_path = "data/DEPT34_T15_sansmer.gpkg",
  smooth_method = "ksmooth",
  max_hole_area = 0  # Seuil en mètres carrés
)

plot(polys_rf_2010_2020)


resultats_rf_data_90_00 <- resultats_rf_spatial(
  df_clusters = krige_data90$df_clusters,  # df_clusters préparé pour 1990-2000
  grille = krige_data90$grille,            # grille pour la prédiction
  covariables = c("mnt", "slope", "aspect", "distmer_lowres", "mrvbf"),  # Covariables utilisées pour la régression
  clusters = c("clus1", "clus2", "clus3", "clus4", "clus5")  # Liste des clusters
)

names(resultats_rf_data_90_00)

plot_rf_data90 <- plot_rf(resultats_rf_data_90_00)

names(plot_rf_data90$rf_clusters$data) %>%
  mutate(max_clus = case_when(
    max_clus == 4 ~ 5,
    max_clus == 5 ~ 4,
    TRUE ~ max_clus
  ))

polys_rf_1990_2000 <- smooth_zoning(
  data = plot_rf_data90$rf_clusters$data,
  majority_window = 5,
  mask_path = "data/DEPT34_T15_sansmer.gpkg",
  smooth_method = "ksmooth",
  max_hole_area = 0  # Seuil en mètres carrés
)

plot(polys_rf_1990_2000)

resultats_rf_data_00_10 <- resultats_rf_spatial(
  df_clusters = krige_data0010$df_clusters,  # df_clusters préparé pour 2000-2010
  grille = krige_data0010$grille,            # grille pour la prédiction
  covariables = c("mnt", "slope", "aspect", "distmer_lowres", "mrvbf"),  # Covariables utilisées pour la régression
  clusters = c("clus1", "clus2", "clus3", "clus4", "clus5")  # Liste des clusters
)

plot_rf_data00_10 <- plot_rf(resultats_rf_data_00_10)

polys_rf_2000_2010 <- smooth_zoning(
  data = plot_rf_data00_10$rf_clusters$data,
  majority_window = 5,
  mask_path = "data/DEPT34_T15_sansmer.gpkg",
  smooth_method = "ksmooth",
  max_hole_area = 0  # Seuil en mètres carrés
)

plot(polys_rf_2000_2010)

polys_rf_1990_2000 <- polys_rf_1990_2000 %>%
  mutate(zone = case_when(
    zone == 4 ~ 5,
    zone == 5 ~ 4,
    TRUE ~ zone
  ))

polys_rf_2000_2010 <- polys_rf_2000_2010 %>%
  mutate(zone = case_when(
    zone == 1 ~ 2,
    zone == 2 ~ 5,
    zone == 5 ~ 1,
    TRUE ~ zone
  ))

library(patchwork)
p1 <- ggplot(polys_rf_1990_2000) +
  geom_sf(aes(fill = factor(zone))) +
  scale_fill_manual(values = palette_prune) +
  labs( fill = "Cluster") +
  theme_minimal()

p2 <- ggplot(polys_rf_2000_2010) +
  geom_sf(aes(fill = factor(zone))) +
  scale_fill_manual(values = palette_prune) +
  labs( fill = "Cluster") +
  theme_minimal()

p3 <- ggplot(polys_rf_2010_2020) +
  geom_sf(aes(fill = factor(zone))) +
  scale_fill_manual(values = palette_prune) +
  labs( fill = "Cluster") +
  theme_minimal()

p1 + p2 + p3 +
  plot_annotation(title = "Clusters prédits avec Random Forest",
                  theme = theme(plot.title = element_text(size = 12, hjust = 0.5)))





########### Prep courbes ##########

common_codes <- Reduce(intersect, list(
  krige_data90$df_clusters$CODE_STATION,
  krige_data0010$df_clusters$CODE_STATION,
  krige_data$df_clusters$CODE_STATION
))

# Afficher les codes communs
print(common_codes)

b90 <- read_csv("data/big_1990-2000.csv")
b0010 <- read_csv("data/big_2000-2010.csv")
b1020 <- read_csv("data/big_2010-2020.csv")

b_all <- bind_rows(b90, b0010, b1020)

b_all <- b_all %>%
  filter(CODE_STATION %in% common_codes) %>%
  mutate(periode = case_when(
    CODE_STATION %in% b90$CODE_STATION ~ "1990-2000",
    CODE_STATION %in% b0010$CODE_STATION ~ "2000-2010",
    CODE_STATION %in% b1020$CODE_STATION ~ "2010-2020"
  ))

courbes_annuelles <- b_all %>%
  group_by(CODE_STATION, ANNEE) %>%
  summarise(
    temp_moy_annuelle = mean((TEMP_MAX_MOY + TEMP_MIN_MOY)/2 , na.rm = TRUE),
    cumul_pluv_annuel = sum(CUMUL_PRECIP , na.rm = TRUE),
    periode = first(periode),
    .groups = "drop"
  )

library(ggplot2)
library(dplyr)

# Supposons que ton dataframe s'appelle courbes_annuelles

# Liste des codes stations uniques
stations <- unique(courbes_annuelles$CODE_STATION)

# Créer une liste de plots, un plot par station
plots_list <- lapply(stations, function(station) {
  data_station <- filter(courbes_annuelles, CODE_STATION == station)

  # Calcul de la regression linéaire
  model <- lm(temp_moy_annuelle ~ ANNEE, data = data_station)
  model_info <- broom::tidy(model)

  # Extraire la pente (coefficient ANNEE)
  slope <- model_info$estimate[model_info$term == "ANNEE"]

  # Préparer le texte à afficher : évolution par an (ex: +0.03 °C / an)
  slope_text <- paste0("Évolution: ", round(slope, 3), " °C par an")

  ggplot(data_station, aes(x = ANNEE, y = temp_moy_annuelle)) +
    geom_line(color = "black", linewidth = 1) +
    geom_point(color = "black", size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "#c64646", linewidth = 1.5) +
    theme_minimal(base_size = 14) +
    labs(x = "Année", y = "Température moyenne annuelle (°C)") +
    annotate("text", x = min(data_station$ANNEE) + 1, y = max(data_station$temp_moy_annuelle),
             label = slope_text, color = "black", hjust = 0, size = 5)
})

names(plots_list) <- stations

plots_list[[1]]

plots_list2 <- lapply(stations, function(station) {
  data_station <- filter(courbes_annuelles, CODE_STATION == station)

  # Calcul de la regression linéaire
  model <- lm(cumul_pluv_annuel ~ ANNEE, data = data_station)
  model_info <- broom::tidy(model)

  # Extraire la pente (coefficient ANNEE)
  slope <- model_info$estimate[model_info$term == "ANNEE"]
  slope_text <- paste0("Évolution: ", round(slope, 3), " mm par an")

  ggplot(data_station, aes(x = ANNEE, y = cumul_pluv_annuel)) +
    geom_line(color = "black", linewidth = 1) +
    geom_point(color = "black", size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "#46a5c6", linewidth = 1.5) +
    theme_minimal(base_size = 14) +
    labs(x = "Année", y = "Cumul pluviomètrique annuel (mm)") +
    annotate("text", x = min(data_station$ANNEE) + 1, y = max(data_station$temp_moy_annuelle),
             label = slope_text, color = "black", hjust = 0, size = 5)
})

names(plots_list2) <- stations

plots_list2[[1]]

plots_combined <- lapply(names(plots_list), function(station) {
  list(
    temp = plots_list[[station]],
    pluie = plots_list2[[station]]
  )
})

# Nommer la liste principale par station
names(plots_combined) <- names(plots_list)

library(patchwork)

plots_combined_patchwork <- lapply(plots_combined, function(plots) {
  plots$temp + plots$pluie + plot_layout(ncol = 1)
})

# Afficher le premier plot combiné
plots_combined_patchwork[[1]]

names(plots_combined_patchwork)



################ Pour Leaflet #######
library(sf)
library(leaflet)
library(leaflet.extras)
library(dplyr)
library(colorspace)
library(htmltools)
library(leafpop)
library(ggplot2)
library(patchwork)

carte_leaflet <- function(map,polys,data_clust, periode = "période",legende=TRUE,plots_combined_patchwork=NULL){
  has_patchwork <- !is.null(plots_combined_patchwork)

  # 1. Reprojection
  zones <- st_transform(polys, 4326)
  stations <- st_transform(st_as_sf(data_clust), 4326)

  names(stations)
  # 2. Ajouter noms de zones (selon les valeurs de 1 à 5) et forcer l'ordre
  zones <- zones %>%
    mutate(nom_zone = case_when(
      zone == 1 ~ "Centre",
      zone == 2 ~ "Sud",
      zone == 3 ~ "Nord-Ouest",
      zone == 4 ~ "Nord-Est",
      zone == 5 ~ "Sud-Est"
    ),
    nom_zone = factor(nom_zone, levels = c("Centre", "Sud", "Nord-Ouest", "Nord-Est", "Sud-Est"))
    )

  # 3. Jointure spatiale : associer chaque station à une zone (avec noms)
  stations <- st_join(stations, zones, join = st_within)

  # 4. Reforcer l'ordre dans stations aussi
  stations$nom_zone <- factor(stations$nom_zone, levels = c("Centre", "Sud", "Nord-Ouest", "Nord-Est", "Sud-Est"))

  # 5. Calcul température moyenne stations par zone
  temp_moy_zone <- stations %>%
    st_set_geometry(NULL) %>%
    group_by(zone, nom_zone) %>%
    summarise(temp_moyenne_zone = mean(temp_moy_annuelle, na.rm = TRUE), .groups = "drop")

  # 6. Ajouter cette info dans zones
  zones <- zones %>%
    left_join(temp_moy_zone, by = c("zone", "nom_zone"))

  # 7. Définir la palette personnalisée
  couleurs_perso <- c(
    "Centre" = "#FEC287",
    "Sud" = "#FB8861",
    "Nord-Ouest" = "#E65164",
    "Nord-Est" = "#B63679",
    "Sud-Est" = "#822681"
  )

  ordre_zones <- c("Centre", "Sud", "Nord-Ouest", "Nord-Est", "Sud-Est")

  pal <- colorFactor(
    palette = couleurs_perso,
    domain = names(couleurs_perso),
    ordered = TRUE
  )

  # 8. Palette sombre pour les cercles
  darker_pal <- function(x) {
    colors <- pal(x)
    colorspace::darken(colors, amount = 0.2)
  }

  # 5. Calcul des moyennes par zone
  vars_a_moyenner <- c(
    "temp_moy_annuelle", "cumul_pluv_annuel", "precip_hiver", "precip_ete",
    "temp_min_hiver", "temp_max_hiver", "temp_min_ete", "temp_max_ete",
    "frequence_mois_secs", "frequence_mois_tres_secs", "mnt", "slope",
    "aspect", "distmer_lowres", "mrvbf"
  )

  temp_moy_zone <- stations %>%
    st_set_geometry(NULL) %>%
    group_by(zone, nom_zone) %>%
    summarise(across(all_of(vars_a_moyenner), ~ round(mean(.x, na.rm = TRUE), 2)))

  zones <- zones %>%
    left_join(temp_moy_zone, by = c("zone", "nom_zone"))

  popup_station <-  function(station_row) {
    paste0(
      "<div style='text-align:center; font-size:18px; font-weight:bold; margin-bottom:2px;'>Station ", station_row$CODE_STATION,"</div>",
      "<div style='text-align:center; font-size:12px; margin-bottom:4px;'>Zone : ", station_row$nom_zone, "</div>",
      "<div style='text-align:center; font-size:14px; font-weight:bold; margin-bottom:5px;'> Variables Climatiques </div>",
      "<b>Température moyenne annuelle:</b> ", round(station_row$temp_moy_annuelle, 2), " °C<br>",
      "<b>Pluviométrie annuelle:</b> ", round(station_row$cumul_pluv_annuel, 2), " mm<br>",
      "<b>Précipitations hiver:</b> ", round(station_row$precip_hiver, 2), " mm<br>",
      "<b>Précipitations été:</b> ", round(station_row$precip_ete, 2), " mm<br>",
      "<b>T° min hiver:</b> ", round(station_row$temp_min_hiver, 2), " °C<br>",
      "<b>T° max hiver:</b> ", round(station_row$temp_max_hiver, 2), " °C<br>",
      "<b>T° min été:</b> ", round(station_row$temp_min_ete, 2), " °C<br>",
      "<b>T° max été:</b> ", round(station_row$temp_max_ete, 2), " °C<br>",
      "<b>Fréquence mois secs:</b> ", round(station_row$frequence_mois_secs, 2), "<br>",
      "<b>Fréquence mois très secs:</b> ", round(station_row$frequence_mois_tres_secs, 2), "<br>",
      "<div style='text-align:center; font-size:14px; font-weight:bold; margin-bottom:5px;'> Variables Géomorphologiques </div>",
      "<b>MNT:</b> ", round(station_row$mnt, 2), " m<br>",
      "<b>Pente:</b> ", round(station_row$slope, 2), "°<br>",
      "<b>Orientation:</b> ", round(station_row$aspect, 2), "°<br>",
      "<b>Distance mer:</b> ", round(station_row$distmer_lowres, 2), " m<br>",
      "<b>MRVBF:</b> ", round(station_row$mrvbf, 2))
  }

  # 6. Popups des stations
  stations <- stations %>%
    rowwise()%>%
    mutate(popup = {
      #debug
      #print(cur_data())
      station_i <- as.character(CODE_STATION)
      if (has_patchwork && station_i %in% names(plots_combined_patchwork)) {
        p1 <- plots_combined_patchwork[[station_i]]
        img_file <- tempfile(fileext = ".png")
        ggsave(img_file, plot = p1, width = 6, height = 7)
        # Encodage base64
        img_html <- base64enc::dataURI(file = img_file, mime = "image/png")
        paste0(
          popup_station(across(everything())),
          "<br><img src='", img_html, "' width='300px'>"
        )
      }
      else {
        popup_station(across(everything()))
      }
    }) %>%
    ungroup()

  # 7. Popups des zones
  popup_zone <- function(zone_row) {
    paste0(
      "<div style='text-align:center; font-size:18px; font-weight:bold; margin-bottom:5px;'>Zone : ", zone_row$nom_zone, "</div>",
      "<div style='text-align:center; font-size:14px; font-weight:bold; margin-bottom:5px;'> Variables Climatiques </div>",
      "<b>T° moyenne annuelle : </b>", zone_row$temp_moy_annuelle, " °C<br>",
      "<b>Pluviométrie annuelle : </b>", zone_row$cumul_pluv_annuel, " mm<br>",
      "<b>Précipitations hiver : </b>", zone_row$precip_hiver, " mm<br>",
      "<b>Précipitations été : </b>", zone_row$precip_ete, " mm<br>",
      "<b>T° min hiver : </b>", zone_row$temp_min_hiver, " °C<br>",
      "<b>T° max hiver : </b>", zone_row$temp_max_hiver, " °C<br>",
      "<b>T° min été : </b>", zone_row$temp_min_ete, " °C<br>",
      "<b>T° max été : </b>", zone_row$temp_max_ete, " °C<br>",
      "<b>Fréquence mois secs : </b>", zone_row$frequence_mois_secs, "<br>",
      "<b>Fréquence mois très secs : </b>", zone_row$frequence_mois_tres_secs, "<br>",
      "<div style='text-align:center; font-size:14px; font-weight:bold; margin-bottom:5px;'> Variables Géomorphologiques </div>",
      "<b>MNT : </b>", zone_row$mnt, " m<br>",
      "<b>Pente : </b>", zone_row$slope, "°<br>",
      "<b>Orientation : </b>", zone_row$aspect, "°<br>",
      "<b>Distance mer : </b>", zone_row$distmer_lowres, " m<br>",
      "<b>MRVBF : </b>", zone_row$mrvbf, "<br>"
    )
  }

  zones <- zones %>%
    rowwise() %>%
    mutate(popup_complete = {
      # Générer les graphiques pour la zone
      zone_i <- nom_zone
      df_i <- stations %>%
        filter(nom_zone == zone_i) %>%
        st_set_geometry(NULL)

      p1 <- ggplot(df_i, aes(x = "", y = temp_moy_annuelle)) +
        geom_violin(fill = couleurs_perso[zone_i], alpha = 0.6, color = couleurs_perso[zone_i]) +
        geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
        labs(title = "Température moyenne annuelle", y = "Température (°C)", x = NULL) +
        theme_minimal() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())

      p2 <- ggplot(df_i, aes(x = "", y = cumul_pluv_annuel)) +
        geom_violin(fill = couleurs_perso[zone_i], alpha = 0.6, color = couleurs_perso[zone_i]) +
        geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
        labs(title = "Cumul pluviométrique", y = "Pluie (mm)", x = NULL) +
        theme_minimal() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())

      img_file <- tempfile(fileext = ".png")
      ggsave(img_file, plot = patchwork::wrap_plots(p1, p2, ncol = 2), width = 6, height = 4)

      # Encodage base64
      img_base64 <- base64enc::dataURI(file = img_file, mime = "image/png")

      # HTML combiné
      paste0(
        popup_zone(across(everything())),
        "<br><img src='", img_base64, "' width='300px'>"
      )
    }) %>%
    ungroup()

  stations <- stations %>%
    mutate(rayon = if (!is.null(plots_combined_patchwork)) {
      ifelse(
        as.character(CODE_STATION) %in% names(plots_combined_patchwork),
        5.5,  # rayon plus grand si plot disponible
        4   # sinon rayon standard
      )
    } else {
      4
    },
    couleuran= if (!is.null(plots_combined_patchwork)) {
      ifelse(
        as.character(CODE_STATION) %in% names(plots_combined_patchwork),
        "red",
        "black"
      )
    } else {
      "black"
    })
  addLegendCustom <- function(map, colors, labels, sizes, group = NULL) {
    legend_content <- mapply(function(color, label, size) {
      paste0(
        "<i style='
        width:", size, "px;
        height:", size, "px;
        border: 2px solid ", color, ";
        border-radius: 50%;
        display: inline-block;
        margin-right: 8px;
        vertical-align: middle;
        background: transparent;
      '></i><span style='vertical-align: middle;'>", label, "</span>"
      )
    }, colors, labels, sizes, SIMPLIFY = FALSE) %>%
      paste(collapse = "<br>")

    addControl(map, html = paste0(
      "<div class='legend leaflet-control' style= font-size: 13px;'>",
      "<strong>Stations</strong><br>", legend_content, "</div>"
    ),
    position = "bottomright", layerId = group)
  }
  # 8. Création carte
  map <- map %>%
    addPolygons(data = zones,
                color = "black",
                weight = 1,
                fillColor = ~pal(nom_zone),
                fillOpacity = 0.5,
                popup = ~popup_complete,
                label = ~paste("Zone :", nom_zone),
                highlightOptions = highlightOptions(
                  weight = 1,
                  color = "lightyellow",
                  fillOpacity = 0.7,
                  bringToFront = FALSE
                ),
                group=periode) %>%

    addCircleMarkers(data = stations,
                     radius = ~rayon,
                     color = ~couleuran,
                     weight = 2,
                     fillColor = ~darker_pal(nom_zone),
                     fillOpacity = 1,
                     popup = ~popup,
                     label = ~paste("Station", CODE_STATION),
                     group=periode)

  if (legende) {
    map <- map %>%
      addLegendCustom(
        colors = c("red", "black"),
        labels = c("Stations avec historique", "Autres stations"),
        sizes = c(11, 8),  # tailles visuellement cohérentes avec rayon 5.5 et 4
        group = periode
      ) %>%
      addLegend("bottomright",
                colors = couleurs_perso,
                labels = ordre_zones,
                title = "Nom des zones",
                group = periode)
  }

  return(map)
}

map_base <- leaflet() %>% addProviderTiles("CartoDB.Positron")

carte_leaflet_10_20 <- carte_leaflet(map_base,polys_rf_2010_2020, krige_data$df_clusters,periode= "2010-2020",legend=TRUE, plots_combined_patchwork = plots_combined_patchwork)

carte_leaflet_00_10 <- carte_leaflet(leaflet() %>%addProviderTiles("CartoDB.Positron"),polys_rf_2000_2010, krige_data0010$df_clusters, periode = "2000-2010",legend=FALSE)

carte_leaflet_90_00 <- carte_leaflet(leaflet() %>%addProviderTiles("CartoDB.Positron"),polys_rf_1990_2000, krige_data90$df_clusters, periode = "1990-2000")


carte_combin <- leaflet() %>%
  addProviderTiles("CartoDB.Positron") %>%
  carte_leaflet(polys_rf_2010_2020, krige_data$df_clusters, periode = "2010-2020", plots_combined_patchwork = plots_combined_patchwork) %>%
  carte_leaflet(polys_rf_2000_2010, krige_data0010$df_clusters, periode = "2000-2010",legend=FALSE, plots_combined_patchwork = plots_combined_patchwork) %>%
  carte_leaflet(polys_rf_1990_2000, krige_data90$df_clusters, periode = "1990-2000",legend=FALSE, plots_combined_patchwork = plots_combined_patchwork) %>%
  addLayersControl(
    baseGroups = c("2010-2020", "2000-2010", "1990-2000"),
    options = layersControlOptions(collapsed = FALSE)
  )

carte_combin
