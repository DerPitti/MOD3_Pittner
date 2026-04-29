library(readxl)
library(dplyr)
library(ggplot2)
library(future.apply)
library(tidyr)
library(dplyr)
library(vegan)
library(dbscan)
library(future)
library(future.apply)
library(clue)
library(mclust)
library(plotly)
library(geometry)

plan(multisession, workers = parallel::detectCores() - 5)

source("clusteranalysis_functions.R")

grunddaten <- read_xlsx("data/tbl_grunddaten.xlsx") # load basic data for each plot, e.g. biotopcodes
plants <- read_xlsx("data/tbl_daten_pflanzen.xlsx") # load plant data
life_forms <- read_xlsx("data/Life_form.xlsx") # load life form data https://floraveg.eu/download/

# group biotope codes
grunddaten$BT_Bund_group <- substr(grunddaten$`Biotoptyp-Bund`,1,5)
grunddaten$BT_Land_group <- substr(grunddaten$`Biotoptyp-Land`,1,2)


# Removal of data sets with non-meaningful plant composition --------------

remove_polygons <- filter(grunddaten, substr(`Biotoptyp-Land`,1,1) %in% c("F", "V", "W")) # maybe G and H as well
# check that every plot has plants
remove_polygons2 <- filter(grunddaten, !Polygon %in% plants$Polygon & !substr(`Biotoptyp-Land`,1,1) %in% c("F", "V", "W"))

remove_poly <- rbind(remove_polygons, remove_polygons2)
plants_sub <- filter(plants, !Polygon %in% remove_poly$Polygon)
grunddaten_sub <- filter(grunddaten, !Polygon %in% remove_poly$Polygon)

# transformation of plant data frame

# remove duplicate rows
plants_clean <- plants_sub %>%
  dplyr::distinct(Polygon, `Wissenschaftlicher Name`, Menge, .keep_all = TRUE)

# only keep row with highest abundance
plants_clean2 <- plants_clean %>%
  dplyr::arrange(Polygon, `Wissenschaftlicher Name`, desc(Menge)) %>%
  dplyr::distinct(Polygon, `Wissenschaftlicher Name`, .keep_all = TRUE
  )

plants_clean2$Menge <- as.numeric(plants_clean2$Menge) # transform abundances to numeric


### find outlier plots

# 01728 will be removed, because only one plant found, and in 04567, Ulmus spec. will be specified to Ulmus glabra accordingly to Beschreibung
plants_clean2 <- plants_clean2 %>%
  filter(Polygon != "01728") %>%
  mutate(`Wissenschaftlicher Name` = if_else(`Wissenschaftlicher Name`=="Ulmus spec.", "Ulmus glabra", `Wissenschaftlicher Name`))

plants_wide <- plant_widening(plants_clean2)

grunddaten_sub <- grunddaten_sub %>% # remove 01728 also from grunddaten
  filter(Polygon != "01728")

# Separation of biotope types ---------------------------------------------

grunddaten_grass <- filter(grunddaten_sub, substr(`Biotoptyp-Land`,1,1) %in% c("E"))
plants_grass <- filter(plants_wide, Polygon %in% grunddaten_grass$Polygon)

grunddaten_forests <- filter(grunddaten_sub, substr(`Biotoptyp-Land`,1,1) %in% c("A"))
plants_forests <- filter(plants_wide, Polygon %in% grunddaten_forests$Polygon)

# further reduction of forest data ---------------------------------------------

# remove "AT... and "AU..."
check_biotope_codes <- unique(grunddaten_forests[,c(4,5)])
grunddaten_forests_red <- filter(grunddaten_forests, !substr(`Biotoptyp-Land`,1,2) %in% c("AT", "AU", "AV"))
plants_forests_red <- filter(plants_forests, Polygon %in% grunddaten_forests_red$Polygon)

### store plants and grunddaten in list for looping later, each list entry contains two list entries with first plants, then grunddaten

data_list <- list()
data_list[["forest_complete"]] <-  list(plants = plants_forests_red, grunddaten = grunddaten_forests_red)

# remove all polygons which only have a few plant entries

forest_complete_wone <- remove_plots(plants_forests_red, grunddaten_forests_red, remove_count = c(0,1))
check_empty_plot(forest_complete_wone[[1]])

data_list[["forest_complete_wone"]] <- forest_complete_wone

forest_complete_wtwo <- remove_plots(plants_forests_red, grunddaten_forests_red, remove_count = c(0:2))
data_list[["forest_complete_wtwo"]] <- forest_complete_wtwo

forest_complete_wthree <- remove_plots(plants_forests_red, grunddaten_forests_red, remove_count = c(0:2))
data_list[["forest_complete_wthree"]] <- forest_complete_wthree

# add life forms
plants_occ <- plant_occurences(plant_data = plants_forests_red)

plants_occ$short <- sub("(\\w+\\s+\\w+).*", "\\1", plants_occ$species)
plants_occ$short <- sub("(\\w+).*", "\\1", plants_occ$species)
life_forms$short <- sub("(\\w+).*", "\\1", life_forms$FloraVeg.Taxon)

plants_LF <- left_join(plants_occ,life_forms, by = "short",multiple = "any")

trees <- filter(plants_LF, Tree == 1)
trees <- filter(trees, species %in% colnames(plants_forests))

plants_forest_trees <- plants_forests_red[c("Polygon",trees$species)]

plants_occ_forests <- plant_occurences(plants_forest_trees)

# remove plants which only occured once or twice across all plots
plants_f_trees_w_zero <- filter(plants_occ_forests, !total %in% c(0)) # maybe also 2
plants_f_trees_wzero <- plants_forests_red[c("Polygon",plants_f_trees_w_zero$species)] # take all 

plants_f_trees_w_one <- filter(plants_occ_forests, !total %in% c(0,1)) # maybe also 2
plants_f_trees_wone <- plants_forests_red[c("Polygon",plants_f_trees_w_one$species)] # take all 
# I'm using the plants_f_trees_wone for further reduction of the data


# check number of plant species per plot -------------------------------------------------

check_empty_plot(plants_f_trees_wone)
plants_trees_wzero <- remove_plots(plants_f_trees_wone, grunddaten_forests_red)
check_empty_plot(plants_trees_wzero[[1]])

data_list[["plants_trees_wzero"]] <- plants_trees_wzero

plants_trees_wone <- remove_plots(plants_f_trees_wone, grunddaten_forests_red, remove_count = c(0,1))
check_empty_plot(plants_trees_wone[[1]])
data_list[["plants_trees_wone"]] <- plants_trees_wone

plants_trees_wtwo <- remove_plots(plants_f_trees_wone, grunddaten_forests_red, remove_count = c(0:2))
check_empty_plot(plants_trees_wtwo[[1]])
data_list[["plants_trees_wtwo"]] <- plants_trees_wtwo

# plants_trees_wthree <- remove_plots(plants_f_trees_wone, grunddaten_forests_red, remove_count = c(0,1))
# check_empty_plot(plants_trees_wthree[[1]])
# data_list[["plants_trees_wthree"]] <- plants_trees_wthree


# NMDS loop ---------------------------------------------------------------
weighting <- data.frame(imbalanced = c(0.01,0.01,0.01,1),
                        balanced = c(0.01,0.05,0.25,0.75))
#stress_vals_backup <- stress_vals
stress_vals <- list(imbalanced = list(), balanced = list())

for(i in 1:2){
  for (o in 1:length(data_list)) {
    plants <- data_list[[o]][[1]]
    plant_mat <- plant_weighting(plants,weighting[1,i],weighting[2,i],weighting[3,i],weighting[4,i])
    plant_mat <- decostand(plant_mat, method = "total") # first relative abundance; maybe for comparison hellinger transformation as well
    # find best dimensionality: via NMDS slowly...
    stress_vals[[i]][[o]] <- future_sapply(2:6, function(k){
      median(replicate(2, metaMDS(plant_mat, k = k, trymax = 3, trace = FALSE)$stress))
    }, future.seed = TRUE)
  }
}

#saveRDS(list(data_list,stress_vals), "data_with_stress_vals.RDS")


par(mfrow = c(2,4))
for(i in stress_vals){
  sapply(i, function(x){
    plot(2:6, x, type = "b")
    # approximate "elbow"
    diff1 <- diff(x)
    diff2 <- diff(diff1)
    
    k_opt <- which.min(diff2) + 1
    print(k_opt)
  })
}

stress_vals_df <- data.frame(data = rep(names(data_list),each = 5),
                             dimensions = rep(seq(2:6),7),
                             imbalanced = rep(0,length(rep(names(data_list),each = 5))),
                             balanced = rep(0, length(rep(names(data_list),each = 5))))


for(i in 1:2){
  results_stress <- c()
  for (o in 1:7){
    
    results_stress <- c(results_stress,stress_vals[[i]][[o]])
  }
  stress_vals_df[,i+2] <- results_stress
  results_stress <- c()
}

stress_vals_df_plot <- pivot_longer(stress_vals_df, cols = imbalanced:balanced,
                                    values_to = "stress", names_to = "balance")

ggplot(stress_vals_df_plot)+
  geom_point(aes(x = dimensions, y = stress, colour = balance))+
  theme_minimal()+
  facet_wrap(facets = "data")

# balanced is ignored
dimensions_balanced <- c(6,5,NA, NA,2,4,NA) # for each list one dimension value to use for hdbscan
dimensions_inbalanced <- c(6,5,5,NA,2,6,6)
nmds_objects_imbalanced <- list(imbalanced = list())

dimensions = dimensions_inbalanced


plants <- list()
for (o in 1:length(data_list)) {
  print(o)
  if(is.na(dimensions[o]) == FALSE){
    print(o)
    plants[[names(data_list)[o]]] <- data_list[[o]][[1]]
  }
}

i = 1
w1 =weighting[1,i]
w2 = weighting[2,i]
w3 = weighting[3,i]
w4 = weighting[4,i]

dimensions_short <- dimensions[!is.na(dimensions)]
nmds_objects_imbalanced[[i]] <- future_lapply(seq_along(plants), function(o){
  plant_mat <- plant_weighting(plants[[o]],w1,w2,w3,w4)
  plant_mat <- decostand(plant_mat, method = "total") # first relative abundance; maybe for comparison hellinger transformation as well
  metaMDS(plant_mat, k = dimensions_short[o], trymax = 15)
}, future.seed = TRUE)

names(nmds_objects_imbalanced$imbalanced) <- names(plants)
nmds_objects_combined <- list(imbalanced = nmds_objects_imbalanced$imbalanced, balanced = nmds_objects$balanced)
# saveRDS(nmds_objects_combined, file = "nmds_objects_all_260404.Rds")


counter = 1
metric_vals_bund = list()
for (i in nmds_objects_combined$imbalanced) {
  metric_vals_bund[[counter]] <- hdbscan_complete(i$points, by = 2,grunddat = data_list[[names(plants)[counter]]][[2]])
  counter = counter+1
}

counter = 1
metric_vals_land = list()
for (i in nmds_objects_combined$imbalanced) {
  metric_vals_land[[counter]] <- hdbscan_complete(i$points, by = 2,grunddat = data_list[[names(plants)[counter]]][[2]], bund = FALSE)
  counter = counter+1
}

# data set 1 with minPts 13 and bt_land

#check whether ordination posses outliers
for(i in nmds_objects_combined){
  outliers <- lapply(i, function(x){
    ordination_outlier_func(x)
  })
  print(outliers)
}

# forest_complete_wone imbalanced has no outliers:
hdbscan_complete(nmds_objects_combined$imbalanced[[2]]$points, by = 1,grunddat = data_list[[names(plants)[2]]][[2]])
evaluation_plt_data <- hdbscan_complete(nmds_objects_combined$imbalanced[[2]]$points, by = 1,
                                        grunddat = data_list[[names(plants)[2]]][[2]], bund = FALSE, coarse = FALSE)
# land best
hdbscan_plot(evaluation_plt_data, name = names(nmds_objects_combined$imbalanced)[[2]])

evaluation_plt_data_coarse <- hdbscan_complete(nmds_objects_combined$imbalanced[[2]]$points, by = 1,
                                               grunddat = data_list[[names(plants)[2]]][[2]], bund = FALSE, coarse = TRUE)
hdbscan_plot(evaluation_plt_data_coarse, name = paste0(names(nmds_objects_combined$imbalanced)[[2]], "_coarse"))




test_one <- nmds_objects_combined$imbalanced[[2]]
ordination_outlier_func(test_one)

test_one_hdbscan <- hdbscan(test_one$points, minPts = 13)
# Purity: How dominant is the main class within each cluster:
# > 0.8	very good; 0.6–0.8	reasonable; < 0.6	weak

# Adjusted Rand Index:
# ~0	random; 0.2–0.4	weak structure; 0.4–0.6	moderate; >0.6	strong; >0.8	excellent


# cluster-wise purity
prop.table(tab, margin = 1)



# 
# check_forest_outliers <- filter(plants_forests_red, Polygon %in% plants_forests_red$Polygon[order(plants_forests_rel_ord$points[,1], decreasing = TRUE)[1:3]])
# check_forest_outliers <- plants_forests_red[outlier_hdbscan,]
# filter(grunddaten_forests_red, Polygon %in% plants_forests_red$Polygon[order(plants_forests_rel_ord$points[,1], decreasing = TRUE)[1:3]])


# display just first two NMDS axes
plot_data <- as.data.frame(test_one$points[, 1:2])
plot_data$cluster <- factor(test_one_hdbscan$cluster)

plot(plot_data)

# project multidemensional data into 2-D
pca <- prcomp(plants_forests_rel_ord_5$points)

plot_data <- as.data.frame(pca$x[, 1:2])
plot_data$cluster <- factor(hdbscan_forest_red_5$cluster)

ggplot(plot_data, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2, alpha = 0.8) +
  #scale_color_manual(values = c("0" = "grey70")) +
  theme_minimal()

plot_data$biotope <- grunddaten_forests_red$`Biotoptyp-Bund`
ggplot(plot_data, aes(PC1, PC2, color = biotope)) +
  geom_point(size = 2) +
  theme_minimal()

ggplot(plot_data, aes(PC1, PC2, color = cluster)) +
  geom_point() +
  stat_ellipse() +
  theme_minimal()


# 3D-plotting -------------------------------------------------------------

df <- hdbscan_mismatch_evaluation(plants_nmds = test_one,
                                  plants_hdbscan = test_one_hdbscan,
                                  grunddat = data_list[[names(plants)[2]]][[2]],
                                  coarse = TRUE)

sum(df$mismatch) # not fitting into the cluster
sum(df$mismatch==FALSE)

plot_ly(df,
        x = ~NMDS1, y = ~NMDS2, z = ~NMDS3,
        color = ~mismatch,
        colors = c("FALSE" = "black", "TRUE" = "red"),
        type = "scatter3d",
        mode = "markers")



plot_ly(df,
        x = ~NMDS1,
        y = ~NMDS2,
        z = ~NMDS3,
        color = ~cluster,
        colors = "Set1",
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 3)) %>%
  layout(scene = list(xaxis = list(title = "NMDS1"),
                      yaxis = list(title = "NMDS2"),
                      zaxis = list(title = "NMDS3")))

# for biotope colouring
colourCount = length(unique(df$biotope))
getPalette = colorRampPalette(colors = c("red","green", "blue"))

plot_ly(df,
        x = ~NMDS1,
        y = ~NMDS2,
        z = ~NMDS3,
        color = ~biotope,
        colors = getPalette(colourCount),
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 3)) %>%
  layout(scene = list(xaxis = list(title = "NMDS1"),
                      yaxis = list(title = "NMDS2"),
                      zaxis = list(title = "NMDS3")))


# with hovering
df$label <- paste("Cluster:", df$cluster,
                  "<br>Biotope:", df$biotope)

plot_ly(df,
        x = ~NMDS1, y = ~NMDS2, z = ~NMDS3,
        color = ~cluster,
        colors = getPalette(colourCount),
        text = ~label,
        hoverinfo = "text",
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 3)) %>%
  layout(scene = list(xaxis = list(title = "NMDS1"),
                      yaxis = list(title = "NMDS2"),
                      zaxis = list(title = "NMDS3")))

hover_3D(df=df)

# remove AG-biotope types -------------------------------------------------

forest_complete_wone_wAG <- forest_complete_wone
wone_wAG <- dplyr::filter(forest_complete_wone[[2]], !substr(`Biotoptyp-Land`,1,2) == "AG")
plants_wone_wAG <- dplyr::filter(forest_complete_wone[[1]], Polygon %in% wone_wAG$Polygon)
wone_wAG_list <- list(plants = plants_wone_wAG, grunddaten = wone_wAG)

stress_vals_wone_wAG <- list(imbalanced = list(), balanced = list())

# future::plan(sequential)
# gc() # optional, frees memory

for(i in 1:2){
  plants <- plants_wone_wAG
  plant_mat <- plant_weighting(plants,weighting[1,i],weighting[2,i],weighting[3,i],weighting[4,i])
  plant_mat <- decostand(plant_mat, method = "total") # first relative abundance; maybe for comparison hellinger transformation as well
  # find best dimensionality: via NMDS slowly...
  stress_vals_wone_wAG[[i]][[1]] <- future_sapply(2:6, function(k){
    median(replicate(2, metaMDS(plant_mat, k = k, trymax = 3, trace = 0)$stress))
  }, future.seed = TRUE)
}

stress_vals_df_wAG <- data.frame(dimensions = seq(2,6),
                                 imbalanced = rep(0, 5),
                                 balanced = rep(0, 5))

for(i in 1:2){
  results_stress_wAG <- c()
  results_stress_wAG <- c(results_stress_wAG,stress_vals_wone_wAG[[i]][[1]])
  stress_vals_df_wAG[,i+1] <- results_stress_wAG
  results_stress_wAG <- c()
}

stress_vals_df_wAG_plot <- pivot_longer(stress_vals_df_wAG, cols = imbalanced:balanced,
                                        values_to = "stress", names_to = "balance")

ggplot(stress_vals_df_wAG_plot)+
  geom_point(aes(x = dimensions, y = stress, colour = balance))+
  theme_minimal()

# best dimension


nmds_object_imbalanced_wAG <- list()
plants_list_wAG <- list(plants_wone_wAG,plants_wone_wAG)
nmds_object_imbalanced_wAG[[1]] <- future_lapply(seq_along(plants_list_wAG), function(o){
  w1 =weighting[1,o]
  w2 = weighting[2,o]
  w3 = weighting[3,o]
  w4 = weighting[4,o]
  plant_mat <- plant_weighting(plants_list_wAG[[o]],w1,w2,w3,w4)
  plant_mat <- decostand(plant_mat, method = "total") # first relative abundance; maybe for comparison hellinger transformation as well
  metaMDS(plant_mat, k = 6, trymax = 15)
}, future.seed = TRUE)

names(nmds_object_imbalanced_wAG[[1]]) <- c("imbalanced", "balanced")
#saveRDS(nmds_object_imbalanced_wAG, "nmds_objects_forests_wone_wAG_dim6.RDS")

####
counter = 1
metric_vals_wAG = list(imbalanced = list(), balanced = list())
for (i in nmds_object_imbalanced_wAG[[1]]) {
  metric_vals_wAG[[names(nmds_object_imbalanced_wAG[[1]])[counter]]][[1]] <- hdbscan_complete(i$points, by = 2,grunddat = wone_wAG)
  metric_vals_wAG[[names(nmds_object_imbalanced_wAG[[1]])[counter]]][[2]] <- hdbscan_complete(i$points, by = 2,grunddat = wone_wAG, bund = FALSE)
  counter = counter+1
}

### imbalanced with minPts 13 and bt_land
#load("250405_1257_Rdata")

#check whether ordination posses outliers
for(i in nmds_object_imbalanced_wAG){
  outliers <- lapply(i, function(x){
    ordination_outlier_func(x)
  })
  print(outliers)
}

# forest_complete_wone imbalanced has no outliers:
evaluation_plt_data <- hdbscan_complete(nmds_object_imbalanced_wAG[[1]][[1]]$points, by = 1,grunddat = wone_wAG, bund = FALSE)

# land best
hdbscan_plot(evaluation_plt_data, name = names(nmds_objects_combined$imbalanced)[[2]])

evaluation_plt_data_coarse <- hdbscan_complete(nmds_object_imbalanced_wAG[[1]][[1]]$points, by = 1,grunddat = wone_wAG, bund = FALSE,coarse = TRUE)

hdbscan_plot(evaluation_plt_data_coarse, name = paste0(names(nmds_objects_combined$imbalanced)[[2]], "_coarse"))

### imbalanced with minPts 10 and bt_land and coarse = 10

hdbscan_wAG <- hdbscan(nmds_object_imbalanced_wAG[[1]][[1]]$points, minPts = 10)

df_wAG <- hdbscan_mismatch_evaluation(plants_nmds = nmds_object_imbalanced_wAG[[1]][[1]],
                                      plants_hdbscan = hdbscan_wAG,
                                      grunddat = wone_wAG,
                                      coarse = TRUE)

sum(df_wAG$mismatch) # not fitting into the cluster
sum(df_wAG$mismatch==FALSE)

hover_3D(df_wAG)
hull_3D(df_wAG, op_hull = 0.6, op_points = 0.3)

table(hdbscan_wAG$cluster, wone_wAG$`Biotoptyp-Land`)
table(hdbscan_wAG$cluster, wone_wAG$`BT_Land_group`)

prop.table(table(hdbscan_wAG$cluster, wone_wAG$`BT_Land_group`), margin=1)
prop.table(t(table(hdbscan_wAG$cluster, wone_wAG$`BT_Land_group`)), margin=1)

ggplot(as.data.frame(table(hdbscan_wAG$cluster, wone_wAG$`BT_Land_group`)),
       aes(Var2, Freq, colour = Var1, fill = Var1)) +
  geom_bar(stat = "identity")+
  theme_minimal()


# GMM ---------------------------------------------------------------------

#load("260406_1001.Rdata")

grid <- expand.grid(
  i = 1:2,
  o = seq_along(data_list)
)


results <- future_lapply(seq_len(nrow(grid)), function(idx) {
  i <- grid$i[idx]
  o <- grid$o[idx]
  
  name <- names(data_list)[o]
  plants <- data_list[[o]][[1]]
  
  plant_mat <- plant_weighting(
    plants,
    weighting[1,i], weighting[2,i],
    weighting[3,i], weighting[4,i]
  )
  
  plant_total <- decostand(plant_mat, method = "total")
  plant_hell  <- decostand(plant_mat, method = "hellinger")
  
  # GMM
  model_total <- Mclust(plant_total)
  model_hell  <- Mclust(plant_hell)
  
  # PCA + GMM
  pca_total <- prcomp(plant_total)
  model_pca_total <- Mclust(pca_total$x[, 1:10])
  
  pca_hell <- prcomp(plant_hell)
  model_pca_hell <- Mclust(pca_hell$x[, 1:10])
  
  list(
    name = name,
    i = i,
    models_total = list(
      raw = model_total,
      pca = model_pca_total
    ),
    models_hell = list(
      raw = model_hell,
      pca = model_pca_hell
    )
  )
}, future.seed = TRUE)

#results <- readRDS("GMM_results_all_260406.Rds")

gmm_models_total <- list()
gmm_models_hellinger <- list()

for (res in results) {
  name <- res$name
  i <- res$i
  
  if (is.null(gmm_models_total[[name]])) {
    gmm_models_total[[name]] <- vector("list", 4)
    gmm_models_hellinger[[name]] <- vector("list", 4)
  }
  
  gmm_models_total[[name]][[i]]     <- res$models_total$raw
  gmm_models_total[[name]][[i+2]]   <- res$models_total$pca
  
  gmm_models_hellinger[[name]][[i]]   <- res$models_hell$raw
  gmm_models_hellinger[[name]][[i+2]] <- res$models_hell$pca
}

# visualize/evaluate results

gmm_total_eval <- evaluate_gmm(gmm_models_total)

for(i in gmm_total_eval$tab_plot){
  print(ggplot(data = i, aes(x = Biotope, y = Cluster, fill = Freq)) +
          geom_tile() +
          scale_fill_gradient(low = "white", high = "blue") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))+
          facet_wrap(facet = "run")
        )
}
for(i in gmm_total_eval$ari_values){
  print(i)
}

gmm_hell_eval <- evaluate_gmm(gmm_models_hellinger)
for(i in gmm_hell_eval$ari_values){
  print(i)
}
for(i in gmm_hell_eval$bic){
  print(i)
}
for(i in gmm_hell_eval$uncertainty){
  print(i)
}
# further analysis --------------------------------------------------------


# reveal indicator species
library(indicspecies)

res <- multipatt(plants_forests_rel_short, hdbscan_forest_red_5$cluster, control = how(nperm=999))
summary(res)

#alternative
pcoa <- cmdscale(plants_forests_bray_short, k = 6)
apply(pcoa, 2, var)
pcoa_hdb <- hdbscan(pcoa[, 1:5], minPts = 10)
table(pcoa_hdb$cluster, grunddaten_forests_red$`Biotoptyp-Bund`)

# Need to check later... --------------------------------------------------




# Hierarchical clustering
hc_forests <- hclust(plants_forests_bray, method = "ward.D2")
plot(hc_forests)
clusters <- cutree(hc_forests, k = 5) # ???
plot(clusters)

# Evaluation of cluster number
library(factoextra)
fviz_nbclust(plants_forests_rel_short, FUN = hcut, method = "silhouette", k.max = 25)
#library(cluster)
sil <- silhouette(clusters, dist_bc)
mean(sil[,3])

# Model-based clusterin
library(mclust)
mc <- Mclust(cmdscale(dist_bc, k=10)$points)

# Spectral clustering
library(kernlab)
sc <- specc(as.matrix(plants_forests_bray), centers = 6)




hdbscan_forest <- hdbscan(plants_forests_bray_short, minPts = 7)
codes_forests_hdb <- cbind(hdbscan_forest$cluster, grunddaten_forests[,c(1,2,4)])

t(table(hdbscan_forest$cluster, grunddaten_forests$`Biotoptyp-Bund`))
t(table(hdbscan_forest$cluster, grunddaten_forests$`Biotoptyp-Land`))

prop.table(table(hdbscan_forest$cluster, grunddaten_forests$`Biotoptyp-Land`), margin=1)


# Reduced dataset ---------------------------------------------------------


# Hellinger transformation to subsequently compute Euclidean distance
plants_hell <- hellinger(plants_wide[,-1])
test_dist <- dist(plants_hell)

#Alternative...
comm_hel <- decostand(plants_wide[,-1], method = "hellinger")
dist_hel <- dist(comm_hel, method = "euclidean")

# Bray-Curtis dissimilarity with raw abundance data

bray_curtis_dist <- vegdist(plants_wide[,-1], method = "bray")

# Gower distance
library(cluster)

# convert species abundances into ordered factors
comm_ord <- as.data.frame(plants_wide[,-1])
comm_ord[] <- lapply(plants_wide[,-1], ordered)

dist_gower <- daisy(comm_ord, metric = "gower")

# Community abundance-weighted transformation
plants_weighted <- plants_wide[,-1]
plants_weighted[plants_weighted == 1] <- 0.01
plants_weighted[plants_weighted == 2] <- 0.05
plants_weighted[plants_weighted == 3] <- 0.25
plants_weighted[plants_weighted == 4] <- 0.75

bray_curtis_dist_weighted <- vegdist(plants_weighted, method = "bray")

# https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/common-distance-measures/

# Ward's algorithm with Bray-Curtis distance metric: https://www.davidzeleny.net/anadat-r/doku.php/en:class-eval_examples

# HDBSCAN: https://rdrr.io/cran/dbscan/f/vignettes/hdbscan.Rmd

# https://r.qcbs.ca/workshop09/book-en/clustering.html



test_hdb <- hdbscan(bray_curtis_dist, minPts = 10)
table(test_hdb$cluster)


test_hdb_bcd_weighted <- hdbscan(bray_curtis_dist_weighted, minPts = 6)
table(test_hdb_bcd_weighted$cluster)





# For euclidean distance
test_hdb2 <- hdbscan(dist_hel, minPts = 10)
table(test_hdb2$cluster)

# For Gower distance
test_hdb_gower <- hdbscan(dist_gower, minPts = 10)
table(test_hdb_gower$cluster)

# for plotting, dimensionality reduction of distance matrix (ordination)
ord <- metaMDS(bray_curtis_dist, k = 2, trymax = 100)
plot(ord, type = "t")


coords <- as.data.frame(scores(ord))
coords$cluster <- factor(grass_hdbscan$cluster)

test_cluster_bt <- data.frame(cbind(cluster=test_hdb_bcd_weighted$cluster, number= test_hdb_bcd_weighted$hc$order))
test_cluster_bt <- test_cluster_bt%>%
  arrange(number)
test_cluster_bt <- cbind(test_cluster_bt, Polygon=plants_wide$Polygon)
join_cluster_code <- inner_join(test_cluster_bt, grunddaten_sub[,c(1,2,4)], by = "Polygon")

ggplot(coords, aes(NMDS1, NMDS2, color = cluster)) +
  geom_point(size = 2, alpha = 0.8) +
  #scale_color_manual(
  # values = c("0" = "grey70", scales::hue_pal()(length(unique(coords$cluster)) - 1))
  #) +
  labs(
    title = "HDBSCAN Clustering (PCoA of Bray–Curtis)",
    x = "PCoA 1",
    y = "PCoA 2",
    color = "Cluster"
  ) +
  theme_minimal()