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
library(cluster)
library(factoextra)

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

data_list[["forest_complete_without_plots_only_one"]] <- forest_complete_wone

forest_complete_wtwo <- remove_plots(plants_forests_red, grunddaten_forests_red, remove_count = c(0:2))
data_list[["forest_complete_without_plots_only_up_to_two"]] <- forest_complete_wtwo

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

### create an additional plant list simplified to genus-level

# collapsing to genus level
plant_forest_genus <- plants_forests_red %>%
  pivot_longer(-Polygon, names_to = "species", values_to = "abundance") %>%
  mutate(genus = sub("(\\w+).*", "\\1", species)) %>%
  group_by(Polygon, genus) %>%
  summarise(abundance = max(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = genus, values_from = abundance, values_fill = 0)

forest_complete_wone_genus <- remove_plots(plant_forest_genus, grunddaten_forests_red, remove_count = c(0,1))
check_empty_plot(forest_complete_wone_genus[[1]])

data_list[["forest_genus_without_plots_only_one"]] <- forest_complete_wone_genus

forest_complete_wtwo_genus <- remove_plots(plant_forest_genus, grunddaten_forests_red, remove_count = c(0:2))
data_list[["forest_genus_without_plots_only_up_to_two"]] <- forest_complete_wtwo_genus


# control number of tree species per plot -------------------------------------------------

check_empty_plot(plants_f_trees_wone)
plants_trees_wzero <- remove_plots(plants_f_trees_wone, grunddaten_forests_red)
check_empty_plot(plants_trees_wzero[[1]])

data_list[["plants_without_zero_trees"]] <- plants_trees_wzero

plants_trees_wone <- remove_plots(plants_f_trees_wone, grunddaten_forests_red, remove_count = c(0,1))
check_empty_plot(plants_trees_wone[[1]])
data_list[["plants_without_one_tree"]] <- plants_trees_wone

plants_trees_wtwo <- remove_plots(plants_f_trees_wone, grunddaten_forests_red, remove_count = c(0:2))
check_empty_plot(plants_trees_wtwo[[1]])
data_list[["plants_without_only_up_to_two_trees"]] <- plants_trees_wtwo

# remove AG-biotope types -------------------------------------------------

data_list[["forest_without_AG&plots_only_one"]] <- remove_land_biotope_code(forest_complete_wone[[1]], forest_complete_wone[[2]],c("AG"))
data_list[["trees_without_AG&plots_without_zero"]] <- remove_land_biotope_code(plants_trees_wzero[[1]], plants_trees_wzero[[2]],c("AG"))
data_list[["trees_without_AG&plots_without_up_to_one"]] <- remove_land_biotope_code(plants_trees_wone[[1]], plants_trees_wone[[2]],c("AG"))

data_list[["trees_genus_without_AG&plots_without_up_to_one"]] <- remove_land_biotope_code(forest_complete_wone_genus[[1]], forest_complete_wone_genus[[2]],c("AG"))


# check direct hdbscan application ----------------------------------------

hdbscan_direct_list <- list()
for(i in 1:length(data_list)){
  weighting_plants <- weight_rel_dist(data_list[[i]]$plants)
  hdbscan_run <- hdbscan_complete(weighting_plants,grunddat = data_list[[i]]$grunddaten, by = 1, print = FALSE)
  hdbscan_direct_list[[names(data_list)[i]]][[1]] <- hdbscan_run
  weighting_plants <- weight_rel_dist(data_list[[i]]$plants, w1 = 0.01, w2 = 0.05, w3 = 0.25, w4 = 0.9)
  hdbscan_run <- hdbscan_complete(weighting_plants,grunddat = data_list[[i]]$grunddaten, by = 1, print = FALSE)
  hdbscan_direct_list[[names(data_list)[i]]][[2]] <- hdbscan_run
}

for (i in 1:length(hdbscan_direct_list)) {
  print(hdbscan_plot(hdbscan_direct_list[[i]][[1]], name = paste0(names(hdbscan_direct_list)[i],"-imbalanced")))
  print(hdbscan_plot(hdbscan_direct_list[[i]][[2]], name = paste0(names(hdbscan_direct_list)[i],"-balanced")))
}

hdbscan_result_df <- function(hdbscan_list, main_list){
  result_df <- bind_rows(lapply(names(hdbscan_list), function(main_name) {
    sublist <- hdbscan_list[[main_name]]
    n_total <- nrow(main_list[[main_name]][[1]])  
    bind_rows(lapply(seq_along(sublist), function(sub_id) {
      df <- sublist[[sub_id]]
    df %>%
      mutate(
        list_name = main_name,
        sublist_id = sub_id,
        noise_prop = noise/n_total
      )
  }))
  }))
  return(result_df)}

hdbscan_direct_result <- hdbscan_result_df(hdbscan_direct_list, data_list)

hdbscan_direct_result$composed_metric <- (0.5*hdbscan_direct_result$ari+
                                            0.4*hdbscan_direct_result$purity+0.1*(1-hdbscan_direct_result$noise_pro))
hdbscan_direct_result$composed_metric_balanced <- (0.3*hdbscan_direct_result$ari+
                                            0.3*hdbscan_direct_result$purity+0.3*(1-hdbscan_direct_result$noise_pro))
hdbscan_direct_result %>%
  group_by(list_name, sublist_id) %>%
  slice_max(order_by = composed_metric, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(composed_metric))

hdbscan_direct_result %>%
  filter(k < 16 & noise < 1000)%>%
  group_by(list_name, sublist_id) %>%
  slice_max(order_by = composed_metric, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(composed_metric))

ggplot(hdbscan_direct_result)+
  geom_line(aes(x=k,y= composed_metric, colour = list_name))+
  theme_classic()+
  facet_wrap(facets = "sublist_id")

#png("results/combined_evaluation_hdbscan.png", width = 3000, height = 2100, res = 300)
ggplot(hdbscan_direct_result)+
  geom_line(aes(x=k,y= composed_metric_balanced, colour = list_name))+
  theme_classic()+
  scale_color_discrete(name = "Plant data composition")+
  labs(title = "equal weighting for performance", y = "performance (ARI + purity + noise)", x = "minimum size of cluster")+
  facet_wrap(facets = "sublist_id", labeller = labeller(sublist_id=c("1"= "imbalanced conversion", "2"= "more balanced conversion")))
#dev.off()

max(hdbscan_direct_result$ari)


# bund and coasre

hdbscan_direct_list_coarse <- list()
for(i in 1:length(data_list)){
  weighting_plants <- weight_rel_dist(data_list[[i]]$plants)
  hdbscan_run <- hdbscan_complete(weighting_plants,grunddat = data_list[[i]]$grunddaten, by = 1, print = FALSE, bund = FALSE, coarse = TRUE)
  hdbscan_direct_list_coarse[[names(data_list)[i]]][[1]] <- hdbscan_run
  weighting_plants <- weight_rel_dist(data_list[[i]]$plants, w1 = 0.01, w2 = 0.05, w3 = 0.25, w4 = 0.9)
  hdbscan_run <- hdbscan_complete(weighting_plants,grunddat = data_list[[i]]$grunddaten, by = 1, print = FALSE, bund = FALSE, coarse = TRUE)
  hdbscan_direct_list_coarse[[names(data_list)[i]]][[2]] <- hdbscan_run
  
}
hdbscan_direct_result_coarse <- hdbscan_result_df(hdbscan_direct_list_coarse,data_list)

hdbscan_direct_result_coarse$composed_metric <- (0.5*hdbscan_direct_result_coarse$ari+
                                                   0.4*hdbscan_direct_result_coarse$purity+0.1*(1-hdbscan_direct_result_coarse$noise_pro))
hdbscan_direct_result_coarse$composed_metric_balanced <- (0.3*hdbscan_direct_result_coarse$ari+
                                                   0.3*hdbscan_direct_result_coarse$purity+0.3*(1-hdbscan_direct_result_coarse$noise_pro))


hdbscan_direct_result_coarse %>%
  group_by(list_name, sublist_id) %>%
  slice_max(order_by = composed_metric, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(composed_metric))

#png("results/combined_evaluation_hdbscan_coarse_noise.png", width = 3000, height = 2100, res = 300)
ggplot(hdbscan_direct_result_coarse)+
  geom_line(aes(x=k,y= composed_metric_balanced, colour = list_name))+
  theme_classic()+
  scale_color_discrete(name = "Plant data composition")+
  labs(title = "equal weighting for performance",y = "performance (ARI + purity + noise)", x = "minimum size of cluster")+
  facet_wrap(facets = "sublist_id", labeller = labeller(sublist_id=c("1"= "imbalanced conversion", "2"= "more balanced conversion")))
#dev.off()

#png("results/combined_evaluation_hdbscan_coarse.png", width = 3000, height = 2100, res = 300)
ggplot(hdbscan_direct_result_coarse)+
  geom_line(aes(x=k,y= composed_metric, colour = list_name))+
  theme_classic()+
  scale_color_discrete(name = "Plant data composition")+
  labs(title = "noise less weighted for performance" ,y = "performance (ARI + purity + noise)", x = "minimum size of cluster")+
  facet_wrap(facets = "sublist_id", labeller = labeller(sublist_id=c("1"= "imbalanced conversion", "2"= "more balanced conversion")))
#dev.off()

hdbscan_complete(plants_dist = weight_rel_dist(data_list[["forest_without_AG&plots_only_one"]][[1]]),
                 grunddat = data_list[["forest_without_AG&plots_only_one"]][[2]], coarse = TRUE, bund = FALSE, by = 1)

clusterVScode(plants_dist = weight_rel_dist(data_list[["forest_without_AG&plots_only_one"]][[1]]),
              grunddat = data_list[["forest_without_AG&plots_only_one"]][[2]], pts = 11, bund = FALSE)

### temporarily visualisation

hdbscan_plot <- hdbscan(weight_rel_dist(data_list[["forest_without_AG&plots_only_one"]][[1]]), minPts = 11)

# project multidemensional data into 2-D
plant_mat <-plant_weighting(data_list[["forest_without_AG&plots_only_one"]][[1]])
plant_mat <- decostand(plant_mat, method = "total")
pca_plants <- prcomp(plant_mat)

# two-dimensional
plot_data <- as.data.frame(pca_plants$x[, 1:2])
plot_data$cluster <- factor(hdbscan_plot$cluster)

ggplot(plot_data, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2, alpha = 0.8) +
  #scale_color_manual(values = c("0" = "grey70")) +
  theme_minimal()

plot_data$biotope <- data_list[["forest_without_AG&plots_only_one"]][[2]]$`Biotoptyp-Land`
ggplot(plot_data, aes(PC1, PC2, color = biotope)) +
  geom_point(size = 2) +
  theme_minimal()

# 3-dimensional; wrong naming to use already coded function!
pca_3d <- list()
pca_3d[["points"]] <- pca_plants$x
hdbscan_eval <- hdbscan_mismatch_evaluation(pca_3d,hdbscan_plot,data_list[["forest_without_AG&plots_only_one"]][[2]], coarse = TRUE)
sum(hdbscan_eval$mismatch==TRUE)
sum(hdbscan_eval$cluster==0)

hover_3D(hdbscan_eval)
hull_3D(hdbscan_eval, op_hull = 0.6, op_points = 0.3)

######

selection_scenarios <- hdbscan_direct_result_coarse %>%
  filter(k > 5) %>%
  group_by(list_name, sublist_id) %>%
  summarise(mean = mean(composed_metric)) %>%
  arrange(desc(mean))

# scenarios to keep for NMDS: 


# NMDS loop ---------------------------------------------------------------
weighting <- data.frame(imbalanced = c(0.01,0.01,0.01,1),
                        balanced = c(0.01,0.05,0.25,0.9))
#stress_vals_backup <- stress_vals
stress_vals <- list(imbalanced = list(), balanced = list())

for(i in 1:2){
  for (o in unique(selection_scenarios$list_name[1:12])) { # after the 12 first scenarios, there appears to be a gap
    plants <- data_list[[o]][[1]]
    plant_mat <- plant_weighting(plants,weighting[1,i],weighting[2,i],weighting[3,i],weighting[4,i])
    plant_mat <- decostand(plant_mat, method = "total") # first relative abundance; maybe for comparison hellinger transformation as well
    # find best dimensionality: via NMDS slowly...
    stress_vals[[i]][[o]] <- future_sapply(2:6, function(k){
      median(replicate(2, metaMDS(plant_mat, k = k, trymax = 3, trace = FALSE)$stress))
    }, future.seed = TRUE)
  }
}

saveRDS(list(data_list,stress_vals), "results/data_with_stress_vals.RDS")


par(mfrow = c(2,4))
for(i in stress_vals){
  sapply(i, function(x){
    plot(2:6, x, type = "b", main = names(i))
    # approximate "elbow"
    diff1 <- diff(x)
    diff2 <- diff(diff1)
    
    k_opt <- which.min(diff2) + 1
    print(k_opt)
  })
}

stress_vals_df <- data.frame(data = rep(unique(selection_scenarios$list_name[1:12]),each = 5),
                             dimensions = rep(seq(2,6),7),
                             imbalanced = rep(0,length(rep(unique(selection_scenarios$list_name[1:12]),each = 5))),
                             balanced = rep(0, length(rep(unique(selection_scenarios$list_name[1:12]),each = 5))))


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

# # balanced is ignored
# dimensions_balanced <- c(6,5,NA, NA,2,4,NA) # for each list one dimension value to use for hdbscan
# dimensions_inbalanced <- c(6,5,5,NA,2,6,6)
# nmds_objects_imbalanced <- list(imbalanced = list())
# 
# dimensions = dimensions_inbalanced
# 
# 
# plants <- list()
# for (o in 1:length(data_list)) {
#   print(o)
#   if(is.na(dimensions[o]) == FALSE){
#     print(o)
#     plants[[names(data_list)[o]]] <- data_list[[o]][[1]]
#   }
# }
# 
# i = 1
# w1 =weighting[1,i]
# w2 = weighting[2,i]
# w3 = weighting[3,i]
# w4 = weighting[4,i]
# 
# dimensions_short <- dimensions[!is.na(dimensions)]
# nmds_objects_imbalanced[[i]] <- future_lapply(seq_along(plants), function(o){
#   plant_mat <- plant_weighting(plants[[o]],w1,w2,w3,w4)
#   plant_mat <- decostand(plant_mat, method = "total") # first relative abundance; maybe for comparison hellinger transformation as well
#   metaMDS(plant_mat, k = dimensions_short[o], trymax = 15)
# }, future.seed = TRUE)
# 
# names(nmds_objects_imbalanced$imbalanced) <- names(plants)
# nmds_objects_combined <- list(imbalanced = nmds_objects_imbalanced$imbalanced, balanced = nmds_objects$balanced)
# # saveRDS(nmds_objects_combined, file = "nmds_objects_all_260404.Rds")
# 
# 
# counter = 1
# metric_vals_bund = list()
# for (i in nmds_objects_combined$imbalanced) {
#   metric_vals_bund[[counter]] <- hdbscan_complete(i$points, by = 2,grunddat = data_list[[names(plants)[counter]]][[2]])
#   counter = counter+1
# }
# 
# counter = 1
# metric_vals_land = list()
# for (i in nmds_objects_combined$imbalanced) {
#   metric_vals_land[[counter]] <- hdbscan_complete(i$points, by = 2,grunddat = data_list[[names(plants)[counter]]][[2]], bund = FALSE)
#   counter = counter+1
# }

###continue

plant_mat <- plant_weighting(data_list[["forest_without_AG&plots_only_one"]][[1]])
plant_mat <- decostand(plant_mat, method = "total") # first relative abundance; maybe for comparison hellinger transformation as well
#nmds_5_wag_wone <- metaMDS(plant_mat, k = 5, trymax = 20)

hdbscan_complete(nmds_5_wag_wone$points, by = 1, data_list[["forest_without_AG&plots_only_one"]][[2]], bund = FALSE, coarse = TRUE)
hdbscan_complete(weight_rel_dist(data_list[["forest_without_AG&plots_only_one"]][[1]]),
                 by = 1, data_list[["forest_without_AG&plots_only_one"]][[2]], bund = FALSE, coarse = TRUE)

#visualisation_3_nmds <- metaMDS(plant_mat, k = 3, trymax = 20)
hdbscan_complete(visualisation_3_nmds$points, by = 1, data_list[["forest_without_AG&plots_only_one"]][[2]], bund = FALSE, coarse = TRUE)

#visualisation_2_nmds <- metaMDS(plant_mat, k = 2, trymax = 20)

nmds_wag_wone <- list(nmds_5 = nmds_5_wag_wone$points, nmds_3 = visualisation_3_nmds$points, nmds_2 = visualisation_2_nmds$points,
                      direct = weight_rel_dist(data_list[["forest_without_AG&plots_only_one"]][[1]]))
saveRDS(nmds_wag_wone, "results/forest_wAG_c(0,1)_nmds.RDS")
nmds_wag_wone_metrics <- list()
for (i in 1:4) {
  evaluation_coarse <- hdbscan_complete(nmds_wag_wone[[i]],by = 1, data_list[["forest_without_AG&plots_only_one"]][[2]], bund = FALSE, coarse = TRUE, print = FALSE)
  evaluation <- hdbscan_complete(nmds_wag_wone[[i]],by = 1, data_list[["forest_without_AG&plots_only_one"]][[2]], bund = FALSE, print = FALSE)
  nmds_wag_wone_metrics[[names(nmds_wag_wone)[i]]][[1]] <- evaluation_coarse
  nmds_wag_wone_metrics[[names(nmds_wag_wone)[i]]][[2]] <- evaluation
}
helper_list <- list()
for(i in 1:4){
  helper_list[[names(nmds_wag_wone)[i]]][[1]] <- data_list[["forest_without_AG&plots_only_one"]][[1]]
}
nmds_wag_wone_compare <- hdbscan_result_df(nmds_wag_wone_metrics, main_list = helper_list)

nmds_wag_wone_compare$composed_metric <- (0.5*nmds_wag_wone_compare$ari+
                                                   0.4*nmds_wag_wone_compare$purity+0.1*(1-nmds_wag_wone_compare$noise_pro))
nmds_wag_wone_compare$composed_metric_balanced <- (0.3*nmds_wag_wone_compare$ari+
                                                            0.3*nmds_wag_wone_compare$purity+0.3*(1-nmds_wag_wone_compare$noise_pro))

#png("results/nmds_evaluation_hdbscan_coarse_noise.png", width = 3000, height = 2100, res = 300)
ggplot(nmds_wag_wone_compare)+
  geom_line(aes(x=k,y= composed_metric_balanced, colour = list_name))+
  theme_classic()+
  scale_color_discrete(name = "Plant data composition")+
  scale_y_continuous(limits = c(0.3, 0.8), breaks = c(0.40,0.60,0.8, 1))+
  labs(title = "forest withoutAG & c(0,1) - equal weighting",y = "performance (ARI + purity + noise)", x = "minimum size of cluster")+
  facet_wrap(facets = "sublist_id", labeller = labeller(sublist_id=c("1"= "coarse evaluation", "2"= "code-level evaluation")))
#dev.off()

#png("results/nmds_evaluation_hdbscan_coarse.png", width = 3000, height = 2100, res = 300)
ggplot(nmds_wag_wone_compare)+
  geom_line(aes(x=k,y= composed_metric, colour = list_name))+
  theme_classic()+
  scale_color_discrete(name = "Plant data composition")+
  scale_y_continuous(limits = c(0.15, 1), breaks = c(0.20,0.40,0.60,0.8, 1))+
  labs(title = "forest withoutAG & c(0,1) - noise less weighted" ,y = "performance (ARI + purity + noise)", x = "minimum size of cluster")+
  facet_wrap(facets = "sublist_id", labeller = labeller(sublist_id=c("1"= "coarse evaluation", "2"= "code-level evaluation")))
#dev.off()
### continue

# at the code level, the nmds dimensionality reduction can have advantages, especially if the aim is to reduce the number of points 
# classified as noise

# data set 1 with minPts 13 and bt_land

#check whether ordination posses outliers
for(i in 1:3){
  outliers <- #lapply(nmds_wag_wone[[i]], function(x){
    ordination_outlier_func(nmds_wag_wone[[i]])
  #})
  print(outliers)
}


# display just first two NMDS axes
plot_data_2d <- as.data.frame(visualisation_2_nmds$points)
plot_data_2d$cluster <- factor(hdbscan_plot$cluster)

plot(plot_data_2d)

ggplot(plot_data_2d, aes(x = MDS1, y = MDS2, color = cluster)) +
  geom_point(size = 2, alpha = 0.8) +
  #scale_color_manual(values = c("0" = "grey70")) +
  theme_minimal()

plot_data_2d$biotope <- data_list[["forest_without_AG&plots_only_one"]][[2]]$BT_Land_group
ggplot(plot_data_2d, aes(x = MDS1, y = MDS2, color = biotope)) +
  geom_point(size = 2) +
  theme_minimal()

ggplot(plot_data_2d, aes(x = MDS1, y = MDS2, color = cluster)) +
  geom_point() +
  stat_ellipse() +
  theme_minimal()

grunddat_temp <- data_list[["forest_without_AG&plots_only_one"]][[2]]
df_plot_3d <- hdbscan_mismatch_evaluation(visualisation_3_nmds, hdbscan_plot,
                            grunddat = grunddat_temp, coarse = TRUE)

sum(df_plot_3d$mismatch) # not fitting into the cluster

sum(df_plot_3d$mismatch[!df_plot_3d$cluster %in% c(0,17)])
sum(df_plot_3d$mismatch==FALSE)

hover_3D(df_plot_3d)
hover_3D(df_plot_3d[!df_plot_3d$cluster %in% c(0,17),])
hull_3D(df_plot_3d, op_hull = 0.6, op_points = 0.3)
hull_3D(df_plot_3d[!df_plot_3d$cluster %in% c(0,17),], op_hull = 0.6, op_points = 0.3)

####
legend_labels <- unique(df_plot_3d[,c(4, 7)])
legend_labels$combined <- paste0(legend_labels$cluster, ": ",legend_labels$cluster_main)
legend_labels <- legend_labels[order(legend_labels$cluster),]

#png("results/hdbscan_tree_coarse.png", width = 3000, height = 2100, res = 300)
plot(hdbscan_plot, show_flat = TRUE, main = "HDBSCAN tree")
legend("topright",
       legend = legend_labels$combined[-1],
       ncol = 3,
       title = "Main biotope code per cluster")
# dev.off()

####

### test 5d nmds into 3d
df_plot_3d5 <- hdbscan_mismatch_evaluation(nmds_5_wag_wone, hdbscan_plot,
                                          grunddat = grunddat_temp, coarse = TRUE)
hover_3D(df_plot_3d5)
hull_3D(df_plot_3d5, op_hull = 0.6, op_points = 0.3)  # nicer hull

####

table(hdbscan_plot$cluster, grunddat_temp$`Biotoptyp-Land`)
table(hdbscan_plot$cluster, grunddat_temp$`BT_Land_group`)

check_cluster_pureness <- as.data.frame(prop.table(table(hdbscan_plot$cluster, grunddat_temp$`BT_Land_group`), margin=1))
ggplot(data = check_cluster_pureness, aes(x = Var1, y = Freq, fill = Var2))+
  geom_bar(stat = "identity")
  
prop.table(t(table(hdbscan_plot$cluster, grunddat_temp$`BT_Land_group`)), margin=1)

ggplot(as.data.frame(table(hdbscan_plot$cluster, grunddat_temp$`BT_Land_group`)),
       aes(Var2, Freq, colour = Var1, fill = Var1)) +
  geom_bar(stat = "identity")+
  theme_minimal()


# UMAP --------------------------------------------------------------------

library(uwot)
library(rnndescent)

dist_mat <- plant_weighting(data_list[["forest_without_AG&plots_only_one"]][[1]])
dist_mat <- decostand(dist_mat, method = "total")

embedding <- umap(
  dist_mat,
  metric = "braycurtis",
  nn_method = "nndescent",
  n_neighbors = 10,
  n_components = 10,
  min_dist = 0.15
)

table(hdbscan_wAG_imbalanced$cluster, wone_wAG$`Biotoptyp-Land`)
table(hdbscan_wAG_imbalanced$cluster, wone_wAG$`BT_Land_group`)

ggplot(data = umap_cluster_pureness, aes(x = Var1, y = Freq, fill = Var2))+
  geom_bar(stat = "identity")

ggplot(as.data.frame(table(umap_hdbscan$cluster, grunddat_temp$`BT_Land_group`)),
       aes(Var2, Freq, colour = Var1, fill = Var1)) +
  geom_bar(stat = "identity")+
  theme_minimal()

### umap with the right setting allows to decrease the noise fraction in return for lower ari

# second cluster round -------------------------------------------------

plants_second_round <- data_list[["forest_without_AG&plots_only_one"]][[1]][hdbscan_plot$cluster %in% c(0,17),]
grunddat_second_round <- filter(grunddaten_forests_red, Polygon %in% plants_second_round$Polygon)

plants_second_max <- apply(plants_second_round[,-1], 1, max)
table(as.numeric(plants_second_max)) # I need to check the 1 and 2s


check_plants <- cbind(plants_second_max, plants_second_round)
check_plants <- right_join(grunddat_temp[,c(1,4,5,25)], check_plants, by = "Polygon")

check_plants[check_plants == 0] = NA
check_polygons <- pivot_longer(check_plants, cols = where(is.numeric)& !Polygon & !plants_second_max ,names_to = "plant", values_to = "abundance", values_drop_na = TRUE) 
                         # plants_second_round[order(plants_second_round$Polygon),]



max_weigthing <- function(plants_mat, w1 = 0.01, w2 = 0.01, w3 = 0.01, w4 = 1, method = "bray"){
  plant_return <- plants_mat[-1]
  plants_max <- apply(plants_mat[,-1], 1, max)
  for(i in 1:nrow(plants_mat[,1])){
    if(plants_max[i]== 4){
      plant_return[i,] = plant_weighting(plants_mat[i,], w1 = w1, w2 = w2, w3 = w3, w4 = w4)
    } else if(plants_max[i]== 3){
      plant_return[i,] = plant_weighting(plants_mat[i,], w1 = w2, w2 = w3, w3 = w4, w4 = w4)
    } else if(plants_max[i]== 2){
      plant_return[i,] = plant_weighting(plants_mat[i,], w1 = w3, w2 = w4, w3 = w4, w4 = w4)
    } else{
      plant_return[i,] = plant_weighting(plants_mat[i,])
    }
   }
  
  plants_rel <- decostand(plant_return, method = "total")
  plants_dist <- vegdist(plants_rel, method = method)
  return(plants_dist)
}


# PAM ---------------------------------------------------------------------
# I want to try PAM instead of agglomerative cluster algorithms because in the end I'm interested in assigning new plots to 
# existing clusters, therefore, the distance to the these clusters must be evaluated, therefore, I think medoids are better suited than 
# hierarchy trees

# https://www.statology.org/k-medoids-in-r/

pam_wAG <- weight_rel_dist(data_list$forest_wone_wAG$plants)

#fviz_nbclust(pam_wAG, FUNcluster = pam, method = "silhouette", diss = TRUE)
# as the function as used above does not work, I#m using the following workaround:

ks <- 2:25

sil_width <- sapply(ks, function(k) {
  pam_fit <- pam(pam_wAG, k, diss = TRUE)
  pam_fit$silinfo$avg.width
})
plot(ks, sil_width, type = "b", xlab = "k", ylab = "Avg silhouette")

pam_model <- pam(pam_wAG, 7, diss = TRUE)
pam_model
plot(pam_model)

pca_wAG <- prcomp(pam_wAG)
plot(pca_wAG$x[,1:2],
     col = pam_model$clustering,
     pch = 16,
     xlab = "PC1", ylab = "PC2")

pam_ord_wAG <- metaMDS(pam_wAG, k = 2)

plot(pam_ord_wAG$points,
     col = pam_model$clustering,
     pch = 16)

pam_model_help <- as.data.frame(pam_model$clustering)
names(pam_model_help) <- c("cluster")
pam_eval <- hdbscan_mismatch_evaluation(plants_nmds = nmds_object_imbalanced_wAG[[1]][[2]],
                            plants_hdbscan = pam_model_help,
                            grunddat = wone_wAG,
                            coarse = TRUE)

sum(pam_eval$mismatch) # not fitting into the cluster
sum(pam_eval$mismatch==FALSE)

adjustedRandIndex(pam_model$clustering, wone_wAG$BT_Land_group)
cl_agreement(
  as.cl_partition(pam_model$clustering),
  as.cl_partition(wone_wAG$BT_Land_group),
  method = "purity"
)

# scrutinize distribution of maximum abundance in plant data set ----------

plants_wAG_max <- apply(data_list$forest_wone_wAG$plants[,-1], 1, max)
table(as.numeric(plants_wAG_max)) # I need to check the 1 and 2s

polygons_few_trees <- data_list$forest_wone_wAG$plants$Polygon[order(plants_wAG_max)]
filter(grunddaten, Polygon %in% polygons_few_trees[1:36])
plants_helper <- inner_join(wone_wAG[,c(1,2,4,5,24,25)],data_list$forest_wone_wAG$plants, by = "Polygon")
plants_helper[plants_helper == 0] = NA
check_polygons <- filter(pivot_longer(plants_helper, cols = where(is.numeric)& !Polygon ,names_to = "plant", values_to = "abundance", values_drop_na = TRUE), 
       Polygon %in% polygons_few_trees[1:36])

# check highest abundance in only trees
trees_wzero_max <- apply(plants_trees_wzero$plants[,-1], 1, max)
table(as.numeric(trees_wzero_max)) # I need to check the 1 and 2s

polygons_few_trees_wzero <- plants_trees_wzero$plants$Polygon[order(trees_wzero_max)]
polygons_few_trees_wzero_codes <- filter(grunddaten, Polygon %in% polygons_few_trees_wzero[1:311])
# maybe exclude not only AGs, but also AVs and based on bund: 43.09.x, 39.02.x

trees_wzero_wAG <- dplyr::filter(plants_trees_wzero[[2]], !BT_Land_group %in% c("AG", "AV")) # maybe &!BT_Land_group %in% c("43.09", "39.02"))
plants_wone_wAG <- dplyr::filter(plants_trees_wzero[[1]], Polygon %in% trees_wzero_wAG$Polygon)
wone_wAGAV_list <- list(plants = plants_wone_wAG, grunddaten = trees_wzero_wAG)

plant_mat <- plant_weighting(wone_wAGAV_list$plants,w1 = 0.01,w2 = 0.02,w3 = 0.1,w4 = 0.95)
plant_mat <- decostand(plant_mat, method = "total")

wone_wAGAV_stress <- list()
wone_wAGAV_stress[[1]] <- future_sapply(2:6, function(k){
  median(replicate(2, metaMDS(plant_mat, k = k, trymax = 3, trace = FALSE)$stress))
}, future.seed = TRUE)
plot(x = 2:6, y = wone_wAGAV_stress[[1]])

wone_wAGAV_NMDS <- metaMDS(plant_mat, k = 2, trymax = 20)

wone_wAGAV_list$grunddaten[c(ordination_outlier_func(wone_wAGAV_NMDS)),]
control <- wone_wAGAV_list$plants[c(ordination_outlier_func(wone_wAGAV_NMDS)),]

hdbscan_complete(wone_wAGAV_NMDS$points, grunddat = wone_wAGAV_list$grunddaten,coarse = TRUE, bund = FALSE, by = 1)

# collapsing to genus level
plant_wAGAV_genus <- wone_wAGAV_list$plants %>%
  pivot_longer(-Polygon, names_to = "species", values_to = "abundance") %>%
  mutate(genus = sub("(\\w+).*", "\\1", species)) %>%
  group_by(Polygon, genus) %>%
  summarise(abundance = max(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = genus, values_from = abundance, values_fill = 0)

plant_wAGAV_genus_dist <- weight_rel_dist(plant_wAGAV_genus)
hdbscan_complete(plant_wAGAV_genus_dist, grunddat = wone_wAGAV_list$grunddaten,coarse = TRUE, bund = FALSE, by = 1)
clusterVScode(plant_wAGAV_genus_dist, grunddat = wone_wAGAV_list$grunddaten,bund = FALSE, pts = 9) #coarse = TRUE, 

plant_wAGAV_genus_hdbscan <- hdbscan(plant_wAGAV_genus_dist,minPts = 9)
cluster_control <- cbind(plant_wAGAV_genus_hdbscan$cluster,wone_wAGAV_list$grunddaten[,c(1,2,4,5,24,25)], plant_wAGAV_genus)
# obviously, sometimes "Hainbuchen-Eichenmischwald" and "Eichen-Hainbuchen" are mixed-up

plant_wAGAV_genus_dist_balanced <- weight_rel_dist(plant_wAGAV_genus,w1 = 0.01,w2 = 0.02,w3 = 0.1,w4 = 0.95)
hdbscan_complete(plant_wAGAV_genus_dist_balanced, grunddat = wone_wAGAV_list$grunddaten,coarse = TRUE, bund = FALSE, by = 1)
plant_wAGAV_genus_hdbscan <- hdbscan(plant_wAGAV_genus_dist2,minPts = 9)
cluster_control <- cbind(plant_wAGAV_genus_hdbscan$cluster,wone_wAGAV_list$grunddaten[,c(1,2,4,5,24,25)], plant_wAGAV_genus)

# maybe two cluster runs, one to get the distinct plots (the ones with a max abundance of 4) and one to seperate the remaining, in this 
# case cluster number 13

# check/prove hdbscan without NMDS performance ----------------------------

plants_weight_wAG <- weight_rel_dist(data_list$forest_wone_wAG$plants)
hdbscan_wAG_plants <- hdbscan_complete(plants_weight_wAG,grunddat = wone_wAG, by = 1)
hdbscan_plot(hdbscan_wAG_plants, name = "HDBSCAN without NMDS")

hdbscan_wAG_plants_wdif <- weight_rel_dist(data_list$forest_wone_wAG$plants, w1 = 0.01, w2 = 0.02, w3 = 0.05, w4 = 0.9)
hdbscan_wAG_plants2 <- hdbscan_complete(hdbscan_wAG_plants_wdif,grunddat = wone_wAG, by = 1)
hdbscan_plot(hdbscan_wAG_plants2, name = "HDBSCAN without NMDS")


# further analysis --------------------------------------------------------


# reveal indicator species
library(indicspecies)

plant_dist_second <- max_weigthing(plants_mat=plants_second_round)
hdbscan_complete(plant_dist_second, grunddat = grunddat_second_round,coarse = TRUE, bund = FALSE, kstop = 13)

# maybe two cluster runs, one to get the distinct plots (the ones with a max abundance of 4) and one to seperate the remaining, in this 
# case cluster number 13

hdbscan_complete(plants_dist = max_weigthing(data_list[["forest_without_AG&plots_only_one"]][[1]]),
                 grunddat = data_list[["forest_without_AG&plots_only_one"]][[2]], coarse = TRUE, bund = FALSE, by = 1)

