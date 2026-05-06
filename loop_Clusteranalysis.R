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
library(uwot)
library(rnndescent)
library(ggrepel)
library(RColorBrewer)
library(colorspace)
library(gridExtra)
library(paletteer)

plan(multisession, workers = parallel::detectCores() - 2)

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

# # further reduction of forest data ---------------------------------------------
# 
# # remove "AT... and "AU..."
# check_biotope_codes <- unique(grunddaten_forests[,c(4,5)])
# grunddaten_forests_red <- filter(grunddaten_forests, !substr(`Biotoptyp-Land`,1,2) %in% c("AT", "AU", "AV"))
# plants_forests_red <- filter(plants_forests, Polygon %in% grunddaten_forests_red$Polygon)
# 
# ### store plants and grunddaten in list for looping later, each list entry contains two list entries with first plants, then grunddaten
# 
# data_list <- list()
# data_list[["forest_complete"]] <-  list(plants = plants_forests_red, grunddaten = grunddaten_forests_red)
# 
# # remove all polygons which only have a few plant entries
# 
# forest_complete_wone <- remove_plots(plants_forests_red, grunddaten_forests_red, remove_count = c(0,1))
# check_empty_plot(forest_complete_wone[[1]])
# 
# data_list[["forest_complete_c(0,1)"]] <- forest_complete_wone
# 
# forest_complete_wtwo <- remove_plots(plants_forests_red, grunddaten_forests_red, remove_count = c(0:2))
# data_list[["forest_complete__c(0:2)"]] <- forest_complete_wtwo
# 
# # add life forms
# plants_occ <- plant_occurences(plant_data = plants_forests_red)
# 
# plants_occ$short <- sub("(\\w+\\s+\\w+).*", "\\1", plants_occ$species)
# plants_occ$short <- sub("(\\w+).*", "\\1", plants_occ$species)
# life_forms$short <- sub("(\\w+).*", "\\1", life_forms$FloraVeg.Taxon)
# 
# plants_LF <- left_join(plants_occ,life_forms, by = "short",multiple = "any")
# 
# trees <- filter(plants_LF, Tree == 1)
# trees <- filter(trees, species %in% colnames(plants_forests))
# 
# plants_forest_trees <- plants_forests_red[c("Polygon",trees$species)]
# 
# plants_occ_forests <- plant_occurences(plants_forest_trees)
# 
# # remove plants which only occured once or twice across all plots
# plants_f_trees_w_zero <- filter(plants_occ_forests, !total %in% c(0)) # maybe also 2
# plants_f_trees_wzero <- plants_forests_red[c("Polygon",plants_f_trees_w_zero$species)] # take all 
# 
# plants_f_trees_w_one <- filter(plants_occ_forests, !total %in% c(0,1)) # maybe also 2
# plants_f_trees_wone <- plants_forests_red[c("Polygon",plants_f_trees_w_one$species)] # take all 
# # I'm using the plants_f_trees_wone for further reduction of the data
# 
# ### create an additional plant list simplified to genus-level
# 
# # collapsing to genus level
# plant_forest_genus <- plants_forests_red %>%
#   pivot_longer(-Polygon, names_to = "species", values_to = "abundance") %>%
#   mutate(genus = sub("(\\w+).*", "\\1", species)) %>%
#   group_by(Polygon, genus) %>%
#   summarise(abundance = max(abundance, na.rm = TRUE), .groups = "drop") %>%
#   pivot_wider(names_from = genus, values_from = abundance, values_fill = 0)
# 
# forest_complete_wone_genus <- remove_plots(plant_forest_genus, grunddaten_forests_red, remove_count = c(0,1))
# check_empty_plot(forest_complete_wone_genus[[1]])
# 
# data_list[["forest_genus_c(0,1)"]] <- forest_complete_wone_genus
# 
# forest_complete_wtwo_genus <- remove_plots(plant_forest_genus, grunddaten_forests_red, remove_count = c(0:2))
# data_list[["forest_genus_c(0,2)"]] <- forest_complete_wtwo_genus
# 
# 
# # control number of tree species per plot -------------------------------------------------
# 
# check_empty_plot(plants_f_trees_wone)
# plants_trees_wzero <- remove_plots(plants_f_trees_wone, grunddaten_forests_red)
# check_empty_plot(plants_trees_wzero[[1]])
# 
# data_list[["trees_c(0)"]] <- plants_trees_wzero
# 
# plants_trees_wone <- remove_plots(plants_f_trees_wone, grunddaten_forests_red, remove_count = c(0,1))
# check_empty_plot(plants_trees_wone[[1]])
# data_list[["trees_c(0,1)"]] <- plants_trees_wone
# 
# plants_trees_wtwo <- remove_plots(plants_f_trees_wone, grunddaten_forests_red, remove_count = c(0:2))
# check_empty_plot(plants_trees_wtwo[[1]])
# data_list[["trees_c(0:2)"]] <- plants_trees_wtwo

# # remove AG-biotope types -------------------------------------------------
# 
# data_list[["forest_AG_c(0,1)"]] <- remove_land_biotope_code(forest_complete_wone[[1]], forest_complete_wone[[2]],c("AG"))
# data_list[["trees_AG&_c(0)"]] <- remove_land_biotope_code(plants_trees_wzero[[1]], plants_trees_wzero[[2]],c("AG"))
# data_list[["trees_AG_c(0,1)"]] <- remove_land_biotope_code(plants_trees_wone[[1]], plants_trees_wone[[2]],c("AG"))
# 
# data_list[["trees_genus_AG_c(0,1)"]] <- remove_land_biotope_code(forest_complete_wone_genus[[1]], forest_complete_wone_genus[[2]],c("AG"))
# 

data_list <- readRDS("data_list.RDS")

# check direct hdbscan application ----------------------------------------

hdbscan_list <- list()
for(i in 1:length(data_list)){
  weighting_plants <- weight_rel_dist(data_list[[i]]$plants)
  hdbscan_run <- lapply(3:20, function(k) hdbscan(weighting_plants, minPts = k))
  hdbscan_list[[names(data_list)[i]]][[1]] <- hdbscan_run
  weighting_plants <- weight_rel_dist(data_list[[i]]$plants, w1 = 0.01, w2 = 0.05, w3 = 0.25, w4 = 0.9)
  hdbscan_run <- lapply(3:20, function(k) hdbscan(weighting_plants, minPts = k))
  hdbscan_list[[names(data_list)[i]]][[2]] <- hdbscan_run
}

hdbscan_metrics_all_imbalanced <- list()
hdbscan_metrics_all_balanced <- list()
for(i in 1:12){
  hdbscan_metrics_all_imbalanced[[names(data_list)[i]]] <- hdbscan_metrics(hdbscan_list[[i]][[1]], grunddat = data_list[[i]][[2]])
  hdbscan_metrics_all_balanced[[names(data_list)[i]]] <- hdbscan_metrics(hdbscan_list[[i]][[2]], grunddat = data_list[[i]][[2]])
}


hdbscan_metrics_all_imbalanced_df <- combine_list_to_df(hdbscan_metrics_all_imbalanced)
hdbscan_metrics_all_balanced_df <- combine_list_to_df(hdbscan_metrics_all_balanced)

check_bund_performance <- hdbscan_metrics_all_imbalanced_df %>%
  filter(`Biotoptyp-Bund_ari` > 0.3 &(`Biotoptyp-Bund_ari`>`Biotoptyp-Land_ari`))

hdbscan_metrics_all_balanced_df %>%
  filter(`Biotoptyp-Bund_ari` > 0.3 &(`Biotoptyp-Bund_ari`>`Biotoptyp-Land_ari`))

hdbscan_metrics_all_imbalanced_df$balance <- "wrong"
hdbscan_metrics_all_balanced_df$balance <- "true"

hdbscan_metrics_all <- rbind(hdbscan_metrics_all_imbalanced_df, hdbscan_metrics_all_balanced_df)

hdbscan_metrics_all <- hdbscan_metrics_all %>%
  mutate(composed_metric_land = (0.5*`Biotoptyp-Land_ari`+ 0.4*`Biotoptyp-Land_purity`+ 0.1*(1-noise_prop)),
         composed_metric_land_coarse = (0.5*`BT_Land_group_ari`+ 0.4*`BT_Land_group_purity`+ 0.1*(1-noise_prop)),
         composed_metric_land_equal = (0.3*`Biotoptyp-Land_ari`+ 0.3*`Biotoptyp-Land_purity`+ 0.3*(1-noise_prop)),
         composed_metric_land_coarse_equal = (0.3*`BT_Land_group_ari`+ 0.3*`BT_Land_group_purity`+ 0.3*(1-noise_prop)))

# bund less performative!




# ARI can be higher on code level than on group level since ARI penalizes oversplitting. In my case, purity might be more relveant???

### check later:
#check best result using weighting 2 with new weighting function and respective weights
hdbscan_complete(max_weighting(data_list[["trees_AG&_c(0)"]][[1]], w1 = 0.01, w2 = 0.05, w3 = 0.25, w4 = 0.9),by = 1, grunddat = data_list[["trees_AG&_c(0)"]][[2]], 
                 bund = FALSE, coarse = FALSE, print = FALSE)

###
# names_plots <- c("Land_ARI_most", "Land_group_ARI_most", "Land_equal_weighting", "Land_group_equal_weighing")
# for(i in 1:4){
#   col_name <- names(hdbscan_metrics_all)[i + 22]
#   png(paste0("results/",names_plots[i],".png"), width = 3000, height = 2100, res = 300)
#   p <- ggplot(hdbscan_metrics_all)+
#     geom_line(aes(x=k,y= .data[[names(hdbscan_metrics_all)[i + 22]]], colour = list_name))+
#     theme_classic()+
#     scale_color_discrete(name = "Plant data composition")+
#     labs(title = names_plots[i], y = "performance (ARI + purity + noise)", x = "minimum size of cluster")+
#     facet_wrap(facets = "balance", labeller = labeller(balance=c("wrong"= "imbalanced conversion", "true"= "more balanced conversion")))
#   print(p)
#   dev.off()
#   if(i < 3){
#     names_purity <- c("Biotoptyp-Land_purity","BT_Land_group_purity")
#     png(paste0("results/",names_purity[i],".png"), width = 3000, height = 2100, res = 300)
#     p <- ggplot(hdbscan_metrics_all)+
#       geom_line(aes(x=k,y= .data[[names_purity[i]]], colour = list_name))+
#       theme_classic()+
#       scale_color_discrete(name = "Plant data composition")+
#       labs(title = names_purity[i], y = "performance (purity)", x = "minimum size of cluster")+
#       facet_wrap(facets = "balance", labeller = labeller(balance=c("wrong"= "imbalanced conversion", "true"= "more balanced conversion")))
#     print(p)
#     dev.off()
#   }
#   
# }

control <- hdbscan_metrics_all[,-c(5:8,11:16.19,20)] %>%
  group_by(list_name, balance) %>%
  slice_max(order_by = composed_metric_land, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(composed_metric_land))



hdbscan_metrics_all[,-c(5:8,11:16.19,20)] %>%
  group_by(list_name, balance) %>%
  slice_max(order_by = composed_metric_land_coarse, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(composed_metric_land_coarse))


## selection

# scenarios to keep for comparison between clustering algorithms:
scenario_names <- c("forest_AG_c(0,1)","trees_genus_AG_c(0,1)", "trees_AG&_c(0)", "trees_AG_c(0,1)", "forest_complete__c(0:2)")
scenario_balance <- c(1,1,2,2,1) # the both tree data sets without AGs perform well with the balanced transformation,
# especially for purity at land level!
scenario_df <- data.frame(cbind(scenario_names, as.numeric(scenario_balance), c(rep("wrong",2),rep("true",2), "wrong")))
names(scenario_df) <- c("list_name", "balance_scenario", "balance")
scenario_df$balance_scenario <-  as.numeric(scenario_balance)

data_list_comp <- list()
data_list_comp <- data_list[names(data_list) %in% c(scenario_names)]

ggplot(hdbscan_metrics_all[hdbscan_metrics_all$list_name %in% scenario_names,])+
  geom_line(aes(x=k,y= composed_metric_land_coarse, colour = list_name))+
  theme_classic()+
  scale_color_discrete(name = "Plant data composition")+
  labs(title = "equal weighting for performance",y = "performance (ARI + purity + noise)", x = "minimum size of cluster")+
  facet_wrap(facets = "balance", labeller = labeller(balance=c( "true"= "more balanced conversion","wrong"= "imbalanced conversion")))

best_metrics_hdbscan <- hdbscan_metrics_all[,-c(5:8,11:16.19,20)] %>%
  semi_join(scenario_df, by = c("list_name", "balance"))%>%
  group_by(list_name, balance) %>%
  filter( rank(-composed_metric_land_coarse_equal) <= 2 |
      rank(-`Biotoptyp-Land_purity`) <= 2)
  
scenario_df$best_k <- c(13,7,6,5,20)  

# hdbscan_complete(plants_dist = weight_rel_dist(data_list[["forest_AG_c(0,1)"]][[1]]),
#                  grunddat = data_list[["forest_AG_c(0,1)"]][[2]], coarse = TRUE, bund = FALSE, by = 1)
# 
# clusterVScode(plants_dist = weight_rel_dist(data_list[["forest_AG_c(0,1)"]][[1]]),
#               grunddat = data_list[["forest_AG_c(0,1)"]][[2]], pts = 11, bund = FALSE)
# 
# ### temporarily visualisation
# 
# hdbscan_plot <- hdbscan(weight_rel_dist(data_list[["forest_AG_c(0,1)"]][[1]]), minPts = 11)
# 
# # project multidemensional data into 2-D
# plant_mat <-plant_weighting(data_list[["forest_AG_c(0,1)"]][[1]])
# plant_mat <- decostand(plant_mat, method = "total")
# pca_plants <- prcomp(plant_mat)
# 
# # two-dimensional
# plot_data <- as.data.frame(pca_plants$x[, 1:2])
# plot_data$cluster <- factor(hdbscan_plot$cluster)
# 
# ggplot(plot_data, aes(x = PC1, y = PC2, color = cluster)) +
#   geom_point(size = 2, alpha = 0.8) +
#   #scale_color_manual(values = c("0" = "grey70")) +
#   theme_minimal()
# 
# plot_data$biotope <- data_list[["forest_AG_c(0,1)"]][[2]]$`Biotoptyp-Land`
# ggplot(plot_data, aes(PC1, PC2, color = biotope)) +
#   geom_point(size = 2) +
#   theme_minimal()
# 
# # 3-dimensional; wrong naming to use already coded function!
# pca_3d <- list()
# pca_3d[["points"]] <- pca_plants$x
# hdbscan_eval <- hdbscan_mismatch_evaluation(pca_3d,hdbscan_plot,data_list[["forest_AG_c(0,1)"]][[2]], coarse = TRUE)
# sum(hdbscan_eval$mismatch==TRUE)
# sum(hdbscan_eval$cluster==0)
# 
# hover_3D(hdbscan_eval)
# hull_3D(hdbscan_eval, op_hull = 0.6, op_points = 0.3)
# 
# ######
# 
# selection_scenarios <- hdbscan_direct_result_coarse %>%
#   filter(k > 5) %>%
#   group_by(list_name, sublist_id) %>%
#   summarise(mean = mean(composed_metric)) %>%
#   arrange(desc(mean))



# comparison of algorithms ------------------------------------------------

plant_comp_dist <- list()
for(i in 1: length(scenario_names)){
  if(scenario_balance[i]== 2){
    dist_mat <- weight_rel_dist(data_list_comp[[i]][[1]],w1 = 0.01, w2 = 0.05, w3 = 0.25, w4 = 0.9)
  } else{
    dist_mat <- weight_rel_dist(data_list_comp[[i]][[1]],)
  }
  plant_comp_dist[[names(data_list_comp)[i]]] <- dist_mat
}
# first PAM
ks <- 2:25


# https://j-sephb-lt-n.github.io/exploring_statistics/PAMS_and_SILHOUETTE_by_hand.html

pam_results <- future_lapply(names(plant_comp_dist), function(name) {
  dist_mat <- plant_comp_dist[[name]]
  res <- lapply(ks, function(k) {
    
    pam_fit <- cluster::pam(dist_mat, k, diss = TRUE)
    
    list(
      k = k,
      sil_width = pam_fit$silinfo$avg.width,
      clustering = pam_fit$clustering
    )
  })
  
  names(res) <- ks
  res
}, future.seed = TRUE)

names(pam_results) <- names(plant_comp_dist)


# evaluate with metrics
pam_evaluation <- evaluate_pam_models(pam_results, data_list_comp)
pam_evaluation_plot <- pivot_longer(pam_evaluation, cols = -c(dataset, k), names_to = "metric", values_to = "value")

# png("results/PAM_evaluation.png", width = 3000, height = 2100, res = 350)
ggplot(pam_evaluation_plot)+
  geom_line(aes(x = k, y = value, colour = dataset))+
  theme_classic()+
  labs(title = "PAM evaluation")+
  facet_wrap(facets = "metric")
# dev.off()


pam_sil_df <- purrr::imap_dfr(pam_results, function(res_list, dataset_name) {
  purrr::map_dfr(res_list, function(x) {
    data.frame(
      dataset = dataset_name,
      k = x$k,
      sil_width = x$sil_width
    )
  })
})
pam_best <- pam_sil_df %>%
  group_by(dataset)%>%
  slice_max(order_by = sil_width, n = 1, with_ties = FALSE) %>%
  ungroup()

ggplot(pam_sil_df, aes(x = k, y = sil_width, col = dataset))+
  geom_line()+
  theme_classic()+
  labs(y = "Avg silhouette")

pam_best_models <- lapply(pam_best$dataset, function(x){
  dist_mat <- plant_comp_dist[[x]]
  k = pam_best$k[pam_best$dataset== x]
  pam_model <- cluster::pam(dist_mat, k, diss = TRUE)
})

names(pam_best_models) <- pam_best$dataset

#save.image("all_inclusive_PAM_0305.RData")

# GMM ---------------------------------------------------------------------

plant_comp_gmm <- list()
plant_comp_gmm_bray <- list()

for(i in 1: length(scenario_names)){
  if(scenario_balance[i]== 2){
    dist_mat <- plant_weighting(data_list_comp[[i]][[1]],w1 = 0.01, w2 = 0.05, w3 = 0.25, w4 = 0.9)
  } else{
    dist_mat <- plant_weighting(data_list_comp[[i]][[1]],)
  }
  rel_mat <- decostand(dist_mat, method = "total")
  plant_comp_gmm[[names(data_list_comp)[i]]] <- dist_mat
  rel_mat <- decostand(dist_mat, method = "hellinger")
  plant_comp_gmm_bray[[names(data_list_comp)[i]]] <- dist_mat
}

gmm_results <- future_lapply(names(plant_comp_gmm), function(name) {
  rel_mat <- plant_comp_gmm[[name]]
  gmm_fit <- Mclust(rel_mat)
    }, future.seed = TRUE)

names(gmm_results) <- names(plant_comp_dist)


# visualize/evaluate results

gmm_total_eval <- evaluate_gmm(gmm_results, data_list_comp)
gmm_total_eval$metrics[,c(1,2,9,10)]

gmm_evaluation_plot <- pivot_longer(gmm_total_eval$metrics, cols = -c(dataset,cluster, bic, uncertainty), names_to = "metric", values_to = "value")

# png("results/GMM_evaluation_total.png", width = 3000, height = 2100, res = 350)
ggplot(gmm_evaluation_plot)+
  geom_path(aes(x = factor(dataset, levels = names(plant_comp_dist)), y = value, colour = metric, group = metric))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(title = "GMM evaluation-total", x = "dataset")
#dev.off()

#################



# # NMDS loop ---------------# NMDS loop ---------------# NMDS loop ---------------------------------------------------------------
# weighting <- data.frame(imbalanced = c(0.01,0.01,0.01,1),
#                         balanced = c(0.01,0.05,0.25,0.9))
# #stress_vals_backup <- stress_vals
# stress_vals <- list(imbalanced = list(), balanced = list())
# 
# for (o in scenario_names) { # after the 12 first scenarios, there appears to be a gap
#     plants <- data_list[[o]][[1]]
#     plant_mat <- plant_weighting(plants,weighting[1,i],weighting[2,i],weighting[3,i],weighting[4,i])
#     plant_mat <- decostand(plant_mat, method = "total") # first relative abundance; maybe for comparison hellinger transformation as well
#     # find best dimensionality: via NMDS slowly...
#     stress_vals[[i]][[o]] <- future_sapply(2:6, function(k){
#       median(replicate(2, metaMDS(plant_mat, k = k, trymax = 3, trace = FALSE)$stress))
#     }, future.seed = TRUE)
#   
# }
# 
# #saveRDS(list(data_list,stress_vals), "results/data_with_stress_vals.RDS")
# #stress_vals <- readRDS("results/data_with_stress_vals.RDS")[[2]]
# 
# par(mfrow = c(2,4))
# for(i in stress_vals){
#   sapply(i, function(x){
#     plot(2:6, x, type = "b", main = names(i))
#     # approximate "elbow"
#     diff1 <- diff(x)
#     diff2 <- diff(diff1)
#     
#     k_opt <- which.min(diff2) + 1
#     print(k_opt)
#   })
# }
# 
# stress_vals_df <- data.frame(data = rep(unique(selection_scenarios$list_name[1:12]),each = 5),
#                              dimensions = rep(seq(2,6),7),
#                              imbalanced = rep(0,length(rep(unique(selection_scenarios$list_name[1:12]),each = 5))),
#                              balanced = rep(0, length(rep(unique(selection_scenarios$list_name[1:12]),each = 5))))
# 
# 
# for(i in 1:2){
#   results_stress <- c()
#   for (o in 1:7){
#     
#     results_stress <- c(results_stress,stress_vals[[i]][[o]])
#   }
#   stress_vals_df[,i+2] <- results_stress
#   results_stress <- c()
# }
# 
# stress_vals_df_plot <- pivot_longer(stress_vals_df, cols = imbalanced:balanced,
#                                     values_to = "stress", names_to = "balance")
# 
# ggplot(stress_vals_df_plot)+
#   geom_point(aes(x = dimensions, y = stress, colour = balance))+
#   theme_minimal()+
#   facet_wrap(facets = "data")
# 
# 
# 
# ###continue
# 
# plant_mat <- plant_weighting(data_list[["forest_AG_c(0,1)"]][[1]])
# plant_mat <- decostand(plant_mat, method = "total") # first relative abundance; maybe for comparison hellinger transformation as well
# #nmds_5_wag_wone <- metaMDS(plant_mat, k = 5, trymax = 20)
# 
# hdbscan_complete(nmds_5_wag_wone$points, by = 1, data_list[["forest_AG_c(0,1)"]][[2]], bund = FALSE, coarse = TRUE)
# hdbscan_complete(weight_rel_dist(data_list[["forest_AG_c(0,1)"]][[1]]),
#                  by = 1, data_list[["forest_AG_c(0,1)"]][[2]], bund = FALSE, coarse = TRUE)
# 
# #visualisation_3_nmds <- metaMDS(plant_mat, k = 3, trymax = 20)
# hdbscan_complete(visualisation_3_nmds$points, by = 1, data_list[["forest_AG_c(0,1)"]][[2]], bund = FALSE, coarse = TRUE)
# 
# #visualisation_2_nmds <- metaMDS(plant_mat, k = 2, trymax = 20)
# 
# nmds_wag_wone <- list(nmds_5 = nmds_5_wag_wone$points, nmds_3 = visualisation_3_nmds$points, nmds_2 = visualisation_2_nmds$points,
#                       direct = weight_rel_dist(data_list[["forest_AG_c(0,1)"]][[1]]))
# saveRDS(nmds_wag_wone, "results/forest_wAG_c(0,1)_nmds.RDS")
# nmds_wag_wone_metrics <- list()
# for (i in 1:4) {
#   evaluation_coarse <- hdbscan_complete(nmds_wag_wone[[i]],by = 1, data_list[["forest_AG_c(0,1)"]][[2]], bund = FALSE, coarse = TRUE, print = FALSE)
#   evaluation <- hdbscan_complete(nmds_wag_wone[[i]],by = 1, data_list[["forest_AG_c(0,1)"]][[2]], bund = FALSE, print = FALSE)
#   nmds_wag_wone_metrics[[names(nmds_wag_wone)[i]]][[1]] <- evaluation_coarse
#   nmds_wag_wone_metrics[[names(nmds_wag_wone)[i]]][[2]] <- evaluation
# }
# helper_list <- list()
# for(i in 1:4){
#   helper_list[[names(nmds_wag_wone)[i]]][[1]] <- data_list[["forest_AG_c(0,1)"]][[1]]
# }
# nmds_wag_wone_compare <- hdbscan_result_df(nmds_wag_wone_metrics, main_list = helper_list)
# 
# nmds_wag_wone_compare$composed_metric <- (0.5*nmds_wag_wone_compare$ari+
#                                                    0.4*nmds_wag_wone_compare$purity+0.1*(1-nmds_wag_wone_compare$noise_pro))
# nmds_wag_wone_compare$composed_metric_balanced <- (0.3*nmds_wag_wone_compare$ari+
#                                                             0.3*nmds_wag_wone_compare$purity+0.3*(1-nmds_wag_wone_compare$noise_pro))
# 
# #png("results/nmds_evaluation_hdbscan_coarse_noise.png", width = 3000, height = 2100, res = 300)
# ggplot(nmds_wag_wone_compare)+
#   geom_line(aes(x=k,y= composed_metric_balanced, colour = list_name))+
#   theme_classic()+
#   scale_color_discrete(name = "Plant data composition")+
#   scale_y_continuous(limits = c(0.3, 0.8), breaks = c(0.40,0.60,0.8, 1))+
#   labs(title = "forest withoutAG & c(0,1) - equal weighting",y = "performance (ARI + purity + noise)", x = "minimum size of cluster")+
#   facet_wrap(facets = "sublist_id", labeller = labeller(sublist_id=c("1"= "coarse evaluation", "2"= "code-level evaluation")))
# #dev.off()
# 
# #png("results/nmds_evaluation_hdbscan_coarse.png", width = 3000, height = 2100, res = 300)
# ggplot(nmds_wag_wone_compare)+
#   geom_line(aes(x=k,y= composed_metric, colour = list_name))+
#   theme_classic()+
#   scale_color_discrete(name = "Plant data composition")+
#   scale_y_continuous(limits = c(0.15, 1), breaks = c(0.20,0.40,0.60,0.8, 1))+
#   labs(title = "forest withoutAG & c(0,1) - noise less weighted" ,y = "performance (ARI + purity + noise)", x = "minimum size of cluster")+
#   facet_wrap(facets = "sublist_id", labeller = labeller(sublist_id=c("1"= "coarse evaluation", "2"= "code-level evaluation")))
# #dev.off()
# ### continue
# 
# # at the code level, the nmds dimensionality reduction can have advantages, especially if the aim is to reduce the number of points 
# # classified as noise
# 
# # data set 1 with minPts 13 and bt_land
# 
# #check whether ordination posses outliers
# for(i in 1:3){
#   outliers <- #lapply(nmds_wag_wone[[i]], function(x){
#     ordination_outlier_func(nmds_wag_wone[[i]])
#   #})
#   print(outliers)
# }
# 
# 
# # display just first two NMDS axes
# plot_data_2d <- as.data.frame(visualisation_2_nmds$points)
# plot_data_2d$cluster <- factor(hdbscan_plot$cluster)
# 
# plot(plot_data_2d)
# 
# ggplot(plot_data_2d, aes(x = MDS1, y = MDS2, color = cluster)) +
#   geom_point(size = 2, alpha = 0.8) +
#   #scale_color_manual(values = c("0" = "grey70")) +
#   theme_minimal()
# 
# plot_data_2d$biotope <- data_list[["forest_AG_c(0,1)"]][[2]]$BT_Land_group
# ggplot(plot_data_2d, aes(x = MDS1, y = MDS2, color = biotope)) +
#   geom_point(size = 2) +
#   theme_minimal()
# 
# ggplot(plot_data_2d, aes(x = MDS1, y = MDS2, color = cluster)) +
#   geom_point() +
#   stat_ellipse() +
#   theme_minimal()
# 
# grunddat_temp <- data_list[["forest_AG_c(0,1)"]][[2]]
# df_plot_3d <- hdbscan_mismatch_evaluation(visualisation_3_nmds, hdbscan_plot,
#                             grunddat = grunddat_temp, coarse = TRUE)
# 
# sum(df_plot_3d$mismatch) # not fitting into the cluster
# 
# sum(df_plot_3d$mismatch[!df_plot_3d$cluster %in% c(0,17)])
# sum(df_plot_3d$mismatch==FALSE)
# 
# hover_3D(df_plot_3d)
# hover_3D(df_plot_3d[!df_plot_3d$cluster %in% c(0,17),])
# hull_3D(df_plot_3d, op_hull = 0.6, op_points = 0.3)
# hull_3D(df_plot_3d[!df_plot_3d$cluster %in% c(0,17),], op_hull = 0.6, op_points = 0.3)
# 
# ####
# legend_labels <- unique(df_plot_3d[,c(4, 7)])
# legend_labels$combined <- paste0(legend_labels$cluster, ": ",legend_labels$cluster_main)
# legend_labels <- legend_labels[order(legend_labels$cluster),]
# 
# #png("results/hdbscan_tree_coarse.png", width = 3000, height = 2100, res = 300)
# plot(hdbscan_plot, show_flat = TRUE, main = "HDBSCAN tree")
# legend("topright",
#        legend = legend_labels$combined[-1],
#        ncol = 3,
#        title = "Main biotope code per cluster")
# # dev.off()
# 
# ####
# 
# ### test 5d nmds into 3d
# df_plot_3d5 <- hdbscan_mismatch_evaluation(nmds_5_wag_wone, hdbscan_plot,
#                                           grunddat = grunddat_temp, coarse = TRUE)
# hover_3D(df_plot_3d5)
# hull_3D(df_plot_3d5, op_hull = 0.6, op_points = 0.3)  # nicer hull
# 
# ####
# 
# table(hdbscan_plot$cluster, grunddat_temp$`Biotoptyp-Land`)
# table(hdbscan_plot$cluster, grunddat_temp$`BT_Land_group`)
# 
# check_cluster_pureness <- as.data.frame(prop.table(table(hdbscan_plot$cluster, grunddat_temp$`BT_Land_group`), margin=1))
# ggplot(data = check_cluster_pureness, aes(x = Var1, y = Freq, fill = Var2))+
#   geom_bar(stat = "identity")
#   
# prop.table(t(table(hdbscan_plot$cluster, grunddat_temp$`BT_Land_group`)), margin=1)
# 
# ggplot(as.data.frame(table(hdbscan_plot$cluster, grunddat_temp$`BT_Land_group`)),
#        aes(Var2, Freq, colour = Var1, fill = Var1)) +
#   geom_bar(stat = "identity")+
#   theme_minimal()
# 

# UMAP --------------------------------------------------------------------

# dist_mat <- plant_weighting(data_list[["forest_AG_c(0,1)"]][[1]])
# dist_mat <- decostand(dist_mat, method = "total")
# 
# embedding <- umap(
#   dist_mat,
#   metric = "braycurtis",
#   nn_method = "nndescent",
#   n_neighbors = 10,
#   n_components = 10,
#   min_dist = 0.15
# )
# 
# table(hdbscan_wAG_imbalanced$cluster, wone_wAG$`Biotoptyp-Land`)
# table(hdbscan_wAG_imbalanced$cluster, wone_wAG$`BT_Land_group`)
# 
# ggplot(data = umap_cluster_pureness, aes(x = Var1, y = Freq, fill = Var2))+
#   geom_bar(stat = "identity")
# 
# ggplot(as.data.frame(table(umap_hdbscan$cluster, grunddat_temp$`BT_Land_group`)),
#        aes(Var2, Freq, colour = Var1, fill = Var1)) +
#   geom_bar(stat = "identity")+
#   theme_minimal()

### umap with the right setting allows to decrease the noise fraction in return for lower ari


# UMAP visualisation ------------------------------------------------------

# firstly, forest_AG_c(0,1) since all three algorithms produced somehow usable results
visualisation_3plots <- c("forest_complete__c(0:2)", "forest_AG_c(0,1)", "trees_genus_AG_c(0,1)")

visualisation_3plots_umap <- list()
visualisation_3plots_umap <- future_lapply(names(plant_comp_gmm), function(name){
  plant_mat <- plant_comp_gmm[[name]]
  umap_2d <- umap(
    plant_mat,
    metric = "braycurtis",
    nn_method = "nndescent",
    n_neighbors = 8,
    n_components = 2,
    min_dist = 0.7
  )
  umap_2d_df <- as.data.frame(umap_2d)
  colnames(umap_2d_df) <- c("UMAP1", "UMAP2")
  return(umap_2d_df)
}, future.seed = TRUE)

names(visualisation_3plots_umap) <- names(plant_comp_gmm)

# UMAP - NMDS comparison --------------------------------------------------

dimesion_comparison <- data_list[["forest_AG_c(0,1)"]]
plant_mat <- plant_comp_gmm[["forest_AG_c(0,1)"]]
comparison_umap <- umap(
  plant_mat,
  metric = "braycurtis",
  nn_method = "nndescent",
  n_neighbors = 15,
  n_components = 5,
  min_dist = 1.5,
  seed = 42
)
hdbscan_complete(comparison_umap, by = 1, grunddat = dimesion_comparison[[2]],
                 bund = FALSE, print = FALSE)
hdbscan_umap <- hdbscan(comparison_umap, minPts = 8)

# stress_vals <- readRDS("results/data_with_stress_vals.RDS")

comparison_nmds <-  metaMDS(
  plant_mat,
  k = 5,
  trymax = 20
)
hdbscan_complete(comparison_nmds$points, by = 1, grunddat = dimesion_comparison[[2]], bund = FALSE)
hdbscan_nmds <- hdbscan(comparison_nmds$points, minPts = 8)

comparison_combined <- visualisation_3plots_umap[["forest_AG_c(0,1)"]]
comparison_combined <- cbind(comparison_combined,
                             umap = hdbscan_umap$cluster,
                             nmds = hdbscan_nmds$cluster,
                             NMDS1 = visualisation_3plots_nmds[["forest_AG_c(0,1)"]]$points[,1],
                             NMDS2 = visualisation_3plots_nmds[["forest_AG_c(0,1)"]]$points[,2],
                             label = dimesion_comparison[[2]]$`Biotoptyp-Land`,
                             label_group = dimesion_comparison[[2]]$BT_Land_group)
comparison_plot_long <- pivot_longer(comparison_combined, cols = -c("UMAP1", "UMAP2", "NMDS1",
                                                                    "NMDS2", "label", "label_group"), names_to = "algorithm", values_to = "cluster")
hull_combined_umap <- comparison_plot_long %>%
  group_by(algorithm, cluster) %>%   # cluster must exist in your data
  slice(chull(UMAP1, UMAP2))


png_name <- paste0("results/forest_AG_c(0,1)_comparison_umap_nmds_umap.png")
#png(png_name, width = 3000, height = 2100, res = 300)
g <- ggplot(comparison_plot_long, aes(x = UMAP1, y = UMAP2, colour = label_group)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_polygon(
    data = hull_combined_umap,
    aes(x = UMAP1, y = UMAP2, group = cluster),
    fill = NA,
    colour = "black",
    linewidth = 0.7
  ) +
  theme_classic()+
  labs(title = "forest_AG_c(0,1)")+
  facet_wrap(facets = "algorithm")
print(g)
#dev.off()

hull_combined_nmds <- comparison_plot_long %>%
  group_by(algorithm, cluster) %>%   # cluster must exist in your data
  slice(chull(NMDS1, NMDS2))

png_name <- paste0("results/forest_AG_c(0,1)_comparison_umap_nmds_ndms.png")
#png(png_name, width = 3000, height = 2100, res = 300)
g <- ggplot(comparison_plot_long, aes(x = NMDS1, y = NMDS2, colour = label_group)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_polygon(
    data = hull_combined_nmds,
    aes(x = NMDS1, y = NMDS2, group = cluster),
    fill = NA,
    colour = "black",
    linewidth = 0.7
  ) +
  theme_classic()+
  labs(title = "forest_AG_c(0,1)")+
  facet_wrap(facets = "algorithm")
print(g)
#dev.off()

# comparison clustering algorithms ----------------------------------------

umap_plot_list <- list()
umap_plot_list_short <- list()
for(i in visualisation_3plots){
  umap_plot <- visualisation_3plots_umap[[i]]
  umap_plot$label <- data_list_comp[[i]][[2]]$`Biotoptyp-Land`
  umap_plot$label_group <- data_list_comp[[i]][[2]]$BT_Land_group
  umap_plot$pam <- as.factor(pam_best_models[[i]]$clustering)
  umap_plot$gmm <- as.factor(gmm_results[[i]]$classification)
  umap_plot$hdbscan <- as.factor(hdbscan_list[[i]][[scenario_df[scenario_df[,1]== i,2]]][[(scenario_df[scenario_df[,1]== i,4])-2]]$cluster)
  umap_plot_list_short[[i]] <- umap_plot
  
  umap_plot_long <- pivot_longer(umap_plot, cols = -c("UMAP1", "UMAP2", "label", "label_group"), names_to = "algorithm", values_to = "cluster")
  #umap_plot_list[[i]] <- umap_plot_long
  png_name <- paste0("results/",substr(i,1,14),"_umap_plot.png") # unfortunately, I used ":" ...
  #png(png_name, width = 3000, height = 2100, res = 300)
  # g <- ggplot(umap_plot_long, aes(x = UMAP1, y = UMAP2, colour = cluster)) +
  #   geom_point(size = 2, alpha = 0.8) +
  #   theme_classic()+
  #   labs(title = i)+
  #   facet_wrap(facets = "algorithm")
  # print(g)
  #dev.off()
  
  hulls <- umap_plot_long %>%
    group_by(algorithm, cluster) %>%   # cluster must exist in your data
    slice(chull(UMAP1, UMAP2))
  
  png_name <- paste0("results/",substr(i,1,14),"_umap_plot_label.png") # unfortunately, I used ":" ...
  #png(png_name, width = 3000, height = 2100, res = 300)
  g <- ggplot(umap_plot_long,
              aes(x = UMAP1, y = UMAP2, colour = label_group)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_polygon(
      data = hulls,
      aes(x = UMAP1, y = UMAP2, group = cluster),
      fill = NA,
      colour = "black",
      linewidth = 0.7
    ) +
    theme_classic()+
    labs(title = i)+
    facet_wrap(facets = "algorithm")
  print(g)
  #dev.off()
}

##### comparison all 4
umap_plot_4 <- umap_plot_list_short[["forest_AG_c(0,1)"]]
umap_plot_4$hdbscan_nmds <- as.factor(hdbscan_nmds$cluster)
umap_plot_long_4 <- pivot_longer(umap_plot_4, cols = -c("UMAP1", "UMAP2", "label", "label_group"), names_to = "algorithm", values_to = "cluster")

hulls <- umap_plot_long_4 %>%
  group_by(algorithm, cluster) %>%   # cluster must exist in your data
  slice(chull(UMAP1, UMAP2))

centers <- umap_plot_long_4 %>%
  group_by(algorithm, cluster) %>%
  summarise(
    cx = mean(UMAP1),
    cy = mean(UMAP2),
    .groups = "drop"
  )
label_pos <- hulls %>%
  left_join(centers, by = c("algorithm", "cluster")) %>%
  mutate(
    dist = (UMAP1 - cx)^2 + (UMAP2 - cy)^2
  ) %>%
  group_by(algorithm, cluster) %>%
  slice_max(dist, n = 1) %>%   # pick furthest hull point
  ungroup()


png_name <- paste0("results/forest_AG_c(0,1)_umap_plot_4_algorithms.png")
#png(png_name, height = 3000, width = 2100, res = 300)
g <- ggplot(umap_plot_long_4,
            aes(x = UMAP1, y = UMAP2, colour = label_group)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_polygon(
    data = hulls,
    aes(x = UMAP1, y = UMAP2, group = cluster),
    fill = NA,
    colour = "black",
    linewidth = 0.7
  ) +
  # geom_text_repel(
  #   data = label_pos,
  #   aes(UMAP1, UMAP2, label = cluster),
  #   min.segment.length = 0,
  #   show.legend = FALSE
  # )+

  geom_text(
    data = label_pos,
    aes(UMAP1, UMAP2, label = cluster),
    #bg.r = 0.15,
    size = 5,
    show.legend = FALSE
  )+
  theme_classic()+
  labs(title = "forest_AG_c(0,1)-comparison_algorithms")+
  facet_wrap(facets = "algorithm")
print(g)
#dev.off()

#save.image("Up_to_Figure20_0405.RData")

### visualisation nmds
plan(multisession, workers = 3)
visualisation_3plots_nmds <- list()
visualisation_3plots_nmds <- future_lapply(visualisation_3plots, function(name){
  plant_mat <- plant_comp_gmm[[name]]
  nmds_2d <- metaMDS(
    plant_mat,
    k = 2,
    trymax = 20
  )

  return(nmds_2d)
}, future.seed = TRUE)

# plan(sequential)
# 
# visualisation_3plots_nmds_backup <- visualisation_3plots_nmds
# visualisation_3plots_nmds <- lapply(visualisation_3plots[1], function(name){
#   plant_mat <- plant_comp_gmm[[name]]
#   
#   metaMDS(plant_mat, k = 2, trymax = 20)
# })
visualisation_3plots_nmds <- c(visualisation_3plots_nmds,visualisation_3plots_nmds_23)
names(visualisation_3plots_nmds) <- visualisation_3plots

nmds_plot_list <- list()
for(i in visualisation_3plots){
  nmds_plot <- as.data.frame(visualisation_3plots_nmds[[i]]$points)
  colnames(nmds_plot) <- c("NMDS1", "NMDS2")
  nmds_plot$label <- data_list_comp[[i]][[2]]$`Biotoptyp-Land`
  nmds_plot$label_group <- data_list_comp[[i]][[2]]$BT_Land_group
  nmds_plot$pam <- as.factor(pam_best_models[[i]]$clustering)
  nmds_plot$gmm <- as.factor(gmm_results[[i]]$classification)
  nmds_plot$hdbscan <- as.factor(hdbscan_list[[i]][[scenario_df[scenario_df[,1]== i,2]]][[(scenario_df[scenario_df[,1]== i,4])-2]]$cluster)
  
  nmds_plot_long <- pivot_longer(nmds_plot, cols = -c("NMDS1", "NMDS2", "label", "label_group"), names_to = "algorithm", values_to = "cluster")
  nmds_plot_list[[i]] <- nmds_plot_long
  png_name <- paste0("results/",substr(i,1,14),"_nmds_plot.png") # unfortunately, I used ":" ...
  png(png_name, width = 3000, height = 2100, res = 300)
  g <- ggplot(nmds_plot_long, aes(x = NMDS1, y = NMDS2, colour = cluster)) +
    geom_point(size = 2, alpha = 0.8) +
    theme_classic()+
    labs(title = i)+
    facet_wrap(facets = "algorithm")
  print(g)
  dev.off()
  
  hulls <- nmds_plot_long %>%
    group_by(algorithm, cluster) %>%   # cluster must exist in your data
    slice(chull(NMDS1, NMDS2))
  
  png_name <- paste0("results/",substr(i,1,14),"_nmds_plot_label.png") # unfortunately, I used ":" ...
  png(png_name, width = 3000, height = 2100, res = 300)
  g <- ggplot(nmds_plot_long, aes(x = NMDS1, y = NMDS2, colour = label_group)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_polygon(
      data = hulls,
      aes(x = NMDS1, y = NMDS2, group = cluster),
      fill = NA,
      colour = "black",
      linewidth = 0.7
    ) +
    theme_classic()+
    labs(title = i)+
    facet_wrap(facets = "algorithm")
  print(g)
  dev.off()
}



### metrics
combined_evaluation_df <- data.frame(dataset = scenario_names)
combined_evaluation_df <- transfer_purity(combined_evaluation_df,hdbscan_metrics_all,
                                          scenario_df)
combined_evaluation_df <- transfer_purity(combined_evaluation_df,pam_evaluation,
                                          pam_best, algorithm = "pam")
combined_evaluation_df <- transfer_purity(combined_evaluation_df,gmm_total_eval$metrics,
                                          NULL, algorithm = "gmm")

######## visualisation all
for(i in names(visualisation_3plots_umap)[5]){
  umap_plot <- visualisation_3plots_umap[[i]]
  umap_plot$label <- data_list_comp[[i]][[2]]$`Biotoptyp-Land`
  umap_plot$label_group <- data_list_comp[[i]][[2]]$BT_Land_group
  umap_plot$pam <- as.factor(pam_best_models[[i]]$clustering)
  umap_plot$hdbscan <- as.factor(hdbscan_list[[i]][[scenario_df[scenario_df[,1]== i,2]]][[(scenario_df[scenario_df[,1]== i,4])-2]]$cluster)
  
  umap_plot_long <- pivot_longer(umap_plot, cols = -c("UMAP1", "UMAP2", "label", "label_group"), names_to = "algorithm", values_to = "cluster")
  #umap_plot_list[[i]] <- umap_plot_long
  png_name <- paste0("results/",substr(i,1,14),"_umap_plot.png") # unfortunately, I used ":" ...
  #png(png_name, width = 3000, height = 2100, res = 300)
  # g <- ggplot(umap_plot_long, aes(x = UMAP1, y = UMAP2, colour = cluster)) +
  #   geom_point(size = 2, alpha = 0.8) +
  #   theme_classic()+
  #   labs(title = i)+
  #   facet_wrap(facets = "algorithm")
  # print(g)
  #dev.off()
  
  hulls <- umap_plot_long %>%
    group_by(algorithm, cluster) %>%   # cluster must exist in your data
    slice(chull(UMAP1, UMAP2))
  
  png_name <- paste0("results/",substr(i,1,14),"_umap_plot_label.png") # unfortunately, I used ":" ...
  #png(png_name, width = 3000, height = 2100, res = 300)
  g <- ggplot(umap_plot_long, aes(x = UMAP1, y = UMAP2, colour = label_group)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_polygon(
      data = hulls,
      aes(x = UMAP1, y = UMAP2, group = cluster),
      fill = NA,
      colour = "black",
      linewidth = 0.7
    ) +
    theme_classic()+
    labs(title = i)+
    facet_wrap(facets = "algorithm")
  print(g)
  #dev.off()
}


# objs <- ls()
# non_functions <- objs[!sapply(objs, function(x) is.function(get(x)))]
# save(list = non_functions, file = "all_until_nmds.RData")



##############
plants_second_round <- data_list[["forest_AG_c(0,1)"]][[1]][hdbscan_plot$cluster %in% c(0,17),]
grunddat_second_round <- filter(grunddaten_forests_red, Polygon %in% plants_second_round$Polygon)

plants_second_max <- apply(plants_second_round[,-1], 1, max)
table(as.numeric(plants_second_max)) # I need to check the 1 and 2s


check_plants <- cbind(plants_second_max, plants_second_round)
check_plants <- right_join(grunddat_temp[,c(1,4,5,25)], check_plants, by = "Polygon")

check_plants[check_plants == 0] = NA
check_polygons <- pivot_longer(check_plants, cols = where(is.numeric)& !Polygon & !plants_second_max ,names_to = "plant", values_to = "abundance", values_drop_na = TRUE) 
                         # plants_second_round[order(plants_second_round$Polygon),]


# HDBSCAN prediction probabilities ----------------------------------------



# control noise -----------------------------------------------------------

control_together <- data.frame(cluster = hdbscan_list[["forest_AG_c(0,1)"]][[1]][[11]]$cluster,
                               label= data_list[["forest_AG_c(0,1)"]][[2]]$BT_Land_group)

# PAM ---------------------------------------------------------------------
# I want to try PAM instead of agglomerative cluster algorithms because in the end I'm interested in assigning new plots to 
# existing clusters, therefore, the distance to the these clusters must be evaluated, therefore, I think medoids are better suited than 
# hierarchy trees

# https://www.statology.org/k-medoids-in-r/

pam_wAG <- weight_rel_dist(data_list$forest_wone_wAG$plants)

#fviz_nbclust(pam_wAG, FUNcluster = pam, method = "silhouette", diss = TRUE)
# as the function as used above does not work, I'm using the following workaround:


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


# figure 1 ----------------------------------------------------------------

sil <- silhouette(pam_results[["trees_AG&_c(0)"]][["9"]]$clustering, plant_comp_dist[["trees_AG&_c(0)"]])

sil_df <- as.data.frame(sil)
sil_df <- sil_df %>%
  group_by(cluster) %>%
  arrange(desc(sil_width)) %>%
  mutate(order = row_number())
ggplot(sil_df, aes(x = factor(cluster), y = sil_width, fill = factor(cluster))) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  theme_minimal()

sil_summary <- sil_df %>%
  group_by(cluster) %>%
  summarise(
    mean_sil = mean(sil_width),
    median_sil = median(sil_width),
    n = n()
  )
ggplot(sil_summary, aes(x = factor(cluster), y = mean_sil)) +
  geom_col() +
  theme_minimal()
####

plant_mat <- plant_comp_gmm[["trees_AG&_c(0)"]]
pam_visual <- umap(
  plant_mat,
  metric = "braycurtis",
  nn_method = "nndescent",
  n_neighbors = 15,
  n_components = 2,
  min_dist = 1.5
)
pam_visual <- as.data.frame(pam_visual)
colnames(pam_visual) <- c("UMAP1", "UMAP2")

#pam_visual <- visualisation_3plots_umap[["trees_AG&_c(0)"]]
pam_fit_visual <- pam_best_models[["trees_AG&_c(0)"]]
#pam_results[["trees_AG&_c(0)"]][["9"]]$clustering == pam_fit_visual$clustering
pam_visual$cluster <- pam_fit_visual$clustering
pam_visual$group <- data_list_comp[["trees_AG&_c(0)"]][[2]]$BT_Land_group

hull_pam <- pam_visual %>%
  group_by(cluster) %>%   # cluster must exist in your data
  slice(chull(UMAP1, UMAP2))
legend_df <- pam_visual[pam_fit_visual$id.med,c(3,4)]
legend_df <- left_join(legend_df,round(sil_summary,2), by = "cluster")
legend_df <- rbind(c("cluster", "group", "sil", "median_sil", "n"), legend_df)

legend_df$row <- 1:10
legend_df$x <- max(pam_visual$UMAP1)-10
legend_df$y <- max(pam_visual$UMAP2)+5
legend_df$x2 <- legend_df$x + 5
legend_df$x3 <- legend_df$x + 10

legend_df$y <- legend_df$y - legend_df$row * 2

groups <- unique(pam_visual$group)
group_colors <- setNames( paletteer_d("ggthemes::Hue_Circle", n = length(groups), direction = -1),
                          nm = groups[sample(1:length(groups), length(groups), replace=FALSE)] # prevent that biotope groups which are directly besides
) # each other, get a similar colour group

#png("results/PAM_medoids_legend.png", width = 1500, height = 1050, res = 150)
ggplot(pam_visual)+
  geom_point(aes(x = UMAP1, y = UMAP2, colour = group))+
  geom_polygon(
    data = hull_pam,
    aes(x = UMAP1, y = UMAP2, group = cluster),
    fill = NA,
    colour = "black",
    linewidth = 0.7
  ) +
  geom_point(data=pam_visual[pam_fit_visual$id.med,], aes(x = UMAP1, y = UMAP2),
             colour = "black",size=4.5)+
  geom_point(data=pam_visual[pam_fit_visual$id.med,], aes(x = UMAP1, y = UMAP2, colour = group),
             size=2,show.legend = FALSE  )+
  geom_text(data=pam_visual[pam_fit_visual$id.med,],
            aes(x = UMAP1+1.5, y = UMAP2,
              label = pam_visual$cluster[pam_fit_visual$id.med]),
             size=4.5, show.legend = FALSE)+
  geom_text(
    data = legend_df[-1,],
    aes(x = x, y = y, label = paste0(cluster,": ", group)),
    hjust = 0,
    size = 4
  )+
  geom_text(
    data = legend_df,
    aes(x = x2, y = y, label = mean_sil),
    hjust = 0,
    size = 4
  )+
  geom_text(
    data = legend_df,
    aes(x = x3, y = y, label = n),
    hjust = 0,
    size = 4
  )+
  theme_classic()+
  scale_colour_manual(values = group_colors,name = "Code group")+
  #scale_color_discrete(name = "Code group")+
  labs(title = "PAM-trees without AG")
#dev.off()

pam_sil_plot <- ggplot(legend_df[-1,], aes(x= cluster))+
  geom_col(data = legend_df[-1,],aes(y = as.numeric(mean_sil), fill = group))+
  geom_point(aes(y = as.numeric(n)/1000), shape = 15)+
  scale_y_continuous(name = "mean silhouette",sec.axis = sec_axis(~ .*1000, name = "number of points"))+
  scale_fill_manual(values = group_colors,name = "Code group")+
  theme_classic()

tab_tree <- table(hdbscan_list[["forest_AG_c(0,1)"]][[1]][[11]]$cluster, data_list_comp[["forest_AG_c(0,1)"]][[2]]$BT_Land_group)
dominant <- apply(tab_tree, 1, function(x) names(which.max(x)))

legend_labels <- data.frame(cluster = 0:(length(dominant)-1),
                            label = dominant)
legend_labels$combined <- paste0(legend_labels$cluster, ": ",legend_labels$label)
png("results/HDBSCAN_tree_forest_AG_c(0,1).png", width = 1500, height = 1050, res = 150)
plot(hdbscan_list[["forest_AG_c(0,1)"]][[1]][[11]], show_flat = TRUE,
     main = "HDBSCAN tree - forest_AG_c(0,1)")
legend("topright",
       legend = legend_labels$combined,
       ncol = 3,
       title="Main biotope group per cluster")
dev.off()

### Figure 3 -exploration
# membership probability to cluster against label
for(i in scenario_names){
  cluster_label_plot <- data.frame(
    cluster = hdbscan_list[[i]][[scenario_df$balance_scenario[scenario_df$list_name == i]]][[(scenario_df$best_k[scenario_df$list_name == i])-2]]$cluster,
                                   #label = data_list_comp[[i]][[2]]$BT_Land_group,
                                    label = data_list_comp[[i]][[2]]$`Biotoptyp-Land`,
                                   prop = hdbscan_list[[i]][[1]][[11]]$membership_prob)
  cluster_label_plot_tab <- as.data.frame(table(cluster_label_plot$cluster, cluster_label_plot$label))
  cluster_label_plot_tab <- as.data.frame(table(cluster_label_plot$label,cluster_label_plot$cluster))
  g <- ggplot(cluster_label_plot_tab,
         aes(Var2, Freq, colour = Var1, fill = Var1)) +
    geom_bar(stat = "identity", colour = "black")+
    theme_minimal()
  print(g)
}



cluster_label_plot_tab <- as.data.frame(table(pam_best_models[["trees_AG&_c(0)"]]$clustering,
                                              data_list_comp[["trees_AG&_c(0)"]][[2]]$BT_Land_group))
cluster_label_plot_tab <- as.data.frame(table(data_list_comp[["trees_AG&_c(0)"]][[2]]$BT_Land_group,
                                              pam_best_models[["trees_AG&_c(0)"]]$clustering))

### PAM
for(i in scenario_names){
  cluster_label_plot <- data.frame(
    cluster = pam_best_models[[i]]$clustering,
    label = data_list_comp[[i]][[2]]$BT_Land_group,
    prop = hdbscan_list[[i]][[1]][[11]]$membership_prob)
  cluster_label_plot_tab <- as.data.frame(table(cluster_label_plot$cluster, cluster_label_plot$label))
  cluster_label_plot_tab <- as.data.frame(table(cluster_label_plot$label,cluster_label_plot$cluster))
  g <- ggplot(cluster_label_plot_tab,
              aes(Var2, Freq, colour = Var1, fill = Var1)) +
    geom_bar(stat = "identity")+
    theme_minimal()
  print(g)
}


#Figure 3
cluster_label_plot_tab <- as.data.frame(table(data_list_comp[["trees_AG_c(0,1)"]][[2]]$BT_Land_group,
                                        hdbscan_list[["trees_AG_c(0,1)"]][[2]][[3]]$cluster))
cluster_label_plot_tab_red <- cluster_label_plot_tab%>%
  group_by(Var2)%>%
  summarise(n_points = sum(Freq), .groups = "drop")%>%
  filter(n_points>12)
names(cluster_label_plot_tab) <- c("biotope code", "cluster", "frequency")

# colour


groups <- unique(cluster_label_plot_tab$`biotope code`)
group_colors <- setNames( paletteer_d("ggthemes::Tableau_20", n = length(groups), direction = -1),
  nm = groups[sample(1:length(groups), length(groups), replace=FALSE)] # prevent that biotope groups which are directly besides
) # each other, get a similar colour group


generate_shades <- function(base_color, n) {
  lighten_vals <- seq(0.01, 0.6, length.out = n)
  colorspace::lighten(base_color, lighten_vals)
}

g1 <- ggplot(cluster_label_plot_tab[cluster_label_plot_tab$cluster %in% cluster_label_plot_tab_red$Var2,],
       aes(x = cluster,y = frequency, fill = `biotope code`)) +
 geom_bar(stat = "identity")+
  theme_minimal()+
  theme_minimal()+
  theme(legend.position = "top",legend.text = element_text(margin = margin(l = 0)),
        legend.title = element_text(hjust = 0),axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())+
  guides(fill = guide_legend(nrow = 1))+
  scale_fill_manual(values = group_colors)

# g1 <- g1+theme(axis.title.x = element_blank(),
#                axis.text.x  = element_blank(),
#                axis.ticks.x = element_blank())

cluster_label_plot_code <- as.data.frame(table(data_list_comp[["trees_AG_c(0,1)"]][[2]]$`Biotoptyp-Land`,
                                              hdbscan_list[["trees_AG_c(0,1)"]][[2]][[3]]$cluster))

cluster_label_prop <- cluster_label_plot_code %>%
  group_by(Var2) %>%                     # per cluster
  mutate(prop = Freq / sum(Freq) * 100) %>%
  ungroup()
names(cluster_label_prop) <- c("biotope code", "cluster", "frequency", "proportion")

subgroups <- unique(cluster_label_prop$`biotope code`)
subgroups_substr <- substr(subgroups, 1,2)
subgroup_colors <- c()
for(i in groups){
  sub_cols <- subgroups[subgroups_substr==i]
  col_temp <- setNames(generate_shades(group_colors[i], length(sub_cols)), sub_cols)
  subgroup_colors <- c(subgroup_colors, col_temp)
  
}


g2 <- ggplot(cluster_label_prop[cluster_label_prop$cluster %in% cluster_label_plot_tab_red$Var2,],
       aes(cluster, proportion, fill = `biotope code`)) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE)+
  guides(fill = guide_legend(nrow = 5))+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_fill_manual(values = subgroup_colors)


#png("results/trees_AG_c(0,1)_prop_without_legends_own_col.png", width = 3000, height = 3000, res = 300)
grid.arrange(g1, g2, ncol = 1, heights = c(1, 2))
#dev.off()


# test predict
add_noise_ordinal <- function(X, p = 0.1) {
  X <- as.matrix(X)
  
  dims <- dim(X)
  
  # probability of changing non-zero vs zero
  p_nonzero <- p          # e.g. 0.1
  p_zero    <- p * 0.2    # much smaller (tune this)
  
  # build mask depending on value
  noise_mask <- matrix(FALSE, nrow = dims[1], ncol = dims[2])
  
  is_nonzero <- X > 0
  is_zero    <- X == 0
  
  noise_mask[is_nonzero] <- runif(sum(is_nonzero)) < p_nonzero
  noise_mask[is_zero]    <- runif(sum(is_zero)) < p_zero
  
  # shifts (still mostly no change)
  shifts <- c(-1, 0, 1)
  probs  <- c(0.3, 0.4, 0.3)
  
  shift_matrix <- matrix(
    sample(shifts, length(X), replace = TRUE, prob = probs),
    nrow = dims[1], ncol = dims[2]
  )
  
  X_noisy <- X
  X_noisy[noise_mask] <- X_noisy[noise_mask] + shift_matrix[noise_mask]
  
  X_noisy[X_noisy < 0] <- 0
  X_noisy[X_noisy > 4] <- 4
  
  changes <- data.frame(
    row = row(noise_mask)[noise_mask],
    col = colnames(X)[col(noise_mask)[noise_mask]],
    old = X[noise_mask],
    new = X_noisy[noise_mask]
  )
  changes <- changes[changes$old != changes$new, ]
  print(changes)
  
  return(X_noisy)
}

test_sample <- sample(1:nrow(data_list[["trees_AG_c(0,1)"]][[1]]),size=10)
test_predict <- data_list[["trees_AG_c(0,1)"]][[1]][test_sample,]
#test_predict_noise <- apply(test_predict[,-1],MARGIN =2, FUN=add_noise_ordinal)
test_predict_noise <- add_noise_ordinal(test_predict[,-1])
test_predict_noise <- as.data.frame(test_predict_noise)
test_predict_noise_dist <- weight_rel_dist(test_predict_noise, ,w1 = 0.01, w2 = 0.05, w3 = 0.25, w4 = 0.9)
predict(object=hdbscan_list[["trees_AG_c(0,1)"]][[2]][[3]], newdata = test_predict_noise,
                data = data_list_comp[["trees_AG_c(0,1)"]][[1]][,-1])

hdbscan_list[["trees_AG_c(0,1)"]][[2]][[3]]$cluster[test_sample]
compare_plant <- data_list_comp[["trees_AG_c(0,1)"]][[1]][,-1][test_sample,]
# Appendix figures --------------------------------------------------------


hist_plot_land <- grunddaten_forests %>%
  group_by(`Biotoptyp-Land`)%>%
  summarise(anzahl = n())%>%
  arrange(.,desc(anzahl))

ggplot(hist_plot_land[1:30,])+
  geom_point(aes(x=`Biotoptyp-Land`, y= anzahl))+
  #scale_y_log10()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_x_discrete(limits = hist_plot_land[1:30,1][[1]])

hist_plot_land_group <- grunddaten_forests %>%
  group_by(`BT_Land_group`)%>%
  summarise(anzahl = n())%>%
  arrange(.,desc(anzahl))

ggplot(hist_plot_land_group)+
  geom_point(aes(x=`BT_Land_group`, y= anzahl))+
  #scale_y_log10()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_x_discrete(limits = hist_plot_land_group[[1]])
