# Daniel Pittner

library(readxl)
library(dplyr)
library(ggplot2)
library(future.apply)
library(tidyr)
library(vegan)
library(dbscan)
library(future)
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
library(ggnewscale)

plan(multisession, workers = parallel::detectCores() - 2)

source("clusteranalysis_functions.R")

# load data: a list of plant, basic data dataset pairs with increasing simplification
data_list <- prepare_data()

# check direct hdbscan application ----------------------------------------
# hdbscan_list <- list()
# for(i in 1:length(data_list)){
#   weighting_plants <- weight_rel_dist(data_list[[i]]$plants) # weighting scheme emphasizing the most dominant species
#   hdbscan_run <- lapply(3:20, function(k) hdbscan(weighting_plants, minPts = k))
#   hdbscan_list[[names(data_list)[i]]][[1]] <- hdbscan_run
#   weighting_plants <- weight_rel_dist(data_list[[i]]$plants, w1 = 0.01, w2 = 0.05, w3 = 0.25, w4 = 0.9) # alternative more balanced weighting scheme
#   hdbscan_run <- lapply(3:20, function(k) hdbscan(weighting_plants, minPts = k))
#   hdbscan_list[[names(data_list)[i]]][[2]] <- hdbscan_run
# }
# 
# saveRDS(hdbscan_list, "results/hdbscan_list.RDS")
test_list <- readRDS("results/hdbscan_list.RDS")

# calculation of metrics for all data sets and hdbscan runs
hdbscan_metrics_all_imbalanced <- list()
hdbscan_metrics_all_balanced <- list()
for(i in 1:12){
  hdbscan_metrics_all_imbalanced[[names(data_list)[i]]] <- hdbscan_metrics(hdbscan_list[[i]][[1]], grunddat = data_list[[i]][[2]])
  hdbscan_metrics_all_balanced[[names(data_list)[i]]] <- hdbscan_metrics(hdbscan_list[[i]][[2]], grunddat = data_list[[i]][[2]])
}

# transforming the single metrics lists stored in lists (imbalanced and balanced) into data frames
hdbscan_metrics_all_imbalanced_df <- combine_list_to_df(hdbscan_metrics_all_imbalanced)
hdbscan_metrics_all_balanced_df <- combine_list_to_df(hdbscan_metrics_all_balanced)

# #### comparison whether bund biotope codes or land biotope codes fit better
check_bund_performance <- hdbscan_metrics_all_imbalanced_df %>%
  filter(`Biotoptyp-Bund_ari` > 0.3 &(`Biotoptyp-Bund_ari`>`Biotoptyp-Land_ari`))

hdbscan_metrics_all_balanced_df %>%
  filter(`Biotoptyp-Bund_ari` > 0.3 &(`Biotoptyp-Bund_ari`>`Biotoptyp-Land_ari`))

# combining both dataframes
hdbscan_metrics_all_imbalanced_df$balance <- "wrong"
hdbscan_metrics_all_balanced_df$balance <- "true"
hdbscan_metrics_all <- rbind(hdbscan_metrics_all_imbalanced_df, hdbscan_metrics_all_balanced_df)

# calculation of combined evaluation metrics
hdbscan_metrics_all <- hdbscan_metrics_all %>%
  mutate(composed_metric_land = (0.5*`Biotoptyp-Land_ari`+ 0.4*`Biotoptyp-Land_purity`+ 0.1*(1-noise_prop)),
         composed_metric_land_coarse = (0.5*`BT_Land_group_ari`+ 0.4*`BT_Land_group_purity`+ 0.1*(1-noise_prop)),
         composed_metric_land_equal = (0.3*`Biotoptyp-Land_ari`+ 0.3*`Biotoptyp-Land_purity`+ 0.3*(1-noise_prop)),
         composed_metric_land_coarse_equal = (0.3*`BT_Land_group_ari`+ 0.3*`BT_Land_group_purity`+ 0.3*(1-noise_prop)))

# bund less suitable!

# ARI can be higher on code level than on group level since ARI penalizes oversplitting.
# In my case, purity may be more relevant?!




##
names_plots <- c("Land_ARI_most", "Land_group_ARI_most", "Land_equal_weighting", "Land_group_equal_weighing")
for(i in 1:4){
  col_name <- names(hdbscan_metrics_all)[i + 22]
  #png(paste0("results/",names_plots[i],".png"), width = 3000, height = 2100, res = 300)
  p <- ggplot(hdbscan_metrics_all)+
    geom_line(aes(x=k,y= .data[[names(hdbscan_metrics_all)[i + 22]]], colour = list_name))+
    theme_classic()+
    scale_color_discrete(name = "Plant data composition")+
    labs(title = names_plots[i], y = "performance (ARI + purity + noise)", x = "minimum size of cluster")+
    facet_wrap(facets = "balance", labeller = labeller(balance=c("wrong"= "imbalanced conversion", "true"= "more balanced conversion")))
  print(p)
  #dev.off()
  if(i < 3){
    names_purity <- c("Biotoptyp-Land_purity","BT_Land_group_purity")
    #png(paste0("results/",names_purity[i],".png"), width = 3000, height = 2100, res = 300)
    p <- ggplot(hdbscan_metrics_all)+
      geom_line(aes(x=k,y= .data[[names_purity[i]]], colour = list_name))+
      theme_classic()+
      scale_color_discrete(name = "Plant data composition")+
      labs(title = names_purity[i], y = "performance (purity)", x = "minimum size of cluster")+
      facet_wrap(facets = "balance", labeller = labeller(balance=c("wrong"= "imbalanced conversion", "true"= "more balanced conversion")))
    print(p)
    #dev.off()
  }

}

# best models per metric and dataset
best_metrics_control <- hdbscan_metrics_all[,-c(5:8,11:16.19,20)] %>%
  group_by(list_name, balance) %>%
  slice_max(order_by = composed_metric_land, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(composed_metric_land))

# Datset selection --------------------------------------------------------

# scenarios to keep for comparison between clustering algorithms (with transformation scheme):
scenario_names <- c("forest_AG_c(0,1)","trees_genus_AG_c(0,1)", "trees_AG&_c(0)", "trees_AG_c(0,1)", "forest_complete__c(0:2)")
scenario_balance <- c(1,1,2,2,1) # the both tree data sets without AGs perform well with the balanced transformation,
# especially for purity at land level!

scenario_df <- data.frame(cbind(scenario_names, as.numeric(scenario_balance), c(rep("wrong",2),rep("true",2), "wrong")))
names(scenario_df) <- c("list_name", "balance_scenario", "balance")
scenario_df$balance_scenario <-  as.numeric(scenario_balance)

# list only with selected datasets
data_list_comp <- list()
data_list_comp <- data_list[names(data_list) %in% c(scenario_names)]

# control plot,that the correct scenarios are selected
ggplot(hdbscan_metrics_all[hdbscan_metrics_all$list_name %in% scenario_names,])+
  geom_line(aes(x=k,y= composed_metric_land_coarse, colour = list_name))+
  theme_classic()+
  scale_color_discrete(name = "Plant data composition")+
  labs(title = "equal weighting for performance",y = "performance (ARI + purity + noise)", x = "minimum size of cluster")+
  facet_wrap(facets = "balance", labeller = labeller(balance=c( "true"= "more balanced conversion","wrong"= "imbalanced conversion")))

# extract metrics from the 5 selected datasets
best_metrics_hdbscan <- hdbscan_metrics_all[,-c(5:8,11:16.19,20)] %>%
  semi_join(scenario_df, by = c("list_name", "balance"))%>%
  group_by(list_name, balance) %>%
  filter( rank(-composed_metric_land_coarse_equal) <= 2 |
      rank(-`Biotoptyp-Land_purity`) <= 2)
  
scenario_df$best_k <- c(13,7,6,5,20)  # selected minPts per data set


# PAM ------------------------------------------------
# https://j-sephb-lt-n.github.io/exploring_statistics/PAMS_and_SILHOUETTE_by_hand.html
# https://www.statology.org/k-medoids-in-r/

# list with distance matrices for selected data sets
plant_comp_dist <- list()
for(i in 1: length(scenario_names)){
  if(scenario_balance[i]== 2){
    dist_mat <- weight_rel_dist(data_list_comp[[i]][[1]],w1 = 0.01, w2 = 0.05, w3 = 0.25, w4 = 0.9)
  } else{
    dist_mat <- weight_rel_dist(data_list_comp[[i]][[1]],)
  }
  plant_comp_dist[[names(data_list_comp)[i]]] <- dist_mat
}


ks <- 2:25 #number of medoids to test for PAM

# run PAM over all medoids for each data set
pam_results <- future_lapply(names(plant_comp_dist), function(name) {
  dist_mat <- plant_comp_dist[[name]]
  res <- lapply(ks, function(k) {
    pam_fit <- pam(dist_mat, k, diss = TRUE)
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


# evaluate PAM results with metrics
pam_evaluation <- evaluate_pam_models(pam_results, data_list_comp)
pam_evaluation_plot <- pivot_longer(pam_evaluation, cols = -c(dataset, k), names_to = "metric", values_to = "value")

# png("results/PAM_evaluation.png", width = 3000, height = 2100, res = 350)
ggplot(pam_evaluation_plot)+
  geom_line(aes(x = k, y = value, colour = dataset))+
  theme_classic()+
  labs(title = "PAM evaluation")+
  facet_wrap(facets = "metric")
# dev.off()

# extracting silhouette width from all pam results
pam_sil_df <- purrr::imap_dfr(pam_results, function(res_list, dataset_name) {
  purrr::map_dfr(res_list, function(x) {
    data.frame(
      dataset = dataset_name,
      k = x$k,
      sil_width = x$sil_width
    )
  })
})

# plot average silhouette width against number of medoids
ggplot(pam_sil_df, aes(x = k, y = sil_width, col = dataset))+
  geom_line()+
  theme_classic()+
  labs(y = "Avg silhouette")

# extract best (highest) silhouette width per data set
pam_best <- pam_sil_df %>%
  group_by(dataset)%>%
  slice_max(order_by = sil_width, n = 1, with_ties = FALSE) %>%
  ungroup()

# run PAM again for each data set with optimal number of medoids
pam_best_models <- lapply(pam_best$dataset, function(x){
  dist_mat <- plant_comp_dist[[x]]
  k = pam_best$k[pam_best$dataset== x]
  pam_model <- cluster::pam(dist_mat, k, diss = TRUE)
})

names(pam_best_models) <- pam_best$dataset


# GMM ---------------------------------------------------------------------

# GMM does not need distance matrices, therefore only weigthing and normilization are applied
# plant_comp_gmm <- list()
# 
# for(i in 1: length(scenario_names)){
#   if(scenario_balance[i]== 2){
#     dist_mat <- plant_weighting(data_list_comp[[i]][[1]],w1 = 0.01, w2 = 0.05, w3 = 0.25, w4 = 0.9)
#   } else{
#     dist_mat <- plant_weighting(data_list_comp[[i]][[1]],)
#   }
#   rel_mat <- decostand(dist_mat, method = "total")
#   plant_comp_gmm[[names(data_list_comp)[i]]] <- dist_mat
# }
# 
# # calculation of GMM for each data set
# gmm_results <- future_lapply(names(plant_comp_gmm), function(name) {
#   rel_mat <- plant_comp_gmm[[name]]
#   gmm_fit <- Mclust(rel_mat)
#     }, future.seed = TRUE)
# 
# names(gmm_results) <- names(plant_comp_dist)
# saveRDS(gmm_results, "results/gmm_results.RDS")
gmm_results <- readRDS("results/gmm_results.RDS")

# visualize/evaluate results

gmm_total_eval <- evaluate_gmm(gmm_results, data_list_comp)
gmm_evaluation_plot <- pivot_longer(gmm_total_eval$metrics, cols = -c(dataset,cluster, bic, uncertainty), names_to = "metric", values_to = "value")

# png("results/GMM_evaluation_total.png", width = 3000, height = 2100, res = 350)
ggplot(gmm_evaluation_plot)+
  geom_path(aes(x = factor(dataset, levels = names(plant_comp_dist)), y = value, colour = metric, group = metric))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(title = "GMM evaluation-total", x = "dataset")
#dev.off()






# UMAP visualisation ------------------------------------------------------

# firstly, forest_AG_c(0,1) since all three algorithms produced somehow usable results
visualisation_3plots <- c("forest_complete__c(0:2)", "forest_AG_c(0,1)", "trees_genus_AG_c(0,1)")

visualisation_3plots_umap <- list()
visualisation_3plots_umap <- future_lapply(names(plant_comp_gmm), function(name){
  plant_mat <- plant_comp_gmm[[name]]
  umap_2d <- umap(
    plant_mat,
    metric = "braycurtis", # using bray-curtis as distance meaure
    nn_method = "nndescent",
    n_neighbors = 8,
    n_components = 2, # two dimension
    min_dist = 0.7 # adjusted so that observation appear in clusters but are not too far away from each other
  )
  umap_2d_df <- as.data.frame(umap_2d)
  colnames(umap_2d_df) <- c("UMAP1", "UMAP2")
  return(umap_2d_df)
}, future.seed = TRUE)

names(visualisation_3plots_umap) <- names(plant_comp_gmm)

# UMAP - NMDS comparison --------------------------------------------------

### visualisation nmds 
# plan(multisession, workers = 3)
# visualisation_3plots_nmds <- list()
# visualisation_3plots_nmds <- future_lapply(visualisation_3plots, function(name){
#   plant_mat <- plant_comp_gmm[[name]]
#   nmds_2d <- metaMDS(
#     plant_mat,
#     k = 2,
#     trymax = 20
#   )
# 
#   return(nmds_2d)
# }, future.seed = TRUE)
# 
# 
# visualisation_3plots_nmds <- c(visualisation_3plots_nmds,visualisation_3plots_nmds_23)
# names(visualisation_3plots_nmds) <- visualisation_3plots
# 
# nmds_plot_list <- list()
# for(i in visualisation_3plots){
#   nmds_plot <- as.data.frame(visualisation_3plots_nmds[[i]]$points)
#   colnames(nmds_plot) <- c("NMDS1", "NMDS2")
#   nmds_plot$label <- data_list_comp[[i]][[2]]$`Biotoptyp-Land`
#   nmds_plot$label_group <- data_list_comp[[i]][[2]]$BT_Land_group
#   nmds_plot$pam <- as.factor(pam_best_models[[i]]$clustering)
#   nmds_plot$gmm <- as.factor(gmm_results[[i]]$classification)
#   nmds_plot$hdbscan <- as.factor(hdbscan_list[[i]][[scenario_df[scenario_df[,1]== i,2]]][[(scenario_df[scenario_df[,1]== i,4])-2]]$cluster)
#   
#   nmds_plot_long <- pivot_longer(nmds_plot, cols = -c("NMDS1", "NMDS2", "label", "label_group"), names_to = "algorithm", values_to = "cluster")
#   nmds_plot_list[[i]] <- nmds_plot_long
#   png_name <- paste0("results/",substr(i,1,14),"_nmds_plot.png") # unfortunately, I used ":" ...
#   #png(png_name, width = 3000, height = 2100, res = 300)
#   g <- ggplot(nmds_plot_long, aes(x = NMDS1, y = NMDS2, colour = cluster)) +
#     geom_point(size = 2, alpha = 0.8) +
#     theme_classic()+
#     labs(title = i)+
#     facet_wrap(facets = "algorithm")
#   print(g)
#   #dev.off()
#   
#   hulls <- nmds_plot_long %>%
#     group_by(algorithm, cluster) %>%   # cluster must exist in your data
#     slice(chull(NMDS1, NMDS2))
#   
#   png_name <- paste0("results/",substr(i,1,14),"_nmds_plot_label.png") # unfortunately, I used ":" ...
#   #png(png_name, width = 3000, height = 2100, res = 300)
#   g <- ggplot(nmds_plot_long, aes(x = NMDS1, y = NMDS2, colour = label_group)) +
#     geom_point(size = 2, alpha = 0.8) +
#     geom_polygon(
#       data = hulls,
#       aes(x = NMDS1, y = NMDS2, group = cluster),
#       fill = NA,
#       colour = "black",
#       linewidth = 0.7
#     ) +
#     theme_classic()+
#     labs(title = i)+
#     facet_wrap(facets = "algorithm")
#   #print(g)
#   dev.off()
# }

#saveRDS(visualisation_3plots_nmds, "results/visualisation_3plots_nmds.RDS")
#visualisation_3plots_nmds <- readRDS("results/visualisation_3plots_nmds.RDS") # whole object too large for github
#saveRDS(visualisation_3plots_nmds[["forest_AG_c(0,1)"]], "results/visualisation_3plots_nmds_AG_c(0,1).RDS")

visualisation_3plots_nmds[["forest_AG_c(0,1)"]] <- readRDS("results/visualisation_3plots_nmds_AG_c(0,1).RDS")

dimesion_comparison <- data_list[["forest_AG_c(0,1)"]]
plant_mat <- plant_comp_gmm[["forest_AG_c(0,1)"]]
# used parameters were adjusted by trial and error
comparison_umap <- umap(
  plant_mat,
  metric = "braycurtis",
  nn_method = "nndescent",
  n_neighbors = 15,
  n_components = 5, # now, five dimensions
  min_dist = 1.5,
  seed = 42
)
# using UMAP data for HDBSCAN
hdbscan_complete(comparison_umap, by = 1, grunddat = dimesion_comparison[[2]],
                 bund = FALSE, print = FALSE)
hdbscan_umap <- hdbscan(comparison_umap, minPts = 8)

# NMDS on same data and also 5 dimensions for comparison
# comparison_nmds <-  metaMDS(
#   plant_mat,
#   k = 5,
#   trymax = 20
# )
# saveRDS(comparison_nmds, "results/comparison_nmds.RDS")
comparison_nmds <- readRDS("results/comparison_nmds.RDS")

# using 5-dimensional NMDS data for HDBSCAN
hdbscan_complete(comparison_nmds$points, by = 1, grunddat = dimesion_comparison[[2]], bund = FALSE)
hdbscan_nmds <- hdbscan(comparison_nmds$points, minPts = 8)
# hdbscan_metrics(list(hdbscan_nmds),dimesion_comparison[[2]])

# combining UMAP and NMDS axes 1 and 2 with clustering results and labels into one dataframe
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
# creating hulls around the clusters resulting from HDBSCAN
hull_combined_umap <- comparison_plot_long %>%
  group_by(algorithm, cluster) %>%   
  slice(chull(UMAP1, UMAP2))

# plotting hdbscan results in UMAP space
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

# plotting hdbscan results in NMDS space
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

# UMAP nicer plot -> will be used for subsequent visualisations

# comparison clustering algorithms ----------------------------------------

umap_plot_list <- list()
umap_plot_list_short <- list()
# using the UMAP representations for the 3 data sets with usable results
for(i in visualisation_3plots){
  umap_plot <- visualisation_3plots_umap[[i]] # UMAP axes
  umap_plot$label <- data_list_comp[[i]][[2]]$`Biotoptyp-Land`
  umap_plot$label_group <- data_list_comp[[i]][[2]]$BT_Land_group
  umap_plot$pam <- as.factor(pam_best_models[[i]]$clustering) # clustering results PAM
  umap_plot$gmm <- as.factor(gmm_results[[i]]$classification) # clustering result GMM
  umap_plot$hdbscan <- as.factor(hdbscan_list[[i]][[scenario_df[scenario_df[,1]== i,2]]][[(scenario_df[scenario_df[,1]== i,4])-2]]$cluster) # clustering result HDBSCAN
  umap_plot_list_short[[i]] <- umap_plot
  
  umap_plot_long <- pivot_longer(umap_plot, cols = -c("UMAP1", "UMAP2", "label", "label_group"), names_to = "algorithm", values_to = "cluster")
  umap_plot_list[[i]] <- umap_plot_long
  png_name <- paste0("results/",substr(i,1,14),"_umap_plot.png") # unfortunately, I used ":" ...
  #png(png_name, width = 3000, height = 2100, res = 300)
  # plot with observations coloured according to the cluster
  g <- ggplot(umap_plot_long, aes(x = UMAP1, y = UMAP2, colour = cluster)) +
    geom_point(size = 2, alpha = 0.8) +
    theme_classic()+
    labs(title = i)+
    facet_wrap(facets = "algorithm")
  print(g)
  #dev.off()
  
  hulls <- umap_plot_long %>%
    group_by(algorithm, cluster) %>%
    slice(chull(UMAP1, UMAP2))
  
  png_name <- paste0("results/",substr(i,1,14),"_umap_plot_label.png") # unfortunately, I used ":" ...
  #png(png_name, width = 3000, height = 2100, res = 300)
  # plot with observations coloured according to the label and clusters cisualized as hulls
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

##### comparison of all 4 (including HDBSCAN on NMDS-5)
# using only one data set for this comparison
umap_plot_4 <- umap_plot_list_short[["forest_AG_c(0,1)"]]
umap_plot_4$hdbscan_nmds <- as.factor(hdbscan_nmds$cluster)
umap_plot_long_4 <- pivot_longer(umap_plot_4, cols = -c("UMAP1", "UMAP2", "label", "label_group"), names_to = "algorithm", values_to = "cluster")

hulls <- umap_plot_long_4 %>%
  group_by(algorithm, cluster) %>%  
  slice(chull(UMAP1, UMAP2))

# calculate center in UMAP space per cluster
centers <- umap_plot_long_4 %>%
  group_by(algorithm, cluster) %>%
  summarise(
    cx = mean(UMAP1),
    cy = mean(UMAP2),
    .groups = "drop"
  )

# label position per cluster at hull
label_pos <- hulls %>%
  left_join(centers, by = c("algorithm", "cluster")) %>%
  mutate(
    dist = (UMAP1 - cx)^2 + (UMAP2 - cy)^2
  ) %>%
  group_by(algorithm, cluster) %>%
  slice_max(dist, n = 1) %>%   # pick furthest hull point
  ungroup()

hulls$highlight <- ifelse(hulls$cluster == 0, "noise", "other") # hull variable used for colouring

# colouring: each biotope group gets one colour, consistent across subsequent plots
set.seed(42)
groups <- unique(umap_plot_4$label_group)
group_colors <- setNames( paletteer_d("ggthemes::Hue_Circle", n = length(groups), direction = -1),
                          nm = groups[sample(1:length(groups), length(groups), replace=FALSE)] # prevent that biotope groups which are directly besides
) # each other, get a similar colour group


png_name <- paste0("results/forest_AG_c(0,1)_umap_plot_4_algorithms.png")
#png(png_name, height = 2800, width = 2100, res = 300)
g <- ggplot(umap_plot_long_4,
            aes(x = UMAP1, y = UMAP2, colour = label_group)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_colour_manual(values = group_colors,name = "Biotope\ngroup")+
  new_scale_color() +   
  geom_polygon(
    data = hulls,
    aes(x = UMAP1, y = UMAP2, group = cluster, colour = highlight),
    fill = NA,
    linewidth = 0.7
  ) +
  scale_colour_manual(values = c("noise" = "red", "other" = "black"), name = "cluster")+
  # geom_text_repel( # option to add segments in case too many labels at a place
  #   data = label_pos,
  #   aes(UMAP1, UMAP2, label = cluster),
  #   min.segment.length = 0,
  #   show.legend = FALSE
  # )+

  geom_text(
    data = label_pos,
    aes(UMAP1, UMAP2, label = cluster),
    size = 5,
    show.legend = FALSE
  )+
  theme_classic()+
  labs(title = "forest_AG_c(0,1)-comparison_algorithms")+
  facet_wrap(facets = "algorithm")
  
print(g)
#dev.off()

#save.image("Up_to_Figure20_0405.RData")




# metrics comparison ------------------------------------------------------

combined_evaluation_df <- data.frame(dataset = scenario_names)
combined_evaluation_df <- transfer_purity(combined_evaluation_df,hdbscan_metrics_all,
                                          scenario_df)
combined_evaluation_df <- transfer_purity(combined_evaluation_df,pam_evaluation,
                                          pam_best, algorithm = "pam")
combined_evaluation_df <- transfer_purity(combined_evaluation_df,gmm_total_eval$metrics,
                                          NULL, algorithm = "gmm")
combined_evaluation_df_print <- cbind(combined_evaluation_df[,1],round(combined_evaluation_df[,-1],2))
#write.csv(t(combined_evaluation_df_print), file = "results/algorithm_evaluation.csv")



# Figure 1 ----------------------------------------------------------------

# calculate silhouette width for each observation based on PAM clustering results for data set with highest PAM purity
sil <- silhouette(pam_results[["trees_AG&_c(0)"]][["9"]]$clustering, plant_comp_dist[["trees_AG&_c(0)"]])
sil_df <- as.data.frame(sil)

sil_df <- sil_df %>%
  group_by(cluster) %>%
  arrange(desc(sil_width)) %>%
  mutate(order = row_number())
# plot distribution of silhouette widths per cluster
ggplot(sil_df, aes(x = factor(cluster), y = sil_width, fill = factor(cluster))) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  theme_minimal()

# calculate mean and silhouette width and number of observations per cluster
sil_summary <- sil_df %>%
  group_by(cluster) %>%
  summarise(
    mean_sil = mean(sil_width),
    n = n()
  )

# calculate UMAP representation with individual parameters for most appealing representation
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
pam_visual$cluster <- pam_fit_visual$clustering
pam_visual$group <- data_list_comp[["trees_AG&_c(0)"]][[2]]$BT_Land_group

# check main biotope group per cluster and accordance to respective medoid
check_max <- table(pam_visual$cluster,pam_visual$group)
colnames(check_max)[max.col(as.matrix(check_max), ties.method = "first")]
colnames(check_max)[max.col(as.matrix(check_max), ties.method = "first")] == pam_visual[pam_fit_visual$id.med,c(4)]

# hull
hull_pam <- pam_visual %>%
  group_by(cluster) %>% 
  slice(chull(UMAP1, UMAP2))


legend_df <- pam_visual[pam_fit_visual$id.med,c(3,4)]
legend_df <- left_join(legend_df,round(sil_summary,2), by = "cluster")



#png("results/PAM_medoids_legend.png", width = 1500, height = 1050, res = 150)
pam_umap <- ggplot(pam_visual)+
  geom_point(aes(x = UMAP1, y = UMAP2, colour = group),show.legend = FALSE)+ # umap points
  geom_polygon( # hulls
    data = hull_pam,
    aes(x = UMAP1, y = UMAP2, group = cluster),
    fill = NA,
    colour = "black",
    linewidth = 0.7
  ) +
  geom_point(data=pam_visual[pam_fit_visual$id.med,], aes(x = UMAP1, y = UMAP2), # medoids black background
             colour = "black",size=4.5)+
  geom_point(data=pam_visual[pam_fit_visual$id.med,], aes(x = UMAP1, y = UMAP2, colour = group),
             size=2,show.legend = FALSE)+ # medoids colour
  geom_text(data=pam_visual[pam_fit_visual$id.med,], # respective cluster number of medoid
            aes(x = UMAP1+1.5, y = UMAP2,
              label = pam_visual$cluster[pam_fit_visual$id.med]),
             size=6.5, show.legend = FALSE)+
  theme_classic()+
  scale_colour_manual(values = group_colors,name = "Code group")+
  theme(text = element_text(size = 15))#+
  #labs(title = "PAM-trees without AG")

#dev.off()
# help data frame to create a legend havong entries for each biotope group
colour_df_point <- data.frame(cluster = rep(1,length(unique(pam_visual$group))),
                              vals = rep(0.01, length(unique(pam_visual$group))),
                              group = unique(pam_visual$group))

pam_sil_plot <- ggplot(legend_df[-1,], aes(x= cluster))+
  geom_col(data = colour_df_point, aes(x = as.factor(cluster), y = vals, fill = group))+
  geom_col(data = legend_df[-1,],aes(y = as.numeric(mean_sil), fill = group))+
  geom_point(aes(y = as.numeric(n)/1000, colour = "number of points"), shape = 16, size = 5)+
  scale_y_continuous(name = "mean silhouette",sec.axis = sec_axis(~ .*1000, name = "number of points"))+
  scale_fill_manual(values = group_colors,name = "Biotope group")+
  scale_colour_manual(values = c("number of points" = "black"), name = "")+
  xlab("cluster")+
  guides(fill = guide_legend(ncol = 3))+
  theme_classic()+
  theme(text = element_text(size = 15),legend.spacing.y = unit(0, "cm"))

#png("results/trees_AG&_c(0)_pam_silwidth.png", width = 3000, height = 3000, res = 300)
grid.arrange(pam_sil_plot, pam_umap, ncol = 1, heights = c(0.8, 2))
#dev.off()

### visualisation of HDBSCAN tree (alternative to Figure 1)
# calculate dominant biotope group per cluster
tab_tree <- table(hdbscan_list[["forest_AG_c(0,1)"]][[1]][[11]]$cluster, data_list_comp[["forest_AG_c(0,1)"]][[2]]$BT_Land_group)
dominant <- apply(tab_tree, 1, function(x) names(which.max(x)))

legend_labels <- data.frame(cluster = 0:(length(dominant)-1), # combine cluster with biotope group
                            label = dominant)
legend_labels$combined <- paste0(legend_labels$cluster, ": ",legend_labels$label)
#png("results/HDBSCAN_tree_forest_AG_c(0,1).png", width = 1500, height = 1050, res = 150)
plot(hdbscan_list[["forest_AG_c(0,1)"]][[1]][[11]], show_flat = TRUE,
     main = "HDBSCAN tree - forest_AG_c(0,1)")
legend("topright",
       legend = legend_labels$combined,
       ncol = 3,
       title="Main biotope group per cluster")
#dev.off()


## Figure 3 exploration ----------------------------------------------------------------

# exploration which data set to use for HDBSCAN cluster visualisation
for(i in scenario_names){
  cluster_label_plot <- data.frame(
    cluster = hdbscan_list[[i]][[scenario_df$balance_scenario[scenario_df$list_name == i]]][[(scenario_df$best_k[scenario_df$list_name == i])-2]]$cluster,
                                   #label = data_list_comp[[i]][[2]]$BT_Land_group,
                                    label = data_list_comp[[i]][[2]]$`Biotoptyp-Land`,
                                   prop = hdbscan_list[[i]][[1]][[11]]$membership_prob)
  #cluster_label_plot_tab <- as.data.frame(table(cluster_label_plot$cluster, cluster_label_plot$label))
  cluster_label_plot_tab <- as.data.frame(table(cluster_label_plot$label,cluster_label_plot$cluster))
  g <- ggplot(cluster_label_plot_tab,
         aes(Var2, Freq, colour = Var1, fill = Var1)) +
    geom_bar(stat = "identity", colour = "black")+
    theme_minimal()
  print(g)
}


### PAM for comparison

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


# Figure 3 ----------------------------------------------------------------
# "trees_AG_c(0,1)" with HDBSCAN
cluster_label_plot_tab <- as.data.frame(table(data_list_comp[["trees_AG_c(0,1)"]][[2]]$BT_Land_group,
                                        hdbscan_list[["trees_AG_c(0,1)"]][[2]][[3]]$cluster))
# select clusters with more than 12 observation
cluster_label_plot_tab_red <- cluster_label_plot_tab%>% 
  group_by(Var2)%>%
  summarise(n_points = sum(Freq), .groups = "drop")%>%
  filter(n_points>12)
names(cluster_label_plot_tab) <- c("biotope group", "cluster", "frequency")

# plot clusters with more than 12 observations and their observations coloured by biotope group
g1 <- ggplot(cluster_label_plot_tab[cluster_label_plot_tab$cluster %in% cluster_label_plot_tab_red$Var2,],
       aes(x = cluster,y = frequency, fill = `biotope group`)) +
 geom_bar(stat = "identity")+
  theme_classic()+
  theme(legend.position = "top",legend.text = element_text(margin = margin(l = 0)),
        legend.title = element_text(hjust = 0),axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),text = element_text(size = 16))+
  guides(fill = guide_legend(nrow = 2))+
  scale_fill_manual(values = group_colors)

# number of observations by biotope code and cluster
cluster_label_plot_code <- as.data.frame(table(data_list_comp[["trees_AG_c(0,1)"]][[2]]$`Biotoptyp-Land`,
                                              hdbscan_list[["trees_AG_c(0,1)"]][[2]][[3]]$cluster))
# proportion of biotope code per cluster
cluster_label_prop <- cluster_label_plot_code %>%
  group_by(Var2) %>%                     # per cluster
  mutate(prop = Freq / sum(Freq) * 100) %>%
  ungroup()
names(cluster_label_prop) <- c("biotope code", "cluster", "frequency", "proportion")

# to generate the shading of biotope codes according to biotope group
subgroups <- unique(cluster_label_prop$`biotope code`)
subgroups_substr <- substr(subgroups, 1,2)
subgroup_colors <- c()
for(i in groups){
  sub_cols <- subgroups[subgroups_substr==i]
  col_temp <- setNames(generate_shades(group_colors[i], length(sub_cols)), sub_cols)
  subgroup_colors <- c(subgroup_colors, col_temp)
  
}

# plot proportions of biotope code per cluster (with n > 12) coloured by biotope code as shade of biotope group colour
g2 <- ggplot(cluster_label_prop[cluster_label_prop$cluster %in% cluster_label_plot_tab_red$Var2,],
       aes(cluster, proportion, fill = `biotope code`)) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE)+
  guides(fill = guide_legend(nrow = 5))+
  theme_minimal()+
  theme(legend.position = "bottom", text = element_text(size = 16))+
  scale_fill_manual(values = subgroup_colors)


#png("results/trees_AG_c(0,1)_prop_without_legends_own_col.png", width = 3000, height = 3200, res = 300)
grid.arrange(g1, g2, ncol = 1, heights = c(1, 2))
#dev.off()