# Daniel Pittner

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(dbscan)
library(plotly)
library(geometry)
library(mclust)
library(clue)
library(uwot)
library(vegan)

source("03_functions.R")


# Data exploration --------------------------------------------------------

data <- load_data()

grunddaten <- data$grunddaten
plants <- data$plants
life_forms <- data$life_forms

length(unique(grunddaten$`Biotoptyp-Bund`)) # number of unique state biotope codes
length(unique(grunddaten$`Biotoptyp-Land`)) # number of unique rhineland-palatinatian biotope codes
length(unique(plants$`Wissenschaftlicher Name`)) # number of unique found plants

grunddaten$BT_Bund_group <- substr(grunddaten$`Biotoptyp-Bund`,1,5)
grunddaten$BT_Land_group <- substr(grunddaten$`Biotoptyp-Land`,1,2)
length(unique(grunddaten$BT_Bund_group))
length(unique(grunddaten$BT_Land_group))

# find plots with completely identical attributes
nrow(grunddaten %>%
  group_by(grunddaten[,-1])%>%
  filter(n()>1) %>%
  ungroup())

# plot number of occurences per state biotope code
ggplot(grunddaten)+
  geom_histogram(aes(`Biotoptyp-Bund`), stat = "count")+
  theme_classic()


hist_plot <- grunddaten %>%
  group_by(`Biotoptyp-Bund`)%>%
  summarise(anzahl = n())%>%
  arrange(.,desc(anzahl))

ggplot(hist_plot[1:30,])+
  geom_point(aes(x=`Biotoptyp-Bund`, y= anzahl))+
  scale_y_log10()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))

hist_plot_land <- grunddaten %>%
  group_by(`Biotoptyp-Land`)%>%
  summarise(anzahl = n())%>%
  arrange(.,desc(anzahl))

ggplot(hist_plot_land[1:30,])+
  geom_point(aes(x=`Biotoptyp-Land`, y= anzahl))+
  scale_y_log10()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))


# Removal of data sets with non-meaningful plant composition --------------

data_list <- prepare_data()

# Appendix Figure --------------------------------------------------------

grunddaten_forests_red <- data_list[["forest_complete"]][[2]]

hist_plot_land <- grunddaten_forests_red %>%
  group_by(`Biotoptyp-Land`)%>%
  summarise(number = n())%>%
  arrange(.,desc(number))

code_level_plot <- ggplot(hist_plot_land[1:29,])+
  geom_point(aes(x=`Biotoptyp-Land`, y= number))+
  scale_y_log10(breaks = c(1,3,10,30,100,300,500))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90),, text = element_text(size = 16))+
  scale_x_discrete(limits = hist_plot_land[1:29,1][[1]])

hist_plot_land_group <- grunddaten_forests_red %>%
  group_by(`BT_Land_group`)%>%
  summarise(number = n())%>%
  arrange(.,desc(number))

group_level_plot <- ggplot(hist_plot_land_group)+
  geom_point(aes(x=`BT_Land_group`, y= number))+
  scale_y_log10(breaks = c(1,3,10,30,100,300,1000))+
  theme_minimal()+
  xlab("Biotoptyp-Land group")+
  theme(axis.text.x = element_text(angle = 90), text = element_text(size = 16))+
  scale_x_discrete(limits = hist_plot_land_group[[1]])

#png("results/number_codes.png", width = 2500, height = 2750, res = 300)
grid.arrange(code_level_plot,group_level_plot, ncol = 1)
#dev.off()

# Different distance measures ---------------------------------------------

#check best result using weighting 2 with new weighting function and respective weights
hdbscan_complete(max_weighting(data_list[["trees_AG&_c(0)"]][[1]], w1 = 0.01, w2 = 0.05, w3 = 0.25, w4 = 0.9),by = 1, grunddat = data_list[["trees_AG&_c(0)"]][[2]], 
                 bund = FALSE, coarse = FALSE, print = FALSE)
hdbscan_complete(weight_rel_dist(data_list[["trees_AG&_c(0)"]][[1]], w1 = 0.01, w2 = 0.05, w3 = 0.25, w4 = 0.9),by = 1, grunddat = data_list[["trees_AG&_c(0)"]][[2]], 
                 bund = FALSE, coarse = FALSE, print = FALSE)
# this max_weighintg approach might improve the results slightly. NOT CONSIDERED IN REPORT!

# comparison of Euclidean against Bray-Curtis distance measure
plants_forests <- data_list[["forest_complete"]][[1]]

#Euclidean
plants_forests_euclidean <- weight_rel_dist(plants_forests, method = "euclidean")
hdbscan_complete(plants_forests_euclidean, by = 1, grunddaten_forests_red, print = FALSE)
hdbscan_evaluation(plants_forests_euclidean, k = 5)
hdb_forests_euc <- hdbscan(plants_forests_euclidean, minPts = 5)
clusterVScode(plants_forests_euclidean, 5, grunddaten_forests_red, FALSE)


#Bray-Curtis
plants_forests_bray <-  weight_rel_dist(plants_forests, method = "bray")
hdbscan_complete(plants_forests_bray, by = 1,grunddaten_forests_red, print = FALSE)
hdbscan_evaluation(plants_forests_bray, k = 6)
clusterVScode(plants_forests_bray, 6, grunddaten_forests_red, FALSE)

# PCA visualisation --------------------------------------------------------

# selected dataset
hdbscan_complete(plants_dist = weight_rel_dist(data_list[["forest_AG_c(0,1)"]][[1]]),
                 grunddat = data_list[["forest_AG_c(0,1)"]][[2]], coarse = TRUE, bund = FALSE, by = 1, print = FALSE)

clusterVScode(plants_dist = weight_rel_dist(data_list[["forest_AG_c(0,1)"]][[1]]),
              grunddat = data_list[["forest_AG_c(0,1)"]][[2]], pts = 11, bund = FALSE)

# hdbscan model for visualisation
hdbscan_plot <- hdbscan(weight_rel_dist(data_list[["forest_AG_c(0,1)"]][[1]]), minPts = 11)

# project multidemensional data into 2-D
plant_mat <-plant_weighting(data_list[["forest_AG_c(0,1)"]][[1]])
plant_mat <- decostand(plant_mat, method = "total")
pca_plants <- prcomp(plant_mat)

# two-dimensional
plot_data <- as.data.frame(pca_plants$x[, 1:2]) # extraction of two first pca-axes
plot_data$cluster <- factor(hdbscan_plot$cluster) # assigning the cluster numbers form HDBSCAN

ggplot(plot_data, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic()
# plot not satisfactory

# colouring according to biotope code
plot_data$biotope <- data_list[["forest_AG_c(0,1)"]][[2]]$`Biotoptyp-Land`
ggplot(plot_data, aes(PC1, PC2, color = biotope)) +
  geom_point(size = 2) +
  theme_minimal()

# 3-dimensional; renaming to use function, which was mainly written for NMDS visualisation...
pca_3d <- list()
pca_3d[["points"]] <- pca_plants$x
hdbscan_eval <- hdbscan_mismatch_evaluation(pca_3d,hdbscan_plot,data_list[["forest_AG_c(0,1)"]][[2]], coarse = TRUE)
sum(hdbscan_eval$mismatch==TRUE) # misclassified observations
sum(hdbscan_eval$cluster==0) # as noise classified observation

hover_3D(hdbscan_eval)
hull_3D(hdbscan_eval, op_hull = 0.6, op_points = 0.3)
# as expected. in 3d, visualisation not better

# NMDS loop --------------- 
# was used in the beginning to identify suitable NMDS representation using stress values. Not needed anymore.

# weighting <- data.frame(imbalanced = c(0.01,0.01,0.01,1),
#                         balanced = c(0.01,0.05,0.25,0.9))
# 
# scenario_names <- c("forest_AG_c(0,1)","trees_genus_AG_c(0,1)", "trees_AG&_c(0)", "trees_AG_c(0,1)", "forest_complete__c(0:2)")

# compute stress values for data sets and NMDS with 2 to 6 dimension
# stress_vals <- list(imbalanced = list(), balanced = list())
# for (o in scenario_names) {
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
# stress_vals_df <- data.frame(data = rep(scenario_names,each = 5),
#                              dimensions = rep(seq(2,6),5),
#                              imbalanced = rep(0,25),
#                              balanced = rep(0,25))

# 
# for(i in 1:2){
#   results_stress <- c()
#   for (o in 1:5){
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




# Visualisation using NMDS and effect on clustering performance ----------

# plant_mat <- plant_weighting(data_list[["forest_AG_c(0,1)"]][[1]])
# plant_mat <- decostand(plant_mat, method = "total") # first relative abundance; maybe for comparison hellinger transformation as well
# nmds_5 <- metaMDS(plant_mat, k = 5, trymax = 20) #trymax = 20
# # hdbscan_complete(nmds_5$points, by = 1, data_list[["forest_AG_c(0,1)"]][[2]], bund = FALSE, coarse = TRUE)
# 
# visualisation_3_nmds <- metaMDS(plant_mat, k = 3, trymax = 20) #trymax = 20
# # hdbscan_complete(visualisation_3_nmds$points, by = 1, data_list[["forest_AG_c(0,1)"]][[2]], bund = FALSE, coarse = TRUE)
# 
# visualisation_2_nmds <- metaMDS(plant_mat, k = 2, trymax = 20) #trymax = 20
# 
#

#saveRDS(visualisation_2_nmds, "results/visualisation_2_nmds.RDS")
#saveRDS(visualisation_3_nmds, "results/visualisation_3_nmds.RDS")
#saveRDS(nmds_5, "results/visualisation_5_nmds.RDS")


visualisation_2_nmds <- readRDS("results/visualisation_2_nmds.RDS")
visualisation_3_nmds <- readRDS("results/visualisation_3_nmds.RDS")
nmds_5 <- readRDS("results/visualisation_5_nmds.RDS")
nmds_comp <- list(nmds_5 = nmds_5$points, nmds_3 = visualisation_3_nmds$points, nmds_2 = visualisation_2_nmds$points,
                                       direct = weight_rel_dist(data_list[["forest_AG_c(0,1)"]][[1]]))

# comparison across all 4 options
nmds_metrics <- list()
for (i in 1:4) {
  evaluation_coarse <- hdbscan_complete(nmds_comp[[i]],by = 1, data_list[["forest_AG_c(0,1)"]][[2]], bund = FALSE, coarse = TRUE, print = FALSE)
  evaluation <- hdbscan_complete(nmds_comp[[i]],by = 1, data_list[["forest_AG_c(0,1)"]][[2]], bund = FALSE, print = FALSE)
  nmds_metrics[[names(nmds_comp)[i]]][[1]] <- evaluation_coarse
  nmds_metrics[[names(nmds_comp)[i]]][[2]] <- evaluation
}
helper_list <- list() #contains the plant data; only for using the hdbscan_result_df-function
for(i in 1:4){
  helper_list[[names(nmds_comp)[i]]][[1]] <- data_list[["forest_AG_c(0,1)"]][[1]]
}
nmds_wag_wone_compare <- hdbscan_result_df(nmds_metrics, main_list = helper_list)

# calculate combined metrics
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


# at the code level, the nmds dimensionality reduction can have advantages, especially if the aim is to reduce the number of points
# classified as noise

# data set 1 with minPts 13 and bt_land

#check whether ordination posses outliers
for(i in 1:3){
  outliers <- #lapply(nmds_wag_wone[[i]], function(x){
    ordination_outlier_func(nmds_comp[[i]])
  #})
  print(outliers)
}


# plot 2-dimensional NMDS 
plot_data_2d <- as.data.frame(visualisation_2_nmds$points)
plot_data_2d$cluster <- factor(hdbscan_plot$cluster)

plot(plot_data_2d)

# coloured by cluster
ggplot(plot_data_2d, aes(x = MDS1, y = MDS2, color = cluster)) +
  geom_point(size = 2, alpha = 0.8) +
  #scale_color_manual(values = c("0" = "grey70")) +
  theme_minimal()

# coloured by biotope group
plot_data_2d$biotope <- data_list[["forest_AG_c(0,1)"]][[2]]$BT_Land_group
ggplot(plot_data_2d, aes(x = MDS1, y = MDS2, color = biotope)) +
  geom_point(size = 2) +
  theme_minimal()

# prepare data for 3d plotting
grunddat_temp <- data_list[["forest_AG_c(0,1)"]][[2]]
df_plot_3d <- hdbscan_mismatch_evaluation(visualisation_3_nmds, hdbscan_plot,
                            grunddat = grunddat_temp, coarse = TRUE)

sum(df_plot_3d$mismatch) # number of assigned codes not fitting into the dominant code category per cluster

hover_3D(df_plot_3d) # interactive 3d plot
hover_3D(df_plot_3d[!df_plot_3d$cluster %in% c(0,17),]) # excluding noise observations and last, mixed cluster
# hull_3D(df_plot_3d, op_hull = 0.6, op_points = 0.3) 
hull_3D(df_plot_3d[!df_plot_3d$cluster %in% c(0,17),], op_hull = 0.6, op_points = 0.3) # hulls indicating clusters
# this visualisation is not useful

####

### test 5d nmds into 3d
df_plot_3d5 <- hdbscan_mismatch_evaluation(nmds_5, hdbscan_plot,
                                          grunddat = grunddat_temp, coarse = TRUE)

hull_3D(df_plot_3d5, op_hull = 0.6, op_points = 0.3)  # nicer hull

# Second clustering round on noise and mixed cluster ----------------------

# extract plants and basic data that are in last cluster or classified as noise
plants_second_round <- data_list[["forest_AG_c(0,1)"]][[1]][hdbscan_plot$cluster %in% c(0,17),]
grunddat_second_round <- filter(grunddaten_forests_red, Polygon %in% plants_second_round$Polygon)

# check the highest abundance per observation
plants_second_max <- apply(plants_second_round[,-1], 1, max)
table(as.numeric(plants_second_max)) # I need to check the 1 and 2s

# use a weighting scheme that also emphasizes lower maximal abundances and cluster extracted observations
plant_dist_second <- max_weighting(plants_mat=plants_second_round)
hdbscan_complete(plant_dist_second, grunddat = grunddat_second_round,coarse = TRUE, bund = FALSE, kstop = 13, print = FALSE)
# too large noise cluster

# Test Predictions -------------------------------------------------------------
# just to test how predicting newly mapped plots into clusters could work

# draw sample out of plant data
test_sample <- sample(1:nrow(data_list[["trees_AG_c(0,1)"]][[1]]),size=10)
test_predict <- data_list[["trees_AG_c(0,1)"]][[1]][test_sample,]

# add noise to the plant data and transform into distance matrix
test_predict_noise <- add_noise_ordinal(test_predict[,-1])
test_predict_noise <- as.data.frame(test_predict_noise)
test_predict_noise_dist <- weight_rel_dist(test_predict_noise,w1 = 0.01, w2 = 0.05, w3 = 0.25, w4 = 0.9)

# create hdbscan object on originial plant data
hdbscan_object_dist <- weight_rel_dist(plants_data = data_list[["trees_AG_c(0,1)"]][[1]], w1 = 0.01, w2 = 0.05, w3 = 0.25, w4 = 0.9)
hdbscan_object <- hdbscan(hdbscan_object_dist, minPts = 5)
# predict new plant data on existing hdbscan obejct
predict(object=hdbscan_object, newdata = test_predict_noise,
        data = data_list[["trees_AG_c(0,1)"]][[1]][,-1])

hdbscan_object$cluster[test_sample] # compare which observations are predicted into another cluster
# only observations to which noise was added should potentially be predicted into another cluster
