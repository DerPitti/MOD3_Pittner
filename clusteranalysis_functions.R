# widening of plant data
plant_widening <- function(plant_data){
  plants_wide <- pivot_wider(plant_data[,c(1,2,4)],names_from = `Wissenschaftlicher Name`,values_from = Menge)
  plants_wide[is.na(plants_wide)] <- 0
  return(plants_wide)
}

# Community abundance-weighted transformation
plant_weighting <- function(plant_data, w1 = 0.01, w2 = 0.01, w3 = 0.01, w4 = 1){
  plants_weighted <- plant_data[,-1]
  plants_weighted[plants_weighted == 1] <- w1 # 0.01
  plants_weighted[plants_weighted == 2] <- w2 # 0.05
  plants_weighted[plants_weighted == 3] <- w3 # 0.25
  plants_weighted[plants_weighted == 4] <- w4 # 0.75
  return(plants_weighted)
}

# relative abundances and euclidean distance
weight_rel_dist <- function(plants_data,w1 = 0.01,w2 = 0.01,w3 = 0.01, w4 = 1, method = "bray"){
  plants_weighted <- plant_weighting(plants_data, w1,w2,w3,w4)
  plants_rel <- decostand(plants_weighted, method = "total")
  plants_dist <- vegdist(plants_rel, method = method)
  return(plants_dist)
}

# find all plants with only few entries
plant_occurences <- function(plant_data){
  plants <- plant_data[,-1] %>%
    mutate(across(everything(), ~ if_else(. > 0, 1L, 0L)))%>%
    summarise(across(everything(), sum, na.rm=TRUE))%>%
    pivot_longer(everything(), names_to = "species", values_to = "total") %>%
    arrange(total)
  return(plants)
}

# calculate number of plants per plant entry
check_empty_plot <- function(plant_data){
  plant_red <- plant_data[-1]
  plant_red[plant_red>0] <- 1
  plot <- rowSums(plant_red)
  print(sort(plot))
  as.vector(plant_data[order(plot),1])
}

# remove polygons with only a certain number of plants
remove_plots <- function(plant_data, grunddaten, remove_count = c(0)){
  plants_data_red <- plant_data[-1]
  plants_data_red[plants_data_red > 0] <- 1
  plants_data <-  filter(plant_data, !rowSums(plants_data_red) %in% remove_count)
  grunddaten <- filter(grunddaten, Polygon %in% plants_data$Polygon)
  return(list(plants = plants_data, grunddaten = grunddaten))
}

# hdbscan -----------------------------------------------------------------


# run hdbscan with increasing minPts
hdbscan_minClusSize <- function(distance_data, by = 2){
  for (k in seq(3, 20,by = by)) {
    h <- hdbscan(distance_data, minPts = k)
    cat("minPts =", k, "-> clusters:", length(unique(h$cluster)), 
        " noise:", sum(h$cluster == 0), "\n")
  }
}

# create table showing number of entries per clusters 
hdbscan_evaluation <- function(dist_mat, k){
  dist_hdbscan <- hdbscan(dist_mat, minPts = k)
  table(dist_hdbscan$cluster)
}

# create table with biotope codes vs. cluster
clusterVScode <- function(plants_dist, pts, grunddat, bund = TRUE){
  hdbscan_plants <- hdbscan(plants_dist, minPts = pts)
  if(bund == TRUE)
  {t(table(hdbscan_plants$cluster, grunddat$`Biotoptyp-Bund`))}
  else {
    t(table(hdbscan_plants$cluster, grunddat$`Biotoptyp-Land`))
  }
}

hdbscan_complete <- function(plants_dist, by = 2, grunddat, bund = TRUE, coarse = FALSE){
  values <- data.frame(k = seq(3, 20,by = by),
                       clusters = rep(0, length(seq(3, 20,by = by))),
                       noise = rep(0, length(seq(3, 20,by = by))),
                      ari = rep(0, length(seq(3, 20,by = by))),
                       purity = rep(0, length(seq(3, 20,by = by))))
  
  for (k in seq(3, 20,by = by)) {
    h <- hdbscan(plants_dist, minPts = k)
    cat("minPts =", k, "-> clusters:", length(unique(h$cluster)), 
        " noise:", sum(h$cluster == 0), "\n")
    values[values$k == k,2] <- length(unique(h$cluster))
    values[values$k == k,3] <- sum(h$cluster == 0)
    #table(h$cluster)
    valid <- h$cluster != 0
    clusters_hdb <- h$cluster[valid]
    labels   <- grunddat[valid,] # check structure!!!
    if (bund) {
      labels_bt <- if (coarse) labels$`BT_Bund_group` else labels$`Biotoptyp-Bund`
    } else {
      labels_bt <- if (coarse) labels$`BT_Land_group` else labels$`Biotoptyp-Land`
    }
    ari <- adjustedRandIndex(clusters_hdb, labels_bt)
    values[values$k == k,4] <- ari  # change format
    purity <- cl_agreement(as.cl_partition(clusters_hdb),
                 as.cl_partition(labels_bt),
                 method = "purity")
    values[values$k == k,5] <- purity
  }
  return(values)
}

# check outlier from hdbscan
ordination_outlier_func <- function(ord_data){
  return(which(is.finite(ord_data$points[,1]) & 
          abs(scale(ord_data$points[,1])) > 3 |
          abs(scale(ord_data$points[,2])) > 3))
}

# evaluate hdbscan graphically
hdbscan_plot <- function(data, name){
  ggplot(data, aes(x = k))+
    geom_line(aes(y= clusters, colour = "Clusters"))+
    geom_line(aes(y= noise, colour = "Noise"))+
    geom_line(aes(y= ari*2000, colour = "ARI"))+
    geom_line(aes(y= purity*2000, colour = "Purity"))+
    scale_y_continuous(sec.axis = sec_axis(~ . /2000, name = "ARI / Purity"))+
    labs(title = name)+
    scale_colour_manual(
      name = "Metric",
      values = c(
        "Clusters" = "black",
        "Noise" = "red",
        "ARI" = "blue",
        "Purity" = "green"
      )
    ) +
    theme_minimal()
}
