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
  plants_weighted[plants_weighted == 9] <- w4 # 0.01, 9 means somehow seldom
  return(plants_weighted)
}

###
combine_list_to_df <- function(lst) {
result <- dplyr::bind_rows(
  lapply(names(lst), function(name) {
    df <- lst[[name]]
    df$list_name <- name
    df
    
  })
)
return(result)
}
###

max_weighting <- function(plants_mat, w1 = 0.01, w2 = 0.01, w3 = 0.01, w4 = 1, method = "bray"){
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

# remove biotope codes
remove_land_biotope_code <- function(plant_data, grunddat, codes = c("AG")){
  grunddat_red <- dplyr::filter(grunddat, !substr(`Biotoptyp-Land`,1,2) %in% codes)
  plant_red <- dplyr::filter(plant_data, Polygon %in% grunddat_red$Polygon)
  return(list(plants = plant_red, grunddaten = grunddat_red))
}

# hdbscan -----------------------------------------------------------------

# transform evaluation metrics of hdbscan loop into a dataframe with columns indicating the data set used and the weigthing scheme applied
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

hdbscan_complete <- function(plants_dist, by = 2, grunddat, bund = TRUE, coarse = FALSE, print = TRUE, kstop = 20){
  values <- data.frame(k = seq(3, kstop,by = by),
                       clusters = rep(0, length(seq(3, kstop,by = by))),
                       noise = rep(0, length(seq(3, kstop,by = by))),
                      ari = rep(0, length(seq(3, kstop,by = by))),
                       purity = rep(0, length(seq(3, kstop,by = by))))
  
  for (k in seq(3, kstop,by = by)) {
    h <- hdbscan(plants_dist, minPts = k)
    if(print == TRUE){
      cat("minPts =", k, "-> clusters:", length(unique(h$cluster)), 
          " noise:", sum(h$cluster == 0), "\n")
    }
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


hdbscan_metrics <- function(hdbscan_objects, grunddat){
  grunddat_short <- grunddat[,c(2,4,24,25)]
  kstop = length(hdbscan_objects)+2
  final <- data.frame(k = seq(3, kstop),
                      clusters = rep(0, length(seq(3, kstop))),
                      noise = rep(0, length(seq(3, kstop))),
                      noise_prop = rep(0, length(seq(3, kstop))))
  for (h in 1:length(hdbscan_objects)) {
    final[final$k == h+2,2] <- length(unique(hdbscan_objects[[h]]$cluster))
    final[final$k == h+2,3] <- sum(hdbscan_objects[[h]]$cluster == 0)
    final[final$k == h+2,4] <- sum(hdbscan_objects[[h]]$cluster == 0)/length(hdbscan_objects[[h]]$cluster)}
  for(o in 1:4){
    values <- data.frame(ari = rep(0, length(seq(3, kstop))),
                         purity = rep(0, length(seq(3, kstop))),
                         ari_all = rep(0, length(seq(3, kstop))),
                         purity_all = rep(0, length(seq(3, kstop))))
    
    for (h in 1:length(hdbscan_objects)) {
      #table(h$cluster)
      valid <- hdbscan_objects[[h]]$cluster != 0
      clusters_hdb <- hdbscan_objects[[h]]$cluster[valid]
      labels_bt   <- as.vector(grunddat_short[valid, o][[1]]) # check structure!!!
      
      ari <- adjustedRandIndex(clusters_hdb, labels_bt)
      values[final$k == h+2,1] <- ari  # change format
      purity <- cl_agreement(as.cl_partition(clusters_hdb),
                             as.cl_partition(labels_bt),
                             method = "purity")
      values[final$k == h+2,2] <- purity
      ## all
      clusters_all <- hdbscan_objects[[h]]$cluster
      labels_bt_all   <- as.vector(grunddat_short[, o][[1]]) # check structure!!!
      
      ari <- adjustedRandIndex(clusters_all, labels_bt_all)
      values[final$k == h+2,3] <- ari  # change format
      purity <- cl_agreement(as.cl_partition(clusters_all),
                             as.cl_partition(labels_bt_all),
                             method = "purity")
      values[final$k == h+2,4] <- purity
    }
    names(values) <- paste0(names(grunddat_short)[[o]],"_", names(values))
    final <- cbind(final, values)
  }
  
  return(final)
}

# check outlier from hdbscan
# ordination_outlier_func <- function(ord_data){
#   return(which(is.finite(ord_data$points[,1]) & 
#           abs(scale(ord_data$points[,1])) > 3 |
#           abs(scale(ord_data$points[,2])) > 3))
# }

ordination_outlier_func <- function(ord_data){
  return(which(is.finite(ord_data[,1]) & 
                 abs(scale(ord_data[,1])) > 3 |
                 abs(scale(ord_data[,2])) > 3))
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

# evaluate number of wrongly predicted data sets
hdbscan_mismatch_evaluation <- function(plants_nmds,plants_hdbscan, grunddat,
                                        bund = FALSE,
                                        coarse = FALSE){
  # Use first 3 NMDS dimensions
  df <- as.data.frame(plants_nmds$points[, 1:3])
  colnames(df) <- c("NMDS1", "NMDS2", "NMDS3")
  df$cluster <- factor(plants_hdbscan$cluster)
  
  if (bund) {
    df$biotope <- if (coarse) grunddat$`BT_Bund_group` else grunddat$`Biotoptyp-Bund`
    df$NC <- grunddat$`NC Biotoptyp-Bund`
  } else {
    df$biotope <- if (coarse) grunddat$`BT_Land_group` else grunddat$`Biotoptyp-Land`
    df$NC <- grunddat$`NC Biotoptyp-Land`
  }
  
  # looking to understand the mismatch in cluster and biotope codes
  tab <- table(plants_hdbscan$cluster, df$biotope)
  dominant <- apply(tab, 1, function(x) names(which.max(x)))
  
  df$cluster_main <- dominant[as.character(df$cluster)]
  df$mismatch <- df$biotope != df$cluster_main
  return(df)
}



# Visualisation -----------------------------------------------------------


hover_3D <- function(df){
  colourCount = length(unique(df$biotope))
  getPalette = colorRampPalette(colors = c("red","green", "blue"))
  
  # with hovering
  df$label <- paste("Cluster:", df$cluster,
                    "<br>Biotope:", df$biotope)
  
  plot_ly(df,
          x = ~NMDS1, y = ~NMDS2, z = ~NMDS3,
          color = ~biotope,
          colors = getPalette(colourCount),
          text = ~label,
          hoverinfo = "text",
          type = "scatter3d",
          mode = "markers",
          marker = list(size = 3)) %>%
    layout(scene = list(xaxis = list(title = "NMDS1"),
                        yaxis = list(title = "NMDS2"),
                        zaxis = list(title = "NMDS3")))
}

# hülle around cluster in 3d
hull_3D <- function(df, op_hull = 0.2, op_points = 0.7){
  colourCount = length(unique(df$biotope))
  getPalette = colorRampPalette(colors = c("red","green", "blue"))
  p <- plot_ly()
  
  for (grp in setdiff(unique(df$cluster), 0)) {
    sub <- df[df$cluster == grp, ]
    
    if (nrow(sub) >= 4) {
      hull <- convhulln(as.matrix(sub[, c("NMDS1","NMDS2","NMDS3")]), 
                        output.options = TRUE)
      
      p <- p %>%
        add_trace(
          type = "mesh3d",
          x = sub$NMDS1,
          y = sub$NMDS2,
          z = sub$NMDS3,
          i = hull$hull[,1] - 1,
          j = hull$hull[,2] - 1,
          k = hull$hull[,3] - 1,
          opacity = op_hull,
          colors = getPalette(colourCount),
          name = paste("Cluster", grp),
          showscale = FALSE
        )
    }
  }
  
  
  # add points on top
  p <- p %>%
    add_trace(
      data = df,
      x = ~NMDS1, y = ~NMDS2, z = ~NMDS3,
      type = "scatter3d",
      mode = "markers",
      color = ~biotope,
      colors = getPalette(colourCount),
      opacity = op_points,
      marker = list(size = 3)
    )
  
  p
  
}

### gmm functions
evaluate_gmm <- function(gmm_list){
  mclust_tabs_total <- list()
  ari_values <- list()
  bic_values <- list()
  uncertainty_values <- list()
  for(i in 1:length(gmm_list)){
    bio_data <- data_list[[i]][[2]]$`BT_Land_group`
    tab <- table(Cluster = gmm_list[[i]][[1]]$classification,
                 Biotope = bio_data)
    df_tab <- as.data.frame(tab)
    df_tab$run <- 1
    ari <- adjustedRandIndex(gmm_list[[i]][[1]]$classification, bio_data)
    bic <- gmm_list[[i]][[1]]$bic
    uncertainty <- mean(gmm_list[[i]][[1]]$uncertainty)
    for(o in 2:4){
      tab <- table(Cluster = gmm_list[[i]][[o]]$classification,
                   Biotope = bio_data)
      df_tab_temp <- as.data.frame(tab)
      df_tab_temp$run <- o
      df_tab <- rbind(df_tab,df_tab_temp)
      ari <-  c(ari,adjustedRandIndex(gmm_list[[i]][[o]]$classification, bio_data))
      bic <- c(bic, gmm_list[[i]][[o]]$bic)
      uncertainty <- c(uncertainty, mean(gmm_list[[i]][[o]]$uncertainty))
    }
    mclust_tabs_total[[i]] <- df_tab
    ari_values[[i]] <- ari
    bic_values[[i]] <- bic
    uncertainty_values[[i]] <- uncertainty
  }
  return(list(tab_plot = mclust_tabs_total, ari_values =ari_values,
              bic = bic_values, uncertainty = uncertainty_values))
}


# PAM ---------------------------------------------------------------------

# pam evaluation function
evaluate_pam_models <- function(pam_list, grunddat_list) {
  results <- lapply(names(pam_list), function(dataset_name) {
    pam_models <- pam_list[[dataset_name]]
    labels <- grunddat_list[[dataset_name]][[2]]["Biotoptyp-Land"][[1]]
    labels_group <- grunddat_list[[dataset_name]][[2]]["BT_Land_group"][[1]]
    
    df <- lapply(seq_along(pam_models), function(i) {
      
      pam_model <- pam_models[[i]]
      clusters <- pam_model$clustering
      
      ari <- mclust::adjustedRandIndex(clusters, labels)
      purity <- clue::cl_agreement(
        clue::as.cl_partition(clusters),
        clue::as.cl_partition(labels),
        method = "purity"
      )
      ari_group <- mclust::adjustedRandIndex(clusters, labels_group)
      purity_group <- clue::cl_agreement(
        clue::as.cl_partition(clusters),
        clue::as.cl_partition(labels_group),
        method = "purity"
      )
      
      data.frame(
        dataset = dataset_name,
        k = length(unique(clusters)),
        ari = ari,
        purity = purity,
        ari_group = ari_group,
        purity_group = purity_group,
        combined = 0.5*ari+0.5*purity,
        combined_group = 0.5*ari_group+0.5* purity_group
      )
    })
    
    dplyr::bind_rows(df)
  })
  
  dplyr::bind_rows(results)
}
