# Daniel Pittner


# Data preparation --------------------------------------------------------

# load data
load_data <- function(){
  grunddaten <- read_xlsx("data/tbl_grunddaten.xlsx") # load basic data for each plot, e.g. biotopcodes
  plants <- read_xlsx("data/tbl_daten_pflanzen.xlsx") # load plant data
  life_forms <- read_xlsx("data/Life_form.xlsx") # load life form data https://floraveg.eu/download/
  list(
    grunddaten = grunddaten,
    plants = plants,
    life_forms = life_forms
  )
}

### prepare reduced datasets
prepare_data <- function() {
  #load data
  data <- load_data()
  
  grunddaten <- data$grunddaten
  plants <- data$plants
  life_forms <- data$life_forms
  # group biotope codes
  grunddaten$BT_Bund_group <- substr(grunddaten$`Biotoptyp-Bund`, 1, 5)
  grunddaten$BT_Land_group <- substr(grunddaten$`Biotoptyp-Land`, 1, 2)
  
  
  # Removal of data sets with non-meaningful plant composition --------------
  
  remove_polygons <- filter(grunddaten,
                            substr(`Biotoptyp-Land`, 1, 1) %in% c("F", "V", "W")) # maybe G and H as well
  # check that every plot has plants
  remove_polygons2 <- filter(
    grunddaten,
    !Polygon %in% plants$Polygon &
      !substr(`Biotoptyp-Land`, 1, 1) %in% c("F", "V", "W")
  )
  
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
    dplyr::distinct(Polygon, `Wissenschaftlicher Name`, .keep_all = TRUE)
  
  plants_clean2$Menge <- as.numeric(plants_clean2$Menge) # transform abundances to numeric
  
  
  ### find outlier plots
  
  # 01728 will be removed, because only one plant found, and in 04567, Ulmus spec. will be specified to Ulmus glabra accordingly to Beschreibung
  plants_clean2 <- plants_clean2 %>%
    filter(Polygon != "01728") %>%
    mutate(
      `Wissenschaftlicher Name` = if_else(
        `Wissenschaftlicher Name` == "Ulmus spec.",
        "Ulmus glabra",
        `Wissenschaftlicher Name`
      )
    )
  
  plants_wide <- plant_widening(plants_clean2)
  
  grunddaten_sub <- grunddaten_sub %>% # remove 01728 also from grunddaten
    filter(Polygon != "01728")
  
  # Separation of biotope types ---------------------------------------------
  
  grunddaten_grass <- filter(grunddaten_sub, substr(`Biotoptyp-Land`, 1, 1) %in% c("E"))
  plants_grass <- filter(plants_wide, Polygon %in% grunddaten_grass$Polygon)
  
  grunddaten_forests <- filter(grunddaten_sub, substr(`Biotoptyp-Land`, 1, 1) %in% c("A"))
  plants_forests <- filter(plants_wide, Polygon %in% grunddaten_forests$Polygon)
  
  # # further reduction of forest data ---------------------------------------------
  
  # remove "AT... and "AU..."
  check_biotope_codes <- unique(grunddaten_forests[, c(4, 5)])
  grunddaten_forests_red <- filter(grunddaten_forests,
                                   !substr(`Biotoptyp-Land`, 1, 2) %in% c("AT", "AU", "AV"))
  plants_forests_red <- filter(plants_forests, Polygon %in% grunddaten_forests_red$Polygon)
  
  ### store plants and grunddaten in list for looping later, each list entry contains two list entries with first plants, then grunddaten
  
  data_list <- list()
  data_list[["forest_complete"]] <-  list(plants = plants_forests_red, grunddaten = grunddaten_forests_red)
  
  # # remove all polygons which only have a few plant entries
  #
  forest_complete_wone <- remove_plots(plants_forests_red,
                                       grunddaten_forests_red,
                                       remove_count = c(0, 1))
  # check_empty_plot(forest_complete_wone[[1]])
  
  data_list[["forest_complete_c(0,1)"]] <- forest_complete_wone
  
  forest_complete_wtwo <- remove_plots(plants_forests_red,
                                       grunddaten_forests_red,
                                       remove_count = c(0:2))
  data_list[["forest_complete__c(0:2)"]] <- forest_complete_wtwo
  
  # add life forms
  plants_occ <- plant_occurences(plant_data = plants_forests_red)
  
  plants_occ$short <- sub("(\\w+\\s+\\w+).*", "\\1", plants_occ$species)
  plants_occ$short <- sub("(\\w+).*", "\\1", plants_occ$species)
  life_forms$short <- sub("(\\w+).*", "\\1", life_forms$FloraVeg.Taxon)
  
  plants_LF <- left_join(plants_occ, life_forms, by = "short", multiple = "any")
  
  trees <- filter(plants_LF, Tree == 1)
  trees <- filter(trees, species %in% colnames(plants_forests))
  
  plants_forest_trees <- plants_forests_red[c("Polygon", trees$species)]
  
  plants_occ_forests <- plant_occurences(plants_forest_trees)
  
  # remove plants which only occured once or twice across all plots
  plants_f_trees_w_zero <- filter(plants_occ_forests, !total %in% c(0)) # maybe also 2
  plants_f_trees_wzero <- plants_forests_red[c("Polygon", plants_f_trees_w_zero$species)] # take all
  
  plants_f_trees_w_one <- filter(plants_occ_forests, !total %in% c(0, 1)) # maybe also 2
  plants_f_trees_wone <- plants_forests_red[c("Polygon", plants_f_trees_w_one$species)] # take all
  # I'm using the plants_f_trees_wone for further reduction of the data
  
  ### create an additional plant list simplified to genus-level
  
  # collapsing to genus level
  plant_forest_genus <- plants_forests_red %>%
    pivot_longer(-Polygon, names_to = "species", values_to = "abundance") %>%
    mutate(genus = sub("(\\w+).*", "\\1", species)) %>%
    group_by(Polygon, genus) %>%
    summarise(abundance = max(abundance, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(
      names_from = genus,
      values_from = abundance,
      values_fill = 0
    )
  
  forest_complete_wone_genus <- remove_plots(plant_forest_genus,
                                             grunddaten_forests_red,
                                             remove_count = c(0, 1))
  # check_empty_plot(forest_complete_wone_genus[[1]])
  
  data_list[["forest_genus_c(0,1)"]] <- forest_complete_wone_genus
  
  forest_complete_wtwo_genus <- remove_plots(plant_forest_genus,
                                             grunddaten_forests_red,
                                             remove_count = c(0:2))
  data_list[["forest_genus_c(0,2)"]] <- forest_complete_wtwo_genus
  
  
  # control number of tree species per plot -------------------------------------------------
  
  # check_empty_plot(plants_f_trees_wone)
  plants_trees_wzero <- remove_plots(plants_f_trees_wone, grunddaten_forests_red)
  # check_empty_plot(plants_trees_wzero[[1]])
  
  data_list[["trees_c(0)"]] <- plants_trees_wzero
  
  plants_trees_wone <- remove_plots(plants_f_trees_wone,
                                    grunddaten_forests_red,
                                    remove_count = c(0, 1))
  # check_empty_plot(plants_trees_wone[[1]])
  data_list[["trees_c(0,1)"]] <- plants_trees_wone
  
  plants_trees_wtwo <- remove_plots(plants_f_trees_wone,
                                    grunddaten_forests_red,
                                    remove_count = c(0:2))
  # check_empty_plot(plants_trees_wtwo[[1]])
  data_list[["trees_c(0:2)"]] <- plants_trees_wtwo
  
  # remove AG-biotope types -------------------------------------------------
  
  data_list[["forest_AG_c(0,1)"]] <- remove_land_biotope_code(forest_complete_wone[[1]], forest_complete_wone[[2]], c("AG"))
  data_list[["trees_AG&_c(0)"]] <- remove_land_biotope_code(plants_trees_wzero[[1]], plants_trees_wzero[[2]], c("AG"))
  data_list[["trees_AG_c(0,1)"]] <- remove_land_biotope_code(plants_trees_wone[[1]], plants_trees_wone[[2]], c("AG"))
  
  data_list[["trees_genus_AG_c(0,1)"]] <- remove_land_biotope_code(forest_complete_wone_genus[[1]],
                                                                   forest_complete_wone_genus[[2]],
                                                                   c("AG"))
  
  
  return(data_list)
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
  plants_weighted[plants_weighted == 9] <- w4 # 0.01, 9 means more or less the same as 4
  return(plants_weighted)
}

# combining abundance transformation, normalization and creation of distance matrix in one function
weight_rel_dist <- function(plants_data,w1 = 0.01,w2 = 0.01,w3 = 0.01, w4 = 1, method = "bray"){
  plants_weighted <- plant_weighting(plants_data, w1,w2,w3,w4)
  plants_rel <- decostand(plants_weighted, method = "total")
  plants_dist <- vegdist(plants_rel, method = method)
  return(plants_dist)
}

# weighting approach that always the highest abundance per observation is emphasized much more than the other observations
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


# helpful functions -------------------------------------------------------

# combine objects in lists built the same into a dataframe
combine_list_to_df <- function(lst) { # written by ChatGPT
result <- dplyr::bind_rows(
  lapply(names(lst), function(name) {
    df <- lst[[name]]
    df$list_name <- name
    df
    
  })
)
return(result)
}

# transform evaluation metrics of hdbscan loop into a dataframe with columns indicating the data set used
# and the weigthing scheme applied
hdbscan_result_df <- function(hdbscan_list, main_list){ # written by ChatGPT
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

#### combine metrics
transfer_purity <- function(combined_evaluation_df,
                            results_df,
                            scenario_df,
                            algorithm = "hdbscan") {
  
  if(algorithm == "hdbscan"){
    # Align column names
    scenario_prepared <- scenario_df %>%
      rename(
        k = best_k
      )
    
    # Keep only relevant rows from metrics
    filtered_metrics <- results_df %>%
      inner_join(
        scenario_prepared,
        by = c("list_name", "balance", "k")
      ) %>%
      select(
        list_name,
        k,
        clusters,
        noise,
        `Biotoptyp-Land_purity`,
        `Biotoptyp-Land_purity_all`,
        `BT_Land_group_purity`,
        `BT_Land_group_purity_all`,
        `BT_Land_group_ari_all`
      )%>%
      mutate(combined_group = 0.5*`BT_Land_group_purity_all`+0.5*`BT_Land_group_ari_all`)
    # Join into your evaluation df
    result <- combined_evaluation_df %>%
      left_join(
        filtered_metrics, #filtered_metrics[,-7],
        by = c("dataset" = "list_name")
      )
  } else if (algorithm == "pam"){
    # Align column names
    scenario_prepared <- scenario_df
    
    # Keep only relevant rows from metrics
    filtered_metrics <- results_df %>%
      inner_join(
        scenario_prepared,
        by = c("dataset", "k")
      ) %>%
      select(
        dataset,
        k,
        `purity`,
        `purity_group`,
        combined_group
      )
    # Join into your evaluation df
    result <- combined_evaluation_df %>%
      left_join(
        filtered_metrics,
        by = c("dataset" = "dataset"),
        suffix = c("_hdbscan", "_pam")
      )} else{
        
        # Keep only relevant rows from metrics
        filtered_metrics <- results_df %>%
          select(
            dataset,
            cluster,
            `purity`,
            `purity_group`,
            combined_group
          )
        # Join into your evaluation df
        result <- combined_evaluation_df %>%
          left_join(
            filtered_metrics,
            by = c("dataset" = "dataset")
          )
      }
  
  return(result)
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

# run hdbscan with increasing minPts and return number of clusters and noise and directly evaluate based on
# desired label from the basic data
hdbscan_complete <- function(plants_dist, by = 2, grunddat, bund = TRUE, coarse = FALSE, print = TRUE, kstop = 20){
  values <- data.frame(k = seq(3, kstop,by = by),
                       clusters = rep(0, length(seq(3, kstop,by = by))),
                       noise = rep(0, length(seq(3, kstop,by = by))),
                      ari = rep(0, length(seq(3, kstop,by = by))),
                       purity = rep(0, length(seq(3, kstop,by = by))))
  
  for (k in seq(3, kstop,by = by)) {
    h <- hdbscan(plants_dist, minPts = k)
    if(print == TRUE){ # directly print number of noise and cluster
      cat("minPts =", k, "-> clusters:", length(unique(h$cluster)), 
          " noise:", sum(h$cluster == 0), "\n")
    }
    values[values$k == k,2] <- length(unique(h$cluster))
    values[values$k == k,3] <- sum(h$cluster == 0)
    #table(h$cluster)
    valid <- h$cluster != 0 # calculate metrics excluding observations classified as noise
    clusters_hdb <- h$cluster[valid]
    labels   <- grunddat[valid,] 
    if (bund) {
      labels_bt <- if (coarse) labels$`BT_Bund_group` else labels$`Biotoptyp-Bund`
    } else {
      labels_bt <- if (coarse) labels$`BT_Land_group` else labels$`Biotoptyp-Land`
    }
    ari <- adjustedRandIndex(clusters_hdb, labels_bt)
    values[values$k == k,4] <- ari  
    purity <- cl_agreement(as.cl_partition(clusters_hdb),
                 as.cl_partition(labels_bt),
                 method = "purity")
    values[values$k == k,5] <- purity
  }
  return(values)
}

# extended hdbscan_complete function for lists with hdbscan objects, ARI and purity calculated with and without noise
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
    final[final$k == h+2,4] <- sum(hdbscan_objects[[h]]$cluster == 0)/length(hdbscan_objects[[h]]$cluster)} # calculate proportion of noise of total observations
  for(o in 1:4){
    values <- data.frame(ari = rep(0, length(seq(3, kstop))),
                         purity = rep(0, length(seq(3, kstop))),
                         ari_all = rep(0, length(seq(3, kstop))),
                         purity_all = rep(0, length(seq(3, kstop))))
    
    for (h in 1:length(hdbscan_objects)) {
      #table(h$cluster)
      valid <- hdbscan_objects[[h]]$cluster != 0
      clusters_hdb <- hdbscan_objects[[h]]$cluster[valid]
      labels_bt   <- as.vector(grunddat_short[valid, o][[1]]) 
      
      ari <- adjustedRandIndex(clusters_hdb, labels_bt)
      values[final$k == h+2,1] <- ari  # change format
      purity <- cl_agreement(as.cl_partition(clusters_hdb),
                             as.cl_partition(labels_bt),
                             method = "purity")
      values[final$k == h+2,2] <- purity
      ## all
      clusters_all <- hdbscan_objects[[h]]$cluster
      labels_bt_all   <- as.vector(grunddat_short[, o][[1]]) 
      
      ari <- adjustedRandIndex(clusters_all, labels_bt_all)
      values[final$k == h+2,3] <- ari  
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

# plot first 3 axes of NMDS representation with colouring according to biotope code and the possibility to hover over the plot
hover_3D <- function(df){ # written by ChatGPT
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

# hull representing cluster around respective observations in 3d NMDS representation
hull_3D <- function(df, op_hull = 0.2, op_points = 0.7){ # written by ChatGPT
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

# shades for biotope codes in Figure 3
generate_shades <- function(base_color, n) {
  lighten_vals <- seq(0.01, 0.6, length.out = n)
  colorspace::lighten(base_color, lighten_vals)
}

# Dimensionality reduction check ------------------------------------------

# checks whether after dimensionality reduction, outliers exist (based on first two axes)
ordination_outlier_func <- function(ord_data){
  return(which(is.finite(ord_data[,1]) & 
                 abs(scale(ord_data[,1])) > 3 |
                 abs(scale(ord_data[,2])) > 3))
}


# GMM ---------------------------------------------------------------------

# gmm evaluation function, metrics on biotope code and biotope group level + plus combined metrics of ARI and purity
# also contingency table as well as gmm specific metrics such as BIC and uncertainty. Contingency table and metrics are returned
# as two objects in one list
evaluate_gmm <- function(gmm_list, grunddat_list) {
  
  results <- lapply(names(gmm_list), function(name) {
    
    gmm_model <- gmm_list[[name]]
    
    labels_group <- grunddat_list[[name]][[2]][["BT_Land_group"]]
    labels <- grunddat_list[[name]][[2]][["Biotoptyp-Land"]]
    
    clusters <- gmm_model$classification
    
    # contingency table
    tab <- as.data.frame(table(
      Cluster = clusters,
      Biotope = labels
    ))
    
    # metrics
    ari <- mclust::adjustedRandIndex(clusters, labels)
    ari_group <- mclust::adjustedRandIndex(clusters, labels_group)

    purity <- as.numeric(clue::cl_agreement(
      clue::as.cl_partition(clusters),
      clue::as.cl_partition(labels),
      method = "purity"
    ))
    purity_group <- as.numeric(clue::cl_agreement(
      clue::as.cl_partition(clusters),
      clue::as.cl_partition(labels_group),
      method = "purity"
    ))
    # gmm metrics
    bic <- max(gmm_model$BIC, na.rm = TRUE)
    
    # uncertainty
    uncertainty <- mean(gmm_model$uncertainty)
    # number cluster
    cluster <- gmm_model$G
    
    list(
      dataset = name,
      tab = tab,
      metrics = data.frame(
        dataset = name,
        cluster = cluster,
        ari = ari,
        ari_group = ari_group,
        purity = purity,
        purity_group = purity_group,
        combined = 0.5*ari+0.5*purity,
        combined_group = 0.5*ari_group+0.5* purity_group,
        bic = bic,
        uncertainty = uncertainty
      )
    )
  })
  
  # combine outputs
  tab_plot <- lapply(results, `[[`, "tab")
  metrics <- dplyr::bind_rows(lapply(results, `[[`, "metrics"))
  
  return(list(
    tab_plot = tab_plot,
    metrics = metrics
  ))
}


# PAM ---------------------------------------------------------------------

# pam evaluation function, metrics on biotope code and biotope group level + plus combined metrics of ARI and purity
evaluate_pam_models <- function(pam_list, grunddat_list) {
  results <- lapply(names(pam_list), function(dataset_name) {
    pam_models <- pam_list[[dataset_name]]
    labels <- grunddat_list[[dataset_name]][[2]]["Biotoptyp-Land"][[1]]
    labels_group <- grunddat_list[[dataset_name]][[2]]["BT_Land_group"][[1]]
    
    df <- lapply(seq_along(pam_models), function(i) {
      
      pam_model <- pam_models[[i]]
      clusters <- pam_model$clustering
      
      ari <- mclust::adjustedRandIndex(clusters, labels)
      purity <- as.numeric(clue::cl_agreement(
        clue::as.cl_partition(clusters),
        clue::as.cl_partition(labels),
        method = "purity"
      ))
      ari_group <- mclust::adjustedRandIndex(clusters, labels_group)
      purity_group <- as.numeric(clue::cl_agreement(
        clue::as.cl_partition(clusters),
        clue::as.cl_partition(labels_group),
        method = "purity"
      ))
      
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

# HDBSCAN Prediction ------------------------------------------------------

# function to add noise to the plant data; completely written by ChatGPT
# probability to change plant abundance >0 is higher than for plants with abundance = 0. This way, the composition
# of the plots changes only slightly. Furthermore, these plants are more relevant for code assignment.
add_noise_ordinal <- function(X, p = 0.1) { # probability to add noise can be adjusted
  X <- as.matrix(X)
  
  dims <- dim(X)
  
  # probability of changing non-zero vs zero
  p_nonzero <- p          
  p_zero    <- p * 0.2    
  
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
