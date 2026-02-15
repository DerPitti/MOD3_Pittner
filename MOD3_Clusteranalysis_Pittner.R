library(readxl)
library(dplyr)
library(ggplot2)
library(future.apply)
library(tidyr)
library(dplyr)
library(vegan)
library("dbscan")

grunddaten <- read_xlsx("data/tbl_grunddaten.xlsx") # load basic data for each plot, e.g. biotopcodes
plants <- read_xlsx("data/tbl_daten_pflanzen.xlsx") # load plant data
life_forms <- read_xlsx("data/Life_form.xlsx") # load life form data https://floraveg.eu/download/


length(unique(grunddaten$`Biotoptyp-Bund`)) # number of unique state biotope codes
length(unique(grunddaten$`Biotoptyp-Land`)) # number of unique rhineland-palatinatian biotope codes
length(unique(plants$`Wissenschaftlicher Name`)) # number of unique found plants

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

#length(grunddaten[`Biotoptyp-Land`in ()])

remove_polygons <- filter(grunddaten, substr(`Biotoptyp-Land`,1,1) %in% c("F", "V"))

# check that every plot has plants
remove_polygons2 <- filter(grunddaten, !Polygon %in% plants$Polygon & !substr(`Biotoptyp-Land`,1,1) %in% c("F", "V"))

remove_poly <- rbind(remove_polygons, remove_polygons2)

plants_sub <- filter(plants, !Polygon %in% remove_polygons2$Polygon)
grunddaten_sub <- filter(grunddaten, !Polygon %in% remove_polygons2$Polygon)


# transformation of plant data frame

# check for completely identical data sets (one plant species with identical abundance in one plot)
plants_sub %>%
  dplyr::summarise(n = dplyr::n(), .by = c(Polygon, `Wissenschaftlicher Name`, Menge)) |>
  dplyr::filter(n > 1L)

# remove duplicate rows
plants_clean <- plants_sub %>%
  dplyr::distinct(Polygon, `Wissenschaftlicher Name`, Menge, .keep_all = TRUE)

# find plots with same plants but different abundance
plants_clean %>%
  dplyr::summarise(n = dplyr::n(), .by = c(Polygon, `Wissenschaftlicher Name`)) |>
  dplyr::filter(n > 1L)

# only keep row with highest abundance
plants_clean2 <- plants_clean %>%
  dplyr::arrange(Polygon, `Wissenschaftlicher Name`, desc(Menge)) %>%
  dplyr::distinct(Polygon, `Wissenschaftlicher Name`, .keep_all = TRUE
  )

plants_clean2$Menge <- as.numeric(plants_clean2$Menge) # transform abundances to numeric

# widening of plant data
plant_widening <- function(plant_data){
  plants_wide <- pivot_wider(plant_data[,c(1,2,4)],names_from = `Wissenschaftlicher Name`,values_from = Menge)
  plants_wide[is.na(plants_wide)] <- 0
  return(plants_wide)
}
plants_wide <- plant_widening(plants_clean2)

# check that every plot has plants
filter(grunddaten_sub, !Polygon %in% plants_wide$Polygon)

plants_sum <- cbind(rowSums(plants_wide[,-1]),plants_wide)
order(rowSums(plants_wide[,-1]))

### find outlier plots

# Community abundance-weighted transformation
plant_weighting <- function(plant_data, w1 = 0.01, w2 = 0.01, w3 = 0.01, w4 = 1){
  plants_weighted <- plant_data[,-1]
  plants_weighted[plants_weighted == 1] <- w1 # 0.01
  plants_weighted[plants_weighted == 2] <- w2 # 0.05
  plants_weighted[plants_weighted == 3] <- w3 # 0.25
  plants_weighted[plants_weighted == 4] <- w4 # 0.75
  return(plants_weighted)
  }

# first attempt to see if clustering might be successful

plants_weighted <- plant_weighting(plant_data = plants_wide)
bray_curtis_dist_weighted <- vegdist(plants_weighted, method = "bray")


#options(future.validate = FALSE)
#plan(multisession, workers = parallel::detectCores() - 2) # start multisession, use all cores except two

ord_new <- metaMDS(bray_curtis_dist_weighted, k = 2) # 2-dimensional NMDS
plot(ord_new, type = "t") # plot NMDS to see if pattern exists -> no..., but outliers
stressplot(ord_new)

# find outliers in nmds ordination
which(is.finite(ord$points[,1]) & 
        abs(scale(ord$points[,1])) > 3 |
        abs(scale(ord$points[,2])) > 3)

ord$points[order(ord$points[,1], decreasing = TRUE),] # 711, 1033, 1901 seem to be outlier


filter(plants_clean2, Polygon %in% plants_wide$Polygon[order(ord$points[,1], decreasing = TRUE)[1:3]])
suspicious_polygons <- filter(grunddaten_sub, Polygon %in% plants_wide$Polygon[order(ord$points[,1], decreasing = TRUE)[1:3]])
# 01728 will be removed, because only one plant found, and in 04567, Ulmus spec. will be specified to Ulmus glabra accordingly to Beschreibung
plants_clean2 <- plants_clean2 %>%
  filter(Polygon != "01728") %>%
  mutate(`Wissenschaftlicher Name` = if_else(`Wissenschaftlicher Name`=="Ulmus spec.", "Ulmus glabra", `Wissenschaftlicher Name`))

plants_wide <- plant_widening(plants_clean2)

grunddaten_sub <- grunddaten_sub %>% # remove 01728 also from grunddaten
  filter(Polygon != "01728")

# find all plants with only few entries
plant_occurences <- function(plant_data){
  plants <- plant_data[,-1] %>%
    mutate(across(everything(), ~ if_else(. > 0, 1L, 0L)))%>%
    summarise(across(everything(), sum, na.rm=TRUE))%>%
    pivot_longer(everything(), names_to = "species", values_to = "total") %>%
    arrange(total)
  return(plants)
}

plants_occ <- plant_occurences(plant_data = plants_wide)

filter(plants_occ, species %in% c("Ulmus spec.","Robinia pseudoacacia"))

# Separation of biotope types ---------------------------------------------

grunddaten_grass <- filter(grunddaten_sub, substr(`Biotoptyp-Land`,1,1) %in% c("E"))
plants_grass <- filter(plants_wide, Polygon %in% grunddaten_grass$Polygon)

grunddaten_forests <- filter(grunddaten_sub, substr(`Biotoptyp-Land`,1,1) %in% c("A"))
plants_forests <- filter(plants_wide, Polygon %in% grunddaten_forests$Polygon)


plants_grass_w <- plant_weighting(plants_grass, 0.01, 0.05, 0.25, 0.75)
grass_bc_dist_weighted <- vegdist(plants_grass_w, method = "bray")

hdbscan_minClusSize <- function(distance_data){
  for (k in seq(3, 15,by = 2)) {
    h <- hdbscan(distance_data, minPts = k)
    cat("minPts =", k, "-> clusters:", length(unique(h$cluster)), 
        " noise:", sum(h$cluster == 0), "\n")
  }
  
}

hdbscan_minClusSize(grass_bc_dist_weighted)

hdbscan_evaluation <- function(dist_mat, k){
  dist_hdbscan <- hdbscan(dist_mat, minPts = k)
  table(dist_hdbscan$cluster)
}
hdbscan_evaluation(grass_bc_dist_weighted, k= 7)

clusterVScode <- function(plants_dist, pts, grunddat, bund = TRUE){
  hdbscan_plants <- hdbscan(plants_dist, minPts = pts)
  if(bund == TRUE)
    {t(table(hdbscan_plants$cluster, grunddat$`Biotoptyp-Bund`))}
  else {
    t(table(hdbscan_plants$cluster, grunddat$`Biotoptyp-Land`))
  }
}

clusterVScode(grass_bc_dist_weighted, 7, grunddaten_grass)
clusterVScode(grass_bc_dist_weighted, 7, grunddaten_grass, FALSE)


ord <- metaMDS(grass_bc_dist_weighted, k = 2, trymax = 100)
plot(ord, type = "t")

# Different distance measures ---------------------------------------------

# relative abundances and euclidean distance
weight_rel_dist <- function(plants_data,w1 = 0.01,w2 = 0.01,w3 = 0.01, w4 = 1, method = "euclidean"){
  plants_weighted <- plant_weighting(plants_data, w1,w2,w3,w4)
  plants_rel <- decostand(plants_weighted, method = "total")
  plants_dist <- vegdist(plants_rel, method = method)
  return(plants_dist)
}

plants_forests_euclidean <- weight_rel_dist(plants_forests)
hdbscan_minClusSize(plants_forests_euclidean)
hdbscan_evaluation(plants_forests_euclidean, k = 5)

clusterVScode(plants_forests_euclidean, 5, grunddaten_forests)
forests_check1 <- as.data.frame(clusterVScode(plants_forests_euclidean, 5, grunddaten_forests, FALSE))

# hdb_forests_euc1 <- hdbscan(plants_forests_euclidean, minPts = 5)
# library(caret)
# forests_euc_conf_mat <- confusionMatrix(factor(hdb_forests_euc1$cluster), factor(grunddaten_forests$`Biotoptyp-Land`), mode = "everything", positive="1")

plants_forests_bray <-  weight_rel_dist(plants_forests, method = "bray")
hdbscan_minClusSize(plants_forests_bray)
hdbscan_evaluation(plants_forests_bray, k = 6)

clusterVScode(plants_forests_bray, 6, grunddaten_forests, FALSE)

# comparing of 2 clustering solutions
hdb_forests_bray <- hdbscan(plants_forests_bray, minPts = 5)
library(fpc)
cluster.stats(d, hdb_forests_euc$cluster, bray$cluster)

kNNdistplot(plants_forests_bray, k = 10)
abline(h = quantile(plants_forests_bray, 0.05), col = "red")

# further reduction of plant data
plants_occurences_forests <- plant_occurences(plants_forests)

# remove plants which only occured once or twice
plants_forests_zero <- filter(plants_occurences_forests, !total %in% c(0,1,2))
filter(plants_occurences_forests, total == 2)

plants_forest_short <- plants_forests[c("Polygon",plants_forests_zero$species)] # take all 

# add life forms
plants_occurences$short <- sub("(\\w+\\s+\\w+).*", "\\1", plants_occurences$species)
plants_occurences$short <- sub("(\\w+).*", "\\1", plants_occurences$species)
life_forms$short <- sub("(\\w+).*", "\\1", life_forms$FloraVeg.Taxon)

# add life form to plants -------------------------------------------------

plants_LF <- left_join(plants_occurences,life_forms, by = "short",multiple = "any")

trees <- filter(plants_LF, Tree == 1)
trees_short <- filter(trees, species %in% colnames(plants_forest_short))

check_empty_plot <- function(plant_data){
  plot <- rowSums(plants_data[,-1])
  plot[order(plot)]
}

check_empty_plot(plants_forest_short)

plants_forest_short <- plants_forest_short[c("Polygon",trees_short$species)]

check_empty_plot <- rowSums(plants_forest_short[-1])
check_empty_plot[order(check_empty_plot)]

# remove "AT... and "AU..."
forests_red_help <- substr(grunddaten_forests$`Biotoptyp-Land`,1,2)
forests_red_help <- which(forests_red_help %in% c("AT", "AU"))
grunddaten_forests_red <- grunddaten_forests[-forests_red_help,]
filter(grunddaten_forests_red, substr(`Biotoptyp-Land`,1,2) %in% c("AT", "AU"))
plants_forests_red <- plants_forest_short[-forests_red_help,]

check_empty_plot <- rowSums(plants_forests_red[,-1])
check_empty_plot[order(check_empty_plot)]
plants_forests_red[order(check_empty_plot),]
grunddaten_forests_red[order(check_empty_plot),c(1,4)]

# check clustering
plants_forests_w_short <- plant_weighting(plants_forest_short)
plants_forests_rel_short <- decostand(plants_forests_w_short, method = "total")

plants_forests_bray_short <- vegdist(plants_forests_rel_short, method = "bray")
hdbscan_minClusSize(plants_forests_bray_short)
hdbscan_evaluation(plants_forests_bray_short, k = 7)

# Hierarchical clustering
hc_forests <- hclust(plants_forests_bray, method = "ward.D2")
plot(hc_forests)
clusters <- cutree(hc_forests, k = 5) # ???
plot(clusters)

# Evaluation of cluster number
library(factoextra)
fviz_nbclust(plants_forests_rel, FUN = hcut, method = "silhouette", k.max = 25)
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
