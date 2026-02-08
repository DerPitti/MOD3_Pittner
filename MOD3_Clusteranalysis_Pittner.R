library(readxl)
library(dplyr)
library(ggplot2)
library(future.apply)

grunddaten <- read_xlsx("data/tbl_grunddaten.xlsx") # load basic data for each plot, e.g. biotopcodes
plants <- read_xlsx("data/tbl_daten_pflanzen.xlsx") # load plant data

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
library(dplyr)
remove_polygons <- filter(grunddaten, substr(`Biotoptyp-Land`,1,1) %in% c("F", "V"))

plants_sub <- filter(plants, !Polygon %in% remove_polygons$Polygon)
grunddaten_sub <- filter(grunddaten, !Polygon %in% remove_polygons$Polygon)


# Hellinger transformation ------------------------------------------------

# transformation of plant data frame
library(tidyr)

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

plants_clean2$Menge <- as.numeric(plants_clean2$Menge)

plants_wide <- pivot_wider(plants_clean2[,c(1,2,4)],names_from = `Wissenschaftlicher Name`,values_from = Menge)
plants_wide[is.na(plants_wide)] <- 0

# check that every plot has plants
summary(rowSums(plants_wide[,-1]))
plants_sum <- cbind(rowSums(plants_wide[,-1]),plants_wide)
order(rowSums(plants_wide[,-1]))

### find outlier plots

# Community abundance-weighted transformation
plants_weighted <- plants_wide[,-1]
plants_weighted[plants_weighted == 1] <- 0.01
plants_weighted[plants_weighted == 2] <- 0.05
plants_weighted[plants_weighted == 3] <- 0.25
plants_weighted[plants_weighted == 4] <- 0.75

bray_curtis_dist_weighted <- vegdist(plants_weighted, method = "bray")

options(future.validate = FALSE)
plan(multisession, workers = parallel::detectCores() - 2) # start multisession, use all cores except two

ord_new <- metaMDS(bray_curtis_dist_weighted, k = 2)
plot(ord_new, type = "t")
stressplot(ord_new)


which(is.finite(ord$points[,1]) & 
        abs(scale(ord$points[,1])) > 3 |
        abs(scale(ord$points[,2])) > 3)

ord$points[order(ord$points[,1], decreasing = TRUE),] # 711, 1033, 1901 seem to be outlier


filter(plants_clean2, Polygon %in% plants_wide$Polygon[order(ord$points[,1], decreasing = TRUE)[1:3]])

# find all plants with only few entries
plants_occurences <- plants_wide[,-1] %>%
  mutate(across(everything(), ~ if_else(. > 0, 1L, 0L)))%>%
  summarise(across(everything(), sum, na.rm=TRUE))%>%
  pivot_longer(everything(), names_to = "species", values_to = "total") %>%
  arrange(total)

filter(plants_occurences, species %in% c("Ulmus spec.","Robinia pseudoacacia"))

suspicious_polygons <- filter(grunddaten_sub, Polygon %in% plants_wide$Polygon[order(ord$points[,1], decreasing = TRUE)[1:3]])
# 01728 will be removed, because only one plant found, and in 04567, Ulmus spec. will be specified to Ulmus glabra accordingly to Beschreibung

plants_clean2 <- plants_clean2 %>%
  filter(Polygon != "01728") %>%
  mutate(`Wissenschaftlicher Name` = if_else(`Wissenschaftlicher Name`=="Ulmus spec.", "Ulmus glabra", `Wissenschaftlicher Name`))

plants_wide <- pivot_wider(plants_clean2[,c(1,2,4)],names_from = `Wissenschaftlicher Name`,values_from = Menge)
plants_wide[is.na(plants_wide)] <- 0

grunddaten_sub <- grunddaten_sub %>%
  filter(Polygon != "01728")

# Hellinger transformation to subsequently compute Euclidean distance
plants_hell <- hellinger(plants_wide[,-1])
test_dist <- dist(plants_hell)

#Alternative...
comm_hel <- decostand(plants_wide[,-1], method = "hellinger")
dist_hel <- dist(comm_hel, method = "euclidean")

# Bray-Curtis dissimilarity with raw abundance data
library(vegan)
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


library("dbscan")
test_hdb <- hdbscan(bray_curtis_dist, minPts = 10)
table(test_hdb$cluster)

for (k in seq(5, 15,by = 2)) {
  h <- hdbscan(bray_curtis_dist_weighted, minPts = k)
  cat("minPts =", k, "-> clusters:", length(unique(h$cluster)), 
      " noise:", sum(h$cluster == 0), "\n")
}

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


coords <- as.data.frame(scores(ord_new))
coords$cluster <- factor(test_hdb_bcd_weighted$cluster)

test_cluster_bt <- data.frame(cbind(cluster=test_hdb_bcd_weighted$cluster, number= test_hdb_bcd_weighted$hc$order))
test_cluster_bt <- test_cluster_bt%>%
  arrange(number)
test_cluster_bt <- cbind(test_cluster_bt, Polygon=plants_wide$Polygon)
join_cluster_code <- inner_join(test_cluster_bt, grunddaten_sub[,c(1,2,4)], by = "Polygon")

ggplot(coords, aes(NMDS1, NMDS2, color = cluster)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(
    values = c("0" = "grey70", scales::hue_pal()(length(unique(coords$cluster)) - 1))
  ) +
  labs(
    title = "HDBSCAN Clustering (PCoA of Bray–Curtis)",
    x = "PCoA 1",
    y = "PCoA 2",
    color = "Cluster"
  ) +
  theme_minimal()
