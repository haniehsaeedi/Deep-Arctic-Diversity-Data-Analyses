# extract the data from GBIF using spocc package

library("spocc")
out <- occ(query = c('Animalia'), geometry='POLYGON((180 66, 180 90, -180 90, -180 66, 180 66))', gbifopts=list(hasCoordinate=TRUE, depth='500,11000'), from = c('gbif'),limit = 1000000)
data_GBIF <- occ2df(out$gbif)
data_GBIF <- apply(data_GBIF,2,as.character)
write.csv(data_GBIF, 'data_GBIF.csv')

# extract the data from OBIS using rOBIS package
library("dplyr") 
library("robis") 
df <- robis::occurrence(c('Animalia'), geometry='POLYGON((180 66, 180 90, -180 90, -180 66, 180 66))', startdepth = 500, enddepth = 11000)
write.csv(df,'data_OBIS.csv')

# add the dataset integrated from authours
data_integrated <- read.csv("Integrated_dataset_cleaned.csv")

#Load data and do some trimming
#OBIS
data_OBIS <- read.csv("data_OBIS.csv")
data_OBIS_Trim <- subset(data_OBIS, terrestrial!="TRUE")
data_OBIS_Trim <- subset(data_OBIS_Trim, absence!="TRUE")
data_OBIS_Trim <- subset(data_OBIS_Trim, decimalLatitude!="NA", decimalLongitude!="NA")
data_OBIS_Trim <- subset(data_OBIS_Trim, basisOfRecord!="FOSSIL_SPECIMEN")
data_OBIS_Trim <- subset(data_OBIS_Trim, coordinateUncertaintyInMeters<= 100000 
                          | is.na(coordinateUncertaintyInMeters)) # remove rows with "Coordinate Uncertainty" exceeding 100 km 
# without deleting rows with "NA"
write.csv(data_GBIF, 'data_GBIF.csv')

#GBIF
data_GBIF <- read.csv("data_GBIF.csv")
data_GBIF_Trim <- subset(data_GBIF, occurrenceStatus!="ABSENT")
data_GBIF_Trim <- subset(data_GBIF_Trim, decimalLatitude!="NA", decimalLongitude!="NA")
data_GBIF_Trim <- subset(data_GBIF_Trim, basisOfRecord!="FOSSIL_SPECIMEN")
data_GBIF_Trim <- subset(data_GBIF_Trim, coordinateUncertaintyInMeters<= 100000 
                          | is.na(coordinateUncertaintyInMeters)) 
#integrated dataset
data_integrated_Trim <- subset(data_integrated, occurrenceStatus!="ABSENT")
data_integrated_Trim <- subset(data_integrated, decimalLatitude!="NA", decimalLongitude!="NA")
data_integrated_Trim <- subset(data_integrated, basisOfRecord!="FOSSIL_SPECIMEN")

# Filter OBIS data
data_OBIS_fil <- data_OBIS_Trim %>%
  dplyr::select(scientificName, occurrenceID, taxonRank, eventID, decimalLatitude, decimalLongitude, depth, occurrenceStatus, basisOfRecord) %>%
  mutate(
    decimalLongitude = round(decimalLongitude, 3),
    decimalLatitude = round(decimalLatitude, 3)
  )

data_GBIF_fil <- data_GBIF_Trim %>%
  dplyr::select(scientificName, occurrenceID, taxonRank, eventID, decimalLatitude, decimalLongitude, depth, occurrenceStatus, basisOfRecord) %>%
  mutate(
    longitude = round(decimalLongitude, 3),
    latitude = round(decimalLatitude, 3)
  )

data_integrated_fil <- data_integrated_Trim %>%
  dplyr::select(scientificName, occurrenceID, taxonRank, eventID, decimalLatitude, decimalLongitude, depth, occurrenceStatus, basisOfRecord) %>%
  mutate(
    longitude = round(decimalLongitude, 3),
    latitude = round(decimalLatitude, 3)
  )

# Merge OBIS, GBIF, and integrated datasets and remove the duplicates by function distinct 
data_merged <- bind_rows(data_GBIF_fil, data_OBIS_fil, data_integrated_fil) %>%
  distinct()
write.csv(data_merged, 'data_merged.csv')

# data cleaning
# data cleaning using obistools
#install.packages("devtools")
#devtools::install_github("iobis/obistools")
library("obistools")

names <- (data_merged$scientificName)
match_taxa(names)# all taxa matched and all matched with worms, the one with multiple selections were resolved by choosing the accepted taxa 
# only Aceros was not accepted and was removed from the data_merged and BOLD data records

#load the taxonmatched data
data_merged_taxmatch <- read.csv("data_merged_taxmatch.csv")

#Plot points on a map
plot_map(data_merged_taxmatch, zoom = TRUE)

#Check points on land
check_onland(data_merged_taxmatch)
check_onland(data_merged_taxmatch, report = TRUE)# it shows there are 273 points on land

#remove points on land using scrubr and other packages
library(xlsx)
library(glue)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(CoordinateCleaner)
library(raster)
library(maptools)
library(rgdal)
library(dismo)
library(scrubr)
land <- check_onland(data_merged_taxmatch)
land_buffer <- check_onland(data_merged_taxmatch, buffer = 1000)
world <- map_data("world")

ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "#dddddd") +
  geom_point(data = land, aes(x = decimalLongitude, y = decimalLatitude), color = "#cc3300") +
  geom_point(data = land_buffer, aes(x = decimalLongitude, y = decimalLatitude), color = "#ff9900") + coord_fixed(1)

#Remove the points on land
data_merged_clean <- anti_join(data_merged_taxmatch, land_buffer)

#Check depth using using obistools
plot_map(check_depth(data_merged_clean, depthmargin = 500), zoom = TRUE)
report <- check_depth(data_merged_clean, report=T, depthmargin = 500)
head(report)# as only max and minimumDepthInMeters are missing, no action is needed 

#report() generates a basic data quality report using obistools
report(data_merged_clean)

#data_merged_clean <- read.csv("data_merged_clean.csv")

#assign taxa to benthic/pelagic
taxon_list_original <- read.csv("Benthic_Pelagic_specieslist.csv")

#taxon_list <- taxon_list_original %>%
  #filter(!is.na(ScientificName_accepted) & ScientificName_accepted != "" & !is.na(Species) & Species != "")

#Join the occurrence data with the species list

data_joined <- data_merged_clean %>%
  left_join(taxon_list_original, by = c('scientificName'))

#Some statistics
stats <- data_joined %>%
  group_by(Group) %>%
  summarize(records = n())
stats

Benthic <- filter(data_joined, Group == 'Benthic')
write.csv(Benthic, 'Benthic.csv')

#map the occurrence points on the Arctic map
#devtools::install_github("MikkoVihtakari/PlotSvalbard") ## Run once
library(PlotSvalbard)
library(sf)
library(tidyverse)
basemap("panarctic", limits = 60)
basemap("panarctic", limits = 60, bathymetry = TRUE) 

# here I had to use QGIS to remove the data in depths less than 500m 
# because the checkdepth function in robis did not distinguish the ones in depth
# less than 500m, the new dataset names is Benthic_Final

# read the data
data<-read.csv('Benthic_Final_With_Habitats.csv')

#install.packages("ggOceanMapsData", repos = c("https://mikkovihtakari.github.io/drat", "https://cloud.r-project.org"))
library(ggOceanMapsData)
library(ggOceanMaps)
  
basemap(limits = 60, glaciers = TRUE, bathymetry = TRUE) + 
  geom_spatial_point(data = data, aes(x = decimalLongitude, y = decimalLatitude), color = "darkorange") 

world <- map_data("world")
worldmap <- ggplot(world, aes(long, lat, group = group)) +
  geom_path() +
  scale_y_continuous(NULL, breaks = (-2:3) * 30, labels = NULL) +
  scale_x_continuous(NULL, breaks = (-4:4) * 45, labels = NULL)

# map the sequence data
data<-read.csv('All_Arctic_deep_sequences.csv')

basemap(limits = 60, glaciers = TRUE, bathymetry = TRUE) + 
  geom_spatial_point(data = data, aes(x = lon, y = lat), color = "darkorange") 


# Some crazier projections
worldmap + coord_map("ortho") 

# Biodiversity analyses
# Import and Formatting

data$depth <- data$depth * -1
data <- mutate(data, decimalLatitude_5 = ceiling(decimalLatitude/5) * 5)
data$decimalLatitude_5 <- formatC(data$decimalLatitude_5, width = 2, format = "d", flag = "0")
data <- mutate(data, depth_rnd = floor(depth/100) * 100)

# Plot (records) latitude / longitude

ggplot(data, aes(x = decimalLongitude, y = decimalLatitude)) +
  geom_jitter(alpha = 0.5) +
  xlab("Longitude [°]") + ylab("Latitude [°]")

# Plot (records) latitude against depth

ggplot(data, aes(x = decimalLatitude, y = depth)) +
  geom_jitter(alpha = 0.5) +
  xlab("Latitude [°]") + ylab("Depth [m]")

# Plot records per latitude

ggplot(data, aes(x = decimalLatitude_5)) + 
  geom_bar(position="dodge") +
  xlab("Latitude [°]") + ylab("Number of Records")

# Plot records per depth

ggplot(data, aes(x = depth_rnd)) + 
  geom_bar(position="dodge") +
  xlab("Depth [m]") + ylab("Number of Records") +
  coord_flip()

# Plot species per latitude

distinct(data, scientificName, decimalLatitude_5, .keep_all = TRUE) %>% 
  ggplot(aes(x = decimalLatitude_5)) + 
  geom_bar(position="dodge") +
  xlab("Latitude [°]") + ylab("Number of Species")

# Plot species per depth

distinct(data, scientificName, depth_rnd, .keep_all = TRUE) %>% 
  ggplot(aes(x = depth_rnd)) + 
  geom_bar(position="dodge") +
  xlab("Depth [m]") + ylab("Number of Species") +
  coord_flip()


# Matrix Function
abu_matrix <- function(data, sites.col, sp.col, keep.n = TRUE) {
  
  stopifnot(
    length(sites.col) == 1,
    length(sp.col) == 1,
    sites.col != sp.col,
    sites.col %in% 1 : ncol(data) | sites.col %in% names(data),
    sp.col %in% 1 : ncol(data) | sp.col %in% names(data),
    is.logical(keep.n)
  )
  
  presabs <- table(data[ , c(sites.col, sp.col)])
  presabs <- as.data.frame(unclass(presabs))
  if (!keep.n)  presabs[presabs > 1] <- 1
  presabs <- data.frame(row.names(presabs), presabs)
  names(presabs)[1] <- names(subset(data, select = sites.col))
  rownames(presabs) <- NULL
  return(presabs)
}

# Matrix
Benthic_lat <- abu_matrix(subset(data), "decimalLatitude_5", "scientificName") %>% 
  column_to_rownames(var = "decimalLatitude_5")

library(readxl)
library(openxlsx)
library(worms)
library(tidyverse)
library(sf)
library(vegan)
library(pvclust)

# Rarefaction

rarecurve(Benthic_lat, step = 20, sample = 50, col = "blue", cex = 0.6, label = TRUE)

# dropped stations with 0 records

# Rarefaction plots

Benthic_lat %>% rarefy(sample = 50, MARGIN = 1) %>% as.data.frame() %>% rownames_to_column("Latitude") %>% 
  ggplot(aes(x = Latitude, y = .)) + 
  geom_col(position="dodge") +
  xlab("Latitude [°]") + ylab("ES50") 


# Shapefile Data

land <- st_read("ne_110m_admin_0_countries.shp")
hex3 <- st_read("hexgrid3.shp")

Benthic_sf <- data[, c("scientificName", "decimalLatitude", "decimalLongitude")]
Benthic_sf <- Benthic_sf %>% st_as_sf(coords = c('decimalLongitude','decimalLatitude'))
st_crs(Benthic_sf) = 4326


# Spatial Function

spatial_fun <- function(land, hex3, Benthic_sf) {
  
  spatial_data <- st_join(hex3, Benthic_sf, join = st_contains) %>% 
    group_by(ID) %>% 
    summarise(Effort = length(scientificName[!is.na(scientificName)]), 
              Richness = n_distinct(scientificName[!is.na(scientificName)]))
  spatial_data <- spatial_data %>% as.data.frame()
  row_sub <- apply(spatial_data[,c(2,3)], 1, function(row) all(row != 0))
  spatial_data <- spatial_data[row_sub,]
  
  spatial_data2 <- st_join(hex3, Benthic_sf, join = st_intersects) %>% 
    as.data.frame() %>% 
    abu_matrix("ID", "scientificName")
  spatial_data_r <- rarefy(spatial_data2[,-1], sample = 20, MARGIN = 1) %>% 
    as.data.frame()
  colnames(spatial_data_r) <- "ES15"
  spatial_data_r$ID <- spatial_data2$ID
  spatial_data <- spatial_data %>% 
    merge(spatial_data_r, by = "ID", all.x = TRUE) %>% 
    st_as_sf()
  
  p1 <- ggplot() +
    geom_sf(data = spatial_data, aes(fill = Effort)) +
    scale_fill_gradient(low = "dodgerblue", high = "royalblue4") +
    geom_sf(data = sf_land, fill = "light grey") +
    coord_sf(xlim = c(-180, 180), ylim = c(0, 90)) +
    theme(legend.position = c(0.15, 0.6))
  
  p2 <- ggplot() +
    geom_sf(data = spatial_data, aes(fill = Richness)) +
    scale_fill_gradient(low = "dodgerblue", high = "royalblue4") +
    geom_sf(data = sf_land, fill = "light grey") +
    coord_sf(xlim = c(-180, 180), ylim = c(0, 90)) +
    theme(legend.position = c(0.15, 0.6))
  
  p3 <- ggplot() +
    geom_sf(data = subset(spatial_data, Effort >= 15), aes(fill = ES15)) +
    scale_fill_gradient(low = "dodgerblue", high = "royalblue4") +
    geom_sf(data = sf_land, fill = "light grey") +
    coord_sf(xlim = c(-180, 180), ylim = c(0, 90)) +
    theme(legend.position = c(0.15, 0.6))
  
  plot_grid(p1, p2, p3, align = "hv", nrow = 1, labels = c("(a)", "(b)", "(c)"))
  
  ggsave("Fig. S19.png", units = "cm", width = 30, height = 15, dpi = 300)
  
}

# Hexagons

spatial_fun(land, hex3, Benthic_sf)

# Ecoregions

eco <- st_read("meow_ecos.shp")
eco <- as.data.frame(eco)
colnames(eco)[1] <- "ID"
eco <- st_as_sf(eco)

spatial_fun(land, eco, Benthic_sf)


#Get citations for OBIS data
# get citations

datasetids <- unique(data$dataset_id)

citations <- data.frame(id = datasetids, citation = NA)

for (i in 1:nrow(citations)) {
  message(citations$id[i])
  d <- dataset(datasetid = citations$id[i])
  citations$citation[i] <- d$citation
}

# output
write.csv(citations, 'citations.csv')