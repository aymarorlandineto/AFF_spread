##############################
# Comparing distributions outside the main native tropical bioregions
##############################

# Load libraries
library(sf)
library(rnaturalearth)
library(dplyr)
library(purrr)
library(ggplot2)
library(trend)
library(stringr)
library(stringi)
library(tidyr)

# -----------------------------
# Load and prepare data
# -----------------------------
setwd("C:/Users/aymar/Downloads")
#amazon native sp occ download GBIF https://doi.org/10.15468/dl.g46bn8
dados_amazon<- read.csv("0029395-250525065834625.csv",sep="\t") # amazon native sp list from amazon_list
#mekong native sp occ download GBIF https://doi.org/10.15468/dl.d6k62t
dados_mekong<- read.csv("0029380-250525065834625.csv",sep="\t") # mekong native sp list from tedesco 2017
#congo native sp occ download GBIF https://doi.org/10.15468/dl.56daff
dados_congo<- read.csv("0029374-250525065834625.csv",sep="\t") # congo native sp list from tedesco 2017

dados_congo <- dados_congo %>% mutate(bacia = "Congo")
dados_mekong <- dados_mekong %>% mutate(bacia = "Mekong")
dados_amazon <- dados_amazon %>% mutate(bacia = "Amazon")

dados_todas_bacias <- bind_rows(dados_congo, dados_mekong,dados_amazon)

GBIF_sem_duplicados <- dados_todas_bacias %>% 
  distinct(species, decimalLatitude, decimalLongitude, year, .keep_all = TRUE)
GBIF_sem_duplicados$scientificName <- str_replace_all(GBIF_sem_duplicados$scientificName,fixed(" "), ".")
GBIF_data_filtered <- GBIF_sem_duplicados %>%
  filter(!is.na(species),
         !is.na(decimalLatitude),
         !is.na(decimalLongitude),
         !is.na(year)) # remove NA and geographic anomalies

GBIF_data_filtered <- st_as_sf(GBIF_data_filtered, 
                               coords = c("decimalLongitude", "decimalLatitude"), 
                               crs = 4326,remove = FALSE) 
ocorrencias <- GBIF_data_filtered
ocorrencias
ocorrencias <- ocorrencias %>% select(species, countryCode,
                                      decimalLatitude,decimalLongitude,
                                      year,bacia,geometry)

ocorrencias<-ocorrencias %>% 
  filter(year < 2025)
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  select(continent, subregion, region_wb, name) %>%
  st_transform(4326)
occurrences <- st_transform(ocorrencias, st_crs(world))

# Tropical bioregions shapefile
# (should contain the union of Neotropical, Ethiopian, and Sino-Oriental)
# Load Leroy freshwater regions (Neotropical = South America,Caribbean and Central America (except North of Mexico),
# Sino-Orirental = Southern Asia, South-Eastern Asia  and Eastern Asia) and Ethiopian (African continent plus parts of Middle East)
bioregions_union_tropical <- st_read("tropical_realms.shp")

ggplot() +
  geom_sf(data = world, color = "black", alpha = 0.5) +
  geom_sf(data = bioregions_union_tropical, aes(color = bioregion), fill = NA, linewidth = 0.8) +
  theme_minimal() +
  labs(title = "Tropical Bioregions")
# Split regions of interest
neotropical   <- bioregions_union_tropical %>% filter(bioregion == "Neotropical")
ethiopian     <- bioregions_union_tropical %>% filter(bioregion == "Ethiopian")
sino_oriental <- bioregions_union_tropical %>% filter(bioregion == "Sino-Oriental")

# Split basins
occ_amazon <- occurrences %>% filter(bacia == "Amazon")
occ_congo  <- occurrences %>% filter(bacia == "Congo")
occ_mekong <- occurrences %>% filter(bacia == "Mekong")

# -----------------------------
# Filter OUTSIDE native bioregions
# -----------------------------
# Amazon outside Neotropical
amazon_outside <- occ_amazon %>%
  filter(lengths(st_intersects(occ_amazon, neotropical, sparse = TRUE)) == 0)

# Congo outside Ethiopian
congo_outside <- occ_congo %>%
  filter(lengths(st_intersects(occ_congo, ethiopian, sparse = TRUE)) == 0)

# Mekong outside Sino-Oriental
mekong_outside <- occ_mekong %>%
  filter(lengths(st_intersects(occ_mekong, sino_oriental, sparse = TRUE)) == 0)

# Join results
occurrences_filtered <- bind_rows(
  amazon_outside,
  congo_outside,
  mekong_outside)
st_write(occurrences_filtered, "amazon_congo_mekong_occ.shp")

occurrences_filtered<- st_read("amazon_congo_mekong_occ.shp")
occurrences_filtered %>%
  st_drop_geometry() %>%
  count(bacia)

occurrences_filtered %>%
  st_drop_geometry() %>%
  distinct(bacia, species) %>%
  count(bacia, name = "n_species")

# Temporal trends data
occurrences_temporal <- occurrences_filtered %>%
  st_drop_geometry() %>%
  group_by(bacia, year) %>%
  summarise(total_occurrences = n(), .groups = "drop") %>%
  filter(year < 2025)

# -----------------------------
# Plot temporal trends
# -----------------------------
base_theme <- theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

temporal_plot <- ggplot(occurrences_temporal,
                        aes(x = year, y = total_occurrences,
                            color = bacia, shape = bacia, linetype = bacia)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#66CC99", "#FFCC33", "#FF9966")) +
  labs(x = "Year", y = "Occurrences",
       color = "Origin Basin", shape = "Origin Basin", linetype = "Origin Basin") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  base_theme

print(temporal_plot)

# -----------------------------
# Temporal Trend tests (Mann-Kendall, Sen’s slope)
# -----------------------------
occurrences_temporal$basin <- as.factor(occurrences_temporal$bacia)
occurrences_temporal <- occurrences_temporal %>% filter(year > 1999)

# Amazon
amazon_trend <- occurrences_temporal %>% filter(bacia == "Amazon")
mk.test(amazon_trend$total_occurrences)
sens.slope(amazon_trend$total_occurrences, conf.level = 0.95)
snh.test(amazon_trend$total_occurrences)

# Congo
congo_trend <- occurrences_temporal %>% filter(bacia == "Congo")
mk.test(congo_trend$total_occurrences)
sens.slope(congo_trend$total_occurrences, conf.level = 0.95)
snh.test(congo_trend$total_occurrences)

# Mekong
mekong_trend <- occurrences_temporal %>% filter(bacia == "Mekong")
mk.test(mekong_trend$total_occurrences)
sens.slope(mekong_trend$total_occurrences, conf.level = 0.95)
snh.test(mekong_trend$total_occurrences)

# Summarize trends for all basins
trends_summary <- occurrences_temporal %>%
  group_by(basin) %>%
  summarise(
    statistic = mk.test(total_occurrences)$statistic,
    p_value = mk.test(total_occurrences)$p.value,
    significance = ifelse(p_value < 0.05, "Significant", "Not significant"),
    magnitude = sens.slope(total_occurrences)$estimates
  )

print(trends_summary)


#############################################
# AMAZONIAN SPECIES RICHNESS AND DISTRIBUTION OUTSIDE THEIR BASIN
#############################################

# Load required packages
library(sf)
library(rnaturalearth)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)
library(countrycode)

# -----------------------------
# Load data
# -----------------------------

# Amazon species list
amazon_list<- read.csv("lista_amazon.csv", 
                       header = T,sep=";")

GBIF<- readr::read_delim("0023879-240906103802322.csv", delim = "\t") #19/09/2024 + atual
species_to_remove <- c("Strongylura marina", "Eleotris pisonis", "Dormitator maculatus")
GBIF <- GBIF[!GBIF$species %in% species_to_remove, ] 

# -----------------------------
# Load Amazon Basin shapefile
# -----------------------------
amazon_watersheds<- st_read("amazon_basin.shp")
amazon_watersheds<-st_transform(amazon_watersheds,4326)
st_is_valid(amazon_watersheds)
amazon_watersheds<-st_make_valid(amazon_watersheds)
amazon_watersheds<-amazon_watersheds %>%filter(BL1 == "Amazon basin")

# World shapefile
land <- ne_countries(scale = "medium", returnclass = "sf")
if (st_crs(amazon_watersheds) != st_crs(land)) {
  land <- st_transform(land, st_crs(amazon_watersheds))
}

# -----------------------------
# Clean and filter GBIF data
# -----------------------------
GBIF_data_filtered <- GBIF %>%
  filter(!is.na(year), !is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
  mutate(year = as.integer(year)) %>%
  filter( year >= 1980)

GBIF_data_filtered <- st_as_sf(GBIF_data_filtered, 
                               coords = c("decimalLongitude", "decimalLatitude"), 
                               crs = 4326,remove = FALSE) 
GBIF_data_sf<-GBIF_data_filtered
names(GBIF_data_sf)

# Remove points inside Amazon Basin
intersects_amazon <- st_intersects(GBIF_data_sf, amazon_watersheds, sparse = FALSE)
GBIF_data_filtered_outside_amazon <- GBIF_data_sf[!apply(intersects_amazon, 1, any), ]
unique(GBIF_data_filtered_outside_amazon$basisOfRecord)

intersects_land <- st_intersects(GBIF_data_filtered_outside_amazon, land, sparse = FALSE)
GBIF_data_filtered_outside_amazon <- GBIF_data_filtered_outside_amazon[apply(intersects_land, 1, any), ]

GBIF_data_filtered_outside_amazon<-GBIF_data_filtered_outside_amazon %>% 
  dplyr::select(gbifID,occurrenceID,order,
                family,genus,species,scientificName,verbatimScientificNameAuthorship,
                countryCode,individualCount,decimalLongitude,decimalLatitude,
                eventDate,day,month,year,
                identifiedBy,dateIdentified,
                geometry)
GBIF_data_filtered_outside_amazon<-GBIF_data_filtered_outside_amazon[!duplicated(GBIF_data_filtered_outside_amazon[,c(1,6,11,12,13)]),]
nrow(GBIF_data_filtered_outside_amazon)
sum(is.na(GBIF_data_filtered_outside_amazon$species))
GBIF_data_filtered_outside_amazon <- drop_na(GBIF_data_filtered_outside_amazon, species)
gbif_filtrado<-GBIF_data_filtered_outside_amazon

gbif_filtrado$country_name <- countrycode::countrycode(gbif_filtrado$countryCode, origin = 'iso2c', destination = 'country.name')
sum(is.na(gbif_filtrado$country_name))

# -----------------------------
# Removing occ where sp are native
# Load Tedesco data
# -----------------------------
Tedesco_ref<- read.csv("References_Table.csv",header = T,sep=";")
Tedesco_ref$Year <- str_extract(Tedesco_ref$X3.Reference, "\\d{4}")
Tedesco<- read.csv("Occurrence_Table_12092019.csv",header = T,sep=",")
Tedesco$X2.Species.Name.in.Source
Tedesco$X1.Basin.Name
Tedesco <- Tedesco %>%
  rename(Taxon = X6.Fishbase.Valid.Species.Name)
Tedesco$Taxon<-str_replace_all(Tedesco$Taxon,
                               fixed("."), " ")
amazon_list$Taxon<-str_replace_all(amazon_list$x,
                                   fixed("."), " ")
Tedesco<- Tedesco %>%
  mutate(Origin = if_else(Taxon %in% amazon_list$Taxon, "Amazon", "Other exotics"))

world <- ne_countries(scale = 'medium')
world <- world %>% dplyr::select(iso_a2,name)
world <-st_transform(world,4326)
colnames(world)[1] <- "countryCode"
colnames(world)[2] <- "country_name"

Tedesco <-Tedesco %>% filter(.,X3.Native.Exotic.Status == "native")
Tedesco <-Tedesco %>% filter(.,Tedesco$X7.Occurrence.Status == "valid")
Tedesco <-Tedesco %>% filter(.,X3.Ecoregion == "Neotropic")
Tedesco <-Tedesco %>% filter(.,Origin == "Amazon")
Tedesco <-Tedesco %>% filter(.,Freshwater != "0")
Tedesco <-Tedesco %>% filter(.,Saltwater != "-1")
Tedesco <-Tedesco %>% filter(.,Brackish != "-1")

Tedesco <- Tedesco[!Tedesco$Taxon %in% species_to_remove, ]
colnames(Tedesco$X1.Basin.Name)

colnames(Tedesco)[13] <- "country_name"
unique(Tedesco$country_name)

# -----------------------------
# Add basins to GBIF occurrences
# -----------------------------
#download freshwater bioregions per basin from Leroy et al 2019
bacias_hidro <- st_read("fish_bioregions_per_basin.gpkg")
bacias_hidro<-st_transform(bacias_hidro,4326)
bacias_hidro <- st_make_valid(bacias_hidro)
ocorrencias <- st_transform(gbif_filtrado, st_crs(bacias_hidro))

gbif_filtrado <- st_join(ocorrencias, bacias_hidro["BasinName"])
Tedesco <-Tedesco %>% filter(.,X7.Occurrence.Status != "questionable")
unique(Tedesco$X7.Occurrence.Status) #keep only "valid"
names(Tedesco)
colnames(Tedesco)[1]<-"BasinName"
colnames(Tedesco)[6]<-"species"

gbif_filtrado <- st_transform(gbif_filtrado, st_crs(bacias_hidro))

non_native_occurrences <- gbif_filtrado %>%
  anti_join(Tedesco, by = c("species", "BasinName"))
non_native_occurrences <- st_transform(non_native_occurrences, st_crs(amazon_watersheds))

# Remove occ inside Amazon basin and add BOL occur
##correction for some countries polygon
non_native_occurrences <- st_filter(non_native_occurrences, amazon_watersheds, .predicate = st_disjoint)
bo_occurrences <- gbif_filtrado %>%
  filter(countryCode == "BO")
non_native_occurrences <- bind_rows(non_native_occurrences, bo_occurrences)
nrow(non_native_occurrences)

# -----------------------------
# Temporal aggregation by period
# -----------------------------
non_native_occurrences <- non_native_occurrences %>%
  mutate(period = cut(year,
                      breaks = c(1980, 1991, 2002, 2013, 2025),
                      labels = c("1980–1990", "1991–2001", "2002–2012", "2013–2024"),
                      right = FALSE)) %>%
  group_by(period, decimalLongitude, decimalLatitude) %>%
  summarise(n = n(), .groups = "drop") %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
  mutate(n_class = case_when(
    n <= 10 ~ "1-10",
    n <= 50 ~ "11-50",
    n <= 200 ~ "51-200",
    n <= 1000 ~ "201-1000",
    TRUE ~ "1000+"
  ),
  n_class = factor(n_class, levels = c("1-10", "11-50", "51-200", "201-1000", "1000+"))
  )

# -----------------------------
# Maps
# -----------------------------
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(continent != "Antarctica") %>%
  st_transform(4326)

# Global map
map_global <- ggplot() +
  geom_sf(data = world, fill = "#f0f0f0", color = "grey70", linewidth = 0.2) +
  geom_sf(data = amazon_watersheds, fill = "#306844", alpha = 0.5, color = NA) +
  geom_sf(data = non_native_occurrences, alpha = 0.5, col = "#de7065") +
  facet_wrap(~period) +
  theme_void()


# South America zoom
south_america <- world %>% filter(continent == "South America")

map_SA <- ggplot() +
  geom_sf(data = south_america, fill = "#f0f0f0", color = "grey70", linewidth = 0.2) +
  geom_sf(data = amazon_watersheds, fill = "#306844", alpha = 0.5, color = NA) +
  geom_sf(data = non_native_occurrences, alpha = 0.5, col = "#de7065") +
  facet_wrap(~period, nrow = 1) +
  coord_sf(xlim = c(-100, -30.5), ylim = c(-58, 15), expand = FALSE) +
  theme_void()

# -----------------------------
# Richness and occurrences per country
# -----------------------------
richness_per_country <- non_native_occurrences %>%
  group_by(countryCode) %>%
  summarise(richness = n_distinct(species)) %>%
  arrange(desc(richness))

occurrences_per_country <- non_native_occurrences %>%
  count(countryCode, name = "total_occurrences")

summary_table <- left_join(occurrences_per_country, st_drop_geometry(richness_per_country),
                           by = "countryCode") %>%
  mutate(country_name = countrycode(countryCode, "iso2c", "country.name"))
summary_table$geometry<-NULL
unique(non_native_occurrences$species)
world <- ne_countries(scale = 'medium')
world <- world %>% dplyr::select(iso_a2,name)
world <-st_transform(world,4326)
colnames(world)[1] <- "countryCode"
colnames(world)[2] <- "country_name"
summary_table$country_name <- countrycode::countrycode(summary_table$countryCode, origin = 'iso2c', destination = 'country.name')
summary_table
# -----------------------------
# summary table
# -----------------------------
summary_table




