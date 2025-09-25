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

# World base map
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  select(continent, subregion, region_wb, name) %>%
  st_transform(4326)

# Tropical bioregions shapefile
# (should contain the union of Neotropical, Ethiopian, and Sino-Oriental)
# Load Leroy freshwater regions (Neotropical = South America,Caribbean and Central America (except North of Mexico),
# Sino-Orirental = Southern Asia, South-Eastern Asia  and Eastern Asia) and Ethiopian (African continent plus parts of Middle East)

bioregions_union_tropical <- st_read("data/tropical_bioreg.shp")

# Plot initial map
ggplot() +
  geom_sf(data = world, color = "black", fill = "lightgrey", alpha = 0.5) +
  geom_sf(data = bioregions_union_tropical,
          aes(color = bioregion, fill = bioregion), linewidth = 1) +
  theme_minimal()

# -----------------------------
# Create non-overlapping buffers
# -----------------------------
bioregions_buffer <- bioregions_union_tropical %>%
  st_transform(4087) %>%              
  st_buffer(dist = 100000) %>%  
  st_transform(st_crs(bioregions_union_tropical))

bioregions_buffer <- st_transform(bioregions_buffer, 4326)

priority_order <- c("Neotropical", "Ethiopian", "Sino-Oriental")

final_buffers <- st_sf(bioregion = character(),
                       geometry = st_sfc(),
                       crs = st_crs(bioregions_buffer))

# Apply sequential buffers without overlaps
for (region in priority_order) {
  current_region <- bioregions_buffer %>% filter(bioregion == region)
  
  if (nrow(final_buffers) > 0) {
    current_region <- current_region %>%
      st_difference(st_union(final_buffers))
  }
  
  final_buffers <- bind_rows(final_buffers, current_region)
}

bioregions_buffer <- final_buffers

ggplot() +
  geom_sf(data = world, color = "black", alpha = 0.5) +
  geom_sf(data = bioregions_buffer, aes(color = bioregion), fill = NA, linewidth = 0.8) +
  theme_minimal() +
  labs(title = "Tropical Bioregions")

# -----------------------------
# Load occurrences
# -----------------------------
occurrences <- st_read("data/08_06_2025_filtered.shp") %>%
  filter(year < 2025)

occurrences <- st_transform(occurrences, st_crs(world))

# Split regions of interest
neotropical   <- bioregions_buffer %>% filter(bioregion == "Neotropical")
ethiopian     <- bioregions_buffer %>% filter(bioregion == "Ethiopian")
sino_oriental <- bioregions_buffer %>% filter(bioregion == "Sino-Oriental")

# Split basins
occ_amazon <- occurrences %>% filter(basin == "Amazon")
occ_congo  <- occurrences %>% filter(basin == "Congo")
occ_mekong <- occurrences %>% filter(basin == "Mekong")

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
  mekong_outside
)
st_write(occurences_filtered, "amazon_congo_mekong_occ.shp")
# -----------------------------
# Summary tables
# -----------------------------

occurrences_filtered<-st_read("amazon_congo_mekong_occ.shp")
occurrences_filtered %>%
  st_drop_geometry() %>%
  count(basin)

occurrences_filtered %>%
  st_drop_geometry() %>%
  distinct(basin, species) %>%
  count(basin, name = "n_species")

# Temporal trends data
occurrences_temporal <- occurrences_filtered %>%
  st_drop_geometry() %>%
  group_by(basin, year) %>%
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
                            color = basin, shape = basin, linetype = basin)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#66CC99", "#FFCC33", "#FF9966")) +
  labs(x = "Year", y = "Occurrences",
       color = "Origin Basin", shape = "Origin Basin", linetype = "Origin Basin") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  base_theme

print(temporal_plot)

# Save figure
ggsave("figures/temporal_trends.png", temporal_plot, width = 8, height = 5, dpi = 300)

# -----------------------------
# Trend tests (Mann-Kendall, Sen’s slope)
# -----------------------------
occurrences_temporal$basin <- as.factor(occurrences_temporal$basin)
occurrences_temporal <- occurrences_temporal %>% filter(year > 1999)

# Amazon
amazon_trend <- occurrences_temporal %>% filter(basin == "Amazon")
mk.test(amazon_trend$total_occurrences)
sens.slope(amazon_trend$total_occurrences, conf.level = 0.95)
snh.test(amazon_trend$total_occurrences)

# Congo
congo_trend <- occurrences_temporal %>% filter(basin == "Congo")
mk.test(congo_trend$total_occurrences)
sens.slope(congo_trend$total_occurrences, conf.level = 0.95)
snh.test(congo_trend$total_occurrences)

# Mekong
mekong_trend <- occurrences_temporal %>% filter(basin == "Mekong")
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
amazon_list <- read.csv("data/list_amazon.csv", header = TRUE, sep = ";")

# GBIF data (latest download 19/09/2024)
GBIF <- read_delim("data/GBIF_occurrences.csv", delim = "\t")

# Species to remove (native to Atlantic)
species_to_remove <- c("Strongylura marina", "Eleotris pisonis", "Dormitator maculatus")

# Filter GBIF
GBIF <- GBIF[!GBIF$species %in% species_to_remove, ]

# -----------------------------
# Load Amazon Basin shapefile
# -----------------------------
amazon_watersheds <- st_read("data/amazon_ecoregion.shp") %>%
  st_transform(4326) %>%
  st_make_valid() %>%
  filter(BL1 == "Amazon basin")

# World shapefile
land <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(st_crs(amazon_watersheds))

# -----------------------------
# Clean and filter GBIF data
# -----------------------------
GBIF_data_filtered <- GBIF %>%
  filter(!is.na(year), !is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
  mutate(year = as.integer(year)) %>%
  filter(year >= 1980)

GBIF_data_sf <- st_as_sf(GBIF_data_filtered,
                         coords = c("decimalLongitude", "decimalLatitude"),
                         crs = 4326, remove = FALSE)

# Remove points inside Amazon Basin
intersects_amazon <- st_intersects(GBIF_data_sf, amazon_watersheds, sparse = FALSE)
GBIF_outside_amazon <- GBIF_data_sf[!apply(intersects_amazon, 1, any), ]

# Keep only points on land
intersects_land <- st_intersects(GBIF_outside_amazon, land, sparse = FALSE)
GBIF_outside_amazon <- GBIF_outside_amazon[apply(intersects_land, 1, any), ]

# Select relevant columns and remove duplicates
GBIF_outside_amazon <- GBIF_outside_amazon %>%
  dplyr::select(gbifID, occurrenceID, order, family, genus, species, scientificName,
                verbatimScientificNameAuthorship, countryCode, individualCount,
                decimalLongitude, decimalLatitude, eventDate, day, month, year,
                identifiedBy, dateIdentified, geometry) %>%
  distinct(gbifID, species, decimalLongitude, decimalLatitude, year, .keep_all = TRUE) %>%
  drop_na(species)

# Add country names
GBIF_outside_amazon$country_name <- countrycode(GBIF_outside_amazon$countryCode,
                                                origin = "iso2c",
                                                destination = "country.name")

# -----------------------------
# Load Tedesco reference data
# -----------------------------
Tedesco_ref <- read.csv("data/Tedesco_2017/References_Table.csv", header = TRUE, sep = ";")
Tedesco_ref$Year <- str_extract(Tedesco_ref$X3.Reference, "\\d{4}")

Tedesco <- read.csv("data/Tedesco_2017/Occurrence_Table_12092019.csv", header = TRUE, sep = ",") %>%
  rename(Taxon = X6.Fishbase.Valid.Species.Name) %>%
  mutate(Taxon = str_replace_all(Taxon, fixed("."), " "))

amazon_list$Taxon <- str_replace_all(amazon_list$x, fixed("."), " ")

Tedesco <- Tedesco %>%
  mutate(Origin = if_else(Taxon %in% amazon_list$Taxon, "Amazon", "Other exotics")) %>%
  filter(X3.Native.Exotic.Status == "native",
         X7.Occurrence.Status == "valid",
         X3.Ecoregion == "Neotropic",
         Origin == "Amazon",
         Freshwater != "0", Saltwater != "-1", Brackish != "-1") %>%
  filter(!Taxon %in% species_to_remove)

colnames(Tedesco)[c(1, 6, 13)] <- c("BasinName", "species", "country_name")

# -----------------------------
# Add basins to GBIF occurrences
# -----------------------------
basins <- st_read("data/Tedesco_2017/fish_bioregions_per_basin.gpkg") %>%
  st_transform(4326) %>%
  st_make_valid()

GBIF_outside_amazon <- st_transform(GBIF_outside_amazon, st_crs(basins))
GBIF_outside_amazon <- st_join(GBIF_outside_amazon, basins["BasinName"])

# Non-native occurrences (outside Amazon and not native in Tedesco)
non_native_occurrences <- GBIF_outside_amazon %>%
  anti_join(Tedesco, by = c("species", "BasinName")) %>%
  st_transform(st_crs(amazon_watersheds)) %>%
  st_filter(amazon_watersheds, .predicate = st_disjoint)

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

ggsave("figures/global_map.png", map_global, width = 10, height = 6, dpi = 300)

# South America zoom
south_america <- world %>% filter(continent == "South America")

map_SA <- ggplot() +
  geom_sf(data = south_america, fill = "#f0f0f0", color = "grey70", linewidth = 0.2) +
  geom_sf(data = amazon_watersheds, fill = "#306844", alpha = 0.5, color = NA) +
  geom_sf(data = non_native_occurrences, alpha = 0.5, col = "#de7065") +
  facet_wrap(~period, nrow = 1) +
  coord_sf(xlim = c(-100, -30.5), ylim = c(-58, 15), expand = FALSE) +
  theme_void()

ggsave("figures/south_america_map.png", map_SA, width = 12, height = 4, dpi = 300)

# -----------------------------
# Save summary table
# -----------------------------
write.csv(summary_table, "results/non_native_summary.csv", row.names = FALSE)


