.libPaths(c(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "work/Rlib_ascii"), .libPaths()))

suppressPackageStartupMessages(library(sf))

wkt <- "POLYGON((111.15 32.45,111.85 32.45,111.85 33.15,111.15 33.15,111.15 32.45))"

water <- st_read("figures/study_area/gis_source/water.shp", wkt_filter = wkt, quiet = TRUE)

stations <- read.csv("figures/study_area/gis_source/stations.csv")

station_sf <- st_as_sf(stations, coords = c("longitude", "latitude"), crs = 4326)

water_m <- st_transform(water, 32649)

station_m <- st_transform(station_sf, 32649)

nearest_rows <- lapply(seq_len(nrow(station_m)), function(i) {
  distances <- as.numeric(st_distance(station_m[i, ], water_m))
  nearest <- order(distances)[1:3]
  data.frame(station = stations$sample[i], rank = seq_along(nearest), osm_id = water$osm_id[nearest],
    distance_m = round(distances[nearest], 1))
})

nearest_table <- do.call(rbind, nearest_rows)

write.csv(nearest_table, "figures/study_area/station_to_water_distance.csv", row.names = FALSE)

print(nearest_table)

area_table <- data.frame(osm_id = water$osm_id, name = water$name, area_km2 = round(as.numeric(st_area(water_m))/1e+06,
  2))

area_table <- area_table[order(-area_table$area_km2), ]

print(head(area_table, 15))
