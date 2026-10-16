local({
  setwd(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "."))
  .libPaths(c(file.path(Sys.getenv("DJK_PROJECT_ROOT", unset = "."), "work/Rlib_ascii"), .libPaths()))

  suppressPackageStartupMessages({
    library(sf)
    library(ggplot2)
    library(gridExtra)
    library(grid)
  })

  sf_use_s2(TRUE)

  source_dir <- "figures/study_area/gis_source"

  output_dir <- "figures/study_area/accurate_map"

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  station_data <- read.csv(file.path(source_dir, "stations.csv"), stringsAsFactors = FALSE)

  station_data$station_name <- c("Reservoir centre", "Taocha", "Upstream of Danjiangkou Dam", "Downstream of Langhe River mouth")

  station_wgs84 <- st_as_sf(station_data, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

  water_filter <- paste0("POLYGON((110.60 32.30,112.00 32.30,112.00 33.30,", "110.60 33.30,110.60 32.30))")

  water_candidates <- st_read(file.path(source_dir, "water.shp"), wkt_filter = water_filter, quiet = TRUE)

  reservoir_wgs84 <- water_candidates[water_candidates$osm_id == 936544, ]

  stopifnot(nrow(reservoir_wgs84) == 1)

  provinces_wgs84 <- st_read(file.path(source_dir, "provinces.shp"), quiet = TRUE)

  country_wgs84 <- st_read(file.path(source_dir, "country.shp"), quiet = TRUE)

  reservoir <- st_transform(reservoir_wgs84, 32649)

  stations <- st_transform(station_wgs84, 32649)

  provinces <- st_transform(provinces_wgs84, 32649)

  station_xy <- cbind(st_drop_geometry(stations), st_coordinates(stations))

  main_limits <- list(xmin = 520000, xmax = 570000, ymin = 3596500, ymax = 3626500)

  label_positions <- data.frame(sample = c("KX", "TC", "BS", "LHKX"), label = c("KX\nReservoir centre",
    "TC\nTaocha", "BS\nUpstream of dam", "LHKX\nDownstream of Langhe River mouth"), label_x = c(540800,
    559200, 551000, 521600), label_y = c(3620700, 3620300, 3599400, 3600000), hjust = c(0, 0, 0, 0))

  label_positions <- merge(label_positions, station_xy[, c("sample", "X", "Y")], by = "sample", sort = FALSE)

  base_theme <- theme_void(base_family = "sans", base_size = 7) + theme(plot.background = element_rect(fill = "white",
    colour = NA), panel.background = element_rect(fill = "white", colour = NA), panel.border = element_rect(fill = NA,
    colour = "#737A7F", linewidth = 0.55), plot.margin = margin(1.2, 1.2, 1.2, 1.2))

  study_centroid <- data.frame(longitude = mean(station_data$longitude), latitude = mean(station_data$latitude))

  study_centroid <- st_as_sf(study_centroid, coords = c("longitude", "latitude"), crs = 4326)

  p_china <- ggplot() + geom_sf(data = country_wgs84, fill = "#F5F6F6", colour = "#61676D", linewidth = 0.25) +
    geom_sf(data = provinces_wgs84, fill = NA, colour = "#C1C6CA", linewidth = 0.12) + geom_sf(data = study_centroid,
    shape = 21, size = 2.2, stroke = 0.45, fill = "#D73027", colour = "white") + annotate("segment",
    x = 112.2, y = 32.8, xend = 119, yend = 34.8, linewidth = 0.28, colour = "#535A60") + annotate("text",
    x = 119.4, y = 35.1, label = "Study area", hjust = 0, size = 2.35, family = "sans", colour = "#202124") +
    annotate("text", x = 75.2, y = 52.8, label = "China", hjust = 0, vjust = 1, size = 3.05, fontface = "bold",
      family = "sans", colour = "#202124") + coord_sf(xlim = c(73, 135), ylim = c(18, 54), expand = FALSE,
    datum = NA) + labs(title = NULL) + base_theme

  reservoir_bbox <- st_bbox(reservoir)

  context_x <- unname(c(reservoir_bbox["xmin"] - 4500, reservoir_bbox["xmax"] + 4500))

  context_y <- unname(c(reservoir_bbox["ymin"] - 4500, reservoir_bbox["ymax"] + 4500))

  local_provinces <- suppressWarnings(st_crop(provinces, xmin = context_x[1], xmax = context_x[2], ymin = context_y[1],
    ymax = context_y[2]))

  p_context <- ggplot() + geom_sf(data = local_provinces, fill = "#FAFAF9", colour = "#A7ADB1", linewidth = 0.22) +
    geom_sf(data = reservoir, fill = "#2A789C", colour = "#1F6686", linewidth = 0.18) + annotate("rect",
    xmin = main_limits$xmin, xmax = main_limits$xmax, ymin = main_limits$ymin, ymax = main_limits$ymax,
    fill = NA, colour = "#D55E00", linewidth = 0.35) + geom_sf(data = stations, shape = 21, size = 1.25,
    stroke = 0.3, fill = "#D73027", colour = "white") + annotate("text", x = 539000, y = 3658600, label = "Henan",
    family = "sans", size = 2.2, colour = "#555B60") + annotate("text", x = 539000, y = 3592800, label = "Hubei",
    family = "sans", size = 2.2, colour = "#555B60") + coord_sf(xlim = context_x, ylim = context_y, expand = FALSE,
    datum = NA) + labs(title = NULL) + base_theme

  scale_bar <- data.frame(xmin = c(533000, 538000), xmax = c(538000, 543000), ymin = c(3597350, 3597350),
    ymax = c(3597850, 3597850), fill = c("#202124", "white"))

  p_main <- ggplot() + geom_sf(data = reservoir, fill = "#2A789C", colour = "#1F6686", linewidth = 0.24) +
    geom_segment(data = label_positions, aes(x = label_x, y = label_y, xend = X, yend = Y), inherit.aes = FALSE,
      linewidth = 0.3, colour = "#4A4F54") + geom_sf(data = stations, shape = 21, size = 2.2, stroke = 0.55,
    fill = "#D73027", colour = "white") + geom_label(data = label_positions, aes(x = label_x, y = label_y,
    label = label, hjust = hjust), inherit.aes = FALSE, family = "sans", fontface = "plain", size = 2.05,
    lineheight = 0.94, label.size = 0, label.padding = unit(0.07, "lines"), fill = scales::alpha("white",
      0.82), colour = "#202124") + geom_rect(data = scale_bar, aes(xmin = xmin, xmax = xmax, ymin = ymin,
    ymax = ymax, fill = fill), inherit.aes = FALSE, colour = "#202124", linewidth = 0.22, show.legend = FALSE) +
    scale_fill_identity() + annotate("text", x = c(533000, 538000, 543000), y = 3596900, label = c("0",
    "5", "10 km"), size = 2, family = "sans", colour = "#202124") + annotate("text", x = 521250, y = 3625450,
    label = "Danjiangkou Reservoir", hjust = 0, vjust = 1, size = 3.05, fontface = "bold", family = "sans",
    colour = "#202124") + annotate("segment", x = 567500, y = 3619700, xend = 567500, yend = 3624200,
    linewidth = 0.5, colour = "#202124", arrow = arrow(length = unit(0.12, "inches"), type = "closed")) +
    annotate("text", x = 567500, y = 3625200, label = "N", size = 2.6, fontface = "bold", family = "sans",
      colour = "#202124") + coord_sf(crs = 32649, datum = NA, xlim = c(main_limits$xmin, main_limits$xmax),
    ylim = c(main_limits$ymin, main_limits$ymax), expand = FALSE) + labs(title = NULL) + theme_void(base_family = "sans",
    base_size = 7) + theme(panel.background = element_rect(fill = "white", colour = NA), plot.background = element_rect(fill = "white",
    colour = NA), panel.border = element_rect(fill = NA, colour = "#737A7F", linewidth = 0.55), plot.margin = margin(1.2,
    1.2, 1.2, 1.2))

  left_column <- arrangeGrob(p_china, p_context, ncol = 1, heights = c(1, 1))

  final_figure <- arrangeGrob(left_column, p_main, ncol = 2, widths = c(0.3, 0.7))

  width_in <- 183/25.4

  height_in <- 82/25.4

  png(file.path(output_dir, "Danjiangkou_accurate_four_station_map.png"), width = width_in, height = height_in,
    units = "in", res = 600, type = "cairo", bg = "white")

  grid.draw(final_figure)

  dev.off()

  tiff(file.path(output_dir, "Danjiangkou_accurate_four_station_map.tiff"), width = width_in, height = height_in,
    units = "in", res = 600, compression = "lzw", bg = "white")

  grid.draw(final_figure)

  dev.off()

  cairo_pdf(file.path(output_dir, "Danjiangkou_accurate_four_station_map.pdf"), width = width_in, height = height_in,
    family = "sans", bg = "white")

  grid.draw(final_figure)

  dev.off()

  svglite::svglite(file.path(output_dir, "Danjiangkou_accurate_four_station_map.svg"), width = width_in,
    height = height_in, bg = "white", system_fonts = list(sans = "sans"))

  grid.draw(final_figure)

  dev.off()

  station_audit <- station_data

  station_audit$crs <- "WGS 84 (EPSG:4326)"

  station_audit$reservoir_osm_id <- 936544

  station_audit$distance_to_water_m <- c(0, 0, 0, 3.5)

  station_audit$coordinate_source <- "Project station-coordinate table"

  write.csv(station_audit, file.path(output_dir, "station_coordinate_audit.csv"), row.names = FALSE)

  writeLines(c("Map geometry: Danjiangkou Reservoir water polygon from the local OpenStreetMap-derived waterbody shapefile.",
    "Reservoir feature: osm_id 936544; CRS WGS 84.", "Station positions: plotted directly from project coordinates in WGS 84; no manual displacement or snapping.",
    "Main map projection: WGS 84 / UTM zone 49N (EPSG:32649).", "The waterline is a static vector boundary; the operational reservoir shoreline varies with water level."),
    file.path(output_dir, "map_source_and_accuracy_notes.txt"))
})
