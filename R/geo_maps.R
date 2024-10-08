# devtools::install_github("UrbanInstitute/urbnmapr")
# devtools::install_github("UrbanInstitute/urbnmapr")
# install.packages("remotes")
# remotes::install_github("UrbanInstitute/urbnthemes")
# devtools::install_github("tidyverse/ggplot2")
#
# suppressPackageStartupMessages({
#   library(ggplot2)
#   library(maps)
#   library(usmap)
#   library(data.table)
#   library(ggsn) # for scale bar `scalebar`
#   library(ggrepel) # if need to repel labels
#   library(deepath)
#   library(urbnmapr)
#   library(dplyr)
#   library(tidyverse)
#   library(urbnmapr)
#   library(urbnthemes)
#   library(stringr)
#   library(cowplot)
# })
suppervised_sil_score <- function(D, metadata, cluster_metadata, cluster_name1, cluster_name2) {
  cluster_a <- colnames(D)[metadata[colnames(D), cluster_metadata] == cluster_name1]
  if (is.na(cluster_name2)) {
    cluster_b <- colnames(D)[metadata[colnames(D), cluster_metadata] != cluster_name1]
  } else {
    cluster_b <- colnames(D)[metadata[colnames(D), cluster_metadata] == cluster_name2]
  }

  diag(D) <- NA
  sil_scores <- vector("numeric", length(cluster_a))
  for (i in 1:length(cluster_a)) {
    if (length(cluster_a) == 1) {
      a <- 0.0
    } else {
      a <- mean(as.numeric(D[cluster_a[i], cluster_a[-i]]), na.rm = TRUE)
      b <- mean(as.numeric(D[cluster_a[i], cluster_b]), na.rm = TRUE)
      s <- (b - a) / max(a, b)
      sil_scores[i] <- s
    }
  }
  print(sil_scores)
  cluster_sil_score <- mean(sil_scores, na.rm = TRUE)
  return(cluster_sil_score)
}

set_geo_coord <- function() {
  dt1 <- as.data.table(copy(state.x77))
  dt1$state <- tolower(rownames(state.x77))
  dt1 <- dt1[, .(state, Population)]
  # only need state name and variable to plot in the input file:
  str(dt1)
  dt2 <- as.data.table(copy(state.x77))
  dt2$state <- tolower(rownames(state.x77))
  dt2 <- dt2[, .(state, Population)]
  dt2 <- as.data.table(urbnmapr::states)
  spatial_data <- left_join(urbnmapr::states,
    get_urbn_map(map = "states", sf = TRUE),
    by = "state_name"
  )
  setkey(dt2, state_name)

  states <- setDT(ggplot2::map_data("state"))
  setkey(states, region)
  return(states)
}
########## within states variations ######
set_cohernace <- function(Distance_file, metadata_file) {
  metadata <- read.delim(
    metadata_file,
    sep = "\t",
    header = TRUE,
    fill = T,
    comment.char = "",
    check.names = F,
    row.names = 1
  )
  rownames(metadata) <- metadata$gisaid_epi_isl

  D <- read.delim(
    Distance_file,
    sep = "\t",
    header = TRUE,
    fill = F,
    comment.char = "",
    check.names = F,
    row.names = 1
  )

  temp <- matrix(NA, ncol = 2, nrow = dim(D)[2])
  D2 <- D[grepl("EPI_ISL_", rownames(D), fixed = TRUE), grepl("EPI_ISL_", colnames(D), fixed = TRUE)]

  rownames(D2) <- str_split_fixed(rownames(D2), "\\|", 3)[, 2]
  colnames(D2) <- str_split_fixed(colnames(D2), "\\|", 3)[, 2]

  intersect_samples <- intersect(rownames(D2), rownames(metadata))
  metadata2 <- metadata[intersect_samples, ]
  D2 <- D2[intersect_samples, intersect_samples]
  # D[is.na(D)] <- 0
  cluster_metadata <- "country_exposure"
  countries <- unique(metadata[, cluster_metadata])
  dt3 <- matrix(NA, nrow = length(countries), ncol = 1)
  colnames(dt3) <- "Coherence"
  rownames(dt3) <- countries
  ### world map
  world <- map_data("world")
  world$Coherence <- NA
  for (country in countries) {
    world[world$region == country, "Coherence"] <- suppervised_sil_score(D2, metadata2,
      cluster_metadata = cluster_metadata,
      cluster_name1 = country, cluster_name2 = NA
    )
  }


  # write.table(world, 'world.txt' , sep = "\t", eol = "\n", na = "", col.names = NA, quote= F, row.names = T)
  # world <- read.delim(
  #  'world.txt',
  #  sep = '\t',
  #  header = TRUE,
  #  fill = F,
  #  comment.char = "" ,
  #  check.names = F,
  #  row.names = 1
  # )
  metadata3 <- metadata2[metadata2$country == "USA", ]
  intersect_samples <- intersect(rownames(D2), rownames(metadata3))
  metadata3 <- metadata[intersect_samples, ]
  D3 <- D2[intersect_samples, intersect_samples]
  usa_covid <- world[world$region == "USA", ]
  cluster_metadata <- "division"
  us_states <- unique(metadata2[metadata2$country == "USA", cluster_metadata])
  dt2$Coherence <- NA


  # create states location and abbreviations for label
  # incl `Population` (the value to plot) in the label dataset, if want to fill with color.
  state_label_dt <- unique(dt2[, .(Coherence, x = mean(range(long)), y = mean(range(lat))), by = state_name])
  snames <- data.table(state = tolower(state.name), abb = state.abb) # these are dataset within R
  setkey(state_label_dt, state_name)
  setkey(snames, state)
  state_label_dt <- snames[state_label_dt]
  state_label_dt$Coherence <- NA
  dt2$Coherence <- as.numeric(dt2$Coherence)
  # All labels for states to the right of lon = -77 will be on the right of lon = -50.
  x_boundary <- -77
  x_limits <- c(-50, NA) # optional, for label repelling

  for (us_state in us_states) {
    print(us_state)
    temp_sil <- suppervised_sil_score(D3, metadata3,
      cluster_metadata = cluster_metadata,
      cluster_name1 = us_state, cluster_name2 = NA
    )
    us_state <- tolower(us_state)
    dt2[dt2$state_name == us_state, "Coherence"] <- temp_sil
    state_label_dt[state_label_dt$state == us_state, "Coherence"] <- temp_sil
  }
  dt2$Coherence <- as.numeric(dt2$Coherence)
  return(world, d2)
}

geo_map <- function(world, d2, output = NA) {
  us_variation_map <- ggplot(data = dt2, aes(x = long, y = lat, group = group)) +
    geom_polygon(aes(fill = Coherence)) +
    geom_path() +
    scale_fill_gradient(low = "#56B1F7", high = "#1313ff", na.value = "#cfcfcf") +
    geom_label_repel(
      data = state_label_dt[x >= x_boundary, ],
      aes(x = x, y = y, label = abb, fill = Coherence),
      arrow = arrow(length = unit(0.02, "npc"), ends = "first"),
      force = 5, hjust = 1, size = 1.75,
      xlim = x_limits, inherit.aes = F
    ) +
    # the normal labels:
    geom_text(
      data = state_label_dt[x < x_boundary, ], aes(x = x, y = y, label = abb),
      size = 1.75, inherit.aes = F
    ) +
    coord_map() +
    theme_classic() +
    labs(fill = "Coherence", x = "Longitude", y = "Latitude") +
    # map scale
    # ggsn::scalebar(data = dt2, dist = 500, dist_unit = "km",
    #                border.size = 0.4, st.size = 4,
    #                box.fill = c('black','white'),
    #                transform = TRUE, model = "WGS84") +
    # put legend at the bottom, adjust legend title and text font sizes
    theme( # legend.position = "bottom",
      legend.title = element_text(size = 12), # font size of the legend
      legend.text = element_text(size = 10),
      axis.title.x = element_blank(), # remove axis, title, ticks
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line = element_blank(),
      axis.title.y = element_blank(), # remove axis, title, ticks
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  us_variation_map
  if (!is.na(output)) {
    ggsave(
      filename = paste(output, "/us_variation.png", sep = ""),
      plot = us_variation_map, width = 7.2, height = 3.75, units = "in", dpi = 350
    )
  }
  world_var_map <- ggplot() +
    geom_map(
      data = world, map = world,
      aes(long, lat, map_id = region, fill = Coherence),
      color = "white", size = 0.1
    ) +
    scale_fill_gradient(low = "#56B1F7", high = "#132B43") +
    # scale_fill_brewer(palette="Blues")+
    theme_classic() +
    theme( # legend.title="Town Name" , #legend.position = "bottom",
      legend.title = element_text(size = 9), # font size of the legend
      legend.text = element_text(size = 6),
      axis.title.x = element_blank(), # remove axis, title, ticks
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line = element_blank(),
      axis.title.y = element_blank(), # remove axis, title, ticks
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) #+guides(fill=guide_legend(title="Starins' homogeneity"))
  if (!is.na(output)) {
    ggsave(
      filename = paste(output, "/world_variation_country.png", sep = ""),
      plot = world_var_map, width = 7.2, height = 3.75, units = "in", dpi = 350
    )
  }

  fig5 <- plot_grid(world_var_map, us_variation_map, #  + theme(legend.position="none"),#final_niche_adaptation, Reference_ratio,
    labels = c("a", "b"),
    ncol = 1, nrow = 2,
    hjust = -0, label_size = 8
  )
  if (!is.na(output)) {
    ggsave(
      filename = paste(output, ",fig5_location_variation.png", sep = ""),
      plot = fig5, width = 7.2, height = 5, units = "in", dpi = 350
    )
  }
  return(fig5)
}
