library(tidyverse)
library(sf)
library(aws.s3)
library(tidyFIA)

bucket <- 'silviaterra-brian'
folder <- 'percent-canopy'

download.file("https://apps.fs.usda.gov/fia/datamart/CSV/REF_SPECIES.csv",
  destfile = "/tmp/ref_species.csv"
)

ref_species <- read_csv('/tmp/ref_species.csv')

s3load(
  object = glue::glue('{folder}/lcw_model.Rda'),
  bucket = bucket
)

s3load(
  object = glue::glue('{folder}/supersection_shape.rda'),
  bucket = bucket
)

bounds <- supersectionShape %>%
  filter(SSection == 'Sierra Nevada') %>%
  st_transform(4326)

# download FIA plots for CA
fia_data <- tidyFIA::tidy_fia(
  aoi = bounds
)

plots <- fia_data[['plot']] %>%
  filter(
    invyr >= 2014,
    invyr != 9999,
    plot_status_cd == 1
  ) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_join(bounds)

trees <- fia_data[['tree']]  %>%
  filter(
    plt_cn %in% plots$cn,
    statuscd == 1,
    dia >= 4.5
  ) %>%
  transmute(
    plt_cn,
    subp,
    statuscd,
    spcd,
    common = ref_species$COMMON_NAME[match(spcd, ref_species$SPCD)],
    dia = round(dia),
    ht,
    crown_width_ft = predict_lcw(
      common = common,
      diameter = dia,
      lcwTable = lcwTable
    ),
    crown_area_ft2 = pi * (crown_width_ft / 2) ^ 2,
    tpa_unadj,
    azimuth,
    dist
  ) %>%
  group_by(
    dia
  ) %>%
  mutate(
    common = factor(common),
    crown_width_ft = case_when(
      is.na(crown_width_ft) ~ mean(crown_width_ft, na.rm = TRUE),
      !is.na(crown_width_ft) ~ crown_width_ft
    )
  ) %>%
  filter(crown_width_ft > 0) %>%
  filter(!is.nan(crown_width_ft))

map_crown_cover <- function(plot_trees) {
  plot_cn <- unique(plot_trees$plt_cn)

  message(plot_cn)
  plot_area <- data.frame(subp = c(1, 2, 3, 4)) %>%
    mutate(
      dist = c(0, 120, 120, 120),
      azimuth = c(0, 0, 120, 240),
      x = dist * sin(azimuth * 2 * pi / 360),
      y = dist * cos(azimuth * 2 * pi / 360)
    ) %>%
    st_as_sf(coords = c("x", "y"), remove = FALSE) %>%
    st_buffer(dist = 24)

  plot_perimeter <- plot_area %>%
    mutate(perimeter = "in") %>%
    group_by(perimeter) %>%
    summarize(.groups = "drop")

  tree_locs <- suppressWarnings(
    plot_trees %>%
      ungroup() %>%
      mutate(
        x0 = plot_area$x[match(subp, plot_area$subp)],
        y0 = plot_area$y[match(subp, plot_area$subp)],
        x1 = dist * sin(azimuth * 2 * pi / 360),
        y1 = dist * cos(azimuth * 2 * pi / 360),
        x = x0 + x1,
        y = y0 + y1
      ) %>%
      st_as_sf(coords = c("x", "y")) %>%
      st_buffer(dist = .$crown_width_ft / 2) %>%
      st_intersection(plot_perimeter) %>%
      mutate(area = st_area(.))
    )

  # get total area covered
  trees_dissolved <- tree_locs %>%
    mutate(crown = "crown") %>%
    group_by(crown) %>%
    summarize(.groups = "drop") %>%
    mutate(area = st_area(.))

  crown_area_dissolved <- trees_dissolved  %>%
    pull(area)

  # overlap %
  plot_stats <- tibble(
    PLT_CN = plot_cn,
    total_crown_area = sum(tree_locs$area),
    area_covered = crown_area_dissolved,
    overlap_prop = (total_crown_area - area_covered) / total_crown_area,
    crown_cover_prop_no_overlap = area_covered / (4 * pi * 24 ^ 2),
    crown_cover_prop_with_overlap = total_crown_area / (4 * pi * 24 ^ 2)
  )

  p <- ggplot() +
    geom_sf(data = tree_locs, aes(fill = common)) +
    geom_sf(data = trees_dissolved, color = "black", alpha = 0) +
    geom_sf(data = plot_perimeter, alpha = 0, color = "black") +
    labs(
      title = paste("PLOT CN:", plot_cn),
      caption = paste(
        "crown area (sq. ft):", round(plot_stats$area_covered),
        "\n",
        "crown cover % (no overlap):", 100 * round(plot_stats$crown_cover_prop_no_overlap, 2),
        "\n",
        "crown cover % (with overlap):", 100 * round(plot_stats$crown_cover_prop_with_overlap, 2)
      )
    )
  p

  return(list(plot_stats = plot_stats, plot_plot = p))
}

trees_list <- trees %>%
  split(f = .$plt_cn)

crown_cover <- purrr::map(
  .x = trees_list,
  .f = ~ map_crown_cover(plot_trees = .x)
)

crown_stats <- purrr::map_dfr(crown_cover, "plot_stats")
crown_plots <- purrr::map(crown_cover, "plot_plot")


plot_sample <- sample(1:length(crown_plots), 4)

for(i in 1:length(plot_sample)){
  tmp_plot <- crown_plots[[plot_sample[i]]]
  print(tmp_plot)
}

cw_quants <- quantile(crown_stats$crown_cover_prop_no_overlap, c(0.5, 0.05, 0.9))

ggplot(crown_stats, aes(x = crown_cover_prop_no_overlap)) +
  geom_histogram(color = "#30363e", fill = "#347b3c") +
  geom_vline(xintercept = cw_quants, color = "black", linetype = "dotted") +
  xlim(0, 1) +
  xlab("Percent canopy cover, corrected for overlap (%)") +
  theme_bw()

plot_stats <- trees %>%
  group_by(plt_cn) %>%
  mutate(stem_ba = 0.005454 * (dia^2)) %>%
  summarise(
    bapa = sum(stem_ba * tpa_unadj),
    tpa = sum(tpa_unadj),
    qmd = sqrt((bapa / tpa) / 0.005454154),
    .groups = "drop"
  ) %>%
  left_join(fia_cond) %>%
  mutate(for_type = factor(fortypcdcalc))

glimpse(crown_stats)

training_data <- crown_stats %>%
  transmute(
    plt_cn = PLT_CN,
    percent_cc = crown_cover_prop_no_overlap
  ) %>%
  left_join(plot_stats) %>%
  filter(!is.na(for_type))

canopy_cover_table <- training_data %>%
  select(
    bapa, tpa, cc_prop = percent_cc
  )

usethis::use_data(canopy_cover_table)
