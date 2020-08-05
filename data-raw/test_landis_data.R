test_landis_training_file <- system.file(
  "extdata", "landis", "10_10.csv.gz",
  package = "backcaster"
)

test_landis_prediction_file <- system.file(
  "extdata", "landis", "10_11.csv.gz",
  package = "backcaster"
)

test_landis_treelist_dir <- system.file(
  "extdata", "trees",
  package = "backcaster"
)

# add a landis data object
test_landis <- process_landis(data.table::fread(test_landis_training_file))

# test clusters object
test_clusters <- cluster_data(
  data = test_landis,
  vars = attr(test_landis, "spp_vars"),
  n = 20
)

# add cluster ids to test_landis
test_landis$cluster <- test_clusters$kmeans$cluster

# add test mapcode match frame
test_match_frame <- get_matching_map_codes(
  new_data = process_landis(data.table::fread(test_landis_prediction_file)),
  lookup = test_landis,
  clusters = test_clusters
)

# read map_code_crosswalk table from s3
map_code_crosswalk <- readr::read_csv(
  "https://silviaterra-delivery.s3.amazonaws.com/TCSI/mapcode_crosswalk_20200803.csv"
)

# add to package data
usethis::use_data(
  test_landis_training_file,
  test_landis_prediction_file,
  test_landis_treelist_dir,
  test_landis,
  test_clusters,
  test_match_frame,
  map_code_crosswalk,
  overwrite = TRUE
)
