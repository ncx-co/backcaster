#' A compressed csv file containing observations of biomass by age class and
#' species to be used as a lookup table for testing backcasting.
#'
"test_landis_training_file"

#' A compressed csv file containing observations of biomass by age class and
#' species to be used as a target dataset for testing backcasting.
#'
"test_landis_prediction_file"

#' The file path to a directory containing the compressed treelist csv that
#' corresponds to \code{\link{test_landis_training_file}}.
#'
"test_landis_treelist_dir"

#' The dataframe created by \code{\link{process_landis}} after reading
#' \code{\link{test_landis_training_file}}.
#'
"test_landis"

#' The output of \code{\link{cluster_data}} when applied to
#' \code{\link{test_landis}} to generate 20 clusters based on biomass by
#' species.
#'
"test_clusters"

#' The output of \code{\link{get_matching_map_codes}} when applied to
#' \code{\link{test_landis_prediction_file}} to identify nearest neighbor
#' matches from \code{\link{test_landis}}.
#'
"test_match_frame"

#' The output of \code{\link{backcast_landis_to_treelists}} when applied to
#' the test data inputs.
#'
"test_backcasted"

#' A dataframe that relates \code{map_code} values to 90 meter pixel
#' coordinates in EPSG:2163 projection.
#'
"map_code_crosswalk"

#' A dataframe with estimates of canopy cover proportion with bapa and tpa
#' observations from FIA data in the Sierra Nevada region
#'
"canopy_cover_table"
