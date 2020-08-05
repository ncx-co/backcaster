#' Process Landis summary data into wide dataframe for knn matching
#'
#' @param data dataframe containing observations of biomass (grams per sq. m) by
#' species (\code{landis_species}), age class (\code{age_class}), and pixel
#' (\code{map_code})
#'
#' @return wide dataframe with one row per pixel, a column with total biomass,
#' one column per species, and one column per age class
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' process_landis(data.table::fread(test_landis_training_file))

process_landis <- function(data) {
  data <- dtplyr::lazy_dt(data)

  # total biomass
  biomass_total <- data %>%
    dplyr::group_by(
      .data[["map_code"]]
    ) %>%
    dplyr::summarize(
      aboveground_biomass_g_per_m2 = sum(
        .data[["aboveground_biomass_g_per_m2"]]
      )
    ) %>%
    data.table::as.data.table()

  # biomass by species
  biomass_x_spp <- data %>%
    dplyr::group_by(
      .data[["map_code"]],
      .data[["landis_species"]]
    ) %>%
    dplyr::summarize(
      aboveground_biomass_g_per_m2 = sum(
        .data[["aboveground_biomass_g_per_m2"]]
      )
    ) %>%
    data.table::as.data.table()

  biomass_x_spp_wide <- biomass_x_spp %>%
    tidyr::pivot_wider(
      names_from = "landis_species",
      values_from = "aboveground_biomass_g_per_m2",
      values_fill = 0
    )

  # biomass by age class
  biomass_x_age <- data %>%
    dplyr::mutate(
      age_class = dplyr::case_when(
        .data[["age_class"]] < 0 ~ as.numeric(0),
        TRUE ~ as.numeric(.data[["age_class"]])
      )
    ) %>%
    dplyr::group_by(
      .data[["map_code"]],
      .data[["age_class"]]
    ) %>%
    dplyr::summarize(
      aboveground_biomass_g_per_m2 = sum(
        .data[["aboveground_biomass_g_per_m2"]]
      )
    ) %>%
    data.table::as.data.table()

  biomass_x_age_wide <- biomass_x_age %>%
    tidyr::pivot_wider(
      names_from = "age_class",
      values_from = "aboveground_biomass_g_per_m2",
      values_fill = 0,
      names_prefix = "age_"
    )

  # combine
  out <- data.frame(biomass_total) %>%
    dplyr::left_join(
      tibble::as_tibble(biomass_x_spp_wide),
      by = "map_code"
    ) %>%
    dplyr::left_join(
      tibble::as_tibble(biomass_x_age_wide),
      by = "map_code"
    )

  # store species and age variables as attributes
  spp_vars <- setdiff(names(biomass_x_spp_wide), "map_code")
  age_vars <- setdiff(names(biomass_x_age_wide), "map_code")

  attr(out, "spp_vars") <- spp_vars
  attr(out, "age_vars") <- age_vars

  return(out)
}

#' Cluster Landis summary data using kmeans clustering algorithm. See
#' \code{stats::\link[stats]{kmeans}} for details.
#'
#' @param data dataframe output of \code{\link{process_landis}}
#' @param vars character vector of variables to be used for grouping pixels
#' into clusters
#' @param n_clusters integer number of clusters
#'
#' @return list containing an object of class \code{kmeans} and the scale
#' parameters of the data passed to the kmeans function.
#' @export
#'
#' @examples
#' cluster_data(
#'   data = test_landis,
#'   vars = attr(test_landis, "spp_vars"),
#'   n = 30
#' )

cluster_data <- function(data, vars, n_clusters) {
  # drop vars that have 0 variance
  zero_vars <- purrr::map_lgl(
    .x = vars,
    .f = ~ sd(data[[.x]]) == 0
  )

  vars <- vars[!zero_vars]

  # center the data on 0 but do not standardize variance
  data_scaled <- scale(
    data[, vars],
    scale = FALSE
  )

  clusters <- stats::kmeans(
    data_scaled,
    centers = n_clusters,
    iter.max = 25
  )

  # return essential pieces
  list(
    kmeans = clusters,
    data_centers = attr(data_scaled, "scaled:center"),
    data_scales = attr(data_scaled, "scaled:scale")
  )

}

#' Assign new data to existing clusters
#'
#' @param data dataframe output of \code{\link{process_landis}}
#' @param data_centers center parameters for scaled data
#' @param data_scales scale parameters for scaled data
#' @param cluster_centers matrix of center values for clusters
#'
#' @return vector of integer cluster assignments
#' @export
#'
#' @examples
#' assign_to_cluster(
#'   data = process_landis(data.table::fread(test_landis_prediction_file)),
#'   data_centers = test_clusters$data_centers,
#'   data_scales = test_clusters$data_scales,
#'   cluster_centers = test_clusters$kmeans$center
#' )

assign_to_cluster <- function(data, data_centers, data_scales,
                              cluster_centers) {

  if (is.null(data_scales)) {
    data_scales = FALSE
  }

  scaled <- scale(
    data[, names(data_centers)],
    center = data_centers,
    scale = data_scales
  )

  matches <- FNN::get.knnx(
    cluster_centers,
    scaled[, names(data_centers)],
    1
  )

  matches$nn.index[, 1] # index == cluster
}

#' Get nearest neighbor matches for new data from pool of pixels
#'
#' @param new_data dataframe object that needs nearest neighbor matches, output
#' of \code{\link{process_landis}}
#' @param lookup dataframe output of \code{\link{process_landis}} from which
#' matches will be drawn
#' @param vars character vector of variables to use for k nearest neighbor
#' matching
#'
#' @return dataframe of map_codes from \code{new_data} and matching map_codes
#' from \code{lookup}
#' @export
#'
#' @examples
#' pull_matches(
#'   new_data = process_landis(data.table::fread(test_landis_prediction_file)),
#'   lookup = test_landis,
#'   vars = c(
#'     "aboveground_biomass_g_per_m2",
#'     attr(test_landis, "spp_vars"),
#'     attr(test_landis, "age_vars")
#'   )
#' )

pull_matches <- function(new_data, lookup, vars) {

  # drop vars that are all 0
  zero_vars <- purrr::map_lgl(
    .x = vars,
    .f = ~ length(unique(lookup[[.x]])) == 1
  )

  vars <- vars[!zero_vars]

  # add 0s to new_data for missing columns
  for (var in vars) {
    if (!var %in% names(new_data)) {
      new_data[[var]] <- 0
    }
  }

  # center the lookup data
  lookup_scaled <- scale(
    lookup[, vars],
    scale = FALSE
  )

  # apply same transformation to the new data
  new_data_scaled <- scale(
    x = new_data[, vars],
    center = attr(lookup_scaled, "scaled:center"),
    scale = FALSE
  )

  # run k nearest neighbors matching, where k = 1
  matches <- FNN::get.knnx(
    data = lookup_scaled,
    query = new_data_scaled,
    k = 1
  )

  # pull index
  match_idxs <- matches$nn.index[, 1]

  # corresponding mapcodes
  match_mapcodes <- lookup$map_code[match_idxs]

  # return data.frame with original map_code and matched map_code
  tibble::tibble(original = new_data$map_code, matched = match_mapcodes)

}

#' Get matching map_codes from lookup data
#'
#' @param new_data dataframe object that needs nearest neighbor matches, output
#' of \code{\link{process_landis}}.
#' @param lookup dataframe output of \code{\link{process_landis}} from which
#' matches will be drawn.
#' @param clusters list output from \code{\link{cluster_data}}.
#'
#' @return dataframe of map_codes from \code{new_data} and matching map_codes
#' from \code{lookup}.
#' @export
#'
#' @examples
#' # add cluster IDs to test_landis
#' test_landis$cluster <- test_clusters$kmeans$cluster
#' get_matching_map_codes(
#'   new_data = process_landis(data.table::fread(test_landis_prediction_file)),
#'   lookup = test_landis,
#'   clusters = test_clusters
#' )

get_matching_map_codes <- function(new_data, lookup, clusters) {

  # match evaluation data to a cluster
  new_data$cluster <- assign_to_cluster(
    data = new_data,
    data_centers = clusters$data_centers,
    data_scales = clusters$data_scales,
    cluster_centers = clusters$kmeans$center
  )

  # within each cluster, pull nearest neighbor matches
  new_data_split <- split(new_data, f = new_data$cluster)
  lookup_split <- split(lookup, f = lookup$cluster)

  # variables for evaluating similarity in knn search
  match_vars <- c(
    "aboveground_biomass_g_per_m2",
    attr(lookup, "spp_vars"),
    attr(lookup, "age_vars")
  )

  # run knn matching within each cluster
  match_frame <- purrr::map2_dfr(
    .x = new_data_split,
    .y = lookup_split[names(new_data_split)],
    .f = ~ pull_matches(
      new_data = .x,
      lookup = .y,
      vars = match_vars
    )
  )

  return(match_frame)
}

#' Get matching tree records
#'
#' @param match_frame dataframe output from
#' \code{\link{get_matching_map_codes}} with \code{map_code} values for the
#' target pixels and for the matched records from the lookup dataset.
#' @param map_code_crosswalk dataframe that relates mapcode pixel identifiers
#' (\code{map_code}) with pixel coordinates (\code{pix_ctr_wkt}), and map code
#' index values.
#' @param treelist_dir file path to directory containing treelist csv files
#' corresponding to the pool of pixels in the lookup dataset.
#'
#' @return dataframe of tree records for the target pixels
#' @export
#'
#' @examples
#' get_trees(
#'   match_frame = test_match_frame,
#'   map_code_crosswalk = map_code_crosswalk,
#'   treelist_dir = test_landis_treelist_dir
#' )

get_trees <- function(match_frame, map_code_crosswalk, treelist_dir) {
  file_frame <- match_frame %>%
    tidyr::extract(
      col = "matched",
      into = "minigrid_group",
      regex = "^([^_]*_[^_]*)_.*$",
      remove = FALSE
    ) %>%
    dplyr::mutate(
      treelist_fn = file.path(
        treelist_dir,
        paste0(.data[["minigrid_group"]], ".csv.gz")
      ),
      pix_ctr_wkt = map_code_crosswalk$pix_ctr_wkt[
        match(.data[["matched"]], map_code_crosswalk$map_code)
      ]
    )

  # list of all treelist files requried
  need_files <- unique(file_frame$treelist_fn)

  # import with data.table::fread
  trees <- do.call(rbind, lapply(need_files, data.table::fread))

  # filter down to pixels in file_frame
  in_trees <- dtplyr::lazy_dt(trees) %>%
    dplyr::filter(
      .data[["pix_ctr_wkt"]] %in% file_frame$pix_ctr_wkt
    ) %>%
    data.table::as.data.table() %>%
    tibble::as_tibble() %>%
    dplyr::left_join(file_frame, by = "pix_ctr_wkt") %>%
    dplyr::mutate(
      pix_ctr_wkt = map_code_crosswalk$pix_ctr_wkt[
        match(.data[["original"]], map_code_crosswalk$map_code)
      ]
    ) %>%
    dplyr::select(
      .data[["pix_ctr_wkt"]],
      map_code = .data[["original"]],
      .data[["pix_area_ha"]],
      .data[["common"]],
      .data[["diameter"]],
      .data[["nTrees"]]
    )

  return(in_trees)
}

#' Summarize landis summary matching performance
#'
#' @param match_frame dataframe output from
#' \code{\link{get_matching_map_codes}} with \code{map_code} values for the
#' target pixels and for the matched records from the lookup dataset.
#' @param new_data dataframe object that needs nearest neighbor matches, output
#' of \code{\link{process_landis}}.
#' @param lookup dataframe output of \code{\link{process_landis}} from which
#' matches are drawn.
#'
#' @return list a containing a dataframe with matching accuracy statistics
#' (\code{comp_stats}), and a dataframe with original and matched values of
#' biomass by species and age class for \code{new_data} (\code{comp_frame})
#' @export
#'
#' @examples
#' summarize_matching_performance(
#'   match_frame = test_match_frame,
#'   new_data = process_landis(data.table::fread(test_landis_prediction_file)),
#'   lookup = test_landis
#' )

summarize_matching_performance <- function(match_frame, new_data, lookup) {
  match_data <- match_frame %>%
    dplyr::left_join(lookup, by = c("matched" = "map_code")) %>%
    dplyr::rename(
      "map_code" = "original"
    ) %>%
    dplyr::select(
      tidyselect::all_of(names(new_data))
    )

  comp_frame <-
    list(
      "original" = new_data,
      "matched" = match_data
    ) %>%
    dplyr::bind_rows(
      .id = "source"
    ) %>%
    tidyr::pivot_longer(
      -tidyselect::all_of(c("map_code", "source")),
      names_to = "attr",
      values_to = "estimate"
    ) %>%
    tidyr::pivot_wider(
      names_from = "source",
      values_from = "estimate"
    )

  comp_stats <- comp_frame %>%
    dplyr::group_by(
      .data[["attr"]]
    ) %>%
    dplyr::summarize(
      mean_original = mean(.data[["original"]]),
      mean_matched = mean(.data[["matched"]]),
      rmse = sqrt(
        mean((.data[["matched"]] - .data[["original"]]) ^ 2)
      ),
      rmse_pct = round(100 * .data[["rmse"]] / .data[["mean_original"]]),
      .groups = "drop"
    )

  # return list
  list(
    comp_frame = comp_frame,
    comp_stats = comp_stats
  )
}

#' Backcast Landis to treelists
#'
#' @param new_data dataframe object that will get backcasted treelists, output
#' of \code{\link{process_landis}}
#' @param lookup dataframe output of \code{\link{process_landis}} from which
#' matches will be drawn.
#' @param n_clusters integer number of clusters for grouping the pool of input
#' pixels in \code{lookup}.
#' @param treelist_dir file path to directory containing treelist csv files
#' corresponding to the pool of pixels in the lookup dataset.
#'
#' @return list containing backcasted treelists for pixels in \code{new_data},
#' a dataframe with matching accuracy statistics (\code{comp_stats}), and a
#' dataframe with original and matched values of biomass by species and age
#' class for \code{new_data} (\code{comp_frame})
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' backcasted <- backcast_landis_to_treelists(
#'   new_data = process_landis(data.table::fread(test_landis_prediction_file)),
#'   lookup = test_landis,
#'   n_clusters = 20,
#'   treelist_dir = test_landis_treelist_dir
#' )

backcast_landis_to_treelists <- function(new_data, lookup, n_clusters = 50,
                                         treelist_dir) {
  # kmeans clusters based on species biomass
  rlang::inform("clustering lookup table on biomass x species")
  spp_clusters <- cluster_data(
    data = lookup,
    vars = attr(lookup, "spp_vars"),
    n_clusters = n_clusters
  )

  # add cluster ids to lookup
  lookup$cluster <- spp_clusters$kmeans$cluster

  # run cluster-wise knn matching
  rlang::inform("identifying nearest neighbor matches")
  match_frame <- get_matching_map_codes(
    new_data = new_data,
    lookup = lookup,
    clusters = spp_clusters
  )

  # get trees
  rlang::inform("collecting matching tree records")
  trees <- get_trees(
    match_frame = match_frame,
    map_code_crosswalk = backcaster::map_code_crosswalk,
    treelist_dir = treelist_dir
  )

  # build table to compare matched attributes to original
  rlang::inform("summarizing original vs matched Landis attributes")

  match_diagnostics <- summarize_matching_performance(
    match_frame = match_frame,
    new_data = new_data,
    lookup = lookup
  )

  # export trees and accompanying matching diagnostics
  list(
    trees = trees,
    comp_stats = match_diagnostics$comp_stats,
    comp_frame = match_diagnostics$comp_frame
  )
}
