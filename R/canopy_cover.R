
#' Estimate canopy cover
#'
#' @param trees dataframe of tree records for a set of pixels
#'
#' @return dataframe with estimate of canopy cover proportion by pixel
#' @export
#'
#' @examples
#' cc <- estimate_canopy_cover(test_backcasted$trees)

estimate_canopy_cover <- function(trees) {
  trees %>%
    dplyr::filter(.data[["diameter"]] >= 4) %>%
    dplyr::group_by(
      .data[["pix_ctr_wkt"]], .data[["map_code"]]
    ) %>%
    dplyr::summarize(
      bapa = sum(
        .data[["nTrees"]] * 0.005454 * .data[["diameter"]] ^ 2 /
          (.data[["pix_area_ha"]] * 2.47)),
      tpa = sum(.data[["nTrees"]] / (.data[["pix_area_ha"]] * 2.47)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      cc = pull_canopy_cover(bapa = .data[["bapa"]], tpa = .data[["tpa"]])
    )
}


#' Pull canopy cover estimate from FIA data
#'
#' @param bapa numeric vector of basal area (sq. ft per acre) in stems 4" and
#' larger.
#' @param tpa numeric vector of trees per acre in stems 4" and larger.
#' @param k integer number of matching records from which average canopy cover
#' estimate will be drawn
#'
#' @return numeric vector of canopy cover estimates expressed as a proportion
#' (0 - 1).
#' @export
#'
#' @examples
#' pull_canopy_cover(bapa = c(100, 200), tpa = c(100, 200), k = 5)

pull_canopy_cover <- function(bapa, tpa, k = 5) {
  # identify nearest neighbors in canopy_cover_table
  nearest_matches <- FNN::get.knnx(
    data = backcaster::canopy_cover_table[, c("bapa", "tpa")],
    query = data.frame(bapa = bapa, tpa = tpa),
    k = k
  )

  # pull average canopy cover value from closest matches in canopy_cover_table
  apply(
    nearest_matches$nn.index,
    MARGIN = 1,
    FUN = function(x) mean(backcaster::canopy_cover_table$cc_prop[x])
  )
}
