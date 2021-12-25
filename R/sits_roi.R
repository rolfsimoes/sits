#' @title Region of interest (roi)
#'
#' @name roi_functions
#'
#' @keywords internal
#'
#' @description
#' Region of interest functions
#'
#' @param  cube  a sits data cube object (inherited from `tbl`)
#' @param  roi   a vector containing either `("xmin", "xmax", "ymin", "ymax")`
#' or `("lon_min", "lat_min", "lon_max", "lat_max")` variables, representing
#' a spatial region of interest
#'
#' @return a `bbox` object (inherited from `tbl`)
NULL

#' @rdname roi_functions
#' @description `.sits_roi_bbox()`: convert a given `roi` to a bbox with same
#' crs as `cube` (that should have only one row). If no `roi` is given,
#' returns the `cube`'s `bbox`.
.sits_roi_bbox <- function(roi, cube) {

    # set caller to show in errors
    .check_set_caller(".sits_roi_bbox")

    # pre-condition
    .check_cube(cube, is_tile = TRUE)

    if (is.null(roi))
        return(.bbox(cube))

    bbox <- .roi_to_bbox(roi, crs = .crs_1(cube))

    return(bbox)
}

#' @keywords internal
.roi_to_bbox <- function(x, crs = NULL) {

    # pre-condition
    .check_chr(crs, len_min = 1, len_max = 1, allow_null = TRUE,
               msg = "invalid 'crs' parameter")

    UseMethod(".roi_to_bbox", x)
}

#' @keywords internal
#' @export
.roi_to_bbox.sf <- function(x, crs = NULL) {

    .check_package("sf")

    bbox <- tibble::as_tibble_row(c(sf::st_bbox(x)))

    .crs(bbox) <- sf::st_crs(x)[["input"]]

    if (is.null(crs))
        return(bbox)

    if (.same_crs(.crs(bbox), crs))
        return(bbox)

    .bbox_transform(bbox, to_crs = crs)
}

#' @keywords internal
#' @export
.roi_to_bbox.default <- function(x, crs = NULL) {

    if (all(c("xmin", "xmax", "ymin", "ymax") %in% names(x)))
        return(.bbox(x, default_crs = crs))


    if (all(c("lon_min", "lon_max", "lat_min", "lat_max") %in% names(x))) {

        ll_names <- c("xmin", "xmax", "ymin", "ymax")
        names(ll_names) <- c("lon_min", "lon_max", "lat_min", "lat_max")
        names(x) <- unname(ll_names[names(x)])

        bbox <- .bbox(x, default_crs = "EPSG:4326")
        return(bbox)
    }

    .check_that(FALSE,
                local_msg = paste("value should contain either",
                                  "('xmin', 'xmax', 'ymin', 'ymax') or",
                                  "('lon_min', 'lon_max', 'lat_min',",
                                  "'lat_max') variables"),
                msg = "invalid roi value")
}
