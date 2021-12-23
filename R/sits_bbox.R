#' @title Get the bounding box of the data
#'
#' @name sits_bbox
#'
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#' @author Rolf Simoes, \email{rolf.simoes@@inpe.br}
#'
#' @description  Obtain a vector of limits (either on lat/long for time series
#'               or in projection coordinates in the case of cubes)
#'
#' @param data      Valid sits tibble (time series or a cube).
#' @param wgs84     Take effect only for data cubes.
#'                  Reproject bbox to WGS84 (EPSG:4326).
#' @param ...       Additional parameters (not implemented).
#'
#' @return named vector with bounding box in WGS84 for time series and
#'         on the cube projection for a data cube unless wgs84 parameter
#'         is TRUE.
#'
#' @export
#'
sits_bbox <- function(data, wgs84 = FALSE, ...) {

    .check_set_caller("sits_bbox")

    # get the meta-type (sits or cube)
    data <- .config_data_meta_type(data)

    UseMethod("sits_bbox", data)
}

#' @export
#'
sits_bbox.sits <- function(data, ...) {
    # is the data a valid set of time series
    .sits_tibble_test(data)

    # get the max and min longitudes and latitudes
    lon_max <- max(data$longitude)
    lon_min <- min(data$longitude)
    lat_max <- max(data$latitude)
    lat_min <- min(data$latitude)
    # create and return the bounding box
    bbox <- c(lon_min, lon_max, lat_min, lat_max)
    names(bbox) <- c("lon_min", "lon_max", "lat_min", "lat_max")
    return(bbox)
}

#' @export
#'
sits_bbox.sits_cube <- function(data, wgs84 = FALSE, ...) {

    # pre-condition
    .cube_check(data)

    # create and return the bounding box
    if (nrow(data) == 1) {
        bbox <- c(xmin = .xmin(data),
                  xmax = .xmax(data),
                  ymin = .ymin(data),
                  ymax = .ymax(data))
    } else {
        bbox <- c(xmin = min(.xmin(data)),
                  xmax = max(.xmax(data)),
                  ymin = min(.ymin(data)),
                  ymax = max(.ymax(data))
        )
    }

    # convert to WGS84?
    if (wgs84) {

        .coords_to_bbox(xmin = bbox[["xmin"]],
                        xmax = bbox[["xmax"]],
                        ymin = bbox[["ymin"]],
                        ymax = bbox[["ymax"]])

    }

    return(bbox)
}

.coords_to_bbox <- function(xmin, xmax, ymin, ymax, crs) {

    pt1 <- c(xmin, ymax)
    pt2 <- c(xmax, ymax)
    pt3 <- c(xmax, ymin)
    pt4 <- c(xmin, ymin)

    bbox <- sf::st_sfc(
        sf::st_polygon(list(rbind(pt1, pt2, pt3, pt4, pt1))), crs = crs
    )

    # create a polygon and transform the proj
    bbox_latlng <- sf::st_bbox(sf::st_transform(bbox, crs = 4326))

    names(bbox_latlng) <- c("lon_min", "lat_min", "lon_max", "lat_max")

    return(bbox_latlng)
}

#' @title Intersection between a bounding box and a cube
#' @name .sits_bbox_intersect
#' @keywords internal
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @param bbox           bounding box for a region of interest
#' @param cube           data cube
#' @return               vector the bounding box intersection
#'
.sits_bbox_intersect <- function(bbox, cube) {
    bbox_out <- vector("double", length = 4)
    names(bbox_out) <- c("xmin", "xmax", "ymin", "ymax")

    if (.xmin(bbox) > .xmax(cube) | .xmax(bbox) < .xmin(cube) |
        .ymin(bbox) > .ymax(cube) | .ymax(bbox) < .ymin(cube)) {
        return(NULL)
    }

    if (.xmin(bbox) < .xmin(cube)) {
        .xmin(bbox_out) <- .xmin(cube)
    } else {
        .xmin(bbox_out) <- .xmin(bbox)
    }

    if (.xmax(bbox) > .xmax(cube)) {
        .xmax(bbox_out) <- .xmax(cube)
    } else {
        .xmax(bbox_out) <- .xmax(bbox)
    }

    if (.ymin(bbox) < .ymin(cube)) {
        .ymin(bbox_out) <- .ymin(cube)
    } else {
        .ymin(bbox_out) <- .ymin(bbox)
    }

    if (.ymax(bbox) > .ymax(cube)) {
        .ymax(bbox_out) <- .ymax(cube)
    } else {
        .ymax(bbox_out) <- .ymax(bbox)
    }

    return(bbox_out)
}


.roi_to_bbox <- function(x, crs = NULL) {

    if (all(c("xmin", "ymin", "xmax", "ymax") %in% names(x)))
        return(.bbox(x, crs = crs))


    if (all(c("lon_min", "lat_min", "lon_max", "lat_max") %in% names(x))) {
        ll_names <- c("xmin", "ymin", "xmax", "ymax")
        names(ll_names) <- c("lon_min", "lat_min", "lon_max", "lat_max")
        names(x) <- unname(ll_names[names(x)])
        return(.bbox(x, crs = "EPSG:4326"))
    }
}

.bbox_vars <- function(x) {

    bbox_vars <- c("xmin", "xmax", "ymin", "ymax", "crs")

    bbox_vars[bbox_vars %in% names(x)]
}

.bbox <- function(x, crs = NULL) {

    UseMethod(".bbox", x)
}

.bbox.default <- function(x, crs = NULL) {

    bbox <- tibble::tibble(xmin = .xmin(x),
                           xmax = .xmax(x),
                           ymin = .ymin(x),
                           ymax = .ymax(x))

    if ("crs" %in% names(x))
        .crs(bbox) <- .crs(x)

    if (!is.null(crs) && !"crs" %in% names(x))
        .crs(bbox) <- crs

    .check_bbox(bbox)

    bbox
}

.check_bbox <- function(x, same_crs = FALSE) {

    .check_chr_contains(names(x), contains = c("xmin", "xmax", "ymin", "ymax"),
                        msg = "object does not have all bbox variables")

    if ("crs" %in% names(x)) {

        crs <- .crs(x)

        .check_that(
            length(crs) == length(.xmin(x)),
            local_msg = "length of crs differs from bbox",
            msg = "invalid bbox value")
    }

    .check_that(
        length(.xmax(x)) == length(.xmin(x)),
        local_msg = "length of xmax differs from xmin",
        msg = "invalid bbox value")

    .check_that(
        length(.ymin(x)) == length(.xmin(x)),
        local_msg = "length of ymin differs from xmin",
        msg = "invalid bbox value")

    .check_that(
        length(.ymax(x)) == length(.xmin(x)),
        local_msg = "length of ymax differs from xmin",
        msg = "invalid bbox value")

    .check_that(
        all(.xmin(x) <= .xmax(x)),
        local_msg = "xmin > xmax",
        msg = "invalid bbox value")

    .check_that(
        all(.ymin(x) <= .ymax(x)),
        local_msg = "ymin > ymax",
        msg = "invalid bbox value")

    if (same_crs)
        .check_that(
            (!"crs" %in% names(x)) || all(.crs(x) == .crs_1(x)),
            local_msg = "different crs in the bbox value",
            msg = "invalid bbox value"
        )

    x
}

.transform_extent <- function(xmin, xmax, ymin, ymax, crs, to_crs) {

    sfc <- sf::st_sfc(
        sf::st_polygon(list(rbind(c(xmin, ymax),
                                  c(xmax, ymax),
                                  c(xmax, ymin),
                                  c(xmin, ymin),
                                  c(xmin, ymax)))),
        crs = crs)

    bbox <- suppressWarnings(
        c(sf::st_bbox(sf::st_transform(sfc, crs = to_crs))))

    .crs(bbox) <- to_crs

    bbox
}

.bbox_to_longlat <- function(x, crs = NULL) {

    .check_bbox(x)

    if (is.null(crs) && "crs" %in% names(x))
        crs <- .crs(x)

    if (length(crs) == 1)
        crs <- rep(crs, length(.xmin(x)))

    .check_that(
        length(crs) == length(.xmin(x)),
        local_msg = paste("length of from_crs should be 1 or",
                          length(.xmin(x))),
        msg = "invalid from_crs parameter")

    bbox_longlat <- purrr::pmap_dfr(
        .xmin(x), .xmax(x), .ymin(x), .ymax(x), crs,
        function(xmin, xmax, ymin, ymax, crs) {
            tibble::tibble_row(
                .transform_extent(xmin = xmin,
                                  xmax = xmax,
                                  ymin = ymin,
                                  ymax = ymax,
                                  crs = crs,
                                  to_crs = "EPSG:4326"))
        })

    x[.bbox_vars(x)] <- c(bbox_longlat)[.bbox_vars(x)]

    x
}

.bbox_from_longlat <- function(x, crs = NULL) {

    .check_bbox(x)

    if ("crs" %in% names(x))
        crs <- .crs(x)

    .check_that(
        length(crs) == 1 || length(crs) == length(.xmin(x)),
        local_msg = paste("length of from_crs should be 1 or",
                          length(.xmin(x))),
        msg = "invalid from_crs parameter")

    if (length(crs) == 1)
        crs <- rep(crs, length(.xmin(x)))

    bbox_longlat <- purrr::pmap_dfr(
        .xmin(x), .xmax(x), .ymin(x), .ymax(x), crs,
        function(xmin, xmax, ymin, ymax, crs) {
            tibble::as_tibble_row(
                .transform_extent(xmin = xmin,
                                  xmax = xmax,
                                  ymin = ymin,
                                  ymax = ymax,
                                  crs = "EPSG:4326",
                                  to_crs = crs))
        })

    x[.bbox_vars(x)] <- c(bbox_longlat)[.bbox_vars(x)]

    x
}

.bbox_inner <- function(x) {

    .check_bbox(x, same_crs = TRUE)

    bbox <- c(xmin = max(.xmin(x)),
              xmax = min(.xmax(x)),
              ymin = max(.ymin(x)),
              ymax = min(.ymax(x)))

    if ("crs" %in% names(x))
        .crs(bbox) <- .crs(x)

    bbox
}

.bbox_outter <- function(x) {

    .check_bbox(x, same_crs = TRUE)

    bbox <- c(xmin = min(.xmin(x)),
              xmax = max(.xmax(x)),
              ymin = min(.ymin(x)),
              ymax = max(.ymax(x)))

    if ("crs" %in% names(x))
        .crs(bbox) <- .crs(x)

    bbox
}


