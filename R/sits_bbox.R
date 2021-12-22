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
