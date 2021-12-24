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
    lon_max <- max(.lon(data))
    lon_min <- min(.lon(data))
    lat_max <- max(.lat(data))
    lat_min <- min(.lat(data))
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

    if (!wgs84)
        .check_that(length(unique(data[["crs"]])) == 1,
                    local_msg = "use `wgs84 = TRUE` for a global bbox",
                    msg = "cube has more than one projection")

    # create and return the bounding box
    bbox <- c(xmin = min(data[["xmin"]]),
              xmax = max(data[["xmax"]]),
              ymin = min(data[["ymin"]]),
              ymax = max(data[["ymax"]]))

    # convert to WGS84?
    if (wgs84)
        bbox <- .sits_coords_to_bbox(
            xmin = bbox[["xmin"]],
            xmax = bbox[["xmax"]],
            ymin = bbox[["ymin"]],
            ymax = bbox[["ymax"]],
            crs = data[["crs"]][[1]])

    return(bbox)
}

.sits_coords_to_bbox <- function(xmin, xmax, ymin, ymax, crs) {

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

    if (all(names(x) %in% c("xmin", "ymin", "xmax", "ymax")))
        return(.bbox(x, default_crs = crs))


    if (all(names(x) %in% c("lon_min", "lat_min", "lon_max", "lat_max"))) {

        ll_names <- c("xmin", "ymin", "xmax", "ymax")
        names(ll_names) <- c("lon_min", "lat_min", "lon_max", "lat_max")
        names(x) <- unname(ll_names[names(x)])

        return(.bbox(x, default_crs = "EPSG:4326"))
    }

    stop("invalid roi value", call. = FALSE)
}

.check_bbox <- function(x, same_crs = FALSE) {

    .check_chr_contains(
        names(x), contains = c("xmin", "xmax", "ymin", "ymax", "crs"),
        msg = "object does not have all bbox variables")

    .check_that(
        length(.crs(x)) == length(.xmin(x)),
        local_msg = "length of crs differs from xmin",
        msg = "invalid bbox value")

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
            all(.crs(x) == .crs_1(x)),
            local_msg = "different crs in the bbox value",
            msg = "invalid bbox value"
        )

    x
}

.bbox <- function(x, default_crs = NULL) {

    UseMethod(".bbox", x)
}

.bbox.default <- function(x, default_crs = NULL) {

    if ("crs" %in% names(x))
        default_crs <- .crs(x)

    .check_that(
        length(default_crs) == length(.xmin(x)) ||
            length(default_crs) == 1,
        local_msg = paste("length of crs should be 1 or", length(.xmin(x))),
        msg = "invalid crs parameter"
    )

    bbox <- tibble::tibble(
        xmin = .xmin(x),
        xmax = .xmax(x),
        ymin = .ymin(x),
        ymax = .ymax(x),
        crs = default_crs)

    .check_bbox(bbox)

    bbox
}

`.bbox<-` <- function(x, value) {

    UseMethod(".bbox<-", x)
}

`.bbox<-.default` <- function(x, value) {

    .check_bbox(value)

    .check_that(is.list(x),
                local_msg = "object does not support all bbox variables",
                msg = "invalid object")

    .xmin(x) <- .xmin(value)
    .xmax(x) <- .xmax(value)
    .ymin(x) <- .ymin(value)
    .ymax(x) <- .ymax(value)
    .crs(x) <- .crs(value)

    x
}

.bbox_inner <- function(x) {

    .check_bbox(x, same_crs = TRUE)

    bbox <- tibble::tibble(
        xmin = max(.xmin(x)),
        xmax = min(.xmax(x)),
        ymin = max(.ymin(x)),
        ymax = min(.ymax(x)),
        crs = .crs_1(x))

    .check_bbox(bbox)

    bbox
}

.bbox_intersection <- function(x, bbox) {

    .check_bbox(x)
    .check_bbox(bbox)

    .check_that(length(.xmin(bbox)) == 1,
                local_msg = "bbox should have length 1",
                msg = "invalid bbox value")

    bbox <- tibble::tibble(
        xmin = pmax(.xmin(x), .xmin(bbox)),
        xmax = pmin(.xmax(x), .xmax(bbox)),
        ymin = pmax(.ymin(x), .ymin(bbox)),
        ymax = pmin(.ymax(x), .ymax(bbox)),
        crs = .crs_1(x)
    )

    .check_bbox(bbox)

    bbox
}

.bbox_intersects <- function(x, bbox) {

    .check_bbox(x)
    .check_bbox(bbox)

    .check_that(length(.xmin(bbox)) == 1,
                local_msg = "bbox should have length 1",
                msg = "invalid bbox value")

    intersects <- pmin(.xmax(x), .xmax(bbox)) <= pmax(.xmin(x), .xmin(bbox)) &
        pmax(.ymin(x), .ymin(bbox)) <= pmin(.ymax(x), .ymax(bbox))

    intersects
}
.bbox_outter <- function(x) {

    .check_bbox(x)

    y <- x
    if (!.bbox_same_crs(x))
        y <- .bbox_transform(x, to_crs = "EPSG:4326")

    .check_bbox(y, same_crs = TRUE)

    bbox <- tibble::tibble(
        xmin = min(.xmin(y)),
        xmax = max(.xmax(y)),
        ymin = min(.ymin(y)),
        ymax = max(.ymax(y)),
        crs = .crs_1(y))

    .check_bbox(bbox)

    if (!.bbox_same_crs(x))
        warning(paste("bbox have different crs",
                      "result transformed to EPSG:4326", sep = "\\n"),
                call. = FALSE)

    bbox
}

.bbox_same_crs <- function(x) {

    all(.crs(x) == .crs_1(x))
}


.bbox_transform <- function(x, to_crs) {

    .check_bbox(x)

    .transform_extent <- function(xmin, xmax, ymin, ymax, crs, to_crs) {

        sfc <- sf::st_sfc(
            sf::st_polygon(list(rbind(c(xmin, ymax),
                                      c(xmax, ymax),
                                      c(xmax, ymin),
                                      c(xmin, ymin),
                                      c(xmin, ymax)))),
            crs = crs)

        bbox <- suppressWarnings(
            tibble::as_tibble_row(c(sf::st_bbox(
                sf::st_transform(sfc, crs = to_crs)))))

        .crs(bbox) <- to_crs

        bbox
    }

    bbox <- purrr::pmap_dfr(.bbox(x), .transform_extent, to_crs = to_crs)

    bbox
}

