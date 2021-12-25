#' @title Bounding box (bbox)
#'
#' @name sits_bbox
#'
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#' @author Rolf Simoes, \email{rolf.simoes@@inpe.br}
#'
#' @description
#' Get the bounding box of a given data.
#'
#' @param data   valid sits tibble or sits cube objects (inherited from `tbl`).
#' @param wgs84  reproject resulting bbox to EPSG:4326
#' (take effect only for sits data cubes).
#' @param ...    Additional parameters (not implemented).
#'
#' @return return a `bbox` object (inherited from `tbl`) containing
#' `("xmin", "xmax", "ymin", "ymax", "crs")` variables
#'
#' @export
sits_bbox <- function(data, wgs84 = FALSE, ...) {

    .check_set_caller("sits_bbox")

    # get the meta-type (sits or cube)
    data <- .config_data_meta_type(data)

    UseMethod("sits_bbox", data)
}

#' @rdname sits_bbox
#' @description
#' `sits_bbox.sits()`: obtain a bbox of time series limits (in EPSG:4326).
#' @export
sits_bbox.sits <- function(data, ...) {

    # is the data a valid set of time series
    .sits_tibble_test(data)

    # create and return the bounding box
    bbox <- .bbox(data)

    return(bbox)
}

#' @rdname sits_bbox
#' @description
#' `sits_bbox.sits_cube()`: obtain a bbox of a data cube. If data cube have
#' more than one crs the result will be converted to EPSG:4326
#' @export
sits_bbox.sits_cube <- function(data, wgs84 = FALSE, ...) {

    # pre-condition
    .check_cube(data)

    # get bbox
    bbox <- .bbox(data)

    # convert to EPSG:4326?
    if (wgs84)
        bbox <- .bbox_transform(bbox, to_crs = "EPSG:4326")

    # get cube bbox
    bbox <- .bbox_outter(bbox)

    return(bbox)
}

#' @title Bounding box (bbox)
#'
#' @name bbox_functions
#'
#' @keywords internal
#'
#' @description
#' Bounding box functions
#'
#' @param cube  a sits data cube object (inherited from `tbl`)
#' @param bbox  a `bbox` object (inherited from `tbl`) containing
#' `("xmin", "xmax", "ymin", "ymax", "crs")` variables
#'
NULL

#' @rdname bbox_functions
#' @description
#' `.sits_bbox_intersection()`: compute intersection between a bbox and a
#' `cube` object (inherited from `tbl`).
#' @return
#' `.sits_bbox_intersection():` a `bbox` object (inherited from `tbl`)
.sits_bbox_intersection <- function(bbox, cube) {

    # set caller to show in errors
    .check_set_caller(".sits_bbox_intersection")

    # this function should be called per tile
    .check_cube(cube, is_tile = TRUE)

    # compute intersection
    bbox <- .bbox_intersection(bbox, bbox = .bbox(cube))

    return(bbox)
}

#' @keywords internal
.bbox <- function(x, ...) {

    UseMethod(".bbox", x)
}

#' @keywords internal
#' @export
.bbox.default <- function(x, ..., default_crs = NULL) {

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

#' @keywords internal
#' @export
.bbox.sits <- function(x, ...) {

    bbox <- tibble::tibble(
        xmin = min(.lon(x)),
        xmax = max(.lon(x)),
        ymin = min(.lat(x)),
        ymax = max(.lat(x)),
        crs = "EPSG:4326")

    .check_bbox(bbox)

    bbox
}

#' @keywords internal
`.bbox<-` <- function(x, value) {

    UseMethod(".bbox<-", x)
}

#' @keywords internal
#' @export
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

#' @keywords internal
.bbox_inner <- function(x) {

    UseMethod(".bbox_inner", x)
}

#' @keywords internal
#' @export
.bbox_inner.default <- function(x) {

    .check_bbox(x, same_crs = TRUE)

    bbox <- .bbox(list(
        xmin = max(.xmin(x)),
        xmax = min(.xmax(x)),
        ymin = max(.ymin(x)),
        ymax = min(.ymax(x)),
        crs = .crs_1(x)))

    bbox
}

#' @keywords internal
.bbox_intersection <- function(x, bbox) {

    UseMethod(".bbox_intersection", x)
}

#' @keywords internal
#' @export
.bbox_intersection.default <- function(x, bbox) {

    .check_bbox(x)
    .check_bbox(bbox)

    .check_that(length(.xmin(bbox)) == 1,
                local_msg = "bbox should have length 1",
                msg = "invalid bbox value")

    # have all bbox the same crs?
    if (.bbox_same_crs(x, y = bbox)) {

        bbox <- .bbox(list(
            xmin = pmax(.xmin(x), .xmin(bbox)),
            xmax = pmin(.xmax(x), .xmax(bbox)),
            ymin = pmax(.ymin(x), .ymin(bbox)),
            ymax = pmin(.ymax(x), .ymax(bbox)),
            crs = .crs_1(x)))

        return(bbox)
    }

    # convert bbox to sfc
    bbox_sfc <- .ext_as_sfc(.xmin(bbox), .xmax(bbox), .ymin(bbox),
                            .ymax(bbox), .crs(bbox))

    .intersection_ext <- function(xmin, xmax, ymin, ymax, crs) {

        sfc <- .ext_as_sfc(xmin, xmax, ymin, ymax, crs)

        bbox <- suppressWarnings(
            tibble::as_tibble_row(c(sf::st_bbox(sf::st_intersection(
                sf::st_transform(bbox_sfc, crs = crs), sfc)))))

        .crs(bbox) <- crs

        bbox
    }

    bbox <- purrr::pmap_dfr(.bbox(x), .intersection_ext)

    .check_bbox(bbox)

    bbox
}

#' @keywords internal
.bbox_intersects <- function(x, bbox) {

    UseMethod(".bbox_intersects", x)
}


#' @keywords internal
#' @export
.bbox_intersects.default <- function(x, bbox) {

    .check_bbox(x)
    .check_bbox(bbox)

    .check_that(length(.xmin(bbox)) == 1,
                local_msg = "bbox should have length 1",
                msg = "invalid bbox value")

    # have all bbox the same crs?
    if (.bbox_same_crs(x, y = bbox)) {

        intersects <-
            pmax(.xmin(x), .xmin(bbox)) <= pmin(.xmax(x), .xmax(bbox)) &
            pmax(.ymin(x), .ymin(bbox)) <= pmin(.ymax(x), .ymax(bbox))

        return(intersects)
    }

    # convert bbox to sfc
    bbox_sfc <- .ext_as_sfc(.xmin(bbox), .xmax(bbox), .ymin(bbox),
                            .ymax(bbox), .crs(bbox))

    .intersect_ext <- function(xmin, xmax, ymin, ymax, crs) {

        sfc <- .ext_as_sfc(xmin, xmax, ymin, ymax, crs)

        intersects <- suppressWarnings(
            apply(sf::st_intersects(
                sf::st_transform(bbox_sfc, crs = crs), sfc), 1, any))

        intersects
    }

    intersects <- purrr::pmap_lgl(.bbox(x), .intersect_ext)

    intersects
}

#' @keywords internal
.bbox_outter <- function(x) {

    UseMethod(".bbox_outter", x)
}

#' @keywords internal
#' @export
.bbox_outter.default <- function(x) {

    .check_bbox(x)

    y <- x
    if (!.bbox_same_crs(x))
        y <- .bbox_transform(x, to_crs = "EPSG:4326")

    .check_bbox(y, same_crs = TRUE)

    bbox <- .bbox(list(
        xmin = min(.xmin(y)),
        xmax = max(.xmax(y)),
        ymin = min(.ymin(y)),
        ymax = max(.ymax(y)),
        crs = .crs_1(y)))

    if (!.bbox_same_crs(x))
        warning(paste("bbox have different crs",
                      "result transformed to EPSG:4326", sep = "\\n"),
                call. = FALSE)

    bbox
}

#' @keywords internal
.bbox_same_crs <- function(x, y = NULL) {

    UseMethod(".bbox_same_crs", x)
}

#' @keywords internal
#' @export
.bbox_same_crs.default <- function(x, y = NULL) {

    .same_crs(.crs(x)) &&
        (is.null(y) || (!is.null(y) && .same_crs(.crs(y)) &&
                            .same_crs(.crs_1(x), .crs_1(y))))
}

#' @keywords internal
.bbox_transform <- function(x, to_crs) {

    UseMethod(".bbox_transform", x)
}

#' @keywords internal
#' @export
.bbox_transform.default <- function(x, to_crs) {

    .check_bbox(x)

    # pre-condition
    .check_length(to_crs, len_min = 1, len_max = 1,
                  "invalid 'to_crs' parameter")

    .transform_ext <- function(xmin, xmax, ymin, ymax, crs, to_crs) {

        sfc <- .ext_as_sfc(xmin, xmax, ymin, ymax, crs)

        bbox <- suppressWarnings(
            tibble::as_tibble_row(c(sf::st_bbox(
                sf::st_transform(sfc, crs = to_crs)))))

        .crs(bbox) <- to_crs

        bbox
    }

    bbox <- purrr::pmap_dfr(.bbox(x), .transform_ext, to_crs = to_crs)

    .check_bbox(bbox)

    bbox
}

#' @keywords internal
.ext_as_sfc <- function(xmin, xmax, ymin, ymax, crs) {

    .check_package("sf")

    sf::st_sfc(
        sf::st_polygon(list(rbind(c(xmin[[1]], ymax[[1]]),
                                  c(xmax[[1]], ymax[[1]]),
                                  c(xmax[[1]], ymin[[1]]),
                                  c(xmin[[1]], ymin[[1]]),
                                  c(xmin[[1]], ymax[[1]])))),
        crs = crs[[1]])
}

#' @keywords internal
.same_crs <- function(crs1, crs2 = NULL) {
    if (is.null(crs2))
        return(all(crs1 == crs1[[1]]))
    all(crs1 == crs2)
}
