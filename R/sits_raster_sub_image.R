#' @title Informs if a spatial ROI intersects a data cube
#' @name .sits_raster_sub_image_intersects
#' @keywords internal

#' @param  cube            input data cube.
#' @param  roi             spatial region of interest
#' @return                 logical
#'
.sits_raster_sub_image_intersects <- function(cube, roi) {

    # if roi is null, returns TRUE
    if (purrr::is_null(roi)) return(TRUE)

    # check if roi is a sf object
    if (inherits(roi, "sf")) {

        # check for roi crs
        if (is.null(sf::st_crs(roi))) {
            stop(".sits_raster_sub_image_intersects: invalid roi crs",
                 call. = FALSE)
        }

        # reproject roi to cube crs
        roi <- suppressWarnings(sf::st_transform(roi, cube$crs[[1]]))

        # region of cube tile
        df <- data.frame(
            X = c(.xmin(cube)[[1]], .xmax(cube)[[1]],
                  .xmax(cube)[[1]], .xmin(cube)[[1]]),
            Y = c(.ymin(cube)[[1]], .ymin(cube)[[1]],
                  .ymax(cube)[[1]], .ymax(cube)[[1]])
        )

        # compute tile polygon
        sf_region <-
            sf::st_as_sf(df, coords = c("X", "Y"), crs = cube$crs[[1]]) %>%
            dplyr::summarise(geometry = sf::st_combine(geometry)) %>%
            sf::st_cast("POLYGON") %>%
            suppressWarnings()

        # check for intersection
        return(apply(sf::st_intersects(sf_region, roi), 1, any))
    }

    # if the ROI is defined, calculate the bounding box
    bbox_roi <- .sits_roi_bbox(roi, cube)

    # calculate the intersection between the bbox of the ROI and the cube
    bbox_in <- .sits_bbox_intersect(bbox_roi, cube)

    return(!purrr::is_null(bbox_in))
}
#' @title Find the dimensions and location of a spatial ROI in a data cube
#' @name .sits_raster_sub_image
#' @keywords internal

#' @param  cube            input data cube.
#' @param  roi             spatial region of interest
#' @return                 vector with information on the subimage
#'
.sits_raster_sub_image <- function(cube, roi) {

    # if the ROI is defined, calculate the bounding box
    bbox_roi <- .sits_roi_bbox(roi, cube)

    # calculate the intersection between the bbox of the ROI and the cube
    bbox_in <- .sits_bbox_intersect(bbox_roi, cube)

    # return the sub_image
    sub_image <- .sits_raster_sub_image_from_bbox(bbox_in, cube)

    return(sub_image)
}
#' @title Find the dimensions of the sub image without ROI
#' @name .sits_raster_sub_image_default
#' @keywords internal

#' @param  cube            input data cube.
#' @param  sf_region       spatial region of interest (sf_object)
#' @return                 vector with information on the subimage
.sits_raster_sub_image_default <- function(cube) {

    # by default, the sub_image has the same dimension as the main cube

    size <- .cube_size(cube)
    bbox <- .cube_tile_bbox(cube)

    sub_image <- c(first_row = 1,
                   first_col = 1,
                   nrows = .nrows(size),
                   ncols = .ncols(size),
                   xmin = .xmin(bbox),
                   xmax = .xmax(bbox),
                   ymin = .ymin(bbox),
                   ymax = .ymax(bbox))

    return(sub_image)
}

#' @title Extract a sub_image from a bounding box and a cube
#' @name .sits_raster_sub_image_from_bbox
#' @keywords internal
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @param bbox           bounding box for a region of interest
#' @param cube           data cube
#' @return               sub_image with additional info on first row,
#'                       first col, nrows, ncols
#'
.sits_raster_sub_image_from_bbox <- function(bbox, cube) {

    # pre-conditions
    .check_num(.xmin(bbox), max = .xmax(bbox),
               msg = "invalid bbox value")

    .check_num(.ymin(bbox), max = .ymax(bbox),
               msg = "invalid bbox value")

    .check_num(.xmin(bbox), min = .xmin(cube), max = .xmax(cube),
               msg = "bbox value is outside the cube")

    .check_num(.xmax(bbox), min = .xmin(cube), max = .xmax(cube),
               msg = "bbox value is outside the cube")

    .check_num(.ymin(bbox), min = .ymin(cube), max = .ymax(cube),
               msg = "bbox value is outside the cube")

    .check_num(.ymax(bbox), min = .ymin(cube), max = .ymax(cube),
               msg = "bbox value is outside the cube")

    # get the resolution
    # throw an error if resolution are not the same
    # for all bands of the cube
    res   <- .cube_resolution(cube)

    # get ncols and nrows
    # throw an error if size are not the same
    size  <- .cube_size(cube)

    # set initial values
    si <- c(first_row = 1, first_col = 1,
            nrows = .nrows(size), ncols = .ncols(size),
            xmin = .xmin(cube), xmax = .xmax(cube),
            ymin = .ymin(cube), ymax = .ymax(cube))

    # find the first row (remember that rows runs from top to bottom and
    # Y coordinates increase from bottom to top)
    si[["first_row"]] <- floor((.ymax(cube) - .ymax(bbox)) / .yres(res)) + 1

    # adjust to fit bbox in cube resolution
    .ymax(si) <- .ymax(cube) - .yres(res) * (si[["first_row"]] - 1)

    # find the first col (remember that rows runs from left to right and
    # X coordinates increase from left to right)
    si[["first_col"]] <- floor((.xmin(bbox) - .xmin(cube)) / .xres(res)) + 1

    # adjust to fit bbox in cube resolution
    .xmin(si) <- .xmin(cube) + .xres(res) * (si[["first_col"]] - 1)

    # find the number of rows (remember that rows runs from top to bottom and
    # Y coordinates increase from bottom to top)
    .nrows(si) <- floor((.ymax(bbox) - .ymin(bbox)) / .yres(res)) + 1

    # adjust to fit bbox in cube resolution
    .ymin(si) <- .ymax(si) - .yres(res) * .nrows(si)

    .ncols(si) <- floor((.xmax(bbox) - .xmin(bbox)) / .xres(res)) + 1

    # adjust to fit bbox in cube resolution
    .xmax(si) <- .xmin(si) + .xres(res) * .ncols(si)

    # pre-conditions
    .check_num(.xmin(si), max = .xmax(si),
               msg = "invalid subimage value")

    .check_num(.ymin(si), max = .ymax(si),
               msg = "invalid subimage value")

    .check_num(.xmin(si), min = .xmin(cube), max = .xmax(cube),
               msg = "invalid subimage value")

    .check_num(.xmax(si), min = .xmin(cube), max = .xmax(cube),
               msg = "invalid subimage value")

    .check_num(.ymin(si), min = .ymin(cube), max = .ymax(cube),
               msg = "invalid subimage value")

    .check_num(.ymax(si), min = .ymin(cube), max = .ymax(cube),
               msg = "invalid subimage value")

    return(si)
}
