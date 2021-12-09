#
#
# data <- list(EVI = terra::rast(cbers_cube$file_info[[1]]$path[[1]]),
#              NDVI = terra::rast(cbers_cube$file_info[[1]]$path[[2]]))
#
#
# test <- function(bands, ...) {
#
#     used_bands <- character(0)
#
#     detect_bands <- lapply(bands, function(band) {
#         quote({
#             used_bands <<- c(used_bands, band)
#             1
#         })
#     })
#
#     names(detect_bands) <- bands
#
#     expr <- substitute(list(...), env = environment())
#
#     return(eval(expr, envir = detect_bands))
#
#     return(used_bands)
# }
#
# test(c("NDVI", "EVI", "B01"), TEST = NDVI * 2)


sits_apply <- function(data, ...) {

    expr <- substitute(list(...), env = environment())

    .apply(data, expr = expr)
}


.apply <- function(data, ...) {

    .check_set_caller(".apply")

    UseMethod(".apply", data)
}

#' @export
.apply.list <- function(data, ...) {

    # pre-condition
    .check_lst(data, min_len = 1, msg = "invalid 'data' value")

    # capture expression
    expr <- substitute(list(...), env = environment())

    # evaluate
    value <- .check_error({
        eval(expr, envir = data)
    }, msg = "evaluation error")

    return(value)
}

#' @export
.apply.sits <- function(data, ...) {

    # pre-condition
    .check_that(inherits(data, "sits"),
                local_msg = "value is not a valid sits tibble",
                msg = "invalid 'data' value")

    .local_mutate <- function(x)
        tibble::as_tibble(c(x, .apply.list(x, ...)))

    value <- .sits_fast_apply(data, col = "time_series", .local_mutate)

    return(value)
}

.apply.raster_cube <- function(data,
                               bbox,
                               res, ...,
                               output_dir = ".",
                               res_method = "bilinear") {

    # pre-condition
    .check_that(inherits(data, "raster_cube"),
                local_msg = "value is not a valid sits tibble",
                msg = "invalid 'data' value")

    .check_that(.cube_is_regular(data),
                msg = "cube is not regular")

    # slide through tiles
    if (nrow(data) > 1)
        slider::slide_dfr(data, .apply.raster_cube, bbox = bbox,
                          resolution = res, output_dir = output_dir)

    # traverse each scene
    # parallel processing
    .sits_parallel_map(sits_timeline(data), function(date) {

        # pull file_info // data has one row at this point
        file_info <- .cube_file_info(data)

        # change this to "fid"
        item <- dplyr::filter(file_info, .data[["date"]] == !!date)

        # create blocks to be processed
        # TODO: get block size from config file
        blocks <- .bbox_split(bbox, ylen = 1024 * res)

        # capture expression
        expr <- substitute(list(...), env = environment())

        # for each block
        output_blocks <- purrr::map(blocks, function(block) {

            # download rasters of each band
            local_block_bands <- purrr::map(item[["bands"]], function(band) {

                # filter band
                item_band <- dplyr::filter(item, .data[["band"]] == !!band)

                # create a local file path
                block_path <- .local_block_path(path = item_band[["path"]],
                                                block = block)

                # skip if file exists
                if (file.exists(block_path))
                    return(block_path)

                # download the image to local path
                gdalUtilities::gdalwarp(
                    srcfile = path,
                    dstfile = block_path,
                    of = "GTiff",
                    te = bbox[c("xmin", "ymin", "xmax", "ymax")],
                    tr = res[c("xres", "yres")],
                    r = res_method)

                # return parallel class error to try again
                # if (!file.exists(block_path))
                #     return()

                block_path
            })

            # read local files
            bands <- purrr::map(local_paths, function(path) {

                .raster_read_rast(file = path)
            })

            # set band names
            names(bands) <- item[["band"]]

            # process
            results <- .apply(bands, ...)

            # save results

            # output file
            output_block_file
        })

        merged_file <- gdalUtilities::gdalwarp(srcfile = output_block_files,
                                               dstfile = output_files)

        # produce each output file
        output_files <- purrr::map(output_files, function(output_file) {


        })

        item
    }, progress = TRUE)
}


# generate a local file path based on input path
.bbox_local_path <- function(path, bbox, output_dir) {

    path <- gsub("^([^?]+)\\?.*$", "\\1", path)
    file <- gsub("^.*/(.*)$", "\\1", path)
    file_name <- gsub("^(.*)\\..+$", "\\1", file)
    file_ext <- gsub("^.*\\.(.+)$", "\\1", file)
    if (!is.null(bbox))
        file_name <- paste(file_name, paste(bbox, collapse = "_"),
                           sep = "_")
    paste(output_dir, paste(file_name, file_ext, sep = "."), sep = "/")
}

.bbox_split <- function(bbox, ...,
                        xlen = NULL,
                        ylen = NULL) {

    .check_num(bbox, len_min = 4, len_max = 4, is_named = TRUE,
               msg = "invalid 'bbox' parameter")

    .check_chr_contains(names(bbox),
                        contains = c("xmin", "xmax", "ymin", "ymax"),
                        msg = "invalid 'bbox' parameter")

    if (is.null(xlen))
        xlen <- bbox[["xmax"]] - bbox[["xmin"]]

    if (is.null(ylen))
        ylen <- bbox[["ymax"]] - bbox[["ymin"]]

    xsplits <- unique(c(seq(from = bbox[["xmin"]],
                            to = bbox[["xmax"]], by = xlen),
                        bbox[["xmax"]]))

    ysplits <- unique(c(seq(from = bbox[["ymin"]],
                            to = bbox[["ymax"]], by = ylen),
                        bbox[["ymax"]]))

    splits <- purrr::map_dfr(seq_len(length(xsplits) - 1), function(i) {
        dplyr::tibble(xmin = xsplits[[i]],
                      xmax = xsplits[[i + 1]],
                      ymin = ysplits[seq_len(length(ysplits) - 1)],
                      ymax = ysplits[-1])
    })

    slider::slide(splits, unlist)
}

.raster_resample <- function(path, ...,
                             bbox = NULL,
                             res = NULL,
                             res_method = "bilinear",
                             data_type = "INT2S",
                             output_dir = ".") {

    .check_set_caller(".raster_resample")

    if (!is.null(bbox)) {

        .check_num(bbox, len_min = 4, len_max = 4, is_named = TRUE,
                   msg = "invalid 'bbox' parameter")

        .check_chr_contains(names(bbox),
                            contains = c("xmin", "xmax", "ymin", "ymax"),
                            msg = "invalid 'bbox' parameter")

        bbox <- bbox[c("xmin", "xmax", "ymin", "ymax")]
    }

    if (!is.null(res)) {

        .check_num(res, len_min = 2, len_max = 2, is_named = TRUE,
                   msg = "invalid 'res' parameter")

        .check_chr_contains(names(res),
                            contains = c("xres", "yres"),
                            msg = "invalid 'res' parameter")

        res <- res[c("xres", "yres")]
    }

    .check_chr(res_method, len_min = 1, len_max = 1,
               msg = "invalid 'res_method' parameter")

    .check_chr_within(res_method, within = c("bilinear", "near"))

    .check_chr(output_dir, len_min = 1, len_max = 1,
               msg = "invalid 'output_dir' parameter")

    .check_that(dir.exists(output_dir),
                local_msg = "directory does not exist",
                msg = "invalid 'output_dir' parameter")

    # create a local file path
    block_path <- .bbox_local_path(path = path,
                                   bbox = bbox,
                                   output_dir = output_dir)

    # skip if file exists
    if (file.exists(block_path))
        return(block_path)

    # open original raster
    r_src <- terra::rast(path)

    # open destination raster
    if (is.null(bbox) && is.null(res))
        r_dst <- terra::rast(r_src)
    else if (is.null(bbox))
        r_dst <- terra::rast(crs = terra::crs(r_src),
                             extent = terra::ext(r_src),
                             resolution = res)
    else if (is.null(res))
        r_dst <- terra::rast(crs = terra::crs(r_src),
                             extent = terra::ext(bbox),
                             resolution = terra::res(r_src))
    else
        r_dst <- r_dst <- terra::rast(crs = terra::crs(r_src),
                                      extent = terra::ext(bbox),
                                      resolution = res)

    # resample original into destination raster
    terra::resample(r_src, r_dst, method = res_method,
                    filename = block_path,
                    gdal = .config_gtiff_default_options(),
                    datatype = data_type)

    # return parallel class error to try again
    .check_file(block_path, msg = "post-condition error")

    return(r_dst)
}


.mutate.matrix <- function(data, ...) {

    .check_set_caller(".mutate")

    # pre-condition
    .check_that(is.matrix(data),
                local_msg = "data is not matrix",
                msg = "invalid data value")

    result <- apply(data, MARGIN = 2, FUN = fn, ...)

    # post-condition
    .check_that(identical(dim(result), dim(data)),
                local_msg = paste("result value dimension should be",
                                  paste(dim(data), collapse = "x"),
                                  "instead of",
                                  paste(dim(result), collapse = "x")),
                msg = "invalid result dimension")

    return(result)
}
