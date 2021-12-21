#' @title Apply a function on a set of time series
#'
#' @name sits_apply
#'
#' @author Rolf Simoes, \email{rolf.simoes@@inpe.br}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @description
#' Apply a named expression to a sits cube or a sits tibble
#' to be evaluated and generate new bands (indices). In the case of sits
#' cubes, it materializes a new band in `output_dir` using `gdalcubes`.
#'
#' @param data          Valid sits tibble or cube
#' @param output_dir    Directory where files will be saved.
#' @param ...           Named expressions to be evaluated.
#'
#' @return A sits tibble or a sits cube with new bands.
#'
#' @examples
#' # Get a time series
#' # apply a normalization function
#'
#' point2 <-
#'   sits_select(point_mt_6bands, "NDVI") %>%
#'   sits_apply(NDVI_norm = (NDVI - min(NDVI)) / (max(NDVI) - min(NDVI))
#' )
#'
NULL

#' @rdname sits_apply
#' @export
sits_apply <- function(data, ...) {

    UseMethod("sits_apply", data)
}

#' @rdname sits_apply
#' @export
sits_apply.sits <- function(data, ...) {

    .check_set_caller("sits_apply.sits")

    .sits_fast_apply(data, col = "time_series", fn = dplyr::mutate, ...)
}

#' @rdname sits_apply
#' @export
sits_apply.raster_cube <- function(data, ...,
                                   impute_fn = sits_impute_linear(),
                                   output_dir = getwd()) {

    .check_set_caller("sits_apply.raster_cube")

    # capture dots as a list of quoted expressions
    list_expr <- lapply(substitute(list(...), env = environment()),
                        unlist, recursive = F)[-1]

    # get new band names from expression
    .check_lst(list_expr, min_len = 1, msg = "invalid expression value")

    # get all input bands
    in_bands <- .cube_bands(data)

    # find which bands are in input expressions
    char_exprs <- paste(unlist(lapply(list_expr, deparse)), collapse = " ")
    used_bands <- purrr::map_lgl(in_bands, grepl, char_exprs)

    # pre-condition
    .check_that(any(used_bands),
                local_msg = "no valid band was informed",
                msg = "invalid expression value")

    # select used bands
    in_bands <- in_bands[used_bands]

    # traverse each tile
    result <- slider::slide_dfr(data, function(tile) {

        # get file_info filtered by bands
        in_file_info <- .file_info(tile, band == !!in_bands)

        # traverse each date in file info
        out_file_info <- .foreach_file_info(
            in_file_info, group_by = "fid", order_by = "band",
            fn = function(file_info_fid) {


            # filter file_info by date
            in_files <- dplyr::filter(in_file_info, band == !!in_bands) %>%
                dplyr::select(band, path) %>%
                tidyr::pivot_wider(names_from = "band",
                                   values_from = "path")

            # load bands data
            in_values <- purrr::map(names(in_files), function(band) {

                # get file path
                file <- in_bands[[band]]

                # read the values
                values <- .raster_read_rast(file)

                # get the missing values, minimum values and scale factors
                missing_value <- .cube_band_missing_value(cube = cube,
                                                          band = band)
                minimum_value <- .cube_band_minimum_value(cube = cube,
                                                          band = band)
                maximum_value <- .cube_band_maximum_value(cube = cube,
                                                          band = band)

                # correct NA, minimum, maximum, and missing values
                values[values < minimum_value] <- NA
                values[values > maximum_value] <- NA
                values[values == missing_value] <- NA

                # scale the data set
                scale_factor <- .cube_band_scale_factor(cube, band = band)
                offset_value <- .cube_band_offset_value(cube = cube,
                                                        band = band)

                # compute scale and offset
                values <- scale_factor * values + offset_value

                # remove NA pixels
                if (!purrr::is_null(impute_fn) && any(is.na(values))) {
                    values <- impute_fn(values)
                }

                return(values)
            })

            # save each output value
            output_files <- purrr::map(names(list_expr), function(band) {

                file_prefix <- paste("cube", tile[["tile"]], band, date,
                                     sep = "_")
                file_name <- paste(file_prefix, "tif", sep = ".")
                file_path <- paste(output_dir, file_name, sep = "/")

                if (file.exists(file_path))
                    return(file_path)

                # evaluate expressions
                # TODO: internal var: 'default_scale_factor'
                out_values <- eval(list_expr[[band]], in_values) *
                    .config_get("raster_cube_scale_factor") +
                    .config_get("raster_cube_offset_value")

                # new raster
                r_obj <- .raster_new_rast(
                    nrows = in_file_info[["nrows"]][[1]],
                    ncols = in_file_info[["ncols"]][[1]],
                    xmin = in_file_info[["xmin"]][[1]],
                    xmax = in_file_info[["xmax"]][[1]],
                    ymin = in_file_info[["ymin"]][[1]],
                    ymax = in_file_info[["ymax"]][[1]],
                    nlayers = 1,
                    crs = tile[["crs"]])

                # set values
                r_obj <- .raster_set_values(r_obj, out_values)

                # write values
                .raster_write_rast(
                    r_obj = r_obj,
                    file = file_path,
                    format = "GTiff",
                    data_type = .config_get("raster_cube_data_type"),
                    gdal_options = .config_gtiff_default_options(),
                    overwrite = FALSE)

                return(file_path)
            })

            # clean memory
            gc()

            # prepare output file_info
            out_file_info <- tibble::tibble(
                fid = file_info[["fid"]][[1]],
                date = file_info[["date"]][[1]],
                band = bands,
                xmin = file_info[["xmin"]][[1]],
                xmax = file_info[["xmax"]][[1]],
                ymin = file_info[["ymin"]][[1]],
                ymax = file_info[["ymax"]][[1]],
                xres = file_info[["xres"]][[1]],
                yres = file_info[["yres"]][[1]],
                nrows = file_info[["nrows"]][[1]],
                ncols = file_info[["ncols"]][[1]],
                path = output_files)

            return(out_file_info)
        })


        # retrieve dates
        dates <- purrr::map(output_files, function(x) {
            dplyr::tibble(date = .gc_get_date(x))
        })

            tile_fid[["file_info"]][[1]] <-
                dplyr::filter(.cube_file_info(tile),
                              .data[["fid"]] == !!fid)

        file_info <- tidyr::unnest(tibble::tibble(
            date = dates,
            band = new_bands,
            res = .cube_resolution(tile),
            path = output_files
        ), cols = c("date", "path"))

            ic <- .gc_create_database(tile_fid,
                                      path_db = tempfile(fileext = ".db"))

            bands <- names(list_expr)

            .check_that(length(bands) == length(list_expr),
                        local_msg = "not all expressions have names",
                        msg = "invalid expressions parameters")


            cv <- .gc_create_cube_view(
                tile = tile_fid,
                period = "P1D",
                res = .cube_resolution(tile_fid),
                roi = NULL,
                toi = toi,
                agg_method = "first",
                resampling = "bilinear"
            )

            rc <- gdalcubes::raster_cube(ic, view = cv)

            output_files <- purrr::map(bands, function(band) {

                cc <- gdalcubes::apply_pixel(rc,
                                             expr = deparse(list_expr[[band]]),
                                             names = band)

                # file prefix
                prefix <- paste("cube", tile_fid[["tile"]], band, "", sep = "_")

                gdalcubes::write_tif(
                    cc,
                    dir = output_dir,
                    prefix = prefix,
                    creation_options = list("COMPRESS" = "LZW",
                                            "BIGTIFF" = "YES"),
                    pack = list(type = "int16", nodata = -9999,
                                scale = 1, offset = 0)
                )

                file_name <- paste0(output_dir, "/", prefix,
                                    .cube_file_info(tile_fid)[["date"]][[1]],
                                    ".tif")
                return(file_name)
            })

            file_info <- .cube_file_info(tile_fid)

            file_info <- tidyr::unnest(tibble::tibble(
                fid = file_info[["fid"]][[1]],
                date = file_info[["date"]][[1]],
                band = bands,
                xmin = file_info[["xmin"]][[1]],
                xmax = file_info[["xmax"]][[1]],
                ymin = file_info[["ymin"]][[1]],
                ymax = file_info[["ymax"]][[1]],
                xres = file_info[["xres"]][[1]],
                yres = file_info[["yres"]][[1]],
                nrows = file_info[["nrows"]][[1]],
                ncols = file_info[["ncols"]][[1]],
                path = output_files
            ), cols = c("date", "path"))

            file_info_fid <- dplyr::bind_rows(tile_fid[["file_info"]][[1]],
                                              file_info) %>%
                dplyr::arrange(date, band)

            return(file_info_fid)
        })

        return(tile)
    })

    return(result)
}

#' @rdname sits_apply
#' @keywords internal
.apply_across <- function(data, fn, ...) {

    .check_set_caller(".apply_across")

    fn_across <- fn
    .sits_fast_apply(data, col = "time_series", fn = function(x, ...) {
        dplyr::mutate(x, dplyr::across(dplyr::matches(sits_bands(data)),
                                       fn_across, ...))
    }, ...)
}
