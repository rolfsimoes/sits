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

.raster_preprocess <- function(file,
                               block,
                               missing_value,
                               minimum_value,
                               maximum_value,
                               scale_factor,
                               offset_value, ...,
                               impute_fn = NULL,
                               filter_fn = NULL,
                               normalize_fn = NULL) {

    # pre-conditions
    .check_block(block, msg = "invalid block parameter")

    .check_num(missing_value, len_min = 1, len_max = 1,
               msg = "invalid missing_value parameter")

    .check_num(minimum_value, len_min = 1, len_max = 1,
               msg = "invalid minimum_value parameter")

    .check_num(maximum_value, len_min = 1, len_max = 1,
               msg = "invalid maximum_value parameter")

    .check_num(scale_factor, len_min = 1, len_max = 1,
               msg = "invalid scale_factor parameter")

    .check_num(offset_value, len_min = 1, len_max = 1,
               msg = "invalid offset_value parameter")

    # read the values
    values <- .raster_read_rast(file, block = block)

    # get input dimension
    d <- dim(values)

    # correct NA, minimum, maximum, and missing values
    values[values == missing_value] <- NA
    values[values < minimum_value] <- NA
    values[values > maximum_value] <- NA

    # impute NA pixels
    if (!is.null(impute_fn) && any(is.na(values))) {

        .check_that(inherits(impute_fn, "function"))

        values <- impute_fn(values)
    }

    # compute scale and offset
    values <- scale_factor * values + offset_value

    # filter the data
    if (!(is.null(filter_fn))) {

        .check_that(inherits(filter_fn, "function"))

        values <- filter_fn(values)
    }

    # normalize the data
    if (!is.null(normalize_fn)) {

        .check_that(inherits(normalize_fn, "function"))

        values <- normalize_fn(values)
    }

    # post-conditions
    .check_that(all(d == dim(values)),
                local_msg = "output has different dimension than input",
                msg = "invalid output values")

    return(values)
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

    # TODO: dryrun expression

    # traverse each tile
    result <- slider::slide_dfr(data, function(tile) {

        # get file_info filtered by bands
        in_file_info <- .cube_file_info(tile)

        # get all fids of file_info
        fids <- unique(in_file_info[["fid"]])

        # traverse each date in file info
        out_file_info <- purrr::map_dfr(fids, function(fid) {

            # filter fid
            in_file_info_fid <- in_file_info %>%
                dplyr::filter(.data[["fid"]] == !!fid)

            # filter bands
            in_files <- in_file_info_fid %>%
                dplyr::filter(.data[["band"]] == !!in_bands) %>%
                dplyr::select(dplyr::all_of(c("band", "path"))) %>%
                tidyr::pivot_wider(names_from = "band",
                                   values_from = "path")

            # get input bands
            in_bands <- names(in_files)

            # load bands data
            in_values <- purrr::map(in_bands, function(band) {

                # get the missing values, minimum values and scale factors
                missing_value <- .cube_band_missing_value(data, band = band)
                minimum_value <- .cube_band_minimum_value(data, band = band)
                maximum_value <- .cube_band_maximum_value(data, band = band)

                # scale the data set
                scale_factor <- .cube_band_scale_factor(data, band = band)
                offset_value <- .cube_band_offset_value(data, band = band)

                # preprocess data
                .raster_preprocess(in_files[[band]],
                                   missing_value = missing_value,
                                   minimum_value = minimum_value,
                                   maximum_value = maximum_value,
                                   scale_factor = scale_factor,
                                   offset_value = offset_value)

            })

            # set band names
            names(in_values) <- names(in_files)

            # get new band names
            new_bands <- toupper(unique(names(list_expr)))

            # pre-condition
            .check_length(new_bands, len_min = length(list_expr),
                          len_max = length(list_expr),
                          msg = "invalid new bands name")

            # save each output value
            output_files <- purrr::map_chr(new_bands, function(new_band) {

                file_prefix <- paste("cube", tile[["tile"]], new_band,
                                     in_file_info_fid[["date"]][[1]],
                                     sep = "_")
                file_name <- paste(file_prefix, "tif", sep = ".")
                file_path <- paste(output_dir, file_name, sep = "/")

                if (file.exists(file_path))
                    return(file_path)

                # evaluate expressions, scale and offset values
                out_values <- eval(list_expr[[new_band]], in_values) /
                    .config_get("raster_cube_scale_factor") -
                    .config_get("raster_cube_offset_value")

                # new raster
                r_obj <- .raster_new_rast(
                    nrows = .nrows(in_file_info_fid)[[1]],
                    ncols = .ncols(in_file_info_fid)[[1]],
                    xmin = .xmin(in_file_info_fid)[[1]],
                    xmax = .xmax(in_file_info_fid)[[1]],
                    ymin = .ymin(in_file_info_fid)[[1]],
                    ymax = .ymax(in_file_info_fid)[[1]],
                    nlayers = 1,
                    crs = .crs(tile))

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
            out_file_info_fid <- tibble::tibble(
                fid = in_file_info_fid[["fid"]][[1]],
                date = in_file_info_fid[["date"]][[1]],
                band = new_bands,
                xmin = .xmin(in_file_info_fid)[[1]],
                xmax = .xmax(in_file_info_fid)[[1]],
                ymin = .ymin(in_file_info_fid)[[1]],
                ymax = .ymax(in_file_info_fid)[[1]],
                xres = .xres(in_file_info_fid)[[1]],
                yres = .yres(in_file_info_fid)[[1]],
                nrows = .nrows(in_file_info_fid)[[1]],
                ncols = .ncols(in_file_info_fid)[[1]],
                path = output_files)

            return(out_file_info_fid)
        })

        out_tile <- .cube_create(
            source = .cube_source(tile),
            collection = .cube_collection(tile),
            tile = .cube_tiles(tile),
            xmin = .xmin(tile),
            xmax = .xmax(tile),
            ymin = .ymin(tile),
            ymax = .ymax(tile),
            crs = .crs(tile),
            labels = tile[["labels"]],
            file_info = out_file_info)

        return(out_tile)
    })

    class(result) <- c("raster_cube", class(result))

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
