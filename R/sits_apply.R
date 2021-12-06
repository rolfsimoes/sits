sits_apply <- function(data, ...) {

    expr <- substitute(list(...), env = environment())

    .apply(data, expr = expr)
}


.apply <- function(data, expr) {

    .check_set_caller(".apply")

    UseMethod(".apply", data)
}

#' @export
.apply.list <- function(data, expr) {

    # pre-condition
    .check_lst(data, min_len = 1, msg = "invalid 'data' value")

    # evaluate
    value <- .check_error({
        eval(expr, envir = data)
    }, msg = "evaluation error")

    return(value)
}

#' @export
.apply.sits <- function(data, expr) {

    # pre-condition
    .check_that(inherits(data, "sits"),
                local_msg = "value is not a valid sits tibble",
                msg = "invalid 'data' value")

    .local_mutate <- function(x)
        tibble::as_tibble(c(x, .apply.list(x, expr)))

    value <- .sits_fast_apply(data, col = "time_series", .local_mutate)

    return(value)
}

.apply.raster_cube <- function(data, ..., output_dir = ".") {

    # pre-condition
    .check_that(inherits(data, "raster_cube"),
                local_msg = "value is not a valid sits tibble",
                msg = "invalid 'data' value")

    .check_that(.cube_is_regular(data),
                msg = "cube is not regular")

    # slide through tiles
    if (nrow(data) > 1)
        slider::slide_dfr(data, .apply.raster_cube)

    # pull file_info // data has one row
    file_info <- .cube_file_info(data)

    # traverse each scene
    # parallel processing
    .sits_parallel_map(sits_timeline(data), function(x) {

        scene <- dplyr::filter(file_info, date == x)

        # open raster
        values <- purrr::map(scene[["path"]], function(path) {
            .raster_read_rast(file = path, block = NULL)
        })

        names(values) <- scene$band

        result <- .apply(values, ...)

        # output file
        output_file <-
        .raster_write_rast()

        scene
    }, progress = TRUE)
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
