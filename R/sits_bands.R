#' @title Informs the names of the bands
#'
#' @name sits_bands
#'
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#' @author Rolf Simoes, \email{rolf.simoes@@inpe.br}
#'
#' @description  Finds/replaces the names of the bands of
#'               a set of time series or of a data cube
#'
#' @param x         Valid sits tibble (time series or a cube)
#'
#' @return A string vector with the names of the bands.
#'
#' @examples {
#' # Retrieve the set of samples for Mato Grosso (provided by EMBRAPA)
#' # show the bands
#' sits_bands(samples_mt_6bands)
#' }
#'
#' @export
#'
sits_bands <- function(x) {
    # get the meta-type (sits or cube)
    x <- .sits_config_data_meta_type(x)

    UseMethod("sits_bands", x)
}

#' @export
#'
sits_bands.sits <- function(x) {

    return(names(sits_time_series(x))[-1])
}

#' @export
#'
sits_bands.cube <- function(x) {

    return(x$bands[[1]])
}

#' @export
#'
sits_bands.patterns <- function(x) {

    return(sits_bands.sits(x))
}

#' @export
#'
sits_bands.sits_model <- function(x) {

    assertthat::assert_that(
        inherits(x, "function"),
        msg = "sits_bands: invalid sits model"
    )

    assertthat::assert_that(
        "data" %in% ls(environment(x)),
        msg = "sits_bands: no samples found in the sits model"
    )

    return(sits_bands.sits(environment(x)$data))
}
