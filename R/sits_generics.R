#' @title generic accessors
#'
#' @name sits_generic
#'
#' @description
#' Internal accessors to be used instead of standard ones (`[[` or \code{$}).
#'
#' This approach has advantage to check the data validity before return
#' value or before set value.
#'
#' @param x  any R object from which value will be set or extracted.
#' @param value any value to set in the \code{x} object.
#'
#' @keywords internal
NULL

#' @rdname sits_generic
#' @description \code{.crs}: access \code{x[["crs"]]} value
.crs <- function(x) {

    UseMethod(".crs", x)
}

#' @export
.crs.default <- function(x) {

    .check_chr_contains(names(x), contains = "crs",
                        msg = "object does not have 'crs' variable")

    x[["crs"]]
}

#' @rdname sits_generic
#' @description \code{.crs}: set \code{x[["crs"]] <- value}
`.crs<-` <- function(x, value) {

    UseMethod(".crs<-", x)
}

#' @export
`.crs<-.default` <- function(x, value) {

    .check_chr(value,
               msg = "invalid crs value")

    x[["crs"]] <- value
    x
}

#' @rdname sits_generic
#' @description \code{.xmin}: access \code{x[["xmin"]]} value
.xmin <- function(x) {

    UseMethod(".xmin", x)
}

#' @export
.xmin.default <- function(x) {

    .check_chr_contains(names(x), contains = "xmin",
                        msg = "object does not have 'xmin' variable")

    x[["xmin"]]
}

#' @rdname sits_generic
#' @description \code{.xmin}: access \code{x[["xmin"]]} value
`.xmin<-` <- function(x, value) {

    UseMethod(".xmin<-", x)
}

#' @export
`.xmin<-.default` <- function(x, value) {

    .check_num(value, msg = "invalid xmin value")

    x[["xmin"]] <- value
    x
}

#' @rdname sits_generic
#' @description \code{.xmax}: access \code{x[["xmax"]]} value
.xmax <- function(x) {

    UseMethod(".xmax", x)
}

#' @export
.xmax.default <- function(x) {

    .check_chr_contains(names(x), contains = "xmax",
                        msg = "object does not have 'xmax' variable")

    x[["xmax"]]
}

#' @rdname sits_generic
#' @description \code{.xmax<-()}: set \code{x[["xmax"]] <- value}
`.xmax<-` <- function(x, value) {

    UseMethod(".xmax<-", x)
}

#' @export
`.xmax<-.default` <- function(x, value) {

    .check_num(value, msg = "invalid xmax value")

    x[["xmax"]] <- value
    x
}

#' @rdname sits_generic
#' @description \code{.ymin}: access \code{x[["ymin"]]} value
.ymin <- function(x) {

    UseMethod(".ymin", x)
}

#' @export
.ymin.default <- function(x) {

    .check_chr_contains(names(x), contains = "ymin",
                        msg = "object does not have 'ymin' variable")

    x[["ymin"]]
}

#' @rdname sits_generic
#' @description \code{.ymin<-()}: set \code{x[["ymin"]] <- value}
`.ymin<-` <- function(x, value) {

    UseMethod(".ymin<-", x)
}

#' @export
`.ymin<-.default` <- function(x, value) {

    .check_num(value, msg = "invalid ymin value")

    x[["ymin"]] <- value
    x
}

#' @rdname sits_generic
#' @description \code{.ymax}: access \code{x[["ymax"]]} value
.ymax <- function(x) {

    UseMethod(".ymax", x)
}

#' @export
.ymax.default <- function(x) {

    .check_chr_contains(names(x), contains = "ymax",
                        msg = "object does not have 'ymax' variable")

    x[["ymax"]]
}

#' @rdname sits_generic
#' @description \code{.ymax<-()}: set \code{x[["ymax"]] <- value}
`.ymax<-` <- function(x, value) {

    UseMethod(".ymax<-", x)
}

#' @export
`.ymax<-.default` <- function(x, value) {

    .check_num(value, msg = "invalid ymax value")

    x[["ymax"]] <- value
    x
}

#' @rdname sits_generic
#' @description \code{.nrows}: access \code{x[["nrows"]]} value
.nrows <- function(x) {

    UseMethod(".nrows", x)
}

#' @export
.nrows.default <- function(x) {

    .check_chr_contains(names(x), contains = "nrows",
                        msg = "object does not have 'nrows' variable")

    x[["nrows"]]
}

#' @rdname sits_generic
#' @description \code{.nrows<-()}: set \code{x[["nrows"]] <- value}
`.nrows<-` <- function(x, value) {

    UseMethod(".nrows<-", x)
}

#' @export
`.nrows<-.default` <- function(x, value) {

    .check_num(value, min = 0, is_integer = TRUE,
               msg = "invalid nrows value")

    x[["nrows"]] <- value
    x
}

#' @rdname sits_generic
#' @description \code{.ncols}: access \code{x[["ncols"]]} value
.ncols <- function(x) {

    UseMethod(".ncols", x)
}

#' @export
.ncols.default <- function(x) {

    .check_chr_contains(names(x), contains = "ncols",
                        msg = "object does not have 'ncols' variable")

    x[["ncols"]]
}

#' @rdname sits_generic
#' @description \code{.ncols<-()}: set \code{x[["ncols"]] <- value}
`.ncols<-` <- function(x, value) {

    UseMethod(".ncols<-", x)
}

#' @export
`.ncols<-.default` <- function(x, value) {

    .check_num(value, min = 0, is_integer = TRUE,
               msg = "invalid ncols value")

    x[["ncols"]] <- value
    x
}

#' @rdname sits_generic
#' @description \code{.xres}: access \code{x[["xres"]]} value
.xres <- function(x) {

    UseMethod(".xres", x)
}

#' @export
.xres.default <- function(x) {

    .check_chr_contains(names(x), contains = "xres",
                        msg = "object does not have 'xres' variable")

    x[["xres"]]
}

#' @rdname sits_generic
#' @description \code{.xres<-()}: set \code{x[["xres"]] <- value}
`.xres<-` <- function(x, value) {

    UseMethod(".xres<-", x)
}


#' @export
`.xres<-.default` <- function(x, value) {

    .check_num(value, min = 0, allow_zero = FALSE,
               msg = "invalid xres value")

    x[["xres"]] <- value
    x
}

#' @rdname sits_generic
#' @description \code{.yres}: access \code{x[["yres"]]} value
.yres <- function(x) {

    UseMethod(".yres", x)
}

#' @export
.yres.default <- function(x) {

    .check_chr_contains(names(x), contains = "yres",
                        msg = "object does not have 'yres' variable")

    x[["yres"]]
}

#' @rdname sits_generic
#' @description \code{.yres<-()}: set \code{x[["yres"]] <- value}
`.yres<-` <- function(x, value) {

    UseMethod(".yres<-", x)
}

#' @export
`.yres<-.default` <- function(x, value) {

    .check_num(value, min = 0, allow_zero = FALSE,
               msg = "invalid yres value")

    x[["yres"]] <- value
    x
}

#' @rdname sits_generic
#' @description \code{.lon}: access \code{x[["lon"]]} value
.lon <- function(x) {

    UseMethod(".lon", x)
}

#' @export
.lon.default <- function(x) {

    .check_chr_contains(names(x), contains = "longitude",
                        msg = "object does not have 'longitude' variable")

    x[["longitude"]]
}

#' @rdname sits_generic
#' @description \code{.lon<-()}: set \code{x[["lon"]] <- value}
`.lon<-` <- function(x, value) {

    UseMethod(".lon<-", x)
}

#' @export
`.lon<-.default` <- function(x, value) {

    .check_num(value, min = -180, max = 180,
               msg = "invalid longitude value")

    x[["longitude"]] <- value
    x
}

#' @rdname sits_generic
#' @description \code{.lat}: access \code{x[["lat"]]} value
.lat <- function(x) {

    UseMethod(".lat", x)
}

#' @export
.lat.default <- function(x) {

    .check_chr_contains(names(x), contains = "latitude",
                        msg = "object does not have 'latitude' variable")

    x[["latitude"]]
}

#' @rdname sits_generic
#' @description \code{.lat<-()}: set \code{x[["lat"]] <- value}
`.lat<-` <- function(x, value) {

    UseMethod(".lat<-", x)
}

#' @export
`.lat<-.default` <- function(x, value) {

    .check_num(value, min = -90, max = 90, msg = "invalid latitude value")

    x[["latitude"]] <- value
    x
}

#' @rdname sits_generic
#' @description \code{.start_date}: access \code{x[["start_date"]]} value
.start_date <- function(x) {

    UseMethod(".start_date", x)
}

#' @export
.start_date.default <- function(x) {

    .check_chr_contains(names(x), contains = "start_date",
                        msg = "object does not have 'start_date' variable")

    date <- suppressWarnings(lubridate::as_date(x[["start_date"]]))

    .check_na(date, msg = "invalid date value")

    date
}

#' @rdname sits_generic
#' @description \code{.start_date<-()}: set \code{x[["start_date"]] <- value}
`.start_date<-` <- function(x, value) {

    UseMethod(".start_date<-", x)
}

#' @export
`.start_date<-.default` <- function(x, value) {

    vaue <- suppressWarnings(lubridate::as_date(value))

    .check_na(value, msg = "invalid date value")

    x[["start_date"]] <- value
    x
}

#' @rdname sits_generic
#' @description \code{.end_date}: access \code{x[["end_date"]]} value
.end_date <- function(x) {

    UseMethod(".end_date", x)
}

#' @export
.end_date.default <- function(x) {

    .check_chr_contains(names(x), contains = "end_date",
                        msg = "object does not have 'end_date' variable")

    date <- suppressWarnings(lubridate::as_date(x[["end_date"]]))

    .check_na(date, msg = "invalid date value")

    date
}

#' @rdname sits_generic
#' @description \code{.end_date<-()}: set \code{x[["end_date"]] <- value}
`.end_date<-` <- function(x, value) {

    UseMethod(".end_date<-", x)
}

#' @export
`.end_date<-.default` <- function(x, value) {

    value <- suppressWarnings(lubridate::as_date(value))

    .check_na(value, msg = "invalid date value")

    x[["end_date"]] <- value
    x
}

#' @rdname sits_generic
#' @description \code{.first_row}: access \code{x[["first_row"]]} value
.first_row <- function(x) {

    UseMethod(".first_row", x)
}

#' @export
.first_row.default <- function(x) {

    .check_chr_contains(names(x), contains = "first_row",
                        msg = "object does not have 'first_row' variable")

    x[["first_row"]]
}

#' @rdname sits_generic
#' @description \code{.first_row<-()}: set \code{x[["first_row"]] <- value}
`.first_row<-` <- function(x, value) {

    UseMethod(".first_row<-", x)
}

#' @export
`.first_row<-.default` <- function(x, value) {

    .check_num(value, min = 1, is_integer = TRUE,
               msg = "invalid first_row value")

    x[["first_row"]] <- value
    x
}

#' @rdname sits_generic
#' @description \code{.first_col}: access \code{x[["first_col"]]} value
.first_col <- function(x) {

    UseMethod(".first_col", x)
}

#' @export
.first_col.default <- function(x) {

    .check_chr_contains(names(x), contains = "first_col",
                        msg = "object does not have 'first_col' variable")

    x[["first_col"]]
}

#' @rdname sits_generic
#' @description \code{.first_col<-()}: set \code{x[["first_col"]] <- value}
`.first_col<-` <- function(x, value) {

    UseMethod(".first_col<-", x)
}

#' @export
`.first_col<-.default` <- function(x, value) {

    .check_num(value, min = 1, is_integer = TRUE,
               msg = "invalid first_col value")

    x[["first_col"]] <- value
    x
}

.crs_1 <- function(x) {.crs(x)[[1]]}

.xmin_1 <- function(x) {.xmin(x)[[1]]}

.xmax_1 <- function(x) {.xmax(x)[[1]]}

.ymin_1 <- function(x) {.ymax(x)[[1]]}

.ymax_1 <- function(x) {.ymax(x)[[1]]}

.nrows_1 <- function(x) {.nrows(x)[[1]]}

.ncols_1 <- function(x) {.ncols(x)[[1]]}

.xres_1 <- function(x) {.xres(x)[[1]]}

.yres_1 <- function(x) {.yres(x)[[1]]}

.lon_1 <- function(x) {.lon(x)[[1]]}

.lat_1 <- function(x) {.lat(x)[[1]]}

.start_date_1 <- function(x) {.start_date(x)[[1]]}

.start_date_n <- function(x) {
    val <- .start_date(x)
    val[[length(val)]]
}
