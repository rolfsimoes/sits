.crs <- function(x) {

    UseMethod(".crs", x)
}

.crs.default <- function(x) {

    .check_chr_contains(names(x), contains = "crs",
                        msg = "object does not have 'crs' variable")

    x[["crs"]]
}

.xmin <- function(x) {

    UseMethod(".xmin", x)
}

.xmin.default <- function(x) {

    .check_chr_contains(names(x), contains = "xmin",
                        msg = "object does not have 'xmin' variable")

    x[["xmin"]]
}

`.xmin<-` <- function(x, value) {

    UseMethod(".xmin<-", x)
}

`.xmin<-.default` <- function(x, value) {

    .check_chr_contains(names(x), contains = "xmin",
                        msg = "object does not have 'xmin' variable")

    x[["xmin"]] <- value
    x
}

.xmax <- function(x) {

    UseMethod(".xmax", x)
}

.xmax.default <- function(x) {

    .check_chr_contains(names(x), contains = "xmax",
                        msg = "object does not have 'xmax' variable")

    x[["xmax"]]
}

`.xmax<-` <- function(x, value) {

    UseMethod(".xmax<-", x)
}

`.xmax<-.default` <- function(x, value) {

    .check_chr_contains(names(x), contains = "xmax",
                        msg = "object does not have 'xmax' variable")

    x[["xmax"]] <- value
    x
}

.ymin <- function(x) {

    UseMethod(".ymin", x)
}

.ymin.default <- function(x) {

    .check_chr_contains(names(x), contains = "ymin",
                        msg = "object does not have 'ymin' variable")

    x[["ymin"]]
}

`.ymin<-` <- function(x, value) {

    UseMethod(".ymin<-", x)
}

`.ymin<-.default` <- function(x, value) {

    .check_chr_contains(names(x), contains = "ymin",
                        msg = "object does not have 'ymin' variable")

    x[["ymin"]] <- value
    x
}

.ymax <- function(x) {

    UseMethod(".ymax", x)
}

.ymax.default <- function(x) {

    .check_chr_contains(names(x), contains = "ymax",
                        msg = "object does not have 'ymax' variable")

    x[["ymax"]]
}

`.ymax<-` <- function(x, value) {

    UseMethod(".ymax<-", x)
}

`.ymax<-.default` <- function(x, value) {

    .check_chr_contains(names(x), contains = "ymax",
                        msg = "object does not have 'ymax' variable")

    x[["ymax"]] <- value
    x
}

.nrows <- function(x) {

    UseMethod(".nrows", x)
}

.nrows.default <- function(x) {

    .check_chr_contains(names(x), contains = "nrows",
                        msg = "object does not have 'nrows' variable")

    x[["nrows"]]
}

`.nrows<-` <- function(x, value) {

    UseMethod(".nrows<-", x)
}

`.nrows<-.default` <- function(x, value) {

    .check_chr_contains(names(x), contains = "nrows",
                        msg = "object does not have 'nrows' variable")

    .check_num(value, min = 0, is_integer = TRUE,
               msg = "invalid nrows value")

    x[["nrows"]] <- value
    x
}

.ncols <- function(x) {

    UseMethod(".ncols", x)
}

.ncols.default <- function(x) {

    .check_chr_contains(names(x), contains = "ncols",
                        msg = "object does not have 'ncols' variable")

    x[["ncols"]]
}

`.ncols<-` <- function(x, value) {

    UseMethod(".ncols<-", x)
}

`.ncols<-.default` <- function(x, value) {

    .check_chr_contains(names(x), contains = "ncols",
                        msg = "object does not have 'ncols' variable")

    .check_num(value, min = 0, is_integer = TRUE,
               msg = "invalid ncols value")

    x[["ncols"]] <- value
    x
}

.xres <- function(x) {

    UseMethod(".xres", x)
}

.xres.default <- function(x) {

    .check_chr_contains(names(x), contains = "xres",
                        msg = "object does not have 'xres' variable")

    x[["xres"]]
}

`.xres<-` <- function(x, value) {

    UseMethod(".xres<-", x)
}


`.xres<-.default` <- function(x, value) {

    .check_chr_contains(names(x), contains = "xres",
                        msg = "object does not have 'xres' variable")

    .check_num(value, allow_zero = FALSE,
               msg = "invalid xres value")

    x[["xres"]] <- value
    x
}

.yres <- function(x) {

    UseMethod(".yres", x)
}

.yres.default <- function(x) {

    .check_chr_contains(names(x), contains = "yres",
                        msg = "object does not have 'yres' variable")

    x[["yres"]]
}

`.yres<-` <- function(x, value) {

    UseMethod(".yres<-", x)
}


`.yres<-.default` <- function(x, value) {

    .check_chr_contains(names(x), contains = "yres",
                        msg = "object does not have 'yres' variable")

    .check_num(value, allow_zero = FALSE,
               msg = "invalid yres value")

    x[["yres"]] <- value
    x
}

.first_row <- function(x) {

    UseMethod(".first_row", x)
}

.first_row.default <- function(x) {

    .check_chr_contains(names(x), contains = "first_row",
                        msg = "object does not have 'first_row' variable")

    x[["first_row"]]
}

`.first_row<-` <- function(x, value) {

    UseMethod(".first_row<-", x)
}

`.first_row<-.default` <- function(x, value) {

    .check_chr_contains(names(x), contains = "first_row",
                        msg = "object does not have 'first_row' variable")

    .check_num(value, min = 1, is_integer = TRUE,
               msg = "invalid first_row value")

    x[["first_row"]] <- value
    x
}

.first_col <- function(x) {

    UseMethod(".first_col", x)
}

.first_col.default <- function(x) {

    .check_chr_contains(names(x), contains = "first_col",
                        msg = "object does not have 'first_col' variable")

    x[["first_col"]]
}

`.first_col<-` <- function(x, value) {

    UseMethod(".first_col<-", x)
}

`.first_col<-.default` <- function(x, value) {

    .check_chr_contains(names(x), contains = "first_col",
                        msg = "object does not have 'first_col' variable")

    .check_num(value, min = 1, is_integer = TRUE,
               msg = "invalid first_col value")

    x[["first_col"]] <- value
    x
}

.crs1 <- function(x) .crs(x)[[1]]

.xmin1 <- function(x) .xmin(x)[[1]]

.xmax1 <- function(x) .xmax(x)[[1]]

.ymin1 <- function(x) .ymax(x)[[1]]

.ymax1 <- function(x) .ymax(x)[[1]]

.nrows1 <- function(x) .nrows(x)[[1]]

.ncols1 <- function(x) .ncols(x)[[1]]

.xres1 <- function(x) .xres(x)[[1]]

.yres1 <- function(x) .yres(x)[[1]]




.check_file_info <- function(file_info) {

    .check_chr_contains(names(file_info),
                        contains = c("fid", "date", "band", "xres",
                                     "yres", "xmin", "xmax", "ymin",
                                     "ymax", "nrows", "ncols", "path"),
                        msg = "invalid file_info parameter")

    .check_that(nrow(file_info) > 0,
                local_msg = "file_info is empty",
                msg = "invalid file_info")
}

.file_info <- function(tile, ...) {

    .cube_check(tile)

    file_info <- tile[["file_info"]][[1]] %>%
        dplyr::filter(...)

    .check_file_info(file_info)

    return(file_info)
}

.foreach_file_info <- function(file_info, group_by, order_by, fn, ...) {

    # pre-condition
    .check_file_info(file_info)

    .check_chr_within(group_by, within = names(file_info),
                      msg = "invalid group_by parameter")

    .check_chr_within(order_by, within = names(file_info),
                      msg = "invalid order_by parameter")

    out_fi <- dplyr::group_by(file_info, dplyr::matches(group_by)) %>%
        dplyr::arrange(dplyr::matches(order_by)) %>%
        dplyr::summarise()
    out_fi <- by(out_fi, out_fi)
        tidyr::nest(.grouped = -dplyr::matches(group_by)) %>%
        slider::slide_dfr(fn, ...)

    # post-condition
    .check_chr_contains(names(out_fi), contains = names(file_info),
                        can_repeat = FALSE, msg = "invalid file_info result")

    return(out_fi)
}

.file_info_ba <- function(file_info) {

    fi <- file_info %>%
        dplyr::select(dplyr::all_of(c("fid", "date", "band", "path"))) %>%
        tidyr::pivot_wider(names_from = "band", values_from = "path") %>%
        dplyr::arrange(dplyr::all_of(c("date", "fid")))

    return(fi)
}

.file_info_fid <- function(tile) {

    return(unique(.file_info(tile)[["fid"]]))
}

.file_info_band <- function(tile) {
    return(unique(.file_info(tile)[["band"]]))
}

.fi_fid <- function(fi) {

    return(fi[["fid"]])
}

