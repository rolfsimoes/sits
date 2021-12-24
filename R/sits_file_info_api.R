
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

