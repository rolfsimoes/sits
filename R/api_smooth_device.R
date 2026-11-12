#' @title Job of a worker in the smoothing on the device
#' @name .smooth_device_job
#' @keywords internal
#' @noRd
#' @description The worker side of \code{.smooth_device()}: the role of the
#'   job selects the step. \code{read} reads one block and returns its
#'   values; \code{write} scales and writes one smoothed block, then frees
#'   memory; \code{merge} merges the blocks of one tile, then crops the
#'   exclusion mask.
#' @param  job  List with the role, the tile index, the block index and the
#'   arguments of the step.
#' @return List with the role, the tile index, the block index and the value
#'   of the step.
.smooth_device_job <- function(job) {
    role <- job[["role"]]
    value <- NULL
    if (role == "read") {
        tile <- job[["tile"]]
        value <- .tile_read_block(
            tile = tile, band = .tile_bands(tile),
            block = .block(job[["chunk"]])
        )
    } else if (role == "write") {
        chunk <- job[["chunk"]]
        block_file <- job[["block_file"]]
        band_conf <- .conf_derived_band(
            derived_class = "probs_cube", band = job[["band"]]
        )
        values <- job[["values"]]
        offset <- .offset(band_conf)
        if (.has(offset) && offset != 0.0) {
            values <- values - offset
        }
        scale <- .scale(band_conf)
        if (.has(scale) && scale != 1.0) {
            values <- values / scale
        }
        .debug_log(
            event = "start_block_write", key = "file", value = block_file
        )
        .raster_write_block(
            files = block_file, block = .block(chunk), bbox = .bbox(chunk),
            values = values, data_type = .data_type(band_conf),
            missing_value = .miss_value(band_conf),
            crop_block = .block(.chunks_no_overlap(chunk))
        )
        .debug_log(event = "end_block_write", key = "file", value = block_file)
        rm(values)
        .debug_log(event = "start_block_gc", key = "file", value = block_file)
        gc()
        .debug_log(event = "end_block_gc", key = "file", value = block_file)
    } else {
        value <- .smooth_device_merge(job)
    }
    list(
        role = role, tile_index = job[["tile_index"]],
        block_index = job[["block_index"]], value = value
    )
}
#' @title Merge the smoothed blocks of a tile in a worker
#' @name .smooth_device_merge
#' @keywords internal
#' @noRd
#' @param  job  List with the tile, the band, the tile file, the block files,
#'   the exclusion mask, the merge threads, the output dir and the version.
#' @return Smoothed tile.
.smooth_device_merge <- function(job) {
    tile <- job[["tile"]]
    band <- job[["band"]]
    out_file <- job[["out_file"]]
    exclusion_mask <- job[["exclusion_mask"]]
    output_dir <- job[["output_dir"]]
    # If exclusion_mask exists, blocks are merged to a different directory
    # than output_dir, which is used to save the final cropped version
    merge_out_file <- out_file
    if (.has(exclusion_mask)) {
        merge_out_file <- .file_derived_name(
            tile = tile, band = band, version = job[["version"]],
            output_dir = file.path(output_dir, ".sits")
        )
    }
    .debug_log(event = "start_tile_merge", key = "file", value = out_file)
    probs_tile <- .tile_derived_merge_blocks(
        file = merge_out_file,
        band = band,
        labels = .tile_labels(tile),
        base_tile = tile,
        block_files = job[["block_files"]],
        derived_class = "probs_cube",
        multicores = job[["multicores"]],
        update_bbox = FALSE
    )
    .debug_log(event = "end_tile_merge", key = "file", value = out_file)
    if (.has(exclusion_mask)) {
        # .crop() starts no cluster here, but resets the debug flag
        log_flag <- .debug()
        log_dir <- sits_env[["output_dir"]]
        probs_tile_crop <- .crop(
            cube = probs_tile, roi = exclusion_mask, output_dir = output_dir,
            multicores = 1L, overwrite = TRUE, progress = FALSE
        )
        .debug(flag = log_flag, output_dir = log_dir)
        unlink(.fi_paths(.fi(probs_tile)))
        probs_tile <- probs_tile_crop
    }
    probs_tile
}
#' @title Smooth the tiles of a cube with the kernel on the device
#' @name .smooth_device
#' @keywords internal
#' @noRd
#' @description This process holds the device and runs the kernel of every
#'   block; the workers of the sits cluster read the blocks, write the
#'   smoothed blocks and merge the tiles. One queue of jobs covers every
#'   tile: a few reads are kept ahead of the kernel, so the device does not
#'   wait for a block, each smoothed block goes to a worker as a write, and
#'   the merge of a tile is sent as soon as its last block is written, so it
#'   runs while the next tile is smoothed. Without a cluster every job runs
#'   in this process, in the order of the current flow. Valid tiles and
#'   block files are reused.
#' @param  cube            Probability data cube.
#' @param  block           Block size of the chunks.
#' @param  overlap         Overlap between chunks.
#' @param  exclusion_mask  Areas to exclude from smoothing, or NULL.
#' @param  smooth_fn       Smoothing function of the kernel.
#' @param  output_dir      Output directory for image files.
#' @param  version         Version of resulting image.
#' @param  progress        Show progress bar of the blocks?
#' @return Smoothed data cube.
.smooth_device <- function(cube, block, overlap, exclusion_mask, smooth_fn,
                           output_dir, version, progress) {
    band <- "bayes"
    tiles <- slider::slide(cube, identity)
    n_tiles <- length(tiles)
    nodes <- .jobs_multicores()
    # the state of the queue, changed by the functions below
    st <- new.env(parent = emptyenv())
    st[["result"]] <- vector("list", n_tiles)
    st[["plans"]] <- vector("list", n_tiles)
    st[["reads"]] <- list()
    st[["writes"]] <- list()
    st[["merges"]] <- list()
    # blocks read and waiting for the kernel
    st[["ready"]] <- list()
    for (i in seq_len(n_tiles)) {
        tile <- tiles[[i]]
        out_file <- .file_derived_name(
            tile = tile, band = band, version = version,
            output_dir = output_dir
        )
        # Resume feature
        if (all(.raster_is_valid(out_file, output_dir = output_dir))) {
            .debug_log(event = "start_tile", key = "file", value = out_file)
            .check_recovery()
            st[["result"]][[i]] <- .tile_derived_from_file(
                file = out_file, band = band, base_tile = tile,
                labels = .tile_labels(tile), derived_class = "probs_cube",
                update_bbox = FALSE
            )
            .debug_log(event = "end_tile", key = "file", value = out_file)
            next
        }
        chunks <- .tile_chunks_create(
            tile = tile, overlap = overlap, block = block
        )
        tile_mask <- exclusion_mask
        if (.has(tile_mask)) {
            chunks <- .chunks_filter_mask(chunks = chunks, mask = tile_mask)
            tile_mask <- .chunks_crop_mask(chunks = chunks, mask = tile_mask)
        }
        chunk_list <- slider::slide(chunks, identity)
        block_files <- purrr::map_chr(chunk_list, function(chunk) {
            .file_block_name(
                pattern = .file_pattern(out_file), block = .block(chunk),
                output_dir = output_dir
            )
        })
        st[["plans"]][[i]] <- list(
            out_file = out_file, block_files = block_files,
            chunks = chunk_list, exclusion_mask = tile_mask,
            left = length(chunk_list), started = FALSE
        )
        tile_reads <- purrr::map(seq_along(chunk_list), function(k) {
            list(
                role = "read", tile_index = i, block_index = k,
                chunk = chunk_list[[k]], tile = tile
            )
        })
        st[["reads"]] <- c(st[["reads"]], tile_reads)
    }
    merge_job <- function(i) {
        plan <- st[["plans"]][[i]]
        list(
            role = "merge", tile_index = i, block_index = NA_integer_,
            tile = tiles[[i]], band = band, out_file = plan[["out_file"]],
            block_files = plan[["block_files"]],
            exclusion_mask = plan[["exclusion_mask"]],
            multicores = nodes, output_dir = output_dir, version = version
        )
    }
    no_blocks <- purrr::map_lgl(st[["plans"]], function(plan) {
        .has(plan) && plan[["left"]] == 0L
    })
    st[["merges"]] <- purrr::map(which(no_blocks), merge_job)
    n_blocks <- length(st[["reads"]])
    pb <- NULL
    if (progress && n_blocks > 0L) {
        pb <- utils::txtProgressBar(min = 0L, max = n_blocks, style = 3L)
        on.exit(close(pb), add = TRUE)
    }
    block_file_of <- function(i, k) st[["plans"]][[i]][["block_files"]][[k]]
    log_tile <- function(event, i) {
        .debug_log(
            event = event, key = "file",
            value = st[["plans"]][[i]][["out_file"]]
        )
    }
    pop <- function(name) {
        jobs <- st[[name]]
        st[[name]] <- jobs[-1L]
        jobs[[1L]]
    }
    # the tile and its blocks start with the read of its first block
    read_sent <- function(job) {
        i <- job[["tile_index"]]
        if (!st[["plans"]][[i]][["started"]]) {
            st[["plans"]][[i]][["started"]] <- TRUE
            log_tile("start_tile", i)
            log_tile("start_tile_blocks", i)
        }
        .debug_log(
            event = "start_block", key = "file",
            value = block_file_of(i, job[["block_index"]])
        )
    }
    block_done <- function(i, k) {
        .debug_log(
            event = "end_block", key = "file", value = block_file_of(i, k)
        )
        left <- st[["plans"]][[i]][["left"]] - 1L
        st[["plans"]][[i]][["left"]] <- left
        if (.has(pb)) {
            utils::setTxtProgressBar(pb, utils::getTxtProgressBar(pb) + 1L)
        }
        if (left == 0L) {
            log_tile("end_tile_blocks", i)
            st[["merges"]] <- c(st[["merges"]], list(merge_job(i)))
        }
    }
    # a block whose file is valid is reused without a read
    read_or_reuse <- function(job) {
        i <- job[["tile_index"]]
        k <- job[["block_index"]]
        read_sent(job)
        if (all(.raster_is_valid(block_file_of(i, k)))) {
            block_done(i, k)
            return(FALSE)
        }
        TRUE
    }
    # the kernel of one block read, then its write job
    kernel <- function(res) {
        i <- res[["tile_index"]]
        k <- res[["block_index"]]
        chunk <- st[["plans"]][[i]][["chunks"]][[k]]
        block_file <- block_file_of(i, k)
        .debug_log(
            event = "start_block_kernel", key = "file", value = block_file
        )
        values <- smooth_fn(values = res[["value"]], block = .block(chunk))
        .debug_log(event = "end_block_kernel", key = "file", value = block_file)
        list(
            role = "write", tile_index = i, block_index = k, chunk = chunk,
            block_file = block_file, band = band, values = values
        )
    }
    done <- function(res) {
        i <- res[["tile_index"]]
        if (res[["role"]] == "read") {
            st[["ready"]] <- c(st[["ready"]], list(res))
        } else if (res[["role"]] == "write") {
            block_done(i, res[["block_index"]])
        } else {
            st[["result"]][[i]] <- res[["value"]]
            log_tile("end_tile", i)
        }
        invisible(NULL)
    }
    if (nodes == 0L) {
        # no cluster: every step in this process, block by block
        while (length(st[["reads"]]) || length(st[["merges"]])) {
            if (length(st[["merges"]])) {
                done(.smooth_device_job(pop("merges")))
                next
            }
            job <- pop("reads")
            if (!read_or_reuse(job)) {
                next
            }
            done(.smooth_device_job(job))
            write_job <- kernel(pop("ready"))
            done(.smooth_device_job(write_job))
        }
    } else {
        st[["free"]] <- seq_len(nodes)
        st[["running"]] <- 0L
        st[["reads_out"]] <- 0L
        # reads kept ahead of the kernel
        ahead <- max(2L, nodes %/% 4L)
        send <- function(job) {
            node <- st[["free"]][[1L]]
            .parallel_submit(
                node = node, fn = .smooth_device_job,
                args = list(job), tag = node
            )
            st[["free"]] <- st[["free"]][-1L]
            st[["running"]] <- st[["running"]] + 1L
        }
        take <- function(timeout = NULL) {
            d <- .parallel_poll(timeout = timeout)
            if (is.null(d)) {
                return(FALSE)
            }
            st[["running"]] <- st[["running"]] - 1L
            st[["free"]] <- c(st[["free"]], d[["node"]])
            # a worker that died leaves no result to resend
            if (inherits(d[["value"]], "retry")) {
                stop(.conf("messages", ".parallel_map"), call. = FALSE)
            }
            .parallel_check_remote_errors(list(d[["value"]]))
            if (d[["value"]][["role"]] == "read") {
                st[["reads_out"]] <- st[["reads_out"]] - 1L
            }
            done(d[["value"]])
            TRUE
        }
        repeat {
            # merges and writes first, then reads up to the ones kept ahead
            while (length(st[["free"]])) {
                if (length(st[["merges"]])) {
                    send(pop("merges"))
                } else if (length(st[["writes"]])) {
                    send(pop("writes"))
                } else if (length(st[["reads"]]) &&
                           st[["reads_out"]] + length(st[["ready"]]) < ahead) {
                    job <- pop("reads")
                    if (read_or_reuse(job)) {
                        send(job)
                        st[["reads_out"]] <- st[["reads_out"]] + 1L
                    }
                } else {
                    break
                }
            }
            if (length(st[["ready"]])) {
                st[["writes"]] <- c(st[["writes"]], list(kernel(pop("ready"))))
                # take the results that arrived during the kernel
                while (st[["running"]] > 0L && take(timeout = 0)) NULL
                next
            }
            if (st[["running"]] == 0L) {
                if (!length(st[["reads"]]) && !length(st[["writes"]]) &&
                    !length(st[["merges"]])) {
                    break
                }
                next
            }
            take()
        }
    }
    # the tiles in the order of the cube, bound as .cube_foreach_tile() does
    slider::slide_dfr(seq_len(n_tiles), function(i) st[["result"]][[i]])
}
