#' @title Smooth a tile
#' @name .smooth_tile
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#'
#' @param  tile            Subset of a data cube containing one tile
#' @param  band            Band to be processed
#' @param  block           Individual block that will be processed
#' @param  overlap         Overlap between tiles (if required)
#' @param  smooth_fn       Smoothing function
#' @param  output_dir      Directory where image will be save
#' @param  version         Version of result
#' @param  progress        Check progress bar?
#' @return                 Smoothed tile-band combination
.smooth_tile <- function(tile,
                         band,
                         block,
                         overlap,
                         exclusion_mask,
                         smooth_fn,
                         output_dir,
                         version,
                         progress) {
    # Output file
    out_file <- .file_derived_name(
        tile = tile, band = band, version = version,
        output_dir = output_dir
    )
    # Resume feature
    if (all(.raster_is_valid(out_file, output_dir = output_dir))) {
        return(.smooth_tile_from_file(
            tile = tile, band = band, out_file = out_file
        ))
    }
    # Create chunks as jobs
    tile_chunks <- .smooth_tile_chunks(
        tile = tile,
        overlap = overlap,
        block = block,
        exclusion_mask = exclusion_mask
    )
    chunks <- tile_chunks[["chunks"]]
    exclusion_mask <- tile_chunks[["exclusion_mask"]]
    # Process jobs in parallel
    block_files <- .jobs_map_parallel_chr(chunks, function(chunk) {
        # Job block
        block <- .block(chunk)
        # Block file name
        block_file <- .file_block_name(
            pattern = .file_pattern(out_file),
            block = block,
            output_dir = output_dir
        )
        # Resume processing in case of failure
        if (all(.raster_is_valid(block_file))) {
            return(block_file)
        }
        # Read and preprocess values
        values <- .tile_read_block(
            tile = tile, band = .tile_bands(tile), block = block
        )
        # Apply the probability function to values
        values <- smooth_fn(values = values, block = block)
        # Prepare probability to be saved
        band_conf <- .conf_derived_band(
            derived_class = "probs_cube", band = band
        )
        offset <- .offset(band_conf)
        if (.has(offset) && offset != 0.0) {
            values <- values - offset
        }
        scale <- .scale(band_conf)
        if (.has(scale) && scale != 1.0) {
            values <- values / scale
        }
        # Job crop block
        crop_block <- .block(.chunks_no_overlap(chunk))
        # Prepare and save results as raster
        .raster_write_block(
            files = block_file, block = block, bbox = .bbox(chunk),
            values = values, data_type = .data_type(band_conf),
            missing_value = .miss_value(band_conf),
            crop_block = crop_block
        )
        # Free memory
        gc()
        # Return block file
        block_file
    }, progress = progress)
    # Merge blocks into a new probs_cube tile
    .smooth_merge_tile(
        tile = tile,
        band = band,
        out_file = out_file,
        block_files = block_files,
        exclusion_mask = exclusion_mask,
        multicores = .jobs_multicores(),
        output_dir = output_dir,
        version = version,
        progress = progress
    )
}
#' @title Tile of a smoothed file that already exists
#' @name .smooth_tile_from_file
#' @keywords internal
#' @noRd
#' @param  tile      Tile of the probability cube
#' @param  band      Band of the smoothed tile
#' @param  out_file  Valid smoothed file of the tile
#' @return           Smoothed tile
.smooth_tile_from_file <- function(tile, band, out_file) {
    .check_recovery()
    .tile_derived_from_file(
        file = out_file,
        band = band,
        base_tile = tile,
        labels = .tile_labels(tile),
        derived_class = "probs_cube",
        update_bbox = FALSE
    )
}
#' @title Chunks of a tile to be smoothed
#' @name .smooth_tile_chunks
#' @keywords internal
#' @noRd
#' @param  tile            Tile of the probability cube
#' @param  overlap         Overlap between chunks
#' @param  block           Block size of the chunks
#' @param  exclusion_mask  Areas to exclude from smoothing, or NULL
#' @return  List with the chunks outside the mask and the mask cropped to them
.smooth_tile_chunks <- function(tile, overlap, block, exclusion_mask) {
    chunks <- .tile_chunks_create(tile = tile, overlap = overlap, block = block)
    # Calculate exclusion mask
    if (.has(exclusion_mask)) {
        # Remove chunks within the exclusion mask
        chunks <- .chunks_filter_mask(
            chunks = chunks,
            mask = exclusion_mask
        )
        exclusion_mask <- .chunks_crop_mask(
            chunks = chunks,
            mask = exclusion_mask
        )
    }
    list(chunks = chunks, exclusion_mask = exclusion_mask)
}
#' @title Merge the smoothed blocks of a tile
#' @name .smooth_merge_tile
#' @keywords internal
#' @noRd
#' @param  tile            Tile of the probability cube
#' @param  band            Band of the smoothed tile
#' @param  out_file        Smoothed file of the tile
#' @param  block_files     Files of the smoothed blocks
#' @param  exclusion_mask  Mask of \code{.smooth_tile_chunks()}, or NULL
#' @param  multicores      Number of threads of the merge
#' @param  output_dir      Directory where image will be save
#' @param  version         Version of result
#' @param  progress        Show progress bar of the crop?
#' @return                 Smoothed tile
.smooth_merge_tile <- function(tile,
                               band,
                               out_file,
                               block_files,
                               exclusion_mask,
                               multicores,
                               output_dir,
                               version,
                               progress) {
    # Check if there is a exclusion_mask
    # If exclusion_mask exists, blocks are merged to a different directory
    # than output_dir, which is used to save the final cropped version
    merge_out_file <- out_file
    if (.has(exclusion_mask)) {
        merge_out_file <- .file_derived_name(
            tile = tile,
            band = band,
            version = version,
            output_dir = file.path(output_dir, ".sits")
        )
    }
    # Merge blocks into a new probs_cube tile
    probs_tile <- .tile_derived_merge_blocks(
        file = merge_out_file,
        band = band,
        labels = .tile_labels(tile),
        base_tile = tile,
        block_files = block_files,
        derived_class = "probs_cube",
        multicores = multicores,
        update_bbox = FALSE
    )
    # Exclude masked areas
    if (.has(exclusion_mask)) {
        # crop
        probs_tile_crop <- .crop(
            cube = probs_tile,
            roi = exclusion_mask,
            output_dir = output_dir,
            multicores = 1L,
            overwrite = TRUE,
            progress = progress
        )

        # delete old files
        unlink(.fi_paths(.fi(probs_tile)))

        # assign new cropped value in the old probs variable
        probs_tile <- probs_tile_crop
    }
    # Return probs tile
    probs_tile
}

#---- Bayesian smoothing ----
#' @title Smooth probability cubes with spatial predictors
#' @noRd
#' @param  cube              Probability data cube.
#' @param  block             Individual block that will be processed
#' @param  window_size       Size of the neighborhood.
#' @param  neigh_fraction    Fraction of neighbors with high probabilities
#'                           to be used in Bayesian inference.
#' @param  smoothness        Estimated variance of logit of class probabilities
#'                           (Bayesian smoothing parameter). It can be either
#'                           a vector or a scalar.
#' @param  use_torch         Use torch to run the Bayesian smoother?
#' @param  use_dataloader    Use the flow of \code{.smooth_dataloader()}?
#' @param  readers           Number of dataloader readers of that flow.
#' @param  multicores        Number of cores to run the smoothing function
#' @param  memsize           Maximum overall memory (in GB) to run the
#'                           smoothing.
#' @param  output_dir        Output directory for image files
#' @param  version           Version of resulting image
#'                           (in the case of multiple tests)
#' @param  progress          Check progress bar?
#' @return                   Smoothed data cube
#'
.smooth <- function(cube,
                    block,
                    window_size,
                    neigh_fraction,
                    smoothness,
                    use_torch,
                    use_dataloader,
                    readers,
                    exclusion_mask,
                    multicores,
                    memsize,
                    output_dir,
                    version,
                    progress) {
    # Overlapping pixels
    overlap <- ceiling(window_size / 2L) - 1L
    # a kernel that does not compile gives the C++ smoother below
    if (use_dataloader) {
        kernel <- tryCatch(
            .torch_smooth_kernel(
                window_size = window_size,
                neigh_fraction = neigh_fraction
            ),
            error = function(e) NULL
        )
        if (.has(kernel)) {
            return(.smooth_dataloader(
                cube = cube,
                block = block,
                overlap = overlap,
                exclusion_mask = exclusion_mask,
                kernel = kernel,
                smoothness = smoothness,
                readers = readers,
                output_dir = output_dir,
                version = version,
                progress = progress
            ))
        }
    }
    # Select smoothing implementation
    smooth_fn <- if (use_torch) {
        .smooth_fn_bayes_torch
    } else {
        .smooth_fn_bayes
    }
    # Create smooth function
    smooth_fn <- smooth_fn(
        window_size = window_size,
        neigh_fraction = neigh_fraction,
        smoothness = smoothness
    )
    # Smoothing
    # Process each tile sequentially
    .cube_foreach_tile(cube, function(tile) {
        # Smooth the data
        .smooth_tile(
            tile = tile,
            band = "bayes",
            block = block,
            overlap = overlap,
            exclusion_mask = exclusion_mask,
            smooth_fn = smooth_fn,
            output_dir = output_dir,
            version = version,
            progress = progress
        )
    })
}
#' @title Split the processes of the dataloader flow
#' @name .smooth_processes
#' @keywords internal
#' @noRd
#' @description \code{multicores} is a budget of processes: half of it reads
#'   blocks in the dataloader, the rest are workers of the sits cluster, at
#'   least one of each.
#' @param  multicores  Number of processes of the call.
#' @return Integer vector with \code{readers} and \code{workers}.
.smooth_processes <- function(multicores) {
    readers <- max(1L, as.integer(multicores) %/% 2L)
    c(readers = readers, workers = max(1L, as.integer(multicores) - readers))
}
#' @title Read the raw values of a block of a probability tile
#' @name .smooth_read_raw
#' @keywords internal
#' @noRd
#' @description Runs in a dataloader reader. The values are those of
#'   \code{.tile_read_block()} before its scale, which
#'   \code{.torch_smooth_block()} applies on the device.
#' @param  tile   Tile of the probability cube.
#' @param  block  Block to be read.
#' @return Matrix \code{[nrows * ncols, nbands]} of integer values, NA where
#'   missing.
.smooth_read_raw <- function(tile, block) {
    band <- .tile_bands(tile)
    values <- .tile_read_block(tile = tile, band = band, block = block)
    band_conf <- .tile_band_conf(tile = tile, band = band)
    # integers are exact in the 32-bit floats of the dataloader
    offset <- .offset(band_conf)
    if (.has(offset) && offset != 0.0) {
        values <- values - offset
    }
    scale <- .scale(band_conf)
    if (.has(scale) && scale != 1.0) {
        values <- values / scale
    }
    round(values)
}
#' @title Job of a sits cluster worker in the dataloader flow
#' @name .smooth_job
#' @keywords internal
#' @noRd
#' @description The one function of the cluster in
#'   \code{.smooth_dataloader()}; the role of the job selects the step.
#'   \code{write} writes one smoothed block, as a worker of
#'   \code{.smooth_tile()} does; \code{merge} merges the blocks of one tile
#'   with \code{.smooth_merge_tile()}.
#' @param  job  List with the role, the tile index, the file, the band and
#'   the arguments of the step.
#' @return List with the role, the tile index and, for a merge, the tile.
.smooth_job <- function(job) {
    role <- job[["role"]]
    .debug_log(
        event = paste0("start_smooth_", role),
        key = "file",
        value = basename(job[["file"]])
    )
    band_conf <- .conf_derived_band(
        derived_class = "probs_cube", band = job[["band"]]
    )
    tile <- switch(role,
        write = {
            chunk <- job[["chunk"]]
            .raster_write_block(
                files = job[["file"]],
                block = .block(chunk),
                bbox = .bbox(chunk),
                values = job[["values"]],
                data_type = .data_type(band_conf),
                missing_value = .miss_value(band_conf),
                crop_block = .block(.chunks_no_overlap(chunk))
            )
            NULL
        },
        merge = .smooth_merge_tile(
            tile = job[["tile"]],
            band = job[["band"]],
            out_file = job[["file"]],
            block_files = job[["block_files"]],
            exclusion_mask = job[["exclusion_mask"]],
            multicores = job[["multicores"]],
            output_dir = job[["output_dir"]],
            version = job[["version"]],
            progress = FALSE
        )
    )
    .debug_log(
        event = paste0("end_smooth_", role),
        key = "file",
        value = basename(job[["file"]])
    )
    list(role = role, tile_index = job[["tile_index"]], tile = tile)
}
#' @title Smooth the tiles of a cube with a dataloader and the device
#' @name .smooth_dataloader
#' @keywords internal
#' @noRd
#' @description
#' One dataloader reads the blocks of every tile in its reader processes.
#' This process holds the device, smooths each block there with
#' \code{.torch_smooth_block()} and sends it, without waiting, as a
#' \code{write} job to the sits cluster. When the last block of a tile is
#' written, the \code{merge} job of the tile is sent. Without a cluster the
#' jobs run in this process. Valid tiles and block files are reused, as in
#' \code{.smooth_tile()}.
#'
#' @param  cube            Probability data cube.
#' @param  block           Block size of the chunks.
#' @param  overlap         Overlap between chunks.
#' @param  exclusion_mask  Areas to exclude from smoothing, or NULL.
#' @param  kernel          Result of \code{.torch_smooth_kernel()}.
#' @param  smoothness      Smoothness vector (one value per class).
#' @param  readers         Number of dataloader readers.
#' @param  output_dir      Output directory for image files.
#' @param  version         Version of resulting image.
#' @param  progress        Show progress bar of the blocks?
#' @param  n_retries       Number of times a job of a failed worker is sent.
#' @return Smoothed data cube.
.smooth_dataloader <- function(cube,
                               block,
                               overlap,
                               exclusion_mask,
                               kernel,
                               smoothness,
                               readers,
                               output_dir,
                               version,
                               progress,
                               n_retries = 3L) {
    band <- "bayes"
    band_conf <- .conf_derived_band(derived_class = "probs_cube", band = band)
    device <- .torch_smooth_device()
    tiles <- slider::slide(cube, identity)
    plans <- purrr::imap(tiles, function(tile, i) {
        out_file <- .file_derived_name(
            tile = tile, band = band, version = version,
            output_dir = output_dir
        )
        # Resume feature
        if (all(.raster_is_valid(out_file, output_dir = output_dir))) {
            return(list(tile = .smooth_tile_from_file(
                tile = tile, band = band, out_file = out_file
            )))
        }
        tile_chunks <- .smooth_tile_chunks(
            tile = tile,
            overlap = overlap,
            block = block,
            exclusion_mask = exclusion_mask
        )
        chunks <- tile_chunks[["chunks"]]
        block_files <- slider::slide_chr(chunks, function(chunk) {
            .file_block_name(
                pattern = .file_pattern(out_file),
                block = .block(chunk),
                output_dir = output_dir
            )
        })
        # Resume processing in case of failure
        written <- purrr::map_lgl(block_files, function(file) {
            all(.raster_is_valid(file))
        })
        chunks[["tile_index"]] <- i
        chunks[["block_file"]] <- block_files
        list(
            out_file = out_file,
            block_files = block_files,
            exclusion_mask = tile_chunks[["exclusion_mask"]],
            in_conf = .tile_band_conf(tile = tile, band = .tile_bands(tile)),
            todo = chunks[!written, ]
        )
    })
    smoothed <- purrr::map(plans, "tile")
    todo <- dplyr::bind_rows(purrr::map(plans, "todo"))
    # blocks of each tile not yet written
    left <- purrr::map_int(plans, function(plan) {
        .default(nrow(plan[["todo"]]), 0L)
    })
    pb <- NULL
    if (progress && nrow(todo) > 0L) {
        pb <- utils::txtProgressBar(min = 0L, max = nrow(todo), style = 3L)
        on.exit(close(pb), add = TRUE)
    }
    nodes <- .jobs_multicores()
    free <- seq_len(nodes)
    # jobs sent to the cluster and not yet answered, by tag
    jobs <- list()
    n_sent <- 0L
    merge_job <- function(i) {
        plan <- plans[[i]]
        list(
            role = "merge",
            tile_index = i,
            file = plan[["out_file"]],
            band = band,
            tile = tiles[[i]],
            block_files = plan[["block_files"]],
            exclusion_mask = plan[["exclusion_mask"]],
            multicores = max(1L, nodes),
            output_dir = output_dir,
            version = version
        )
    }
    done <- function(result) {
        i <- result[["tile_index"]]
        if (result[["role"]] == "merge") {
            smoothed[[i]] <<- result[["tile"]]
            return(invisible(NULL))
        }
        left[[i]] <<- left[[i]] - 1L
        if (.has(pb)) {
            utils::setTxtProgressBar(pb, utils::getTxtProgressBar(pb) + 1L)
        }
        if (left[[i]] == 0L) {
            send(merge_job(i))
        }
        invisible(NULL)
    }
    # takes one result; with a timeout it returns FALSE when none arrived
    receive <- function(timeout = NULL) {
        d <- .parallel_poll(timeout = timeout)
        if (is.null(d)) {
            return(FALSE)
        }
        free <<- c(free, d[["node"]])
        tag <- as.character(d[["tag"]])
        job <- jobs[[tag]]
        jobs[[tag]] <<- NULL
        if (inherits(d[["value"]], "retry")) {
            job[["retries"]] <- .default(job[["retries"]], 0L) + 1L
            if (job[["retries"]] > n_retries) {
                stop(.conf("messages", ".parallel_map"), call. = FALSE)
            }
            send(job)
            return(TRUE)
        }
        .parallel_check_remote_errors(list(d[["value"]]))
        done(d[["value"]])
        TRUE
    }
    send <- function(job) {
        if (nodes == 0L) {
            return(done(.smooth_job(job)))
        }
        while (length(free) == 0L) {
            receive()
        }
        node <- free[[1L]]
        free <<- free[-1L]
        sent <- tryCatch(
            {
                .parallel_submit(
                    node = node, fn = .smooth_job, args = list(job),
                    tag = n_sent + 1L
                )
                NULL
            },
            error = function(e) e
        )
        if (.has(sent)) {
            # a broken send leaves the worker inside a message it never ends
            .parallel_reset_node(node)
            free <<- c(free, node)
            job[["retries"]] <- .default(job[["retries"]], 0L) + 1L
            if (job[["retries"]] > n_retries) {
                stop(sprintf(
                    .conf("messages", ".smooth_dataloader"),
                    n_retries + 1L, conditionMessage(sent)
                ), call. = FALSE)
            }
            return(send(job))
        }
        n_sent <<- n_sent + 1L
        jobs[[as.character(n_sent)]] <<- job
        invisible(NULL)
    }
    # on an error, a worker with a job would answer the next user of the
    # cluster, or never; it is replaced instead of waited for
    on.exit(
        purrr::walk(setdiff(seq_len(nodes), free), .parallel_reset_node),
        add = TRUE
    )
    # tiles whose blocks were all written before
    purrr::walk(
        which(purrr::map_lgl(smoothed, is.null) & left == 0L),
        function(i) send(merge_job(i))
    )
    if (nrow(todo) > 0L) {
        iter <- .dataloader_iter(
            chunks = todo,
            tiles = tiles,
            read_fn = .smooth_read_raw,
            readers = readers
        )
        repeat {
            item <- .dataloader_next(iter)
            if (is.null(item)) {
                break
            }
            chunk <- todo[item[["index"]], ]
            block <- .block(chunk)
            values <- .torch_smooth_block(
                input = item[["input"]],
                kernel = kernel,
                smoothness = smoothness,
                device = device,
                in_conf = plans[[chunk[["tile_index"]]]][["in_conf"]],
                out_conf = band_conf
            )
            .check_processed_values(values, .nrows(block) * .ncols(block))
            send(list(
                role = "write",
                tile_index = chunk[["tile_index"]],
                file = chunk[["block_file"]],
                band = band,
                chunk = chunk,
                values = values
            ))
            # answers that arrived free their workers for the next blocks
            while (nodes > 0L && receive(timeout = 0L)) NULL
        }
    }
    while (length(jobs) > 0L) {
        receive()
    }
    slider::slide_dfr(seq_along(smoothed), function(i) smoothed[[i]])
}
#' @title Define smoothing function
#' @noRd
#' @param  window_size       Size of the neighborhood.
#' @param  neigh_fraction    Fraction of neighbors with high probabilities
#'                           to be used in Bayesian inference.
#' @param  smoothness        Estimated variance of logit of class probabilities
#'                           (Bayesian smoothing parameter). It can be either
#'                           a vector or a scalar.
#' @return Function to be applied to smoothen data
.smooth_fn_bayes <- function(window_size,
                             neigh_fraction,
                             smoothness) {
    # Check window size
    .check_int_parameter(window_size, min = 5L, max = 21L, is_odd = TRUE)
    # Check neigh_fraction
    .check_num_parameter(neigh_fraction, min = 0.1, max = 1.0)

    prob_max <- .conf("smooth_prob_max")
    prob_min <- .conf("smooth_prob_min")
    # Define smooth function
    smooth_fn <- function(values, block) {
        # Check values length
        input_pixels <- nrow(values)
        # Compute logit
        # adjust values to avoid -Inf or +Inf in logits
        values[values == 1.0] <- prob_max
        values[values == 0.0] <- prob_min
        # tranform to logits
        values <- log(values / (rowSums(values) - values))
        # Process Bayesian
        values <- bayes_smoother_fraction(
            logits = values,
            nrows = .nrows(block),
            ncols = .ncols(block),
            window_size = window_size,
            smoothness = smoothness,
            neigh_fraction = neigh_fraction
        )
        # Compute inverse logit
        values <- exp(values) / (exp(values) + 1.0)
        # Are the results consistent with the data input?
        .check_processed_values(values, input_pixels)
        # Return values
        values
    }
    # Return a closure
    smooth_fn
}
#' @title Smooth a vector probability cube
#' @name .smooth_vector
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#'
#' @param  cube              Probability vector cube.
#' @param  neigh_fraction    Fraction of neighbors with high probabilities
#'                           to be used in Bayesian inference.
#' @param  smoothness        Estimated variance of logit of class probabilities
#'                           (Bayesian smoothing parameter).
#' @param  output_dir        Output directory for image files.
#' @param  version           Version of resulting image.
#' @param  progress          Show progress bar?
#' @param ...                Additional parameters to exactextractr.
#' @return                   Smoothed data cube.
.smooth_vector <- function(cube,
                           neigh_fraction,
                           smoothness,
                           output_dir,
                           version,
                           progress, ...) {
    # Process each tile sequentially
    smooth_cube <- .cube_foreach_tile(cube, function(tile) {
        .smooth_vector_tile(
            tile = tile,
            neigh_fraction = neigh_fraction,
            smoothness = smoothness,
            output_dir = output_dir,
            version = version,
            progress = progress, ...
        )
    })
    # Set probs_vector_cube class chain
    vector_classes <- c(
        .conf_vector_s3class("probs_vector_cube"),
        class(smooth_cube)
    )
    smooth_cube <- .cube_set_class(smooth_cube, vector_classes)
    return(smooth_cube)
}
#' @title Smooth a vector probability tile
#' @name .smooth_vector_tile
#' @keywords internal
#' @noRd
#' @param tile              Tile of a data cube.
#' @param neigh_fraction    Fraction of neighbors with high probabilities
#'                          to be used in Bayesian inference.
#' @param smoothness        Smoothness vector (one value per class).
#' @param output_dir        Output directory for image files.
#' @param version           Version of resulting image.
#' @param progress          Show progress bar?
#' @param ...               Additional parameters to exactextractr.
#' @return                  Smoothed tile.
.smooth_vector_tile <- function(tile,
                                neigh_fraction,
                                smoothness,
                                output_dir,
                                version,
                                progress, ...) {
    band <- "bayes"
    # Output file
    out_file <- .file_derived_name(
        tile = tile,
        band = band,
        version = version,
        output_dir = output_dir
    )
    # Resume feature
    if (all(.raster_is_valid(out_file, output_dir = output_dir))) {
        .check_recovery()
        smooth_tile <- .tile_derived_from_file(
            file = out_file,
            band = band,
            base_tile = tile,
            labels = .tile_labels(tile),
            derived_class = "probs_cube",
            update_bbox = FALSE
        )
        smooth_tile[["vector_info"]] <- tile[["vector_info"]]
        vector_classes <- c(
            .conf_vector_s3class("probs_vector_cube"),
            class(smooth_tile)
        )
        return(.cube_set_class(smooth_tile, vector_classes))
    }

    # Get tile band
    tile_band <- .tile_bands(tile)
    # Get band configuration
    band_conf <- .tile_band_conf(tile = tile, band = tile_band)
    # Get labels
    labels <- .tile_labels(tile)
    # Read the segments
    segments <- .segments_read_vec(tile)
    # Create ID for segments
    segments[["ID"]] <- seq_len(nrow(segments))
    # Open probability raster (all bands)
    probs_path <- .tile_path(tile)
    probs_rast <- .raster_open_rast(probs_path)

    # Extract values
    values <- .segments_extract_features(
        rast = probs_rast,
        segments = segments,
        fun = NULL,
        seg_id_col = "ID",
        include_cell = TRUE,
        ...
    )
    # Get ID and cell from segments
    values_id <- values[["ID"]]
    values_cell <- values[["cell"]]
    # Probability columns (all bands in the probs raster)
    prob_cols <- setdiff(colnames(values), c("ID", "cell"))
    # Convert values to matrix
    values <- as.matrix(values[, prob_cols])
    # Apply scale and offset to extracted probability values
    values <- .tile_scale(tile = tile, band = tile_band, values = values)
    # Avoid zero or one values to prevent -Inf/Inf/NaN in logit
    values[values <= 0.00001] <- 0.00001
    values[values >= 0.99999] <- 0.99999

    # Calculate denom
    row_sums <- rowSums(values)
    denom <- row_sums - values
    denom[denom <= 0.00001] <- 0.00001
    # Calculate logits
    values <- log(values / denom)

    # Call C++ function to perform segment-based Bayesian smoothing
    values <- segment_bayes(
        logits = values,
        ids = as.integer(values_id),
        n_segments = nrow(segments),
        neigh_fraction = neigh_fraction,
        smoothness = smoothness
    )

    # Convert logits back to probabilities (inverse logit)
    values <- exp(values) / (exp(values) + 1.0)

    # Unscale values
    values <- .tile_unscale(
        tile = tile,
        band = band,
        values = values,
        derived_class = "probs_cube"
    )

    # Create empty raster with same structure
    smooth_rast <- .raster_rast(probs_rast, nlayers = length(prob_cols))
    .raster_init_value(smooth_rast, NA_real_)
    names(smooth_rast) <- prob_cols

    # Assign smoothed values to corresponding cells
    smooth_rast[values_cell] <- values

    # Set missing value
    smooth_rast <- .raster_set_na(smooth_rast, .miss_value(band_conf))
    # Write raster
    .raster_write_rast(
        rast = smooth_rast,
        file = out_file,
        data_type = .data_type(band_conf),
        overwrite = TRUE,
        missing_value = .miss_value(band_conf)
    )
    # Create output tile from file
    smooth_tile <- .tile_derived_from_file(
        file = out_file,
        band = band,
        base_tile = tile,
        derived_class = "probs_cube",
        labels = labels,
        update_bbox = FALSE
    )
    # Preserve vector_info from the input probs tile
    smooth_tile[["vector_info"]] <- tile[["vector_info"]]
    # Set probs_vector_cube + probs_cube chain
    vector_classes <- c(
        .conf_vector_s3class("probs_vector_cube"),
        class(smooth_tile)
    )
    .cube_set_class(smooth_tile, vector_classes)
}
