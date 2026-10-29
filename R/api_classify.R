#' @title Classify a chunk of raster data  using multicores in CPU
#' @name .classify_tile_cpu
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#'
#' @description Classifies a block of data using multicores, breaking
#' the data into blocks and divides them between the available cores.
#' The size of the blocks is optimized to account for COG files and
#' for the balance of multicores and memory size.
#'
#' After all cores process their blocks,
#' it joins the result and then writes it
#' in the classified images for each corresponding year.
#'
#' @param  tile            Single tile of a data cube.
#' @param  out_band        Band to be produced.
#' @param  out_file        Name of the output probabibility file.
#' @param  bands           Bands to extract time series
#' @param  base_bands      Base bands to extract values
#' @param  ml_model        Model trained by \code{\link[sits]{sits_train}}.
#' @param  block           Optimized block to be read into memory.
#' @param  roi             Region of interest.
#' @param  exclusion_mask  Area to be excluded from classification
#' @param  impute_fn       Imputation function.
#' @param  output_dir      Output directory.
#' @param  version         Version of result.
#' @param  verbose         Print processing information?
#' @param  progress        Show progress bar?
#' @return List of the classified raster layers.
.classify_tile_cpu <- function(tile,
                               out_band,
                               out_file,
                               bands,
                               base_bands,
                               ml_model,
                               block,
                               roi,
                               exclusion_mask,
                               impute_fn,
                               output_dir,
                               version,
                               verbose,
                               progress) {

    # Initial time for tile classification
    tile_start_time <- .tile_classif_start(
        tile = tile,
        verbose = verbose
    )
    # Create chunks to be allocated to jobs in parallel
    chunks <- .tile_chunks_create(
        tile = tile,
        overlap = 0L,
        block = block
    )
    # Create a variable to control updating of bounding box
    # by default, update_bbox is FALSE
    update_bbox <- FALSE
    if (.has(exclusion_mask)) {
        # How many chunks there are in tile?
        nchunks <- nrow(chunks)
        # Remove chunks within the exclusion mask
        chunks <- .chunks_filter_mask(
            chunks = chunks,
            mask = exclusion_mask
        )
        # Create crop region
        chunks["mask"] <- .chunks_crop_mask(
            chunks = chunks,
            mask = exclusion_mask
        )
        # Should bbox of resulting tile be updated?
        update_bbox <- nrow(chunks) != nchunks
    }
    if (.has(roi)) {
        # How many chunks do we need to process?
        nchunks <- nrow(chunks)
        # Intersect chunks with ROI
        chunks <- .chunks_filter_spatial(
            chunks = chunks,
            roi = roi
        )
        # Update bbox to account for ROI
        update_bbox <- nrow(chunks) != nchunks
    }
    # Obtain configuration parameters for probability cube
    band_conf <- .conf_derived_band(
        derived_class = "probs_cube",
        band = out_band
    )
    # Process jobs in parallel - one job per chunk
    # In parallel processing the model was exported once to the workers'
    # global environment; in sequential processing it is passed directly
    block_files <- .jobs_map_parallel_chr(
        jobs = chunks,
        fn = .classify_chunk_cpu,
        tile = tile,
        base_bands = base_bands,
        bands = bands,
        band_conf = band_conf,
        impute_fn = impute_fn,
        output_dir = output_dir,
        out_file = out_file,
        ml_model = if (.parallel_is_open()) NULL else ml_model,
        progress = progress
    )
    # Merge blocks into a new probs_cube tile
    # If ROI exists, blocks are merged to a different directory
    # than output_dir, which is used to save the final cropped version
    merge_out_file <- out_file
    if (.has(roi)) {
        merge_out_file <- .file_derived_name(
            tile = tile,
            band = out_band,
            version = version,
            output_dir = file.path(output_dir, ".sits")
        )
    }
    probs_tile <- .tile_derived_merge_blocks(
        file = merge_out_file,
        band = out_band,
        labels = .ml_labels_code(ml_model),
        base_tile = tile,
        block_files = block_files,
        derived_class = "probs_cube",
        multicores = .jobs_multicores(),
        update_bbox = update_bbox
    )
    # if there is a ROI, crop the probability cube
    if (.has(roi)) {
        probs_tile_crop <- .crop(
            cube = probs_tile,
            roi = roi,
            output_dir = output_dir,
            multicores = 1L,
            progress = progress
        )
        unlink(.fi_paths(.fi(probs_tile)))
    }
    # show final time for classification
    .tile_classif_end(
        tile = tile,
        start_time = tile_start_time,
        verbose = verbose
    )
    # Return probs tile (cropped version in case of ROI)
    if (.has(roi)) {
        probs_tile_crop
    } else {
        probs_tile
    }
}
#' @title Classify a tile of raster data in GPU using multicores
#' @name .classify_tile_gpu
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#'
#' @description Classifies one tile with \code{.classify_tiles_gpu}.
#'
#' @param  tile            Single tile of a data cube.
#' @param  out_file        Name of the output probabibility file.
#' @inheritParams .classify_tiles_gpu
#' @return The probs tile.
.classify_tile_gpu <- function(tile,
                               out_band,
                               out_file,
                               bands,
                               base_bands,
                               ml_model,
                               block,
                               roi,
                               exclusion_mask,
                               impute_fn,
                               output_dir,
                               version,
                               verbose,
                               progress) {
    .classify_tiles_gpu(
        tiles = tile,
        out_files = list(out_file),
        out_band = out_band,
        bands = bands,
        base_bands = base_bands,
        ml_model = ml_model,
        block = block,
        roi = roi,
        exclusion_mask = exclusion_mask,
        impute_fn = impute_fn,
        output_dir = output_dir,
        version = version,
        verbose = verbose,
        progress = progress
    )[[1L]]
}
#' @title Classify the tiles of a cube in GPU with one dataloader
#' @name .classify_tiles_gpu
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#'
#' @description Reads the blocks of every tile with the workers of one
#' torch dataloader, so the workers start once per call and not once per
#' tile. Each block carries the index of its tile; the callback writes it
#' to the block files of that tile. A tile is merged into its probs file
#' as soon as its last block is written, in the order of the tiles, while
#' the workers read the next blocks.
#'
#' A tile whose probs file exists is not classified again; the valid
#' block files of a tile are not computed again.
#'
#' @param  tiles           Data cube with the tiles to classify.
#' @param  out_files       List with the probs file of each tile.
#' @param  out_band        Band to be produced.
#' @param  bands           Bands to extract time series
#' @param  base_bands      Base bands to extract values
#' @param  ml_model        Model trained by \code{\link[sits]{sits_train}}.
#' @param  block           Optimized block to be read into memory.
#' @param  roi             Region of interest.
#' @param  exclusion_mask  Area to be excluded from classification
#' @param  impute_fn       Imputation function.
#' @param  output_dir      Output directory.
#' @param  version         Version of result.
#' @param  verbose         Print processing information?
#' @param  progress        Show progress bar?
#' @return List with the probs tile of each tile.
.classify_tiles_gpu <- function(tiles,
                                out_files,
                                out_band,
                                bands,
                                base_bands,
                                ml_model,
                                block,
                                roi,
                                exclusion_mask,
                                impute_fn,
                                output_dir,
                                version,
                                verbose,
                                progress) {
    n_tiles <- nrow(tiles)
    probs_tiles <- vector("list", n_tiles)
    # A tile whose probs file exists is built from it
    done <- vapply(out_files, function(f) all(file.exists(f)), logical(1L))
    for (i in which(done)) {
        .check_recovery()
        probs_tiles[[i]] <- .classify_tile_from_files(
            out_files[[i]], tiles[i, ], ml_model
        )
    }
    # Pending chunks, recovered block files and update_bbox of each tile
    plans <- lapply(seq_len(n_tiles), function(i) {
        if (done[[i]]) {
            return(NULL)
        }
        .classify_tile_gpu_chunks(
            tile = tiles[i, ],
            out_file = out_files[[i]],
            block = block,
            roi = roi,
            exclusion_mask = exclusion_mask
        )
    })
    # Block files of each tile, the recovered ones first
    block_files <- lapply(plans, function(plan) plan[["recovered_files"]])
    # Number of blocks still to be written in each tile
    n_pending <- vapply(plans, function(plan) {
        if (is.null(plan)) 0L else nrow(plan[["chunks"]])
    }, integer(1L))
    # Tiles to be merged, in the order of the cube
    to_merge <- which(!done)
    tile_start_time <- NULL
    # Merge, in order, each tile that has no pending block left
    merge_ready <- function() {
        while (length(to_merge) > 0L) {
            i <- to_merge[[1L]]
            tile <- tiles[i, ]
            if (is.null(tile_start_time)) {
                tile_start_time <<- .tile_classif_start(tile, verbose)
            }
            if (n_pending[[i]] > 0L) {
                return(invisible(NULL))
            }
            # Merge blocks into a new probs_cube tile, cropped when there
            # is a ROI
            probs_tiles[[i]] <<- .classify_merge(
                tile = tile,
                out_file = out_files[[i]],
                block_files = block_files[[i]],
                ml_model = ml_model,
                out_band = out_band,
                version = version,
                roi = roi,
                update_bbox = plans[[i]][["update_bbox"]],
                progress = progress
            )
            # The post loop of the caller starts on this line
            .tile_classif_end(
                tile = tile,
                start_time = tile_start_time,
                verbose = verbose
            )
            tile_start_time <<- NULL
            to_merge <<- to_merge[-1L]
            gc()
        }
        invisible(NULL)
    }
    # Called by the callback after it writes a block of tile i
    on_block <- function(i, files) {
        block_files[[i]] <<- c(block_files[[i]], files)
        n_pending[[i]] <<- n_pending[[i]] - 1L
        merge_ready()
    }
    # Tiles at the head with all blocks recovered are merged now
    merge_ready()
    if (sum(n_pending) > 0L) {
        # Pending chunks of every tile, each with the index of its tile
        # the crop mask column is dropped: the callback writes blocks
        # without it, and its sf geometries differ in crs among tiles
        chunks <- .map_dfr(which(n_pending > 0L), function(i) {
            chunks <- plans[[i]][["chunks"]]
            chunks[["mask"]] <- NULL
            chunks[["tile_index"]] <- i
            chunks
        })
        # Build the chunk dataset
        dataset <- .torch_chunks_dataset(
            chunks = chunks,
            tile = tiles,
            read_fn = .classify_data_read,
            bands = bands,
            base_bands = base_bands,
            stats = .ml_stats(ml_model),
            ml_features_name = .ml_features_name(ml_model),
            ml_temporal_model = .ml_torch_is_temporal(ml_model),
            impute_fn = impute_fn,
            verbose = verbose,
            output_dir = output_dir
        )
        # Get band configuration
        band_conf <- .conf_derived_band(
            derived_class = "probs_cube",
            band = out_band
        )
        # Define post-process callback
        # This callback reconstructs + writes each block in its tile
        callback <- .callback_post_classify(
            output_dir = output_dir,
            out_file = out_files,
            out_band = out_band,
            band_conf = band_conf,
            ml_labels = .ml_labels(ml_model),
            crs = slider::slide_chr(tiles, .tile_crs),
            on_block = on_block
        )
        # Classify!
        new_files <- unlist(ml_model(
            list(dataset = dataset, callback = callback)
        ))
        # Log end of chunk group save
        .debug_log(
            event = "end_chunk_group_save",
            key = "n_blocks",
            value = length(new_files)
        )
    }
    # Every block was written, so every tile was merged
    .check_that(
        length(to_merge) == 0L,
        msg = "the dataloader did not return every block of the tiles"
    )
    probs_tiles
}
#' @title Chunks of a tile to classify in GPU
#' @name .classify_tile_gpu_chunks
#' @keywords internal
#' @noRd
#' @param  tile            Single tile of a data cube.
#' @param  out_file        Probs file of the tile.
#' @param  block           Optimized block to be read into memory.
#' @param  roi             Region of interest.
#' @param  exclusion_mask  Area to be excluded from classification
#' @return List with the pending chunks, the recovered block files and
#'   update_bbox.
.classify_tile_gpu_chunks <- function(tile, out_file, block, roi,
                                      exclusion_mask) {
    # Create chunks to be allocated to jobs in parallel
    chunks <- .tile_chunks_create(
        tile = tile,
        overlap = 0L,
        block = block
    )
    # Create a variable to control updating of bounding box
    # by default, update_bbox is FALSE
    update_bbox <- FALSE
    if (.has(exclusion_mask)) {
        # How many chunks there are in tile?
        nchunks <- nrow(chunks)
        # Remove chunks within the exclusion mask
        chunks <- .chunks_filter_mask(
            chunks = chunks,
            mask = exclusion_mask
        )
        # Create crop region
        chunks["mask"] <- .chunks_crop_mask(
            chunks = chunks,
            mask = exclusion_mask
        )
        # Should bbox of resulting tile be updated?
        update_bbox <- nrow(chunks) != nchunks
    }
    if (.has(roi)) {
        # How many chunks do we need to process?
        nchunks <- nrow(chunks)
        # Intersect chunks with ROI
        chunks <- .chunks_filter_spatial(
            chunks = chunks,
            roi = roi
        )
        # Update bbox to account for ROI
        update_bbox <- nrow(chunks) != nchunks
    }
    # Regenerate all block files
    block_files_all <- slider::slide(chunks, function(chunk) {
        .classify_block_files(out_file = out_file, block = .block(chunk))
    })
    # Define which blocks are already processed
    recovered <- purrr::map_lgl(block_files_all, function(f) {
        all(.raster_is_valid(f))
    })
    list(
        chunks = chunks[!recovered, , drop = FALSE],
        recovered_files = unlist(block_files_all[recovered]),
        update_bbox = update_bbox
    )
}
#' @title Probs file of a tile
#' @name .classify_out_file
#' @keywords internal
#' @noRd
.classify_out_file <- function(tile, ml_model, version, output_dir) {
    .file_derived_name(
        tile = tile, band = "probs", version = version,
        output_dir = output_dir
    )
}
#' @title Block files of the probs files of a tile, in their dirs
#' @name .classify_block_files
#' @keywords internal
#' @noRd
.classify_block_files <- function(out_file, block) {
    stats::setNames(vapply(out_file, function(file) {
        .file_block_name(
            pattern = .file_pattern(file), block = block,
            output_dir = dirname(file)
        )
    }, character(1L), USE.NAMES = FALSE), names(out_file))
}
#' @title Probs tile from its file
#' @name .classify_tile_from_files
#' @keywords internal
#' @noRd
.classify_tile_from_files <- function(out_file, tile, ml_model) {
    .tile_derived_from_file(
        file = out_file,
        band = "probs",
        base_tile = tile,
        labels = .ml_labels_code(ml_model),
        derived_class = "probs_cube",
        update_bbox = TRUE
    )
}
#' @title Merge the blocks of a tile into its probs file
#' @name .classify_merge
#' @keywords internal
#' @noRd
#' @description With a ROI, the merged file goes to .sits and is cropped
#'   into its dir.
.classify_merge <- function(tile, out_file, block_files, ml_model, out_band,
                            version, roi, update_bbox, progress) {
    merge_out_file <- out_file
    if (.has(roi)) {
        merge_out_file <- .file_derived_name(
            tile = tile,
            band = out_band,
            version = version,
            output_dir = file.path(dirname(out_file), ".sits")
        )
    }
    probs_tile <- .tile_derived_merge_blocks(
        file = merge_out_file,
        band = out_band,
        labels = .ml_labels_code(ml_model),
        base_tile = tile,
        block_files = unname(block_files),
        derived_class = "probs_cube",
        multicores = .jobs_multicores(),
        update_bbox = update_bbox
    )
    if (!.has(roi)) {
        return(probs_tile)
    }
    probs_tile_crop <- .crop(
        cube = probs_tile,
        roi = roi,
        output_dir = dirname(out_file),
        multicores = 1L,
        progress = progress
    )
    unlink(.fi_paths(.fi(probs_tile)))
    probs_tile_crop
}
.classify_write_block <- function(values, block, output_dir, out_file, out_band) {
    # Get block geometry
    block_geom <- .block(block)
    # Create a temporary block file name
    block_file <- .file_block_name(
        pattern = .file_pattern(out_file),
        block = block_geom,
        output_dir = output_dir
    )
    # Resume processing in case of failure
    if (all(.raster_is_valid(block_file))) {
        return(block_file)
    }
    # Obtain configuration parameters for probability cube
    band_conf <- .conf_derived_band(
        derived_class = "probs_cube",
        band = out_band
    )
    # Log start of block saving
    .debug_log(
        event = "start_block_data_save",
        key = "file",
        value = block_file
    )
    # Prepare and save results as raster
    .raster_write_block(
        files = block_file,
        block = block_geom,
        bbox = .bbox(block),
        values = values,
        data_type = .data_type(band_conf),
        missing_value = .miss_value(band_conf),
        crop_block = block[["mask"]]
    )
    # Log end of block saving
    .debug_log(
        event = "end_block_data_save",
        key = "file",
        value = block_file
    )
    # Free memory
    gc()
    # Returned block file
    block_file
}
#' @title Read a block of values from a set of raster images
#' @name  .classify_data_read
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#'
#' @param  tile            Input tile to read data.
#' @param  block           Bounding box in (col, row, ncols, nrows).
#' @param  bands           Bands to extract time series
#' @param  base_bands      Base bands to extract values
#' @param  ml_features_name Names of features used in model trained by
#'                         \code{\link[sits]{sits_train}}.
#' @param  impute_fn       Imputation function
#' @return A matrix with values for classification.
.classify_data_read <- function(tile, block, bands, base_bands,
                                ml_features_name, impute_fn) {
    # For cubes that have a time limit to expire (MPC cubes only)
    tile <- .cube_token_generator(tile)
    # Read and preprocess values of cloud
    # Get cloud values (NULL if data does not exist)
    cloud_mask <- .tile_cloud_read_block(
        tile = tile,
        block = block
    )
    # Read and preprocess values of each eo band
    values <- purrr::map(bands, function(band) {
        # Get band values (stops if band not found)
        values <- .tile_read_block(
            tile = tile,
            band = band,
            block = block
        )
        # Log
        .debug_log(
            event = "start_block_data_process",
            key = "process",
            value = "cloud-impute-filter"
        )
        # Remove cloud masked pixels
        if (.has(cloud_mask)) {
            values[cloud_mask] <- NA
        }
        # are there NA values? interpolate them
        if (anyNA(values)) {
            values <- impute_fn(values)
        }
        # Log
        .debug_log(
            event = "end_block_data_process",
            key = "band",
            value = band
        )
        # Return values
        as.data.frame(values)
    })
    # Read and preprocess values of each base band
    values_base <- purrr::map(base_bands, function(band) {
        # Read and preprocess values of each base band
        values_base <- .tile_read_block(
            tile = .tile_base_info(tile),
            band = band,
            block = block
        )
        # Return values
        as.data.frame(values_base)
    })
    # Combine two lists
    values <- c(values, values_base)
    # collapse list to get data.frame
    values <- suppressMessages(
        purrr::list_cbind(values, name_repair = "universal")
    )
    # Compose final values
    values <- as.matrix(values)
    # Set values features name
    if (.has(ml_features_name)) {
        colnames(values) <- ml_features_name
    }
    # Return values
    values
}
#' @title Classify a sits tibble using machine learning models
#' @name .classify_ts
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#'
#' @description Returns a sits tibble with the results of the ML classifier.
#'
#' @param  samples    Tibble with sits samples
#' @param  ml_model   Model trained by \code{\link[sits]{sits_train}}.
#' @param  impute_fn  Imputation function (to remove NA)
#' @param  multicores number of threads to process the time series.
#' @param  gpu_memory Memory available in GPU
#' @param  progress   Show progress bar?
#' @return A tibble with the predicted labels.
.classify_ts <- function(samples,
                         ml_model,
                         impute_fn,
                         multicores,
                         gpu_memory,
                         progress) {
    # Prepare parallel processing
    started <- .parallel_start(workers = multicores)
    on.exit(.parallel_stop(started), add = TRUE)

    # Get bands from model
    bands <- .ml_bands(ml_model)
    # Update samples bands order
    if (any(bands != .samples_bands(samples))) {
        samples <- .samples_select_bands(
            samples = samples,
            bands = bands
        )
    }
    # Apply imputation filter
    if (.has(impute_fn)) {
        samples <- .apply_across(
            data = samples,
            fn = impute_fn
        )
    }
    # Compute the breaks in time for multiyear classification
    class_info <- .timeline_class_info(
        data = samples,
        samples = .ml_samples(ml_model)
    )
    # Split long time series of samples in a set of small time series
    if (length(class_info[["dates_index"]][[1L]]) > 1L) {
        splitted <- .samples_split(
            samples = samples,
            split_intervals = class_info[["dates_index"]][[1L]]
        )
        pred <- .predictors(
            samples = splitted,
            ml_model = ml_model
        )
        # Post condition: is predictor data valid?
        .check_predictors(
            pred = pred,
            samples = splitted
        )
    } else {
        # Convert samples time series in predictors and preprocess data
        pred <- .predictors(
            samples = samples,
            ml_model = ml_model
        )
    }
    # choose between GPU and CPU
    if (.torch_gpu_classification() && .ml_is_torch_model(ml_model)) {
        prediction <- .classify_ts_gpu(
            pred = pred,
            ml_model = ml_model,
            gpu_memory = gpu_memory
        )
    } else {
        prediction <- .classify_ts_cpu(
            pred = pred,
            ml_model = ml_model,
            multicores = multicores,
            progress = progress
        )
    }
    # Store the result in the input data
    if (length(class_info[["dates_index"]][[1L]]) > 1L) {
        prediction <- .tibble_prediction_multiyear(
            data = samples,
            class_info = class_info,
            prediction = prediction
        )
    } else {
        prediction <- .tibble_prediction(
            data = samples,
            prediction = prediction
        )
    }
    # Set result class and return it
    prediction <- .set_class(
        prediction, paste(class(samples)[[1L]], "predicted", sep = "_"),
        class(samples)
    )
    prediction
}
#' @title Classify predictors using CPU
#' @name .classify_ts_cpu
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#'
#' @description Returns a sits tibble with the results of the ML classifier.
#'
#' @param  pred       a tibble with predictors
#' @param  ml_model   model trained by \code{\link[sits]{sits_train}}.
#' @param  multicores number of threads to process the time series.
#' @param  progress   Show progress bar?
#' @return A tibble with the predicted values.
.classify_ts_cpu <- function(pred,
                             ml_model,
                             multicores,
                             progress) {
    # Divide samples predictors in chunks to parallel processing
    parts <- .pred_create_partition(
        pred = pred,
        partitions = multicores
    )
    # Do parallel process
    prediction <- .jobs_map_parallel_dfr(parts, function(part) {
        # Extract predictors of a given partition
        # Extract features to classify
        # Classify, normalize and calibrate value
        values <- part |>
            .pred_part() |>
            .pred_features() |>
            ml_model() |>
            .ml_normalize(ml_model)
        # Extract columns
        values_columns <- colnames(values)
        # Transform classification results
        values <- tibble::tibble(as.data.frame(values))
        # Fix column names to avoid errors with non-standard column name
        # (e.g., with spaces, icons)
        colnames(values) <- values_columns
        # Return classification
        values
    }, progress = progress)
    prediction
}
#' @title Classify predictors using GPU
#' @name .classify_ts_gpu
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#'
#' @description Returns a sits tibble with the results of the ML classifier.
#'
#' @param  pred       a tibble with predictors
#' @param  ml_model   model trained by \code{\link[sits]{sits_train}}.
#' @param  gpu_memory memory available in GPU
#' @return A tibble with the predicted values.
.classify_ts_gpu <- function(pred,
                             ml_model,
                             gpu_memory) {
    # estimate size of GPU memory required (in GB)
    pred_size <- nrow(pred) * ncol(pred) * 8.0 / 1000000000.0
    # estimate how should we partition the predictors
    num_parts <- ceiling(pred_size / gpu_memory)
    # Divide samples predictors in chunks to parallel processing
    parts <- .pred_create_partition(
        pred = pred,
        partitions = num_parts
    )
    prediction <- slider::slide_dfr(parts, function(part) {
        # Get predictors of a given partition
        # Get predictors features to classify
        # Classify
        # normalize and calibrate values
        values <- part |>
            .pred_part() |>
            .pred_features() |>
            ml_model() |>
            .ml_normalize(ml_model)
        # Extract columns
        values_columns <- colnames(values)
        # Transform classification results
        values <- tibble::tibble(as.data.frame(values))
        # Fix column names to avoid errors with non-standard column name
        # (e.g., with spaces, icons)
        colnames(values) <- values_columns
        # Clean GPU memory
        .ml_gpu_clean(ml_model)
        values
    })
    prediction
}
#' @title Start recording processing time
#' @name .classify_verbose_start
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @description Prints the block size and computes
#' start time for processing
#'
#' @param  verbose TRUE/FALSE
#' @param  block   block size
#' @return start time for processing
.classify_verbose_start <- function(verbose, block) {
    if (verbose) {
        msg <- paste0(
            .conf("messages", ".verbose_block_size"), " ",
            .nrows(block), " x ", .ncols(block)
        )
        message(msg)
    }
    Sys.time()
}
#' @title End recording processing time
#' @name .classify_verbose_end
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @description Prints the processing time
#' @param  verbose TRUE/FALSE
#' @param  start_time   initial processing time
#' @return elapsed processing time for processing
.classify_verbose_end <- function(verbose, start_time) {
    if (verbose) {
        end_time <- Sys.time()
        message("")
        message(.conf("messages", ".verbose_task_end"), end_time)
        message(
            .conf("messages", ".verbose_task_elapsed"),
            format(round(end_time - start_time, digits = 2L))
        )
    }
}
#' @title Classify a chunk of data on CPU
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @keywords internal
#' @noRd
#' @param  chunk        Chunk to be processed
#' @param  tile         Input data tile
#' @param  base_bands  Base bands
#' @param  bands        Bands to be used
#' @param  band_conf    Band configuration
#' @param  impute_fn    Imputation function
#' @param  output_dir   Output directory
#' @param  out_file     Output file
#' @param  ml_model     Model closure (sequential processing only); in
#'                      parallel processing it is NULL and the model is
#'                      resolved from the worker global environment, where
#'                      it was exported once by .parallel_start()
#' @return              Block file path
.classify_chunk_cpu <- function(chunk,
                                tile,
                                base_bands,
                                bands,
                                band_conf,
                                impute_fn,
                                output_dir,
                                out_file,
                                ml_model = NULL) {
    # Get exported ml_model (parallel processing)
    if (is.null(ml_model)) {
        ml_model <- get("ml_model", envir = globalenv())
    }
    # Retrive block to be processed
    block <- .block(chunk)
    # Create a temporary block file name
    block_file <- .file_block_name(
        pattern = .file_pattern(out_file),
        block = block,
        output_dir = output_dir
    )
    # Resume processing in case of failure
    if (all(.raster_is_valid(block_file))) {
        return(block_file)
    }
    # Read and preprocess values from files
    values <- .classify_data_read(
        tile = tile,
        block = block,
        bands = bands,
        base_bands = base_bands,
        ml_features_name = .ml_features_name(ml_model),
        impute_fn = impute_fn
    )
    # Get mask of NA pixels
    na_mask <- C_mask_na(values)
    # Filter out NA pixels - only classify valid pixels
    values <- values[!na_mask, , drop = FALSE]
    # Define control variable to check for correct termination
    input_pixels <- nrow(values)
    # Start log file
    .debug_log(
        event = "start_block_data_classification",
        key = "model",
        value = .ml_class(ml_model)
    )
    # Apply the classification model only to valid (non-NA) pixels
    if (input_pixels > 0L) {
        # Apply the classification model to values
        # Uses the closure created by sits_train
        values <- ml_model(values)
        # Normalize and calibrate the values
        # Perform softmax for torch models
        values <- .ml_normalize(values, ml_model)
        # Are the results consistent with the data input?
        .check_processed_values(
            values = values,
            input_pixels = input_pixels
        )
        # apply scale and offset
        offset <- .offset(band_conf)
        if (.has(offset) && offset != 0.0) {
            values <- values - offset
        }
        scale <- .scale(band_conf)
        if (.has(scale) && scale != 1.0) {
            values <- values / scale
        }
    }
    # Log end of block
    .debug_log(
        event = "end_block_data_classification",
        key = "model",
        value = .ml_class(ml_model)
    )
    # Reconstruct full output matrix with NA for masked pixels
    n_labels <- length(.ml_labels(ml_model))
    full_values <- matrix(
        NA_real_,
        nrow = length(na_mask),
        ncol = n_labels,
        dimnames = list(NULL, .ml_labels(ml_model))
    )
    if (input_pixels > 0L) {
        full_values[!na_mask, ] <- values
    }
    # Log start of block saving
    .debug_log(
        event = "start_block_data_save",
        key = "file",
        value = block_file
    )
    # Prepare and save results as raster
    .raster_write_block(
        files = block_file,
        block = block,
        bbox = .bbox(chunk),
        values = full_values,
        data_type = .data_type(band_conf),
        missing_value = .miss_value(band_conf),
        crop_block = chunk[["mask"]]
    )
    # Log end of block saving
    .debug_log(
        event = "end_block_data_save",
        key = "file",
        value = block_file
    )
    # Free memory
    gc()
    # Returned block file
    block_file
}
