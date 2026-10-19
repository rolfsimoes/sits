#' @title Stack torch models to classify a cube with one read
#' @name sits_stack_model
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#'
#' @description Takes named torch models trained on the same bands and
#' timeline and returns one model. \code{\link[sits]{sits_classify}} with
#' it reads each block of the cube once and runs every model on it, each
#' with its own normalization statistics and its own softmax. Each
#' model writes its own probability file, whose version is the model name,
#' to its own \code{output_dir}; \code{sits_classify} returns a named list
#' of probability cubes, one per model.
#'
#' @param ...  Two or more named torch models
#'             (trained by \code{\link[sits]{sits_train}}).
#' @return A model to be used by \code{\link[sits]{sits_classify}}, whose
#'   \code{output_dir} is one directory, or one per model, named by model.
#'
#' @examples
#' if (sits_run_examples()) {
#'     labels <- sits_labels(samples_modis_ndvi)
#'     a <- sits_train(samples_modis_ndvi, sits_tempcnn(epochs = 2))
#'     b <- sits_train(
#'         samples_modis_ndvi[samples_modis_ndvi$label != labels[[1]], ],
#'         sits_tempcnn(epochs = 2)
#'     )
#'     stack <- sits_stack_model(a = a, b = b)
#'     cube <- sits_cube(
#'         source = "BDC",
#'         collection = "MOD13Q1-6.1",
#'         data_dir = system.file("extdata/raster/mod13q1", package = "sits")
#'     )
#'     probs <- sits_classify(cube, ml_model = stack, output_dir = tempdir())
#'     probs[["a"]]
#' }
#' @export
sits_stack_model <- function(...) {
    models <- list(...)
    .check_that(
        length(models) >= 2L && !is.null(names(models)) &&
            all(nzchar(names(models))) &&
            !anyDuplicated(names(models)),
        msg = "sits_stack_model() needs two or more models with unique names"
    )
    .check_that(
        !any(grepl(":", names(models), fixed = TRUE)),
        msg = "model names prefix the labels as name:label, so they have no colon"
    )
    .check_that(
        all(vapply(models, .ml_is_torch_model, logical(1L))),
        msg = "sits_stack_model() takes torch models only"
    )
    first <- models[[1L]]
    same <- function(fn) {
        all(vapply(models, function(m) identical(fn(m), fn(first)), logical(1L)))
    }
    .check_that(
        same(.ml_bands) && same(function(m) .samples_ntimes(.ml_samples(m))) &&
            same(.ml_features_name) && same(.ml_torch_is_temporal),
        msg = "stacked models need the same bands, timeline and layout"
    )
    # labels in each model's own order, prefixed by the model name
    groups <- lapply(models, .ml_labels)
    labels <- unlist(Map(paste0, paste0(names(models), ":"), groups),
                     use.names = FALSE)
    samples <- .ml_samples(first)
    # only the names of the predictors are read (.ml_features_name)
    train_samples <- environment(first)[["train_samples"]][0L, ]
    ml_stats <- NULL
    if (!.ml_torch_is_temporal(first)) {
        # .ml_torch_is_temporal() looks for this name
        mlp_model <- TRUE
    }
    predict_fun <- function(values) {
        if (.torch_use_dataloader(values)) {
            # one module that normalizes the raw block for each model
            vessel <- .torch_model_restore(models[[1L]])
            vessel[["model"]] <- .stack_module(
                modules = lapply(models, function(m) {
                    .torch_model_restore(m)[["model"]]
                }),
                stats = lapply(models, .ml_stats),
                temporal = .ml_torch_is_temporal(models[[1L]]),
                n_bands = length(.ml_bands(models[[1L]])),
                n_times = .samples_ntimes(.ml_samples(models[[1L]]))
            )
            return(unlist(.torch_predict_chunks(
                torch_model = vessel,
                dataset = values[["dataset"]],
                callback = values[["callback"]]
            )))
        }
        # predict_fun normalizes in place, so each model gets a copy
        values <- do.call(cbind, lapply(models, function(m) m(values + 0.0)))
        colnames(values) <- labels
        values
    }
    .set_class(predict_fun, "stack_model", "torch_model", "sits_model",
               class(predict_fun))
}

#' @title The models of a stack, named
#' @keywords internal
#' @noRd
.stack_models <- function(ml_model) {
    environment(ml_model)[["models"]]
}

#' @title One output dir per model of a stack, named by model
#' @keywords internal
#' @noRd
#' @param output_dir One dir, or one per model, named by model.
.stack_output_dir <- function(output_dir, ml_model) {
    models <- names(.stack_models(ml_model))
    if (length(output_dir) == 1L && is.null(names(output_dir))) {
        output_dir <- stats::setNames(rep(output_dir, length(models)), models)
    }
    .check_that(
        length(output_dir) == length(models) &&
            setequal(names(output_dir), models),
        msg = "output_dir is one dir, or one per model, named by model"
    )
    for (dir in output_dir) .check_output_dir(dir)
    output_dir[models]
}

#' @title Probs file of a tile: one per model of a stack, named by model
#' @keywords internal
#' @noRd
#' @description A stack names each file with its model name as version.
.classify_out_file <- function(tile, ml_model, version, output_dir) {
    if (!inherits(ml_model, "stack_model")) {
        return(.file_derived_name(
            tile = tile, band = "probs", version = version,
            output_dir = output_dir
        ))
    }
    models <- names(.stack_models(ml_model))
    stats::setNames(vapply(models, function(name) {
        .file_derived_name(
            tile = tile, band = "probs", version = name,
            output_dir = output_dir[[name]]
        )
    }, character(1L), USE.NAMES = FALSE), models)
}

#' @title Block files of the probs files of a tile, in their dirs
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

#' @title Write the probs of a block, each model of a stack to its file
#' @keywords internal
#' @noRd
.classify_block_write <- function(files, values, ml_model, block, bbox,
                                  band_conf, crop_block) {
    columns <- list(seq_len(ncol(values)))
    if (inherits(ml_model, "stack_model")) {
        model_of <- sub(":.*$", "", .ml_labels(ml_model))
        columns <- lapply(names(files), function(name) which(model_of == name))
    }
    for (i in seq_along(files)) {
        .raster_write_block(
            files = files[[i]],
            block = block,
            bbox = bbox,
            values = values[, columns[[i]], drop = FALSE],
            data_type = .data_type(band_conf),
            missing_value = .miss_value(band_conf),
            crop_block = crop_block
        )
    }
    files
}

#' @title Probs tile from its files: a named list of tiles for a stack
#' @keywords internal
#' @noRd
.classify_tile_from_files <- function(out_file, tile, ml_model) {
    one <- function(file, model) {
        .tile_derived_from_file(
            file = file,
            band = "probs",
            base_tile = tile,
            labels = .ml_labels_code(model),
            derived_class = "probs_cube",
            update_bbox = TRUE
        )
    }
    if (!inherits(ml_model, "stack_model")) {
        return(one(out_file, ml_model))
    }
    Map(one, out_file, .stack_models(ml_model)[names(out_file)])
}

#' @title Merge the blocks of a tile into its probs files
#' @keywords internal
#' @noRd
#' @description One probs tile, or a named list of tiles for a stack. With
#'   a ROI, the merged file goes to .sits and is cropped into its dir.
.classify_merge <- function(tile, out_file, block_files, ml_model, out_band,
                            version, roi, update_bbox, progress) {
    one <- function(file, blocks, model, file_version) {
        merge_out_file <- file
        if (.has(roi)) {
            merge_out_file <- .file_derived_name(
                tile = tile,
                band = out_band,
                version = file_version,
                output_dir = file.path(dirname(file), ".sits")
            )
        }
        probs_tile <- .tile_derived_merge_blocks(
            file = merge_out_file,
            band = out_band,
            labels = .ml_labels_code(model),
            base_tile = tile,
            block_files = blocks,
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
            output_dir = dirname(file),
            multicores = 1L,
            progress = progress
        )
        unlink(.fi_paths(.fi(probs_tile)))
        probs_tile_crop
    }
    if (!inherits(ml_model, "stack_model")) {
        return(one(out_file, unname(block_files), ml_model, version))
    }
    models <- .stack_models(ml_model)
    stats::setNames(lapply(names(out_file), function(name) {
        one(out_file[[name]], unname(block_files[names(block_files) == name]),
            models[[name]], name)
    }), names(out_file))
}

#' @title A torch module that runs the models of a stack on one block
#' @name .stack_module
#' @keywords internal
#' @noRd
#' @description The block comes raw (the stack has no statistics); each
#' model gets it normalized with its own q02 and q98, as
#' .torch_chunks_dataset does for one model.
.stack_module <- torch::nn_module(
    ".stack_module",
    initialize = function(modules, stats, temporal, n_bands, n_times) {
        self$stack <- torch::nn_module_list(unname(modules))
        self$n_models <- length(modules)
        for (k in seq_len(self$n_models)) {
            # stats are band-major, one value per band and date
            shape_stats <- function(x) {
                x <- torch::torch_tensor(as.numeric(x), dtype = torch::torch_float())
                if (temporal) x <- x$reshape(c(n_bands, n_times))$t()
                x
            }
            q02 <- shape_stats(.stats_q02(stats[[k]]))
            q98 <- shape_stats(.stats_q98(stats[[k]]))
            self[[paste0("min_", k)]] <- torch::nn_buffer(q02)
            self[[paste0("range_", k)]] <- torch::nn_buffer(q98 - q02)
        }
    },
    forward = function(x) {
        outputs <- lapply(seq_len(self$n_models), function(k) {
            values <- torch::torch_clamp(
                (x - self[[paste0("min_", k)]]) / self[[paste0("range_", k)]],
                min = 0.0001, max = 1.0
            )
            module <- self$stack[[k]]
            if (is.null(module$predict)) module(values) else module$predict(values)
        })
        torch::torch_cat(outputs, dim = 2L)
    }
)
