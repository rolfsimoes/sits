#' @title Stack torch models to classify a cube with one read
#' @name sits_stack_model
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#'
#' @description Takes named torch models trained on the same bands and
#' timeline and returns one model. \code{\link[sits]{sits_classify}} with
#' it reads each block of the cube once and runs every model on it, each
#' with its own normalization statistics and its own softmax. The
#' probability cube has the labels of all models, as \code{name:label};
#' \code{\link[sits]{sits_stack_split}} writes one probability cube per
#' model.
#'
#' @param ...  Two or more named torch models
#'             (trained by \code{\link[sits]{sits_train}}).
#' @return A model to be used by \code{\link[sits]{sits_classify}}.
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
#'     probs_ab <- sits_stack_split(probs, output_dir = tempdir())
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

#' @title Write one probability cube per model of a stack
#' @name sits_stack_split
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#'
#' @description Takes the probability cube made by
#' \code{\link[sits]{sits_classify}} with a model of
#' \code{\link[sits]{sits_stack_model}} and writes, for each model, a
#' probability cube with its labels, whose version is the model name.
#' A file that exists is kept.
#'
#' @param cube        Probability cube with \code{name:label} labels.
#' @param output_dir  Directory of the output files, or one directory per
#'                    model, named by model.
#' @return A named list of probability cubes, one per model.
#' @export
sits_stack_split <- function(cube, output_dir) {
    .check_is_probs_cube(cube)
    labels <- .cube_labels(cube)
    .check_that(
        all(grepl(":", labels, fixed = TRUE)),
        msg = "sits_stack_split() needs the labels of a stack (name:label)"
    )
    names_of <- sub(":.*$", "", labels)
    models <- unique(names_of)
    if (length(output_dir) == 1L && is.null(names(output_dir))) {
        output_dir <- stats::setNames(rep(output_dir, length(models)), models)
    }
    .check_that(
        setequal(names(output_dir), models) &&
            length(output_dir) == length(models),
        msg = "output_dir is one dir, or one per model, named by model"
    )
    for (dir in output_dir) .check_output_dir(dir)
    cubes <- lapply(models, function(name) {
        index <- which(names_of == name)
        .cube_foreach_tile(cube, function(tile) {
            out_file <- .file_derived_name(
                tile = tile, band = "probs", version = name,
                output_dir = output_dir[[name]]
            )
            if (!file.exists(out_file)) {
                tmp_file <- paste0(out_file, ".tmp.tif")
                params <- c(
                    list("-of" = "GTiff"),
                    stats::setNames(as.list(as.character(index)),
                                    rep("-b", length(index))),
                    stats::setNames(as.list(.conf("gdal_presets", "image", "co")),
                                    rep("-co", length(.conf("gdal_presets", "image", "co"))))
                )
                .gdal_translate(
                    file = tmp_file, base_file = .tile_path(tile),
                    params = params, quiet = TRUE
                )
                file.rename(tmp_file, out_file)
            }
            # codes from 1, as .ml_labels_code() of the model alone
            model_labels <- sub("^[^:]*:", "", unname(labels[index]))
            names(model_labels) <- seq_along(model_labels)
            .tile_derived_from_file(
                file = out_file, band = "probs", base_tile = tile,
                derived_class = "probs_cube", labels = model_labels
            )
        })
    })
    stats::setNames(cubes, models)
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
