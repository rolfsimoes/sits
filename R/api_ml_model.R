#' @title Return machine learning model inside a closure
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#'
#' @param  ml_model  Closure that contains ML model and its environment
#' @return           ML model as specified by the original ML function
.ml_model <- function(ml_model) {
    env <- environment(ml_model)

    if (.ml_is_torch_model(ml_model)) {
        return(.torch_model_restore(ml_model))
    }
    if ("model" %in% ls(env, all.names = TRUE)) {
        return(env[["model"]])
    }


    stop(.conf("messages", ".ml_model"))
}
#' @title Return statistics of ML model inside a closure (new version)
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#'
#' @param  ml_model  Closure that contains ML model and its environment
#' @return           Data statistics contained in the model closure
.ml_stats <- function(ml_model) {
    # New stats variable
    environment(ml_model)[["ml_stats"]]
}
#' @title Symbols used by code
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @description Walks calls, pairlists, lists and the formals and body of
#' functions, and returns the names of the symbols found. Constants are
#' ignored. The right side of \code{$} and \code{@@} and the package of
#' \code{::} are names, not variables, and are not returned.
#' @param  x  A function, a call or a list of them
#' @return    Character vector with symbol names
.ml_symbols <- function(x) {
    if (is.symbol(x)) {
        return(setdiff(as.character(x), ""))
    }
    if (is.function(x)) {
        if (is.primitive(x)) {
            return(character(0L))
        }
        return(union(.ml_symbols(formals(x)), .ml_symbols(body(x))))
    }
    if (is.call(x)) {
        head <- x[[1L]]
        if (identical(head, as.name("::")) || identical(head, as.name(":::"))) {
            return(character(0L))
        }
        if (identical(head, as.name("$")) || identical(head, as.name("@"))) {
            return(.ml_symbols(x[[2L]]))
        }
    }
    if (is.call(x) || is.pairlist(x) || is.expression(x) || is.list(x)) {
        x <- as.list(x)
        symbols <- character(0L)
        for (i in seq_along(x)) {
            # an argument without default is the empty symbol
            if (!identical(x[[i]], quote(expr = ))) {
                symbols <- union(symbols, .ml_symbols(x[[i]]))
            }
        }
        return(symbols)
    }
    character(0L)
}
#' @title Environment with the objects that code reads
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @description Looks up the symbols of \code{code} in \code{env} and its
#' parents, up to the namespace, and copies the objects found into a new
#' environment whose parent is the namespace. A function copied from those
#' frames is walked too and moved to the new environment.
#' @param  code  A function, a call or a list of them
#' @param  env   Environment where \code{code} was defined
#' @return       New environment
.ml_env <- function(code, env) {
    namespace <- topenv(env)
    frames <- list()
    while (!identical(env, namespace) && !identical(env, emptyenv())) {
        frames[[length(frames) + 1L]] <- env
        env <- parent.env(env)
    }
    is_frame <- function(e) any(purrr::map_lgl(frames, identical, e))
    result <- new.env(parent = namespace)
    pending <- .ml_symbols(code)
    while (.has(pending)) {
        name <- pending[[1L]]
        pending <- pending[-1L]
        frame <- purrr::detect(frames, function(e) {
            exists(name, envir = e, inherits = FALSE)
        })
        if (is.null(frame) || exists(name, envir = result, inherits = FALSE)) {
            next
        }
        # a binding to a missing argument has no value
        value <- tryCatch(
            get(name, envir = frame, inherits = FALSE),
            error = function(e) NULL
        )
        if (is.function(value) && is_frame(environment(value))) {
            pending <- union(pending, .ml_symbols(value))
            environment(value) <- result
        }
        assign(name, value, envir = result)
    }
    result
}
#' @title Detach a function from the frame where it was defined
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @description For a function defined in a training function and kept by
#' the model (e.g. a loss kept by a fitted luz model): the function keeps
#' only the objects its code reads, not the training frame.
#' @param  fun  Function
#' @return      The function, with a new environment
.ml_function <- function(fun) {
    environment(fun) <- .ml_env(fun, environment(fun))
    fun
}
#' @title Keep in a model only what prediction reads
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @description Gives the predict function a new environment with the
#' objects its code reads and the metadata of the training samples (labels,
#' bands, timeline, reference dates, feature names), so the training data
#' and the frames of the constructor are not kept.
#' @param  predict_fun  Predict function of a model or encoder
#' @param  samples      Training samples
#' @param  ...          Other metadata of the model (e.g. embedding_dim)
#' @param  keep         Names of other objects of the training frame that
#'                      users read from the model environment (torch_seed)
#' @return              The predict function, with its new environment
.ml_compact <- function(predict_fun, samples, ..., keep = character(0L)) {
    code <- c(list(predict_fun), purrr::map(keep, as.name))
    env <- .ml_env(code, environment(predict_fun))
    env[["ml_metadata"]] <- list(
        labels = .samples_labels(samples),
        bands = .samples_bands(samples),
        timeline = .samples_timeline(samples),
        start_date = samples[["start_date"]][[1L]],
        end_date = samples[["end_date"]][[1L]],
        features_name = names(.predictors(samples[1L, ]))[-2L:0L],
        ...
    )
    environment(predict_fun) <- env
    predict_fun
}
#' @title Return the metadata of the training samples of a model
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @param  ml_model  Closure that contains ML model and its environment
#' @return           List with labels, bands, timeline, reference dates,
#'                   feature names and other metadata
.ml_metadata <- function(ml_model) {
    env <- environment(ml_model)
    if (.has(env[["ml_metadata"]])) {
        return(env[["ml_metadata"]])
    }
    # models trained before the metadata was kept hold their samples
    samples <- env[["samples"]]
    list(
        labels = .samples_labels(samples),
        bands = .samples_bands(samples),
        timeline = .samples_timeline(samples),
        start_date = samples[["start_date"]][[1L]],
        end_date = samples[["end_date"]][[1L]],
        features_name = names(env[["train_samples"]])[-2L:0L],
        embedding_dim = env[["embedding_dim"]]
    )
}
#' @title Return class of ML model inside a closure
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @param  ml_model  Closure that contains ML model and its environment
#' @return           ML model class
.ml_class <- function(ml_model) {
    class(ml_model)[[1L]]
}
#' @title Return names of features used to train ML model
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @param  ml_model  Closure that contains ML model and its environment
#' @return           Features used to build the model
.ml_features_name <- function(ml_model) {
    .ml_metadata(ml_model)[["features_name"]]
}
#' @title Return names of bands used to train ML model
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @param  ml_model  Closure that contains ML model and its environment
#' @return           Bands used to build the model
.ml_bands <- function(ml_model) {
    .ml_metadata(ml_model)[["bands"]]
}
#' @title Return timeline of the samples used to train ML model
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @param  ml_model  Closure that contains ML model and its environment
#' @return           Timeline of the first training sample
.ml_timeline <- function(ml_model) {
    as.Date(.ml_metadata(ml_model)[["timeline"]])
}
#' @title Return labels of samples of used to train ML model
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @param  ml_model  Closure that contains ML model and its environment
#' @return           Sample labels used to build the model
.ml_labels <- function(ml_model) {
    .ml_metadata(ml_model)[["labels"]]
}
#' @title Return codes of sample labels of used to train ML model
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @param  ml_model  Closure that contains ML model and its environment
#' @return           Codes of sample labels used to build the model
.ml_labels_code <- function(ml_model) {
    labels <- .ml_labels(ml_model)
    names(labels) <- seq_along(labels)
    labels
}
#' @title Clean GPU memory allocation
#' @keywords internal
#' @noRd
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @param  ml_model  Closure that contains ML model and its environment
#' @return           Called for side effects
.ml_gpu_clean <- function(ml_model) {
    # Clean torch allocations
    if (.torch_cuda_enabled() && .ml_is_torch_model(ml_model)) {
        torch::cuda_empty_cache()
    }
    invisible(NULL)
}

#' @title normalize the probability results
#' @keywords internal
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @noRd
#' @param  values    Values to be normalized
#' @param  ml_model  Closure that contains ML model and its environment
#' @return           Normalized values
#'
.ml_normalize <- function(values, ml_model) {
    UseMethod(".ml_normalize", ml_model)
}
#' @export
#'
.ml_normalize.torch_model <- function(values, ml_model) {
    # Correct the default behaviour of softmax in torch models
    # Run softmax here instead of inside a torch model
    column_names <- colnames(values)
    values[is.na(values)] <- 0
    values <- softmax(values)
    colnames(values) <- column_names
    values
}
#' @export
#'
.ml_normalize.default <- function(values, ml_model) {
    values[is.na(values)] <- 0
    values
}
#' @title Update multicores for models that do internal multiprocessing
#' @keywords internal
#' @noRd
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @param  ml_model   Closure that contains ML model and its environment
#' @param  multicores Current multicores setting
#' @return            Updated multicores
#'
.ml_update_multicores <- function(ml_model, multicores) {
    # xgboost model has internal multiprocessing
    if ("xgb_model" %in% .ml_class(ml_model)) {
        multicores <- 1L
    }
    # torch in GPU has internal multiprocessing
    if (.torch_gpu_classification() && .ml_is_torch_model(ml_model)) {
        multicores <- 1L
    }

    multicores
}
#' @title Is the ML model a torch model?
#' @keywords internal
#' @noRd
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @param  ml_model   Closure that contains ML model and its environment
#' @return            TRUE/FALSE
.ml_is_torch_model <- function(ml_model) {
    inherits(ml_model, "torch_model")
}

#' @title Is the ML torch model temporal?
#' @keywords internal
#' @noRd
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @param  ml_model   Closure that contains ML model and its environment
#' @return            TRUE/FALSE
.ml_torch_is_temporal <- function(ml_model) {
    # models trained before the mlp_model class keep the module generator
    mlp <- inherits(ml_model, "mlp_model") ||
        "mlp_model" %in% ls(environment(ml_model))
    .ml_is_torch_model(ml_model) && !mlp
}
