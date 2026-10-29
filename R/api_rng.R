#' @title Enter the random number stream of a call
#' @name .rng_enter
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @description Switches the R generator to the L'Ecuyer-CMRG stream of a
#' call and returns the function that leaves it, to be passed to `on.exit()`.
#' The caller separates its own independent draws by stream; its parallel
#' tasks get substreams (`.rng_substreams()`).
#' @param seed `NULL`, an integer, or a L'Ecuyer-CMRG state (the integer
#' vector of length 7 kept in `.Random.seed`).
#' @return A function with no arguments that restores the generator kind
#' and state. With `seed`, it is the state as received; without `seed`, the
#' state includes the draws spent.
.rng_enter <- function(seed = NULL) {
    .check_set_caller(".rng_enter")
    if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        set.seed(NULL)
    }
    received <- .rng_state()
    if (.has(seed)) {
        .check_that(
            is.numeric(seed) && !anyNA(seed) && all(seed == round(seed)) &&
                (length(seed) == 1L ||
                    (length(seed) == 7L && seed[[1L]] %% 100L == 7L))
        )
        if (length(seed) == 7L) {
            .rng_set_state(as.integer(seed))
        } else {
            set.seed(seed, kind = "L'Ecuyer-CMRG")
        }
        restore <- received
    } else if (.rng_is_lecuyer(received)) {
        # the call goes on along the global stream
        restore <- NULL
    } else {
        stream <- .rng_state_from_draws(received)
        # the received generator, after the six draws spent
        restore <- .rng_state()
        .rng_set_state(stream)
    }
    function() {
        if (!is.null(restore)) {
            .rng_set_state(restore)
        }
    }
}
#' @title Fill a L'Ecuyer-CMRG state with six draws of another generator
#' @name .rng_state_from_draws
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @description One integer gives 2^31 starting points; six draws give about
#' 2^191, so two calls do not share a stream by chance.
#' @param received `.Random.seed` of the global generator.
#' @return A L'Ecuyer-CMRG state, with the normal and sample kinds received.
.rng_state_from_draws <- function(received) {
    moduli <- c(rep(4294967087, 3L), rep(4294944443, 3L))
    values <- floor(stats::runif(6L) * moduli)
    # a component of all zeros is not a valid state
    if (all(values[1L:3L] == 0)) values[[1L]] <- 1
    if (all(values[4L:6L] == 0)) values[[4L]] <- 1
    # R keeps the unsigned 32-bit values as signed integers
    values <- ifelse(values >= 2^31, values - 2^32, values)
    kind <- received[[1L]] - received[[1L]] %% 100L + 7L
    as.integer(c(kind, values))
}
#' @title Substreams for the tasks of a call
#' @name .rng_substreams
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @description Task i gets the i-th substream after the current state, so
#' the draws of a task depend on its index only. The current state moves to
#' the substream after the last one, so the next draws of the caller and the
#' tasks of its next parallel map do not reuse them.
#' @param n Number of tasks.
#' @return A list of n L'Ecuyer-CMRG states.
.rng_substreams <- function(n) {
    .check_set_caller(".rng_substreams")
    state <- .rng_state()
    .check_that(.rng_is_lecuyer(state))
    substreams <- vector("list", n)
    for (i in seq_len(n)) {
        state <- parallel::nextRNGSubStream(state)
        substreams[[i]] <- state
    }
    .rng_set_state(parallel::nextRNGSubStream(state))
    substreams
}
#' @title Wrap a task so that it draws from its substream
#' @name .rng_task
#' @keywords internal
#' @noRd
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#' @description Defined at the top level so that the wrapper carries only
#' `fn` when it is sent to a worker.
#' @param fn Function applied to the value of each task.
#' @return A function of a task, a list with `value` and `state`.
.rng_task <- function(fn) {
    force(fn)
    function(task, ...) {
        assign(".Random.seed", task[["state"]], envir = globalenv())
        fn(task[["value"]], ...)
    }
}
#' @title Is a state of the L'Ecuyer-CMRG generator?
#' @keywords internal
#' @noRd
#' @param state `.Random.seed` value.
#' @return TRUE or FALSE.
.rng_is_lecuyer <- function(state) {
    length(state) == 7L && state[[1L]] %% 100L == 7L
}
#' @title Read the state of the R generator
#' @keywords internal
#' @noRd
#' @return `.Random.seed` of the global environment.
.rng_state <- function() {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
}
#' @title Set the state of the R generator
#' @keywords internal
#' @noRd
#' @description The first element of the state holds the generator kind, so
#' the kind follows the state.
#' @param state `.Random.seed` value.
#' @return Called for side effects.
.rng_set_state <- function(state) {
    assign(".Random.seed", state, envir = globalenv())
    invisible(state)
}
