# Expected values computed outside sits; the first one, 0.127011122046577,
# is the value published with L'Ecuyer's RngStreams
rng_base <- c(10407L, rep(12345L, 6L))
rng_set_state <- function(state) {
    assign(".Random.seed", state, envir = globalenv())
}
rng_get_state <- function() {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
}
test_that("substreams have the expected L'Ecuyer-CMRG states", {
    old <- RNGkind()
    on.exit(do.call(RNGkind, as.list(old)), add = TRUE)
    RNGkind("L'Ecuyer-CMRG")
    rng_set_state(rng_base)
    subs <- sits:::.rng_substreams(3L)
    expect_identical(
        subs[[1L]][-1L],
        c(870504860L, -1653269569L, 884013853L,
          339352413L, -1920660590L, -643363409L)
    )
    expect_identical(
        subs[[2L]][-1L],
        c(460387934L, 1532391390L, 877287553L,
          120103512L, -2141851355L, 335837774L)
    )
    expect_identical(
        subs[[3L]][-1L],
        c(-519857236L, -1086671252L, 1257177538L,
          378684317L, -1427855118L, -2093661213L)
    )
    # the caller moves past the last substream it handed out
    expect_identical(
        rng_get_state()[-1L],
        c(1870130441L, 490396226L, 1081325149L,
          -609881575L, -1946924678L, 1094489500L)
    )
    rng_set_state(subs[[3L]])
    expect_equal(
        stats::runif(2L), c(0.5032122888761005, 0.16517391832456343),
        tolerance = 1e-15
    )
})
test_that("a stream given as seed is used as given", {
    set.seed(1L, kind = "Mersenne-Twister")
    before <- rng_get_state()
    f <- function(seed) {
        rng_exit <- sits:::.rng_enter(seed)
        on.exit(rng_exit(), add = TRUE)
        stats::runif(3L)
    }
    expect_equal(
        f(rng_base),
        c(0.12701112204657714, 0.3185275653967945, 0.3091860155832701),
        tolerance = 1e-15
    )
    # with seed, the state as received, kind included
    expect_identical(rng_get_state(), before)
    expect_identical(RNGkind()[[1L]], "Mersenne-Twister")
})
test_that("an integer seed gives the same draws and leaves the state", {
    set.seed(5L, kind = "Mersenne-Twister")
    before <- rng_get_state()
    f <- function(seed) {
        rng_exit <- sits:::.rng_enter(seed)
        on.exit(rng_exit(), add = TRUE)
        c(RNGkind()[[1L]], stats::runif(2L))
    }
    a <- f(42L)
    expect_identical(a[[1L]], "L'Ecuyer-CMRG")
    expect_identical(f(42L), a)
    expect_false(identical(f(43L), a))
    expect_identical(rng_get_state(), before)
})
test_that("without seed, a L'Ecuyer-CMRG global state goes on", {
    old <- RNGkind()
    on.exit(do.call(RNGkind, as.list(old)), add = TRUE)
    RNGkind("L'Ecuyer-CMRG")
    rng_set_state(rng_base)
    f <- function() {
        rng_exit <- sits:::.rng_enter(NULL)
        on.exit(rng_exit(), add = TRUE)
        stats::runif(3L)
    }
    expect_equal(
        f(),
        c(0.12701112204657714, 0.3185275653967945, 0.3091860155832701),
        tolerance = 1e-15
    )
    # the restored state includes the draws spent
    expect_identical(
        rng_get_state()[-1L],
        c(-1271176443L, -1271176443L, -909607723L,
          -1816685032L, 1655725443L, 2057415812L)
    )
    expect_equal(stats::runif(1L), 0.8258468629271136, tolerance = 1e-15)
})
test_that("without seed, another global generator gives six draws", {
    set.seed(7L, kind = "Mersenne-Twister")
    invisible(stats::runif(6L))
    after_six <- rng_get_state()
    set.seed(7L, kind = "Mersenne-Twister")
    f <- function() {
        rng_exit <- sits:::.rng_enter(NULL)
        on.exit(rng_exit(), add = TRUE)
        c(RNGkind()[[1L]], stats::runif(2L))
    }
    a <- f()
    expect_identical(a[[1L]], "L'Ecuyer-CMRG")
    # the kind comes back, and the state includes the six draws
    expect_identical(RNGkind()[[1L]], "Mersenne-Twister")
    expect_identical(rng_get_state(), after_six)
    # two calls in a row get two streams; set.seed() repeats them
    b <- f()
    expect_false(identical(a, b))
    set.seed(7L, kind = "Mersenne-Twister")
    expect_identical(f(), a)
})
test_that("an invalid seed is refused", {
    f <- function(seed) {
        rng_exit <- sits:::.rng_enter(seed)
        on.exit(rng_exit(), add = TRUE)
        stats::runif(1L)
    }
    expect_error(f("a"))
    expect_error(f(c(1L, 2L)))
    expect_error(f(c(10403L, rep(1L, 6L))))
})
