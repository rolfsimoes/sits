test_that("torch bayes smoother matches C++ output (neigh_fraction = 1.0)", {
    set.seed(42)
    nrows <- 20L
    ncols <- 20L
    nbands <- 4L
    window_size <- 5L
    smoothness <- rep(20.0, nbands)
    neigh_fraction <- 1.0

    # Synthetic logit matrix — values in a range that avoids ±Inf
    logits <- matrix(
        runif(nrows * ncols * nbands, min = -3, max = 3),
        nrow = nrows * ncols,
        ncol = nbands
    )

    cpp_result <- bayes_smoother_fraction(
        logits         = logits,
        nrows          = nrows,
        ncols          = ncols,
        window_size    = window_size,
        smoothness     = smoothness,
        neigh_fraction = neigh_fraction
    )

    torch_result <- .torch_smooth_bayes_fraction(
        logits     = logits,
        nrows      = nrows,
        ncols      = ncols,
        smoothness = smoothness,
        kernel     = .torch_smooth_kernel(window_size, neigh_fraction)
    )

    # Interior pixels (excluding the leg-width border) must match closely
    leg      <- window_size %/% 2L
    row_ids  <- (leg + 1L):(nrows - leg)
    col_ids  <- (leg + 1L):(ncols - leg)
    interior <- as.vector(
        outer(row_ids - 1L, col_ids, FUN = function(r, c) r * ncols + c)
    )

    max_diff <- max(abs(cpp_result[interior, ] - torch_result[interior, ]))
    expect_lt(max_diff, 1e-3)
})

test_that("torch bayes smoother matches C++ output (neigh_fraction = 0.5)", {
    set.seed(123)
    nrows <- 20L
    ncols <- 20L
    nbands <- 4L
    window_size <- 5L
    smoothness <- c(10.0, 20.0, 15.0, 25.0)
    neigh_fraction <- 0.5

    logits <- matrix(
        runif(nrows * ncols * nbands, min = -3, max = 3),
        nrow = nrows * ncols,
        ncol = nbands
    )

    cpp_result <- bayes_smoother_fraction(
        logits         = logits,
        nrows          = nrows,
        ncols          = ncols,
        window_size    = window_size,
        smoothness     = smoothness,
        neigh_fraction = neigh_fraction
    )

    torch_result <- .torch_smooth_bayes_fraction(
        logits     = logits,
        nrows      = nrows,
        ncols      = ncols,
        smoothness = smoothness,
        kernel     = .torch_smooth_kernel(window_size, neigh_fraction)
    )

    leg      <- window_size %/% 2L
    row_ids  <- (leg + 1L):(nrows - leg)
    col_ids  <- (leg + 1L):(ncols - leg)
    interior <- as.vector(
        outer(row_ids - 1L, col_ids, FUN = function(r, c) r * ncols + c)
    )

    max_diff <- max(abs(cpp_result[interior, ] - torch_result[interior, ]))
    expect_lt(max_diff, 1e-3)
})

test_that("compiled torch smoother matches the R tensor kernel", {
    skip_if_not_installed("torch")
    # R tensor kernel of sits dev (c8f8d744), the reference of the compiled one
    r_kernel <- function(logits, nrows, ncols, window_size, smoothness,
                         neigh_fraction) {
        device <- "cpu"
        nbands <- ncol(logits)
        leg <- window_size %/% 2L
        win_sq <- as.integer(window_size * window_size)
        x <- torch::torch_tensor(
            t(logits),
            dtype  = torch::torch_float32(),
            device = device
        )$view(c(nbands, nrows, ncols))$unsqueeze(1L)
        x_pad <- torch::nnf_pad(
            x,
            pad  = c(leg, leg, leg, leg),
            mode = "reflect"
        )
        zero_t <- torch::torch_tensor(
            0.0,
            dtype = torch::torch_float32(),
            device = device
        )
        neginf_t <- torch::torch_tensor(
            -Inf,
            dtype = torch::torch_float32(),
            device = device
        )
        pos <- torch::torch_tensor(
            seq_len(win_sq),
            dtype  = torch::torch_float32(),
            device = device
        )$view(c(1L, win_sq))
        x0_all <- torch::torch_tensor(
            logits,
            dtype = torch::torch_float32(),
            device = device
        )
        band_results <- purrr::map(seq_len(nbands), function(b) {
            x_b <- x_pad[, b, , , drop = FALSE]
            wins_b <- torch::nnf_unfold(
                x_b,
                kernel_size = window_size,
                stride      = 1L,
                padding     = 0L
            )$squeeze(1L)$t()
            nan_mask_b <- torch::torch_isnan(wins_b)
            if (neigh_fraction == 1.0) {
                selected_b <- torch::torch_where(nan_mask_b, zero_t, wins_b)
                sel_mask_b <- !nan_mask_b
                n_b <- (!nan_mask_b)$sum(dim = -1L)$to(
                    dtype = torch::torch_float32()
                )
            } else {
                wins_sort_b <- torch::torch_where(nan_mask_b, neginf_t, wins_b)
                sorted_b <- torch::torch_sort(
                    wins_sort_b,
                    dim = -1L,
                    descending = TRUE
                )[[1L]]
                valid_b <- (!nan_mask_b)$sum(dim = -1L)$to(
                    dtype = torch::torch_float32()
                )
                neigh_hi_b <- torch::torch_ceil(
                    neigh_fraction * valid_b
                )$clamp_min(1L)
                sel_mask_b <- pos$le(neigh_hi_b$unsqueeze(-1L))
                selected_b <- torch::torch_where(sel_mask_b, sorted_b, zero_t)
                n_b        <- neigh_hi_b
            }
            m0_b <- selected_b$sum(dim = -1L) / n_b
            diff_sq_b <- torch::torch_where(
                sel_mask_b,
                (selected_b - m0_b$unsqueeze(-1L))$pow(2L),
                zero_t
            )
            s0_b <- torch::torch_where(
                n_b$gt(1.0),
                diff_sq_b$sum(dim = -1L) / (n_b - 1.0),
                zero_t
            )
            x0_b    <- x0_all[, b]
            w_b     <- s0_b / (s0_b + smoothness[b])
            bayes_b <- w_b * x0_b + (1.0 - w_b) * m0_b
            use_m0  <- torch::torch_isnan(x0_b) | s0_b$lt(1e-4)
            torch::torch_where(use_m0, m0_b, bayes_b)
        })
        result <- torch::torch_stack(band_results, dim = 2L)
        as.matrix(
            result$cpu()$to(dtype = torch::torch_float64())
        )
    }
    set.seed(7)
    nrows <- 40L
    ncols <- 30L
    nbands <- 4L
    window_size <- 7L
    smoothness <- c(10.0, 20.0, 15.0, 25.0)
    logits <- matrix(
        runif(nrows * ncols * nbands, min = -3, max = 3),
        nrow = nrows * ncols,
        ncol = nbands
    )
    # a patch of NA smaller than the window and scattered NA pixels
    na_patch <- as.vector(
        outer(9:11, 15:17, FUN = function(r, c) r * ncols + c)
    )
    logits[c(na_patch, sample(nrows * ncols, 20L)), ] <- NA
    purrr::walk(c(0.5, 1.0), function(neigh_fraction) {
        r_result <- r_kernel(
            logits         = logits,
            nrows          = nrows,
            ncols          = ncols,
            window_size    = window_size,
            smoothness     = smoothness,
            neigh_fraction = neigh_fraction
        )
        compiled_result <- .torch_smooth_bayes_fraction(
            logits     = logits,
            nrows      = nrows,
            ncols      = ncols,
            smoothness = smoothness,
            kernel     = .torch_smooth_kernel(window_size, neigh_fraction)
        )
        expect_equal(is.na(compiled_result), is.na(r_result))
        expect_lt(max(abs(compiled_result - r_result), na.rm = TRUE), 1e-4)
    })
})

test_that("torch smoother availability supports fallback and override", {
    old_force_cpp <- Sys.getenv("SITS_SMOOTH_FORCE_CPP", unset = NA_character_)
    on.exit(
        if (is.na(old_force_cpp)) {
            Sys.unsetenv("SITS_SMOOTH_FORCE_CPP")
        } else {
            Sys.setenv(SITS_SMOOTH_FORCE_CPP = old_force_cpp)
        }
    )
    smooth_available <- function(functional, force_cpp = "FALSE") {
        testthat::local_mocked_bindings(
            .torch_is_functional = function() functional
        )
        Sys.setenv(SITS_SMOOTH_FORCE_CPP = force_cpp)
        .torch_smooth_available()
    }
    expect_true(smooth_available(functional = TRUE))
    expect_false(smooth_available(functional = FALSE))
    expect_false(smooth_available(functional = TRUE, force_cpp = "TRUE"))
    # a kernel that does not compile gives the C++ smoother
    old_doc_mode <- Sys.getenv("SITS_DOCUMENTATION_MODE")
    Sys.setenv(SITS_DOCUMENTATION_MODE = "FALSE")
    on.exit(Sys.setenv(SITS_DOCUMENTATION_MODE = old_doc_mode), add = TRUE)
    testthat::local_mocked_bindings(
        .torch_smooth_kernel = function(...) stop("no reflection_pad2d")
    )
    expect_warning(
        fallback_fn <- .smooth_fn_bayes_torch(5L, 0.5, rep(20.0, 3L)),
        .conf("messages", ".smooth_fn_bayes_torch"),
        fixed = TRUE
    )
    cpp_fn <- .smooth_fn_bayes(5L, 0.5, rep(20.0, 3L))
    set.seed(5)
    values <- matrix(runif(300L), ncol = 3L)
    block <- list(nrows = 10L, ncols = 10L)
    expect_equal(fallback_fn(values, block), cpp_fn(values, block))
})

test_that("torch smoother accounts for materialized neighborhoods", {
    expect_equal(
        .torch_smooth_block_memsize(
            block_size = 100,
            window_size = 9L
        ),
        0.000324
    )
})

test_that("Bayesian smoothers validate neighborhood parameters", {
    probs_cube <- structure(list(), class = "probs_cube")
    probs_vector_cube <- structure(list(), class = "probs_vector_cube")
    testthat::local_mocked_bindings(
        .check_raster_cube_files = function(cube) invisible(cube)
    )

    expect_error(
        sits_smooth(
            probs_cube,
            window_size = 3L,
            output_dir = tempdir()
        )
    )
    expect_error(
        sits_smooth(
            probs_cube,
            window_size = 23L,
            output_dir = tempdir()
        )
    )
    expect_error(
        sits_smooth(
            probs_cube,
            neigh_fraction = 0.0,
            output_dir = tempdir()
        )
    )
    expect_error(
        sits_smooth(
            probs_vector_cube,
            neigh_fraction = 0.0,
            output_dir = tempdir()
        )
    )
    expect_error(
        sits_smooth(
            probs_cube,
            progress = "TRUE",
            output_dir = tempdir()
        )
    )
})

test_that("torch functional check disables automatic installation", {
    old_torch_install <- Sys.getenv("TORCH_INSTALL", unset = NA_character_)
    on.exit(
        if (is.na(old_torch_install)) {
            Sys.unsetenv("TORCH_INSTALL")
        } else {
            Sys.setenv(TORCH_INSTALL = old_torch_install)
        }
    )
    Sys.setenv(TORCH_INSTALL = "1")
    testthat::local_mocked_bindings(
        .torch_package_is_installed = function() TRUE,
        .torch_native_is_installed = function() {
            expect_equal(Sys.getenv("TORCH_INSTALL"), "0")
            TRUE
        }
    )
    expect_true(.torch_is_functional())
    expect_equal(Sys.getenv("TORCH_INSTALL"), "1")
})

test_that("sits_smooth selects CPU and torch implementations", {
    skip_on_cran()
    skip_if_not_installed("torch")

    data_dir <- system.file("extdata/raster/mod13q1", package = "sits")
    raster_cube <- sits_cube(
        source     = "BDC",
        collection = "MOD13Q1-6.1",
        data_dir   = data_dir,
        tiles      = "012010",
        bands      = "NDVI",
        start_date = "2013-09-14",
        end_date   = "2014-08-29",
        multicores = 1L,
        progress   = FALSE
    )
    rfor_model <- sits_train(
        samples_modis_ndvi,
        sits_rfor(num_trees = 40L)
    )

    out_dir <- file.path(tempdir(), "smooth-torch-comparison")
    # Purge stale cache files from any previous run
    unlink(out_dir, recursive = TRUE)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    probs_cube <- suppressWarnings(
        sits_classify(
            data = raster_cube,
            ml_model = rfor_model,
            multicores = 1L,
            output_dir = out_dir,
            progress = FALSE
        )
    )

    # Exercise both implementations through the public API. Torch falls back
    # to CPU tensors when no CUDA or MPS device is available.
    smooth_backend <- function(use_torch, version) {
        testthat::local_mocked_bindings(
            .torch_smooth_available = function() use_torch,
            .torch_gpu_available = function() FALSE
        )
        sits_smooth(
            cube = probs_cube,
            output_dir = out_dir,
            version = version,
            multicores = 1L,
            progress = FALSE
        )
    }
    smooth_cpp <- smooth_backend(use_torch = FALSE, version = "v1-cpp")
    smooth_torch <- smooth_backend(use_torch = TRUE, version = "v1-torch")

    # Compare interior pixels only.
    # The absolute tile boundary (leg = window_size %/% 2 = 4 pixels) uses
    # reflect padding in torch vs. a boundary-inclusive mirror in C++, so
    # those edge pixels are expected to differ — they are cropped by the
    # overlap mechanism for all interior chunks anyway.
    cpp_rast <- .raster_open_rast(.fi_paths(.fi(smooth_cpp)))
    torch_rast <- .raster_open_rast(.fi_paths(.fi(smooth_torch)))
    nrow_r <- .raster_nrows(cpp_rast)
    ncol_r <- .raster_ncols(cpp_rast)

    leg <- 9L %/% 2L
    row_ids <- (leg + 1L):(nrow_r - leg)
    col_ids <- (leg + 1L):(ncol_r - leg)
    interior <- as.vector(
        outer(row_ids - 1L, col_ids, FUN = function(r, c) r * ncol_r + c)
    )

    cpp_vals <- .raster_get_values(cpp_rast)
    torch_vals <- .raster_get_values(torch_rast)

    valid <- is.finite(cpp_vals[interior, 1L]) &
        is.finite(torch_vals[interior, 1L])
    pix <- interior[valid]

    max_diff <- max(abs(cpp_vals[pix, ] - torch_vals[pix, ]))
    expect_lt(max_diff, 5.0)
})
