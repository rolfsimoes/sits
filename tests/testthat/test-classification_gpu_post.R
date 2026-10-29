# One small TempCNN model of samples_modis_ndvi
gpu_post_model <- function() {
    set.seed(101L)
    torch::torch_manual_seed(101L)
    sits_train(samples_modis_ndvi, sits_tempcnn(epochs = 2L, verbose = FALSE))
}

gpu_post_cube <- function() {
    sits_cube(
        source = "BDC",
        collection = "MOD13Q1-6.1",
        data_dir = system.file("extdata/raster/mod13q1", package = "sits"),
        progress = FALSE
    )
}

# "cpu": the block path in the R session (one core, so the code under test
# runs, not the installed package of the workers); "dataloader": the GPU
# path, mocked on, with CPU tensors
gpu_post_classify <- function(path, cube, ml_model, output_dir, version,
                              ...) {
    run <- function() {
        sits_classify(
            data = cube, ml_model = ml_model, memsize = 4L,
            multicores = 1L, output_dir = output_dir, version = version,
            progress = FALSE, ...
        )
    }
    testthat::with_mocked_bindings(
        run(), .torch_gpu_classification = function() path == "dataloader"
    )
}

gpu_post_values <- function(cube) {
    .raster_open_rast(.fi(cube)[["path"]])[]
}

# Four blocks of the 147 x 255 test cube
gpu_post_blocks <- c(nrows = 40L, ncols = 255L)

test_that("the dataloader path takes softmax and scaling from the device", {
    skip_on_cran()
    skip_if_not_installed("torch")
    skip_if_not_installed("luz")
    model <- gpu_post_model()
    cube <- gpu_post_cube()
    root <- file.path(tempdir(), "device_probs_one")
    unlink(root, recursive = TRUE)
    dir.create(root)
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    dirs <- file.path(root, c("cpu", "gpu"))
    for (d in dirs) dir.create(d)
    cpu <- gpu_post_classify("cpu", cube, model, dirs[[1L]], "v",
                             block_size = gpu_post_blocks)
    # the R softmax is not called on the dataloader path
    gpu <- testthat::with_mocked_bindings(
        gpu_post_classify("dataloader", cube, model, dirs[[2L]], "v",
                          block_size = gpu_post_blocks),
        softmax = function(values) stop("softmax ran in R")
    )
    ref <- gpu_post_values(cpu)
    got <- gpu_post_values(gpu)
    # within one unit of the INT2U scale; NA pixels in the same place
    expect_equal(is.na(got), is.na(ref))
    expect_lte(max(abs(got - ref), na.rm = TRUE), 1L)
})

test_that("the dataloader path writes a block after the next is predicted", {
    skip_on_cran()
    skip_if_not_installed("torch")
    skip_if_not_installed("luz")
    model <- gpu_post_model()
    cube <- gpu_post_cube()
    dir <- file.path(tempdir(), "deferred_write")
    unlink(dir, recursive = TRUE)
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    events <- character(0L)
    read <- .classify_data_read
    write <- .classify_write_block
    # one core: the dataloader reads each block in the R session, just
    # before the model runs on it
    probs <- testthat::with_mocked_bindings(
        gpu_post_classify("dataloader", cube, model, dir, "v",
                          block_size = gpu_post_blocks),
        .classify_data_read = function(...) {
            events <<- c(events, "read")
            read(...)
        },
        .classify_write_block = function(...) {
            events <<- c(events, "write")
            write(...)
        }
    )
    n_blocks <- ceiling(147L / gpu_post_blocks[["nrows"]])
    expect_equal(
        events,
        c("read", rep(c("read", "write"), n_blocks - 1L), "write")
    )
    expect_s3_class(probs, "probs_cube")
    expect_false(all(is.na(gpu_post_values(probs))))
})

test_that("the dataloader path frees the CUDA cache after each tile", {
    skip_on_cran()
    skip_if_not_installed("torch")
    skip_if_not_installed("luz")
    model <- gpu_post_model()
    cube <- gpu_post_cube()
    # two tiles, the second a copy of the first under another name
    cube <- .set_class(rbind(cube, cube), class(cube))
    cube[["tile"]][[2L]] <- "copy"
    dir <- file.path(tempdir(), "gpu_clean")
    unlink(dir, recursive = TRUE)
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    events <- character(0L)
    merge <- .classify_merge
    clean <- .ml_gpu_clean
    testthat::with_mocked_bindings(
        gpu_post_classify("dataloader", cube, model, dir, "v"),
        .classify_merge = function(...) {
            events <<- c(events, "merge")
            merge(...)
        },
        .ml_gpu_clean = function(...) {
            events <<- c(events, "clean")
            clean(...)
        }
    )
    # one clean per tile, then the one at the end of sits_classify()
    expect_equal(events, c("merge", "clean", "merge", "clean", "clean"))
    # without CUDA the clean does not touch the cache
    testthat::with_mocked_bindings(
        expect_null(.ml_gpu_clean(model)),
        .torch_cuda_enabled = function() FALSE
    )
    emptied <- 0L
    testthat::with_mocked_bindings(
        testthat::with_mocked_bindings(
            .ml_gpu_clean(model),
            .torch_cuda_enabled = function() TRUE
        ),
        cuda_empty_cache = function() emptied <<- emptied + 1L,
        .package = "torch"
    )
    expect_equal(emptied, 1L)
})
