# Two small TempCNN models on different labels of samples_modis_ndvi
stack_test_models <- function() {
    labels <- .samples_labels(samples_modis_ndvi)
    train <- function(keep, seed, factor) {
        set.seed(seed)
        torch::torch_manual_seed(seed)
        samples <- samples_modis_ndvi[samples_modis_ndvi[["label"]] %in% keep, ]
        # another NDVI scale gives the model other normalization statistics
        samples[["time_series"]] <- lapply(samples[["time_series"]], function(ts) {
            ts[["NDVI"]] <- ts[["NDVI"]] * factor
            ts
        })
        sits_train(samples, sits_tempcnn(epochs = 2L, verbose = FALSE))
    }
    list(
        a = train(labels, 101L, 1.0),
        b = train(labels[-1L], 202L, 0.7)
    )
}

stack_test_cube <- function() {
    sits_cube(
        source = "BDC",
        collection = "MOD13Q1-6.1",
        data_dir = system.file("extdata/raster/mod13q1", package = "sits"),
        progress = FALSE
    )
}

# "cpu": the block path in the R session (one core, so the code under test
# runs, not the installed package of the workers); "workers": the same path
# on two cores; "dataloader": the GPU path, mocked on, with CPU tensors
classify_with <- function(path, cube, ml_model, output_dir, version, ...) {
    run <- function() {
        sits_classify(
            data = cube, ml_model = ml_model, memsize = 4L,
            multicores = if (path == "workers") 2L else 1L,
            output_dir = output_dir, version = version, progress = FALSE, ...
        )
    }
    if (path == "dataloader") {
        testthat::with_mocked_bindings(
            run(), .torch_gpu_classification = function() TRUE
        )
    } else {
        testthat::with_mocked_bindings(
            run(), .torch_gpu_classification = function() FALSE
        )
    }
}

probs_values <- function(cube) {
    .raster_open_rast(.fi(cube)[["path"]])[]
}

test_that("sits_stack_model checks the models and names the labels", {
    skip_on_cran()
    skip_if_not_installed("torch")
    models <- stack_test_models()
    stack <- sits_stack_model(a = models[["a"]], b = models[["b"]])
    expect_s3_class(stack, "stack_model")
    expect_s3_class(stack, "torch_model")
    expect_equal(
        .ml_labels(stack),
        c(paste0("a:", .ml_labels(models[["a"]])),
          paste0("b:", .ml_labels(models[["b"]])))
    )
    expect_equal(.ml_bands(stack), .ml_bands(models[["a"]]))
    expect_null(.ml_stats(stack))
    # unnamed, a single model, a non-torch model
    expect_error(sits_stack_model(models[["a"]], models[["b"]]))
    expect_error(sits_stack_model(a = models[["a"]]))
    # the name prefixes the labels, so it has no colon
    expect_error(sits_stack_model(`a:1` = models[["a"]], b = models[["b"]]))
    rfor <- sits_train(samples_modis_ndvi, sits_rfor(num_trees = 20L))
    expect_error(sits_stack_model(a = models[["a"]], r = rfor))
    # other bands: NDVI and a copy of it as EVI
    samples_2b <- samples_modis_ndvi
    samples_2b[["time_series"]] <- lapply(
        samples_2b[["time_series"]], function(ts) {
            ts[["EVI"]] <- ts[["NDVI"]]
            ts
        }
    )
    set.seed(3L)
    other <- sits_train(samples_2b, sits_tempcnn(epochs = 1L, verbose = FALSE))
    expect_error(sits_stack_model(a = models[["a"]], o = other))
})

# Every .tif directly in a dir, and its number of bands
tifs_in <- function(dir) {
    files <- list.files(dir, pattern = "\\.tif$", full.names = TRUE)
    vapply(files, function(f) as.integer(.raster_nlayers(.raster_open_rast(f))), integer(1L))
}

test_that("a stack writes each model's probs to its own dir, on every path", {
    skip_on_cran()
    skip_if_not_installed("torch")
    skip_if_not_installed("luz")
    models <- stack_test_models()
    stack <- sits_stack_model(a = models[["a"]], b = models[["b"]])
    cube <- stack_test_cube()
    for (path in c("cpu", "workers", "dataloader")) {
        root <- file.path(tempdir(), paste0("stack_", path))
        unlink(root, recursive = TRUE)
        dirs <- c(a = file.path(root, "a"), b = file.path(root, "b"),
                  ref = file.path(root, "ref"))
        for (d in c(root, dirs)) dir.create(d)
        probs_a <- classify_with(path, cube, models[["a"]], dirs[["ref"]], "a")
        probs_b <- classify_with(path, cube, models[["b"]], dirs[["ref"]], "b")
        probs <- classify_with(path, cube, stack, dirs[c("a", "b")], "ignored")
        expect_named(probs, c("a", "b"))
        expect_s3_class(probs[["a"]], "probs_cube")
        expect_equal(.cube_labels(probs[["a"]]), .cube_labels(probs_a))
        expect_equal(.cube_labels(probs[["b"]]), .cube_labels(probs_b))
        # one file per model in its own dir, with the model's labels only
        expect_equal(unname(tifs_in(dirs[["a"]])), length(.ml_labels(models[["a"]])))
        expect_equal(unname(tifs_in(dirs[["b"]])), length(.ml_labels(models[["b"]])))
        expect_equal(dirname(.fi(probs[["a"]])[["path"]]), unname(dirs[["a"]]))
        expect_match(.fi(probs[["b"]])[["path"]], "_probs_b\\.tif$")
        # no block file is left
        expect_length(list.files(file.path(dirs[c("a", "b")], ".sits"),
                                 pattern = "block"), 0L)
        # within one unit of the INT2U scale, as the CPU/GPU test
        expect_lte(max(abs(probs_values(probs[["a"]]) - probs_values(probs_a))), 1L)
        expect_lte(max(abs(probs_values(probs[["b"]]) - probs_values(probs_b))), 1L)
        unlink(root, recursive = TRUE)
    }
})

test_that("a stack writes to one dir, named by model, and keeps done files", {
    skip_on_cran()
    skip_if_not_installed("torch")
    skip_if_not_installed("luz")
    models <- stack_test_models()
    stack <- sits_stack_model(a = models[["a"]], b = models[["b"]])
    cube <- stack_test_cube()
    dir <- file.path(tempdir(), "stack_one_dir")
    unlink(dir, recursive = TRUE)
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    probs <- classify_with("cpu", cube, stack, dir, "ignored")
    expect_setequal(basename(names(tifs_in(dir))),
                    basename(c(.fi(probs[["a"]])[["path"]], .fi(probs[["b"]])[["path"]])))
    expect_match(.fi(probs[["a"]])[["path"]], "_probs_a\\.tif$")
    # a second call keeps the files
    path <- .fi(probs[["b"]])[["path"]]
    mtime <- file.mtime(path)
    Sys.sleep(1.1)
    again <- classify_with("cpu", cube, stack, dir, "ignored")
    expect_equal(file.mtime(path), mtime)
    expect_equal(.cube_labels(again[["b"]]), .cube_labels(probs[["b"]]))
    # dirs not named by the models of the stack
    expect_error(classify_with("cpu", cube, stack, c(a = dir, c = dir), "x"))
})

# Four blocks of the 147 x 255 test cube
test_blocks <- c(nrows = 40L, ncols = 255L)

test_that("a stack takes the softmax of each model from the device", {
    skip_on_cran()
    skip_if_not_installed("torch")
    skip_if_not_installed("luz")
    models <- stack_test_models()
    stack <- sits_stack_model(a = models[["a"]], b = models[["b"]])
    cube <- stack_test_cube()
    root <- file.path(tempdir(), "device_probs")
    unlink(root, recursive = TRUE)
    dir.create(root)
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    dirs <- file.path(root, c("stack_cpu", "stack_gpu"))
    for (d in dirs) dir.create(d)
    cpu <- classify_with("cpu", cube, stack, dirs[[1L]], "v",
                         block_size = test_blocks)
    # the R softmax is not called on the dataloader path
    gpu <- testthat::with_mocked_bindings(
        classify_with("dataloader", cube, stack, dirs[[2L]], "v",
                      block_size = test_blocks),
        softmax = function(values) stop("softmax ran in R")
    )
    for (m in names(cpu)) {
        ref <- probs_values(cpu[[m]])
        got <- probs_values(gpu[[m]])
        # within one unit of the INT2U scale; NA pixels in the same place
        expect_equal(is.na(got), is.na(ref))
        expect_lte(max(abs(got - ref), na.rm = TRUE), 1L)
    }
})
