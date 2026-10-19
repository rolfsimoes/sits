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
classify_with <- function(path, cube, ml_model, output_dir, version) {
    run <- function() {
        sits_classify(
            data = cube, ml_model = ml_model, memsize = 4L,
            multicores = if (path == "workers") 2L else 1L,
            output_dir = output_dir, version = version, progress = FALSE
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

test_that("a stack gives each model its own probs, on both paths", {
    skip_on_cran()
    skip_if_not_installed("torch")
    skip_if_not_installed("luz")
    models <- stack_test_models()
    stack <- sits_stack_model(a = models[["a"]], b = models[["b"]])
    cube <- stack_test_cube()
    n_a <- length(.ml_labels(models[["a"]]))
    for (path in c("cpu", "workers", "dataloader")) {
        dir <- file.path(tempdir(), paste0("stack_", path))
        unlink(dir, recursive = TRUE)
        dir.create(dir)
        probs_a <- classify_with(path, cube, models[["a"]], dir, "a")
        probs_b <- classify_with(path, cube, models[["b"]], dir, "b")
        probs_s <- classify_with(path, cube, stack, dir, "ab")
        expect_equal(unname(.cube_labels(probs_s)), .ml_labels(stack))
        s <- probs_values(probs_s)
        # within one unit of the INT2U scale, as the CPU/GPU test
        expect_lte(max(abs(s[, seq_len(n_a)] - probs_values(probs_a))), 1L)
        expect_lte(max(abs(s[, -seq_len(n_a)] - probs_values(probs_b))), 1L)
        unlink(dir, recursive = TRUE)
    }
})

test_that("sits_stack_split writes one probs cube per model", {
    skip_on_cran()
    skip_if_not_installed("torch")
    skip_if_not_installed("luz")
    models <- stack_test_models()
    stack <- sits_stack_model(a = models[["a"]], b = models[["b"]])
    cube <- stack_test_cube()
    dir <- file.path(tempdir(), "stack_split")
    unlink(dir, recursive = TRUE)
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    probs_a <- classify_with("cpu", cube, models[["a"]], dir, "a")
    probs_s <- classify_with("cpu", cube, stack, dir, "ab")
    split <- sits_stack_split(probs_s, output_dir = dir)
    expect_named(split, c("a", "b"))
    expect_s3_class(split[["a"]], "probs_cube")
    expect_equal(.cube_labels(split[["a"]]), .cube_labels(probs_a))
    expect_equal(.cube_labels(split[["b"]]), .ml_labels_code(models[["b"]]))
    n_a <- length(.ml_labels(models[["a"]]))
    s <- probs_values(probs_s)
    expect_equal(probs_values(split[["a"]]), s[, seq_len(n_a)], ignore_attr = TRUE)
    expect_equal(probs_values(split[["b"]]), s[, -seq_len(n_a)], ignore_attr = TRUE)
    # a second call keeps the files
    path <- .fi(split[["b"]])[["path"]]
    mtime <- file.mtime(path)
    Sys.sleep(1.1)
    sits_stack_split(probs_s, output_dir = dir)
    expect_equal(file.mtime(path), mtime)
})

test_that("sits_stack_split writes each model to its own dir", {
    skip_on_cran()
    skip_if_not_installed("torch")
    skip_if_not_installed("luz")
    models <- stack_test_models()
    stack <- sits_stack_model(a = models[["a"]], b = models[["b"]])
    cube <- stack_test_cube()
    root <- file.path(tempdir(), "stack_dirs")
    unlink(root, recursive = TRUE)
    dirs <- c(a = file.path(root, "a"), b = file.path(root, "b"))
    for (d in c(root, dirs)) dir.create(d)
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    probs_s <- classify_with("cpu", cube, stack, root, "ab")
    split <- sits_stack_split(probs_s, output_dir = dirs)
    expect_equal(dirname(.fi(split[["a"]])[["path"]]), unname(dirs[["a"]]))
    expect_equal(dirname(.fi(split[["b"]])[["path"]]), unname(dirs[["b"]]))
    expect_length(list.files(dirs[["a"]]), 1L)
    expect_length(list.files(dirs[["b"]]), 1L)
    # dirs not named by the models of the cube
    expect_error(sits_stack_split(probs_s, output_dir = c(a = root, c = root)))
})
