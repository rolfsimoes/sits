# The dataloader path (GPU) reads the blocks of every tile with one
# dataloader. The tiles are cut from the mod13q1 sample cube: tile A is its
# west 127 columns, tile B the next 127, so a block written to the wrong
# tile has another bbox and other values. The tiles have the same size
# because .check_block_size() takes one size per cube. The reference is
# the CPU path, which classifies each tile on its own without a dataloader.
gpu_tiles_cube <- function(root) {
    cube <- sits_cube(
        source = "BDC",
        collection = "MOD13Q1-6.1",
        data_dir = system.file("extdata/raster/mod13q1", package = "sits"),
        progress = FALSE
    )
    split_x <- cube[["xmin"]] + 127L * .tile_xres(cube)
    crop <- function(xmin, xmax, name) {
        dir <- file.path(root, paste0("input_", name))
        dir.create(dir, recursive = TRUE)
        roi <- sf::st_as_sfc(sf::st_bbox(
            c(xmin = xmin, xmax = xmax,
              ymin = cube[["ymin"]], ymax = cube[["ymax"]]),
            crs = sf::st_crs(cube[["crs"]])
        ))
        # gdalwarp warns that it moves -3000 to -3001, the same in both
        tile <- suppressWarnings(sits_cube_copy(
            cube, roi = roi, multicores = 1L, output_dir = dir,
            progress = FALSE
        ))
        tile[["tile"]] <- name
        tile
    }
    dplyr::bind_rows(
        crop(cube[["xmin"]], split_x, "A"),
        crop(split_x, split_x + 127L * .tile_xres(cube), "B")
    )
}

gpu_tiles_model <- function(seed) {
    set.seed(seed)
    torch::torch_manual_seed(seed)
    sits_train(
        samples_modis_ndvi,
        sits_mlp(
            layers = c(32L, 32L), dropout_rates = c(0.2, 0.2), epochs = 2L,
            verbose = FALSE
        )
    )
}

# Classify on the dataloader path (TRUE) or on the CPU path (FALSE), in
# the R session, with blocks of 64 x 64 pixels (several per tile).
# Returns the probs and how many dataloaders were created.
classify_gpu_tiles <- function(gpu, cube, ml_model, output_dir,
                               version = "v1", verbose = FALSE) {
    n_loaders <- 0L
    trace(
        "dataloader",
        tracer = function() n_loaders <<- n_loaders + 1L,
        where = asNamespace("torch"),
        print = FALSE
    )
    on.exit(untrace("dataloader", where = asNamespace("torch")), add = TRUE)
    probs <- testthat::with_mocked_bindings(
        sits_classify(
            data = cube, ml_model = ml_model, memsize = 4L, multicores = 1L,
            block_size = c(ncols = 64L, nrows = 64L),
            output_dir = output_dir, version = version, verbose = verbose,
            progress = FALSE
        ),
        .torch_gpu_classification = function() gpu
    )
    list(probs = probs, n_loaders = n_loaders)
}

tile_values <- function(probs, name) {
    .raster_open_rast(.fi(probs[probs[["tile"]] == name, ])[["path"]])[]
}

test_that("one dataloader classifies all tiles, each as on its own", {
    skip_on_cran()
    skip_if_not_installed("torch")
    skip_if_not_installed("luz")
    root <- file.path(tempdir(), "gpu_tiles")
    unlink(root, recursive = TRUE)
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    cube <- gpu_tiles_cube(root)
    model <- gpu_tiles_model(2971L)
    dirs <- file.path(root, c("ref", "gpu"))
    for (d in dirs) dir.create(d)
    # the CPU path counts luz's dataloaders of each block, not checked here
    ref <- classify_gpu_tiles(FALSE, cube, model, dirs[[1L]])
    # record, when each tile reports its end, which probs files exist
    probs_file <- function(name) {
        .file_derived_name(
            tile = cube[cube[["tile"]] == name, ], band = "probs",
            version = "v1", output_dir = dirs[[2L]]
        )
    }
    at_end <- list()
    gpu <- withCallingHandlers(
        classify_gpu_tiles(TRUE, cube, model, dirs[[2L]], verbose = TRUE),
        message = function(m) {
            tile <- regmatches(
                conditionMessage(m),
                regexec("^Tile '(.*)' finished", conditionMessage(m))
            )[[1L]]
            if (length(tile) == 2L) {
                at_end[[tile[[2L]]]] <<- c(
                    A = file.exists(probs_file("A")),
                    B = file.exists(probs_file("B"))
                )
            }
            invokeRestart("muffleMessage")
        }
    )
    expect_equal(gpu[["n_loaders"]], 1L)
    expect_equal(gpu[["probs"]][["tile"]], c("A", "B"))
    # tile A is merged before tile B ends, tile B when it ends
    expect_equal(at_end[["A"]], c(A = TRUE, B = FALSE))
    expect_equal(at_end[["B"]], c(A = TRUE, B = TRUE))
    for (name in c("A", "B")) {
        # within one unit of the INT2U scale, as the CPU/GPU test
        expect_lte(
            max(abs(tile_values(gpu[["probs"]], name) -
                    tile_values(ref[["probs"]], name)), na.rm = TRUE),
            1L
        )
        expect_equal(
            is.na(tile_values(gpu[["probs"]], name)),
            is.na(tile_values(ref[["probs"]], name))
        )
    }
    # no block file is left
    expect_length(
        list.files(dirs[[2L]], pattern = "block", recursive = TRUE,
                   all.files = TRUE),
        0L
    )
})

test_that("the dataloader path skips done tiles and recovered blocks", {
    skip_on_cran()
    skip_if_not_installed("torch")
    skip_if_not_installed("luz")
    root <- file.path(tempdir(), "gpu_tiles_recovery")
    unlink(root, recursive = TRUE)
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    cube <- gpu_tiles_cube(root)
    model <- gpu_tiles_model(2971L)
    dir <- file.path(root, "gpu")
    dir.create(dir)
    first <- classify_gpu_tiles(TRUE, cube, model, dir)
    path_a <- .fi(first[["probs"]][1L, ])[["path"]]
    path_b <- .fi(first[["probs"]][2L, ])[["path"]]
    values_b <- tile_values(first[["probs"]], "B")
    mtime_a <- file.mtime(path_a)
    Sys.sleep(1.1)
    # a done tile is not classified again
    unlink(path_b)
    again <- classify_gpu_tiles(TRUE, cube, model, dir)
    expect_equal(again[["n_loaders"]], 1L)
    expect_equal(file.mtime(path_a), mtime_a)
    expect_true(file.exists(path_b))
    expect_lte(max(abs(tile_values(again[["probs"]], "B") - values_b),
                   na.rm = TRUE), 1L)
    # tile B fails after all its blocks are written: they are kept
    unlink(path_b)
    merge <- .classify_merge
    expect_error(testthat::with_mocked_bindings(
        classify_gpu_tiles(TRUE, cube, model, dir),
        .classify_merge = function(tile, ...) {
            if (tile[["tile"]] == "B") stop("merge of tile B failed")
            merge(tile, ...)
        }
    ), "merge of tile B failed")
    expect_false(file.exists(path_b))
    expect_gt(length(list.files(dir, pattern = "block", recursive = TRUE,
                                all.files = TRUE)), 0L)
    # every block of tile B is recovered, so no dataloader is created
    last <- classify_gpu_tiles(TRUE, cube, model, dir)
    expect_equal(last[["n_loaders"]], 0L)
    expect_equal(file.mtime(path_a), mtime_a)
    expect_lte(max(abs(tile_values(last[["probs"]], "B") - values_b),
                   na.rm = TRUE), 1L)
})

test_that("a stack classifies all tiles with one dataloader", {
    skip_on_cran()
    skip_if_not_installed("torch")
    skip_if_not_installed("luz")
    root <- file.path(tempdir(), "gpu_tiles_stack")
    unlink(root, recursive = TRUE)
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    cube <- gpu_tiles_cube(root)
    models <- list(a = gpu_tiles_model(101L), b = gpu_tiles_model(202L))
    stack <- sits_stack_model(a = models[["a"]], b = models[["b"]])
    dirs <- c(a = file.path(root, "a"), b = file.path(root, "b"),
              ref = file.path(root, "ref"))
    for (d in dirs) dir.create(d)
    gpu <- classify_gpu_tiles(TRUE, cube, stack, dirs[c("a", "b")])
    expect_equal(gpu[["n_loaders"]], 1L)
    expect_named(gpu[["probs"]], c("a", "b"))
    for (name in c("a", "b")) {
        ref <- classify_gpu_tiles(FALSE, cube, models[[name]], dirs[["ref"]],
                                  version = name)
        probs <- gpu[["probs"]][[name]]
        expect_equal(probs[["tile"]], c("A", "B"))
        expect_equal(unique(dirname(.fi(probs)[["path"]])),
                     unname(dirs[[name]]))
        for (tile in c("A", "B")) {
            expect_lte(max(abs(tile_values(probs, tile) -
                               tile_values(ref[["probs"]], tile)),
                           na.rm = TRUE), 1L)
        }
    }
})
