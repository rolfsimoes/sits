# two tiles of 512 x 512 with 3 labels, file blocks of 128 x 128
smooth_dataloader_cube <- function(dir) {
    set.seed(11)
    labels <- c("1" = "A", "2" = "B", "3" = "C")
    purrr::walk(1:2, function(i) {
        xmin <- 500000 + (i - 1) * 5120
        r <- terra::rast(
            nrows = 512, ncols = 512, nlyrs = 3,
            xmin = xmin, xmax = xmin + 5120, ymin = 8000000, ymax = 8005120,
            crs = "EPSG:32723"
        )
        raw <- matrix(runif(512 * 512 * 3, 0.05, 1.0), ncol = 3)
        raw[] <- pmax(1, round(raw / rowSums(raw) * 10000))
        # scattered missing pixels and a patch smaller than the window
        patch <- as.vector(outer(200:202, 300:302, function(r, c) r * 512 + c))
        raw[c(patch, sample(512 * 512, 50)), ] <- 65535
        terra::values(r) <- raw
        terra::writeRaster(
            r,
            file.path(dir, sprintf(
                "SYN_X_00000%d_2020-01-01_2020-12-31_probs_v1.tif", i
            )),
            datatype = "INT2U", NAflag = 65535,
            gdal = c("TILED=YES", "BLOCKXSIZE=128", "BLOCKYSIZE=128")
        )
    })
    sits_cube(
        source = "BDC",
        collection = "MOD13Q1-6.1",
        data_dir = dir,
        bands = "probs",
        labels = labels,
        progress = FALSE
    )
}
test_that("sits_smooth tile flow with a dataloader on the CPU device", {
    skip_on_cran()
    skip_if_not(.torch_is_functional())
    old_env <- Sys.getenv(
        c("SITS_SMOOTH_FORCE_DATALOADER", "SITS_SMOOTH_FORCE_CPP"),
        unset = NA_character_
    )
    on.exit(
        purrr::iwalk(old_env, function(value, name) {
            if (is.na(value)) Sys.unsetenv(name) else
                do.call(Sys.setenv, stats::setNames(list(value), name))
        }),
        add = TRUE
    )
    Sys.unsetenv(c("SITS_SMOOTH_FORCE_DATALOADER", "SITS_SMOOTH_FORCE_CPP"))
    base_dir <- file.path(tempdir(), "smooth-dataloader")
    unlink(base_dir, recursive = TRUE)
    on.exit(unlink(base_dir, recursive = TRUE), add = TRUE)
    dirs <- purrr::map_chr(
        c(probs = "probs", ref = "ref", ref_mask = "ref_mask", flow = "flow",
          flow_mask = "flow_mask", log = "log"),
        function(d) {
            path <- file.path(base_dir, d)
            dir.create(path, recursive = TRUE)
            path
        }
    )
    probs_cube <- smooth_dataloader_cube(dirs[["probs"]])
    expect_equal(nrow(probs_cube), 2L)
    exclusion_mask <- sf::st_as_sfc(
        sf::st_bbox(
            c(xmin = 501000, ymin = 8001000, xmax = 502000, ymax = 8002000),
            crs = sf::st_crs("EPSG:32723")
        )
    )
    smooth <- function(output_dir, memsize, multicores, mask = NULL) {
        sits_smooth(
            cube = probs_cube,
            exclusion_mask = mask,
            memsize = memsize,
            multicores = multicores,
            output_dir = output_dir,
            progress = FALSE
        )
    }
    values_of <- function(cube) {
        purrr::map(seq_len(nrow(cube)), function(i) {
            .raster_get_values(
                .raster_open_rast(.fi_paths(.fi(.tile(cube[i, ]))))
            )
        })
    }
    expect_same_values <- function(flow_cube, ref_cube) {
        purrr::walk2(values_of(flow_cube), values_of(ref_cube), function(a, b) {
            expect_equal(is.na(a), is.na(b))
            expect_lte(max(abs(a - b), na.rm = TRUE), 1)
        })
    }
    # reference: the CPU torch path, one block per tile
    ref <- smooth(dirs[["ref"]], memsize = 4L, multicores = 1L)
    ref_mask <- smooth(
        dirs[["ref_mask"]], memsize = 4L, multicores = 1L, mask = exclusion_mask
    )
    Sys.setenv(SITS_SMOOTH_FORCE_DATALOADER = "TRUE")
    # the flow with blocks of 256 x 128: 2 readers and the 2 cluster workers
    # of sits_parallel, whose workers log their events
    sits_parallel(workers = 2L, log = TRUE, output_dir = dirs[["log"]])
    worker_pids <- unlist(parallel::clusterCall(sits_parallel(), Sys.getpid))
    flow <- smooth(dirs[["flow"]], memsize = 1L, multicores = 4L)
    sits_parallel(workers = 0L)
    expect_same_values(flow, ref)
    events <- purrr::map_dfr(
        list.files(file.path(dirs[["log"]], ".sits"), pattern = "\\.log$",
                   full.names = TRUE),
        utils::read.csv, strip.white = TRUE
    )
    writes <- events[events[["event"]] == "end_smooth_write", ]
    merges <- events[events[["event"]] == "end_smooth_merge", ]
    chunks <- .tile_chunks_create(
        tile = .tile(probs_cube), overlap = 4L,
        block = c(ncols = 256L, nrows = 128L)
    )
    expect_equal(nrow(writes), 2L * nrow(chunks))
    expect_equal(length(unique(writes[["value"]])), nrow(writes))
    expect_true(all(writes[["pid"]] %in% worker_pids))
    expect_equal(nrow(merges), 2L)
    expect_true(all(merges[["pid"]] %in% worker_pids))
    expect_false(Sys.getpid() %in% c(writes[["pid"]], merges[["pid"]]))
    expect_length(
        list.files(file.path(dirs[["flow"]], ".sits"), pattern = "_block_"), 0L
    )
    # a second call finds the valid tiles and reads no block
    again <- testthat::with_mocked_bindings(
        suppressMessages(smooth(dirs[["flow"]], memsize = 1L, multicores = 4L)),
        .dataloader_iter = function(...) stop("blocks read again")
    )
    paths_of <- function(cube) {
        purrr::map(seq_len(nrow(cube)), function(i) .fi_paths(.fi(cube[i, ])))
    }
    expect_equal(paths_of(again), paths_of(flow))
    # the flow with its own cluster and an exclusion mask
    flow_mask <- smooth(
        dirs[["flow_mask"]], memsize = 1L, multicores = 4L, mask = exclusion_mask
    )
    expect_null(sits_parallel())
    expect_same_values(flow_mask, ref_mask)
    expect_length(
        list.files(file.path(dirs[["flow_mask"]], ".sits"), pattern = "\\.tif$"),
        0L
    )
})
test_that("sits_smooth dataloader flow survives or reports a failed send", {
    skip_on_cran()
    skip_if_not(.torch_is_functional())
    skip_if_not_installed("callr")
    old_env <- Sys.getenv(
        c("SITS_SMOOTH_FORCE_DATALOADER", "SITS_SMOOTH_FORCE_CPP"),
        unset = NA_character_
    )
    on.exit(
        purrr::iwalk(old_env, function(value, name) {
            if (is.na(value)) Sys.unsetenv(name) else
                do.call(Sys.setenv, stats::setNames(list(value), name))
        }),
        add = TRUE
    )
    Sys.unsetenv(c("SITS_SMOOTH_FORCE_DATALOADER", "SITS_SMOOTH_FORCE_CPP"))
    base_dir <- file.path(tempdir(), "smooth-dataloader-send")
    unlink(base_dir, recursive = TRUE)
    on.exit(unlink(base_dir, recursive = TRUE), add = TRUE)
    dirs <- purrr::map_chr(
        c(probs = "probs", ref = "ref", once = "once", always = "always"),
        function(d) {
            path <- file.path(base_dir, d)
            dir.create(path, recursive = TRUE)
            path
        }
    )
    probs_cube <- smooth_dataloader_cube(dirs[["probs"]])
    ref <- sits_smooth(
        probs_cube, memsize = 4L, multicores = 1L,
        output_dir = dirs[["ref"]], progress = FALSE
    )
    # a send broken by a signal, as Rprof gives, fails in the main process;
    # a separate process bounds the wait of a flow that does not end
    flow <- function(output_dir, fail) {
        callr::r(
            function(cube, output_dir, fail) {
                Sys.setenv(SITS_SMOOTH_FORCE_DATALOADER = "TRUE")
                submit <- sits:::.parallel_submit
                calls <- 0L
                utils::assignInNamespace(".parallel_submit", function(...) {
                    calls <<- calls + 1L
                    if (calls %in% fail) stop("error writing to connection")
                    submit(...)
                }, ns = "sits")
                tryCatch({
                    cube <- sits::sits_smooth(
                        cube, memsize = 1L, multicores = 4L,
                        output_dir = output_dir, progress = FALSE
                    )
                    purrr::map_chr(seq_len(nrow(cube)), function(i) {
                        sits:::.fi_paths(sits:::.fi(cube[i, ]))
                    })
                }, error = function(e) paste("error:", conditionMessage(e)))
            },
            args = list(cube = probs_cube, output_dir = output_dir,
                        fail = fail),
            timeout = 180
        )
    }
    # one failed send: the block is sent again and the output is complete
    once <- flow(dirs[["once"]], fail = 2L)
    ref_paths <- purrr::map_chr(seq_len(nrow(ref)), function(i) {
        .fi_paths(.fi(ref[i, ]))
    })
    expect_length(once, 2L)
    purrr::walk2(once, ref_paths, function(a, b) {
        a <- .raster_get_values(.raster_open_rast(a))
        b <- .raster_get_values(.raster_open_rast(b))
        expect_equal(is.na(a), is.na(b))
        expect_lte(max(abs(a - b), na.rm = TRUE), 1)
    })
    # every send fails: an error, not a wait
    always <- flow(dirs[["always"]], fail = seq_len(1000L))
    expect_length(always, 1L)
    expect_true(startsWith(always, "error:"))
})
