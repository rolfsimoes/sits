test_that("onLoad", {
    .onLoad()
    expect(TRUE, "onLoad failed to load")
})
test_that("onLoad makes GDAL take the folder of a raster as empty", {
    # the no-data value lives only in a sidecar file, which GDAL reads only
    # when it looks at the folder: the NA count shows the option in effect
    dir <- tempfile("readdir_")
    dir.create(dir)
    file <- file.path(dir, "pam.tif")
    rast <- terra::rast(
        nrows = 2L, ncols = 2L, xmin = 0.0, xmax = 2.0, ymin = 0.0,
        ymax = 2.0, crs = "EPSG:32723", vals = 1:4
    )
    terra::writeRaster(rast, file, datatype = "INT2S")
    writeLines(
        c(
            "<PAMDataset><PAMRasterBand band=\"1\">",
            "<NoDataValue>1</NoDataValue></PAMRasterBand></PAMDataset>"
        ),
        paste0(file, ".aux.xml")
    )
    na_count <- function(file) sum(is.na(terra::values(terra::rast(file))))
    # global environment, so a worker does not receive the test frame
    environment(na_count) <- globalenv()
    # NA count in the main process (terra and sf) and in two sits workers
    na_counts <- function() {
        info <- sf::gdal_utils("info", file, quiet = TRUE)
        started <- .parallel_start(workers = 2L)
        on.exit(.parallel_stop(started), add = TRUE)
        workers <- parallel::clusterCall(sits_env[["cluster"]], na_count, file)
        c(
            terra = na_count(file),
            sf = as.integer(grepl("NoData Value=1", info, fixed = TRUE)),
            workers = unlist(workers)
        )
    }
    old <- Sys.getenv("GDAL_DISABLE_READDIR_ON_OPEN", unset = NA)
    on.exit(
        if (is.na(old)) {
            Sys.unsetenv("GDAL_DISABLE_READDIR_ON_OPEN")
        } else {
            Sys.setenv(GDAL_DISABLE_READDIR_ON_OPEN = old)
        },
        add = TRUE
    )
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    read <- c(terra = 1L, sf = 1L, workers1 = 1L, workers2 = 1L)
    ignored <- c(terra = 0L, sf = 0L, workers1 = 0L, workers2 = 0L)
    # stand-in check: without the option every reader finds the sidecar
    Sys.unsetenv("GDAL_DISABLE_READDIR_ON_OPEN")
    expect_equal(na_counts(), read)
    # after loading sits, GDAL reports the option and ignores the sidecar
    .onLoad()
    expect_equal(
        unname(terra::getGDALconfig("GDAL_DISABLE_READDIR_ON_OPEN")),
        "EMPTY_DIR"
    )
    expect_equal(na_counts(), ignored)
    # a value set by the user after loading sits wins at the next open
    Sys.setenv(GDAL_DISABLE_READDIR_ON_OPEN = "NO")
    expect_equal(na_counts(), read)
    # a value set by the user before loading sits is kept
    .onLoad()
    expect_equal(Sys.getenv("GDAL_DISABLE_READDIR_ON_OPEN"), "NO")
    expect_equal(na_counts(), read)
})
