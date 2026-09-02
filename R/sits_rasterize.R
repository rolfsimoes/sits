#' @title Rasterize a vector layer onto a data cube's grid
#' @name sits_rasterize
#'
#' @author Rolf Simoes, \email{rolfsimoes@@gmail.com}
#'
#' @description
#' Burns a vector layer (polygons or lines, any format OGR reads) onto the
#' grid of an existing raster data cube, tile by tile and block by block,
#' producing a new \code{class_cube}. Useful to bring a vector product with
#' no matching raster release into the cube framework, on the cube's own
#' grid.
#'
#' @details
#' Pixels outside the vector's coverage are set to \code{init}.
#'
#' @param cube        A \code{raster_cube}. Only its geometry (tiles, block
#'                     layout, CRS, extent) is used, not its pixel values.
#' @param vector       Path to a vector file readable by GDAL/OGR.
#' @param field        Attribute field to burn. A literal field on
#'                      \code{vector}, or (when \code{sql} is given) a
#'                      column \code{sql} projects.
#' @param sql          Optional OGR SQL \code{SELECT} statement selecting
#'                      or deriving the features and fields \code{field}
#'                      is read from.
#' @param init         Fill value for pixels the vector does not cover
#'                      (default 0).
#' @param band         Output band name (default "class").
#' @param labels       Named vector mapping the burned integer codes to
#'                      class names (as in \code{\link[sits]{sits_reclassify}}).
#' @param memsize      Memory available for processing (in GB).
#' @param multicores   Number of cores to be used for processing.
#' @param output_dir   Directory where files will be saved.
#' @param version      Version of resulting image (character).
#' @param progress     Show progress bar?
#'
#' @return A \code{class_cube} with the vector burned onto \code{cube}'s grid.
#'
#' @examples
#' if (sits_run_examples()) {
#'     data_dir <- system.file("extdata/raster/mod13q1", package = "sits")
#'     cube <- sits_cube(
#'         source = "BDC", collection = "MOD13Q1-6.1", data_dir = data_dir
#'     )
#'     shp_file <- system.file(
#'         "extdata/shapefiles/mato_grosso/mt.shp", package = "sits"
#'     )
#'     class_cube <- sits_rasterize(
#'         cube = cube, vector = shp_file, field = "NM_ESTADO",
#'         labels = c("1" = "Mato Grosso"), output_dir = tempdir()
#'     )
#' }
#' @export
sits_rasterize <- function(cube,
                           vector,
                           ...,
                           field,
                           sql = NULL,
                           init = 0L,
                           band = "class",
                           labels,
                           memsize = 4L,
                           multicores = 2L,
                           output_dir,
                           version = "v1",
                           progress = TRUE) {
    .check_set_caller("sits_rasterize")
    # Force lazily-evaluated arguments now: parallel workers are separate
    # processes and cannot resolve a promise against the caller's environment.
    force(vector); force(field); force(sql); force(init)
    force(band); force(labels)
    # Preconditions
    .check_is_raster_cube(cube)
    .check_that(.has(field))
    .check_int_parameter(memsize, min = 1L, max = 16384L)
    .check_int_parameter(multicores, min = 1L, max = 2048L)
    .check_output_dir(output_dir)
    version <- .message_version(version)
    progress <- .message_progress(progress)

    # The following functions define optimal parameters for parallel processing
    block <- .raster_file_blocksize(.raster_open_rast(.tile_path(cube)))
    job_block_memsize <- .jobs_block_memsize(
        block_size = .block_size(block = block, overlap = 0L),
        npaths = 1L, nbytes = 4L,
        proc_bloat = .conf("processing_bloat")
    )
    multicores <- .jobs_max_multicores(
        job_block_memsize = job_block_memsize,
        memsize = memsize, multicores = multicores
    )
    block <- .jobs_optimal_block(
        job_block_memsize = job_block_memsize, block = block,
        image_size = .tile_size(.tile(cube)), memsize = memsize,
        multicores = multicores
    )
    # Prepare parallel processing
    started <- .parallel_start(workers = multicores)
    on.exit(.parallel_stop(started), add = TRUE)
    # Process each tile sequentially (blocks within a tile run in parallel)
    class_cube <- .cube_foreach_tile(cube, function(tile) {
        .rasterize_tile(
            tile = tile, vector = vector, field = field, sql = sql,
            init = init, band = band, labels = labels,
            block = block, multicores = multicores, memsize = memsize,
            output_dir = output_dir, version = version, progress = progress
        )
    })
    class(class_cube) <- c("class_cube", class(class_cube))
    class_cube
}
