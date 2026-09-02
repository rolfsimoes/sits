#' @title Rasterize a vector layer onto one tile's grid
#' @keywords internal
#' @noRd
#' @param  tile        Subset of a data cube (defines the destination grid)
#' @param  vector      Path to a vector file readable by OGR
#' @param  field       Attribute field to burn. A literal field on `vector`,
#'                      or (when `sql` is given) a column `sql` projects
#' @param  sql         Optional OGR SQL SELECT statement selecting or
#'                      deriving the features and fields `field` is read from
#' @param  init        Fill value for pixels the vector does not cover
#' @param  band        Output band name
#' @param  labels      Output labels
#' @param  block       Image block to be processed
#' @param  multicores  Number of cores to run the function
#' @param  memsize     Maximum overall memory (in GB) to run the function
#' @param  output_dir  Directory where the image will be saved
#' @param  version     Version of result
#' @param  progress    Show progress bar?
#' @return Rasterized class_cube tile
.rasterize_tile <- function(tile,
                            vector,
                            field,
                            sql,
                            init,
                            band,
                            labels,
                            block,
                            multicores,
                            memsize,
                            output_dir,
                            version,
                            progress) {
    # Output file
    out_file <- .file_derived_name(
        tile = tile, band = band, version = version, output_dir = output_dir
    )
    # Resume feature
    if (all(.raster_is_valid(out_file, output_dir = output_dir))) {
        .check_recovery()
        class_tile <- .tile_derived_from_file(
            file = out_file,
            band = band,
            base_tile = tile,
            derived_class = "class_cube",
            labels = labels,
            update_bbox = FALSE
        )
        return(.tile_update_label(
            tile = class_tile, labels = labels,
            multicores = multicores, memsize = memsize
        ))
    }
    band_conf <- .conf_derived_band(derived_class = "class_cube", band = band)
    .check_that(init != .miss_value(band_conf))
    # Create chunks as jobs
    chunks <- .tile_chunks_create(tile = tile, overlap = 0L, block = block)
    # start parallel process
    block_files <- .jobs_map_parallel_chr(chunks, function(chunk) {
        # Get job block
        block <- .block(chunk)
        # Output file name
        block_file <- .file_block_name(
            pattern = .file_pattern(out_file),
            block = block,
            output_dir = output_dir
        )
        # Resume processing in case of failure
        if (all(.raster_is_valid(block_file))) {
            return(block_file)
        }
        # Template fill is the cube's missing value, not `init` (a real
        # class value), so uncovered pixels come back tagged as NoData.
        miss_value <- .miss_value(band_conf)
        tmp_file <- .file_block_name(
            pattern = .file_pattern(out_file, suffix = "_tmp"),
            block = block, output_dir = output_dir
        )
        .gdal_template_block(
            block = block, bbox = .bbox(chunk), file = tmp_file,
            nlayers = 1L, miss_value = miss_value,
            data_type = .data_type(band_conf)
        )
        # Burn the vector into the template
        .gdal_rasterize_into(
            file = tmp_file, vector = vector, field = field, sql = sql
        )
        # terra reads NoData back as NA, so uncovered pixels surface here
        # as NA. Replace them with `init` before writing the block.
        values <- .raster_get_values(.raster_open_rast(tmp_file))
        values[is.na(values)] <- init
        .raster_write_block(
            files = block_file, block = block, bbox = .bbox(chunk),
            values = values, data_type = .data_type(band_conf),
            missing_value = miss_value, crop_block = NULL
        )
        unlink(tmp_file)
        # Free memory
        gc()
        # Returned value
        block_file
    }, progress = progress)
    # Merge blocks into a new class_cube tile
    class_tile <- .tile_derived_merge_blocks(
        file = out_file,
        band = band,
        labels = labels,
        base_tile = tile,
        block_files = block_files,
        derived_class = "class_cube",
        multicores = .jobs_multicores(),
        update_bbox = FALSE
    )
    # Update tile labels
    .tile_update_label(
        tile = class_tile, labels = labels,
        multicores = multicores, memsize = memsize
    )
}
