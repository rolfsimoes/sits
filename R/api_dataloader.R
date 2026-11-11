#' @title Torch dataset of the blocks of a set of tiles
#' @name .dataloader_blocks_dataset
#' @keywords internal
#' @noRd
#' @description
#' Each item is one block of one tile, read by \code{read_fn} in a reader
#' process of the dataloader. Reader processes are \code{callr} sessions
#' started with the \code{.libPaths()} of the main process, so they load the
#' same \code{sits} as the main process.
#'
#' @param chunks   Chunks of every tile, with a \code{tile_index} column.
#' @param tiles    List of the tiles of the chunks.
#' @param read_fn  Function \code{(tile, block)} that returns a matrix
#'                 \code{[nrows * ncols, nbands]}.
#' @return A torch dataset whose item \code{i} is a list with \code{input}, a
#'   float tensor \code{[nbands, nrows, ncols]}, and \code{index}, which is
#'   \code{i}.
.dataloader_blocks_dataset <- torch::dataset(
    "sits_blocks_dataset",
    initialize = function(chunks, tiles, read_fn) {
        self$chunks <- chunks
        self$tiles <- tiles
        self$read_fn <- read_fn
    },
    .getitem = function(i) {
        chunk <- self$chunks[i, ]
        tile <- self$tiles[[chunk[["tile_index"]]]]
        block <- .block(chunk)
        values <- self$read_fn(tile = tile, block = block)
        input <- torch::torch_tensor(t(values), dtype = torch::torch_float32())
        list(
            input = input$view(c(ncol(values), .nrows(block), .ncols(block))),
            index = i
        )
    },
    .length = function() {
        nrow(self$chunks)
    }
)

#' @title Start the iteration of a dataloader over the blocks of tiles
#' @name .dataloader_iter
#' @keywords internal
#' @noRd
#' @param chunks   Chunks of every tile, with a \code{tile_index} column.
#' @param tiles    List of the tiles of the chunks.
#' @param read_fn  Function \code{(tile, block)} that reads one block.
#' @param readers  Number of reader processes.
#' @return A dataloader iterator; the blocks come in the order of
#'   \code{chunks}, each with a batch dimension of one.
.dataloader_iter <- function(chunks, tiles, read_fn, readers) {
    dataset <- .dataloader_blocks_dataset(
        chunks = chunks,
        tiles = tiles,
        read_fn = read_fn
    )
    loader <- torch::dataloader(
        dataset = dataset,
        batch_size = 1L,
        num_workers = readers
    )
    torch::dataloader_make_iter(loader)
}

#' @title Next block of a dataloader iterator
#' @name .dataloader_next
#' @keywords internal
#' @noRd
#' @param iter  Result of \code{.dataloader_iter()}.
#' @return List with \code{input}, the tensor of the next block, and
#'   \code{index}, its row in the chunks, or \code{NULL} after the last block.
.dataloader_next <- function(iter) {
    batch <- torch::dataloader_next(iter, completed = NULL)
    if (is.null(batch)) {
        return(NULL)
    }
    list(input = batch[["input"]], index = batch[["index"]]$item())
}
