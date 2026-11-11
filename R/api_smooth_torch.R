#' @title Check if the torch smoother is available
#' @name .torch_smooth_available
#' @keywords internal
#' @noRd
#' @description Use torch when it is functional unless
#'   SITS_SMOOTH_FORCE_CPP is set to TRUE.
#'
#' @return A logical value
.torch_smooth_available <- function() {
    force_cpp <- Sys.getenv("SITS_SMOOTH_FORCE_CPP", unset = "FALSE")
    if (toupper(force_cpp) == "TRUE") {
        return(FALSE)
    }
    .torch_is_functional()
}

#' @title Check if the torch smoother runs in the dataloader flow
#' @name .torch_smooth_dataloader
#' @keywords internal
#' @noRd
#' @description The flow of \code{.smooth_dataloader()} runs when torch
#'   smooths on a GPU, or on the CPU device when SITS_SMOOTH_FORCE_DATALOADER
#'   is set to TRUE.
#'
#' @param use_torch  Result of \code{.torch_smooth_available()}.
#' @return A logical value
.torch_smooth_dataloader <- function(use_torch) {
    if (!use_torch) {
        return(FALSE)
    }
    force_dataloader <- Sys.getenv(
        "SITS_SMOOTH_FORCE_DATALOADER", unset = "FALSE"
    )
    toupper(force_dataloader) == "TRUE" || .torch_gpu_available()
}

#' @title Device of the torch smoother
#' @name .torch_smooth_device
#' @keywords internal
#' @noRd
#' @return "cuda" or "mps" when a GPU is available, otherwise "cpu".
.torch_smooth_device <- function() {
    if (!.torch_gpu_available()) {
        return("cpu")
    }
    if (torch::cuda_is_available()) "cuda" else "mps"
}

#' @title Estimate memory used by torch smoothing windows
#' @name .torch_smooth_block_memsize
#' @keywords internal
#' @noRd
#' @description Estimate the memory used to materialize the neighborhood
#'   tensors for one raster block. Bands are processed sequentially and their
#'   intermediate tensors are released between iterations, so the dominant
#'   window allocation does not grow with the number of bands.
#'
#' @param block_size Number of pixels in the block, including overlap.
#' @param window_size Side length of the square neighborhood.
#'
#' @return Estimated memory in GB.
.torch_smooth_block_memsize <- function(block_size, window_size) {
    .jobs_block_memsize(
        block_size = block_size,
        npaths = window_size^2,
        nbytes = 4L,
        proc_bloat = .conf("processing_bloat_smooth_torch")
    )
}

#' @title TorchScript source of the Bayesian smoother
#' @name .torch_smooth_script
#' @keywords internal
#' @noRd
#'
#' @description
#' TorchScript function \code{smooth(x, smoothness)} with the calculation of
#' \code{bayes_smoother_fraction()}. Its tensors are freed by reference count
#' in C++, so a block leaves no garbage for the R collector. \code{topk}
#' takes the highest values, since at most \code{top} values are selected.
#' The variance has the denominator of R \code{var()}.
#'
#' @param  window_size     Side length of the square neighbourhood.
#' @param  neigh_fraction  Fraction of neighbourhood values (highest logits)
#'                         used to estimate the prior.
#' @return A character string with the TorchScript source.
.torch_smooth_script <- function(window_size, neigh_fraction) {
    top <- ceiling(neigh_fraction * window_size * window_size)
    sprintf("
def smooth(x: Tensor, smoothness: Tensor):
    window = %1$d
    leg = window // 2
    fraction = float(%2$.17g)
    top = %3$d
    win_sq = window * window
    x_pad = torch.reflection_pad2d(x, [leg, leg, leg, leg])
    pos = torch.arange(1, top + 1).float().to(x.device).view(1, top)
    results: List[Tensor] = []
    for b in range(x.size(1)):
        wins = x_pad[0, b].unfold(0, window, 1).unfold(1, window, 1)
        wins = wins.reshape(-1, win_sq)
        nan_mask = torch.isnan(wins)
        valid = torch.logical_not(nan_mask).sum(dim=-1).float()
        if fraction == 1.0:
            selected = wins.masked_fill(nan_mask, 0.0)
            sel_mask = torch.logical_not(nan_mask)
            n = valid
        else:
            ordered = torch.topk(
                wins.masked_fill(nan_mask, float('-inf')), top, dim=-1,
                largest=True, sorted=True
            )[0]
            n = torch.ceil(fraction * valid).clamp(min=1.0)
            sel_mask = pos.le(n.unsqueeze(-1))
            selected = ordered.masked_fill(torch.logical_not(sel_mask), 0.0)
        not_sel = torch.logical_not(sel_mask)
        m0 = selected.sum(dim=-1) / n
        diff_sq = (selected - m0.unsqueeze(-1)).pow(2).masked_fill(not_sel, 0.0)
        s0 = diff_sq.sum(dim=-1) / (n - 1.0)
        s0 = s0.masked_fill(torch.logical_not(n.gt(1.0)), 0.0)
        x0 = x[0, b].reshape(-1)
        w = s0 / (s0 + smoothness[b])
        bayes = w * x0 + (1.0 - w) * m0
        results.append(torch.where(torch.isnan(x0) | s0.lt(1e-4), m0, bayes))
    return torch.stack(results, dim=1)
", as.integer(window_size), neigh_fraction, as.integer(top))
}

#' @title Compile the TorchScript Bayesian smoother
#' @name .torch_smooth_kernel
#' @keywords internal
#' @noRd
#'
#' @description
#' Compiles the source of \code{.torch_smooth_script()}. The result keeps the
#' source and the process id, since a compiled unit sent to a worker process
#' is an invalid pointer there (see \code{.torch_smooth_kernel_unit()}).
#' A torch without \code{reflection_pad2d} or \code{topk} fails here.
#'
#' @param  window_size     Side length of the square neighbourhood.
#' @param  neigh_fraction  Fraction of neighbourhood values (highest logits)
#'                         used to estimate the prior.
#' @return A list with the source, the compiled unit and the process id.
.torch_smooth_kernel <- function(window_size, neigh_fraction) {
    source <- .torch_smooth_script(window_size, neigh_fraction)
    list(
        source = source,
        unit = torch::jit_compile(source),
        pid = Sys.getpid()
    )
}

#' @title Compiled unit of the Bayesian smoother in this process
#' @name .torch_smooth_kernel_unit
#' @keywords internal
#' @noRd
#'
#' @description
#' Returns the compiled unit of \code{kernel} in the process that compiled it.
#' A worker process compiles the source once and keeps the unit in
#' \code{sits_env} for the next blocks of the same source.
#'
#' @param  kernel  Result of \code{.torch_smooth_kernel()}.
#' @return A TorchScript compilation unit.
.torch_smooth_kernel_unit <- function(kernel) {
    if (identical(kernel[["pid"]], Sys.getpid())) {
        return(kernel[["unit"]])
    }
    cached <- sits_env[["smooth_torch_kernel"]]
    if (!identical(cached[["source"]], kernel[["source"]]) ||
        !identical(cached[["pid"]], Sys.getpid())) {
        cached <- list(
            source = kernel[["source"]],
            unit = torch::jit_compile(kernel[["source"]]),
            pid = Sys.getpid()
        )
        sits_env[["smooth_torch_kernel"]] <- cached
    }
    cached[["unit"]]
}

#' @title Torch-backed Bayesian smoother for probability cubes
#' @name .torch_smooth_bayes_fraction
#' @keywords internal
#' @noRd
#' @author Alexandre Assuncao, \email{alexcarssuncao@@gmail.com}
#'
#' @description
#' Drop-in replacement for the C++ \code{bayes_smoother_fraction()} function.
#' This is the default raster smoothing backend whenever torch and its native
#' dependencies are functional. It runs the compiled kernel of
#' \code{.torch_smooth_kernel()} on CPU, or on a GPU when CUDA or MPS is
#' available.
#'
#' Padding note: PyTorch \code{"reflect"} mode mirrors from the element
#' adjacent to the boundary (boundary pixel is not repeated), whereas the
#' C++ implementation repeats the boundary pixel once.  The mismatch only
#' affects the outermost \code{leg = floor(window_size / 2)} rows/columns of
#' each processed block; those pixels are discarded by the overlap/crop_block
#' mechanism in \code{.smooth_tile()}, so the final output is unaffected for
#' all interior blocks.
#'
#' @param  logits          Numeric matrix \code{[nrows*ncols, nbands]} of
#'                         logit-transformed probabilities.
#' @param  nrows           Number of rows in the block.
#' @param  ncols           Number of columns in the block.
#' @param  smoothness      Numeric vector of length \code{nbands}: the prior
#'                         variance parameter per class.
#' @param  kernel          Result of \code{.torch_smooth_kernel()}, which
#'                         holds the window size and the neighbour fraction.
#' @return Numeric matrix \code{[nrows*ncols, nbands]} of smoothed logits.
#'
.torch_smooth_bayes_fraction <- function(logits,
                                         nrows,
                                         ncols,
                                         smoothness,
                                         kernel) {
    # Select device
    device <- .torch_smooth_device()
    # Number of bands in data
    nbands <- ncol(logits)
    # Build input tensor [1, nbands, nrows, ncols]
    x <- torch::torch_tensor(
        t(logits),
        dtype  = torch::torch_float32(),
        device = device
    )$view(c(nbands, nrows, ncols))$unsqueeze(1L)
    smoothness <- torch::torch_tensor(
        smoothness,
        dtype  = torch::torch_float32(),
        device = device
    )
    unit <- .torch_smooth_kernel_unit(kernel)
    # Result has dim [npix, nbands]
    result <- unit$smooth(x, smoothness)
    as.matrix(
        result$cpu()$to(dtype = torch::torch_float64())
    )
}

#' @title Smooth one block of raw probabilities on the device
#' @name .torch_smooth_block
#' @keywords internal
#' @noRd
#'
#' @description
#' The steps of the smoother of \code{.smooth_fn_bayes_torch()} and of the
#' writer of \code{.smooth_tile()}, in the same order, on the device: scale of
#' the input, logit, compiled kernel, inverse logit, scale of the output.
#' Only the result comes back to the host, as integers.
#'
#' @param  input       Float tensor \code{[1, nbands, nrows, ncols]} of raw
#'                     values of the input band, NaN where missing.
#' @param  kernel      Result of \code{.torch_smooth_kernel()}.
#' @param  smoothness  Numeric vector of length \code{nbands}.
#' @param  device      Device of the calculation.
#' @param  in_conf     Band configuration of the input.
#' @param  out_conf    Band configuration of the output.
#' @return Integer matrix \code{[nrows * ncols, nbands]} with the values to be
#'   written, NA where the smoother gives NaN.
.torch_smooth_block <- function(input, kernel, smoothness, device,
                                in_conf, out_conf) {
    float64 <- torch::torch_float64()
    x <- input$to(device = device, dtype = float64)
    # the operations of .tile_scale(), so the values are those of the reader
    in_scale <- .scale(in_conf)
    if (.has(in_scale) && in_scale != 1.0) {
        x <- x * in_scale
    }
    in_offset <- .offset(in_conf)
    if (.has(in_offset) && in_offset != 0.0) {
        x <- x + in_offset
    }
    x <- x$masked_fill(x$eq(1.0), .conf("smooth_prob_max"))
    x <- x$masked_fill(x$eq(0.0), .conf("smooth_prob_min"))
    # logit in double precision, as R does before the kernel
    logits <- torch::torch_log(x / (x$sum(dim = 2L, keepdim = TRUE) - x))
    smoothness <- torch::torch_tensor(
        smoothness,
        dtype  = torch::torch_float32(),
        device = device
    )
    unit <- .torch_smooth_kernel_unit(kernel)
    values <- unit$smooth(logits$to(dtype = torch::torch_float32()), smoothness)
    values <- values$to(dtype = float64)
    values <- torch::torch_exp(values)
    values <- values / (values + 1.0)
    out_offset <- .offset(out_conf)
    if (.has(out_offset) && out_offset != 0.0) {
        values <- values - out_offset
    }
    out_scale <- .scale(out_conf)
    if (.has(out_scale) && out_scale != 1.0) {
        values <- values / out_scale
    }
    na_mask <- torch::torch_isnan(values)
    # terra truncates the doubles of the CPU path written as integers
    values <- torch::torch_floor(values)
    # R reads the lowest 32-bit integer as NA_integer_
    values <- values$masked_fill(na_mask, -.Machine[["integer.max"]] - 1.0)
    torch::as_array(values$to(dtype = torch::torch_int32())$cpu())
}

#' @title Torch Bayesian smoother closure factory
#' @name  .smooth_fn_bayes_torch
#' @keywords internal
#' @noRd
#' @author Alexandre Assuncao, \email{alexcarssuncao@@gmail.com}
#'
#' @description
#' Mirrors \code{.smooth_fn_bayes()} from \file{api_smooth.R} but calls the
#' torch-backed \code{.torch_smooth_bayes_fraction()} instead of the C++
#' function. It applies the same logit and inverse-logit transforms on CPU,
#' CUDA, or MPS. The kernel is compiled once here and used for every block;
#' when it does not compile, the C++ smoother of \code{.smooth_fn_bayes()}
#' is returned.
#'
#' @param  window_size     Size of the neighbourhood
#'                         (odd integer, min 5, max 21).
#' @param  neigh_fraction  Fraction of highest-valued neighbours to use.
#' @param  smoothness      Numeric vector (length = number of classes).
#' @return A closure \code{function(values, block)} compatible with
#'         \code{.smooth_tile()}.
#'
.smooth_fn_bayes_torch <- function(window_size,
                                   neigh_fraction,
                                   smoothness) {
    # Check window size
    .check_int_parameter(window_size, min = 5L, max = 21L, is_odd = TRUE)
    # Check neigh_fraction
    .check_num_parameter(neigh_fraction, min = 0.1, max = 1.0)
    kernel <- tryCatch(
        .torch_smooth_kernel(
            window_size = window_size,
            neigh_fraction = neigh_fraction
        ),
        error = function(e) NULL
    )
    if (is.null(kernel)) {
        if (.message_warnings()) {
            warning(.conf("messages", ".smooth_fn_bayes_torch"),
                call. = FALSE
            )
        }
        return(.smooth_fn_bayes(
            window_size = window_size,
            neigh_fraction = neigh_fraction,
            smoothness = smoothness
        ))
    }
    prob_max <- .conf("smooth_prob_max")
    prob_min <- .conf("smooth_prob_min")
    smooth_fn <- function(values, block) {
        # Record expected number of pixels
        input_pixels <- nrow(values)
        # Clamp to avoid ±Inf in logit
        values[values == 1.0] <- prob_max
        values[values == 0.0] <- prob_min
        # Transform to logits
        values <- log(values / (rowSums(values) - values))
        # Apply torch Bayesian smoother
        values <- .torch_smooth_bayes_fraction(
            logits     = values,
            nrows      = .nrows(block),
            ncols      = .ncols(block),
            smoothness = smoothness,
            kernel     = kernel
        )
        # Inverse logit
        values <- exp(values) / (exp(values) + 1.0)
        # Sanity check
        .check_processed_values(values, input_pixels)
        values
    }
    smooth_fn
}
