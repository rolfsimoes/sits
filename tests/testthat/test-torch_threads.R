# A fake cgroup directory with the given files, for .torch_threads()
fake_cgroup <- function(files) {
    dir <- tempfile("cgroup")
    for (f in names(files)) {
        path <- file.path(dir, f)
        dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
        writeLines(files[[f]], path)
    }
    dir
}

test_that(".torch_threads reads the CPU quota of the cgroup", {
    cores <- parallel::detectCores()
    # cgroup v2: quota and period in microseconds
    expect_equal(.torch_threads(fake_cgroup(c("cpu.max" = "200000 100000"))), 2L)
    # a fractional quota rounds down, never below 1
    expect_equal(.torch_threads(fake_cgroup(c("cpu.max" = "250000 100000"))), 2L)
    expect_equal(.torch_threads(fake_cgroup(c("cpu.max" = "50000 100000"))), 1L)
    # no quota: the cores
    expect_equal(.torch_threads(fake_cgroup(c("cpu.max" = "max 100000"))), cores)
    # cgroup v1
    expect_equal(.torch_threads(fake_cgroup(c(
        "cpu/cpu.cfs_quota_us" = "300000", "cpu/cpu.cfs_period_us" = "100000"
    ))), 3L)
    expect_equal(.torch_threads(fake_cgroup(c(
        "cpu/cpu.cfs_quota_us" = "-1", "cpu/cpu.cfs_period_us" = "100000"
    ))), cores)
    # no cgroup at all
    expect_equal(.torch_threads(tempfile("none")), cores)
    # a quota above the cores: the cores
    big <- paste(as.integer(cores + 8L) * 100000L, 100000L)
    expect_equal(.torch_threads(fake_cgroup(c("cpu.max" = big))), cores)
})

test_that("sits_train trains with the CPU quota of torch threads, then restores", {
    skip_if_not(torch::torch_is_installed())
    before <- torch::torch_get_num_threads()
    on.exit(torch::torch_set_num_threads(before), add = TRUE)
    torch::torch_set_num_threads(3L)
    seen <- NA_integer_
    # a training method that records the threads it was given
    ml_method <- function(samples) {
        seen <<- torch::torch_get_num_threads()
        function(values) values
    }
    local_mocked_bindings(.torch_threads = function(...) 2L)
    sits_train(samples_modis_ndvi, ml_method)
    expect_equal(seen, 2L)
    expect_equal(torch::torch_get_num_threads(), 3L)
    # threads given by the caller win over the quota
    sits_train(samples_modis_ndvi, ml_method, threads = 1L)
    expect_equal(seen, 1L)
    expect_equal(torch::torch_get_num_threads(), 3L)
    # restored when the training fails
    expect_error(sits_train(samples_modis_ndvi, function(samples) stop("fails")))
    expect_equal(torch::torch_get_num_threads(), 3L)
    expect_error(sits_train(samples_modis_ndvi, ml_method, threads = 0L))
})
