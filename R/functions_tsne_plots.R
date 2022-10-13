

bin_values = function (x, n_bins, xrng = range(x)) {
    stopifnot(length(xrng) == 2)
    floor(rescale_capped(x, 0:1, xrng) * (n_bins - 1e-05)) +
        1
}

rescale_capped = function (x, to = c(0, 1), from = range(x, na.rm = TRUE, finite = TRUE)) {
    y = scales::rescale(x, to, from)
    y[y > max(to)] = max(to)
    y[y < min(to)] = min(to)
    y
}


