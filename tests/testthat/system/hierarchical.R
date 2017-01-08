context("\tHierarchical")

# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

# =================================================================================================
# multiple k
# =================================================================================================

test_that("Multiple k works as expected.", {
    hc_k <- dtwclust(data_reinterpolated, type = "h", k = 2L:5L,
                     distance = "L2", seed = 938)

    expect_identical(length(hc_k), 4L)

    skip_on_cran()

    hc_k <- lapply(hc_k, reset_nondeterministic)
    assign("hc_k", hc_k, persistent)
})

# =================================================================================================
# hierarchical algorithms
# =================================================================================================

test_that("Hierarchical clustering works as expected.", {
    ## ---------------------------------------------------------- all
    hc_all <- dtwclust(data, type = "hierarchical", k = 20L,
                       distance = "sbd", method = "all")

    hc_all <- lapply(hc_all, reset_nondeterministic)

    assign("hc_all", hc_all, persistent)

    ## ---------------------------------------------------------- non-symmetric
    hc_lbi <- dtwclust(data_reinterpolated, type = "hierarchical", k = 20L,
                       distance = "lbi", method = "all",
                       control = list(window.size = 17L))

    hc_lbi <- lapply(hc_lbi, reset_nondeterministic)

    assign("hc_lbi", hc_lbi, persistent)

    ## ---------------------------------------------------------- custom centroid
    hc_cent <- dtwclust(data, type = "hierarchical", k = 20L,
                        distance = "sbd", method = "all",
                        preproc = zscore, centroid = shape_extraction,
                        seed = 320)

    hc_cent <- lapply(hc_cent, reset_nondeterministic)

    assign("hc_cent", hc_cent, persistent)
})

# =================================================================================================
# cumstom hierarchical function
# =================================================================================================

test_that("A valid custom hierarchical function works as expected.", {
    require(cluster)

    hc_diana <- dtwclust(data, type = "hierarchical", k = 20L,
                         distance = "sbd", method = diana)

    hc_diana <- reset_nondeterministic(hc_diana)
    hc_diana$call <- NULL

    assign("hc_diana", hc_diana, persistent)
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))