## ----setup, include = FALSE, cache = FALSE-------------------------------
library(knitr)
library(dtwclust)

data(uciCT)

# knitr defaults
opts_chunk$set(fig.width = 8, fig.asp = 0.5625,
               out.width = "\\linewidth",
               fig.align = "center", fig.pos = "htbp",
               cache = TRUE, echo = FALSE, autodep = TRUE)

## ----dtw-intuition, out.width = "0.75\\linewidth", fig.cap = "Sample alignment performed by the DTW algorithm between two series. The dashed blue lines exemplify how some points are mapped to each other, which shows how they can be warped in time. Note that the vertical position of each series was artificially altered for visualization."----
dtw_example <- dtw(CharTraj[[1L]], CharTraj[[2L]], keep.internals = TRUE)
plot(dtw_example, type = "two",
     offset = 1, match.indices = 30,
     match.col = "blue",
     xlab = "Time", ylab = "Series")

## ----dtw-path, fig.cap = "Visual representation of the optimum path found. The big square in the center represents the LCM created for these specific series."----
plot(dtw_example, type = "three")

## ----step-patterns, out.width = "0.45\\linewidth", fig.width = 4, fig.asp = 1, fig.cap = "Two common step patterns used by DTW when traversing the LCM. At each step, the lines denote the allowed directions that can be taken, as well as the weight associated with each one.", fig.subcap = c("\\code{symmetric1} step pattern", "\\code{symmetric2} step pattern")----
plot(symmetric1)
plot(symmetric2)

## ----dtw-window-plot, out.width = "0.6\\linewidth", fig.asp = 1, fig.cap = "Visual representation of the Sakoe-Chiba constraint for DTW. The red elements will not be considered by the algorithm when traversing the LCM."----
dtwWindow.plot(sakoeChibaWindow, window.size = 2, reference = 10, query = 10)

## ----envelop-plot, out.width = "0.75\\linewidth", fig.cap = "Visual representation of a time-series (shown as a solid black line) and its corresponding envelops based on the Sakoe-Chiba window. The green dashed line represents the upper envelop, while the red dashed line represents the lower envelop."----
lbk <- lb_keogh(CharTraj[[1L]], CharTraj[[2L]], window.size = 15)
matplot(cbind(lbk$lower.env, lbk$upper.env),
        type = "l", lty = 2, col = 2:3,
        xlab = "Time", ylab = "Series")
lines(CharTraj[[2L]])

## ----sbd-alignment, out.width = "0.45\\linewidth", fig.width = 6, fig.cap = "Visualization of the NCCc-based alignment performed on two sample series. After alignment, the second (red) series is either truncated and/or prepended/appended with zeros so that its length matches the first(black) series.", fig.subcap = c("Series before alignment", "Series after alignment")----
matplot(cbind(CharTraj[[61L]], CharTraj[[65L]]),
        type = "l", lty = 1L,
        xlab = "Time", ylab = "Series")

sbd_align <- SBD(CharTraj[[61L]], CharTraj[[65L]])

matplot(cbind(CharTraj[[61L]], sbd_align$yshift),
        type = "l", lty = 1L,
        xlab = "Time", ylab = "Series")

## ----dendrogram, fig.width = 12, fig.cap = "Sample dendrogram created by using SBD and average linkage."----
hc <- dtwclust(CharTraj, type = "h", method = "average", distance = "sbd")
plot(hc)

## ----example-register-proxy, echo = TRUE, message = FALSE----------------
require(TSclust)

proxy::pr_DB$set_entry(FUN = diss.ACF, names=c("ACFD"),
                       loop = TRUE, type = "metric", distance = TRUE,
                       description = "Autocorrelation-based distance")

# Taking just a subset of the data
# Note that subsetting with single brackets preserves the list format
proxy::dist(CharTraj[3:8], method = "ACFD", upper = TRUE)

## ----example-dtw-lb, echo = TRUE-----------------------------------------
# Reinterpolate to same length
data <- lapply(CharTraj, reinterpolate, newLength = 180L)

# Calculate the DTW distances between all elements
system.time(D1 <- proxy::dist(data[1L:5L], data[6L:50L],
                              method = "DTW",
                              window.type = "sakoechiba",
                              window.size = 20L))

# Nearest neighbors
NN1 <- apply(D1, 1L, which.min)

# Calculate the distance matrix with dtw_lb
system.time(D2 <- dtw_lb(data[1L:5L], data[6L:50L],
                         window.size = 20L))

# Nearest neighbors
NN2 <- apply(D2, 1L, which.min)

# Same results?
all(NN1 == NN2)

## ----example-hc, echo = TRUE, fig.cap = c("Resulting dendrogram after hierarchical clustering.", "Obtained clusters and their respective prototypes (centroids) shown as dashed lines.")----
hc_sbd <- dtwclust(CharTraj, type = "h", k = 20L,
                   method = "average", preproc = zscore,
                   distance = "sbd", centroid = shape_extraction,
                   seed = 899, control = list(trace = TRUE))

# Cluster sizes
table(hc_sbd@cluster)

# By default, the dendrogram is plotted in hierarchical clustering
plot(hc_sbd)

# The series and the obtained prototypes can be plotted too
plot(hc_sbd, type = "sc")

## ----example-hc-2, echo = TRUE, out.width = "0.45\\linewidth", fig.cap = "Side by side comparison of the series in the first cluster and the obtained prototype.", fig.subcap = c("Series in cluster 1", "Prototype obtained by applying shape extraction to cluster 1")----
# Focusing on the first cluster
plot(hc_sbd, type = "series", clus = 1L)
plot(hc_sbd, type = "centroids", clus = 1L)

## ----example-pc, echo = TRUE, warning = FALSE----------------------------
# Reinterpolate to same length
data <- lapply(CharTraj, reinterpolate,
               newLength = max(lengths(CharTraj)))

# z-normalization
data <- zscore(data[60L:100L])

pc_dtw <- dtwclust(data, k = 4L,
                   distance = "dtw2", centroid = "dba",
                   seed = 8, control = list(window.size = 20L,
                                            norm = "L2",
                                            trace = TRUE))

pc_dtwlb <- dtwclust(data, k = 4L,
                     distance = "dtw_lb", centroid = "dba",
                     seed = 8, control = list(window.size = 20L,
                                              norm = "L2",
                                              trace = TRUE))

pc_ks <- dtwclust(data, k = 4L,
                  distance = "sbd", centroid = "shape",
                  seed = 8, control = list(trace = TRUE))

pc_tp <-dtwclust(data, k = 4L,
                 type = "tadpole", dc = 1.5,
                 seed = 8, control = list(window.size = 20L,
                                          trace = TRUE))

sapply(list(DTW = pc_dtw, DTW_LB = pc_dtwlb, kShape = pc_ks, TADPole = pc_tp),
       cvi, b = CharTrajLabels[60L:100L], type = "VI")

## ----example-fuzzy, echo = TRUE, fig.cap = "Visualization of the clusters that would result from changing the fuzzy partition to a crisp one. Note that the original time-series are used, so the centroids are not shown, since the centroids are made of autocorrelation coefficients."----
# Calculate autocorrelation up to 50th lag
acf_fun <- function(dat) {
     lapply(dat, function(x) {
          as.numeric(acf(x, lag.max = 50, plot = FALSE)$acf)}
     )
}

# Fuzzy c-means
fc <- dtwclust(CharTraj[1:20], type = "f", k = 4L,
               preproc = acf_fun, distance = "L2",
               seed = 42)

# Fuzzy membership matrix
fc@fcluster

# Are constraints fulfilled?
all.equal(rep(1, 20), rowSums(fc@fcluster), check.attributes = FALSE)

# Plot crisp partition in the original space
plot(fc, data = CharTraj[1:20], type = "series")

## ----example-doParallel, echo = TRUE, message = FALSE--------------------
require(doParallel)

# Create parallel workers
workers <- makeCluster(2L)

# Preload dtwclust in each worker; not necessary but useful
invisible(clusterEvalQ(workers, library(dtwclust)))

# Register the backend; this step MUST be done
registerDoParallel(workers)

# Calling dtwclust
pc_par <- dtwclust(CharTraj[1L:20L], k = 4L,
                   distance = "dtw", centroid = "dba",
                   seed = 938, control = list(trace = TRUE,
                                              window.size = 15L))

# Stop parallel workers
stopCluster(workers)

# Go back to sequential computation
registerDoSEQ()

