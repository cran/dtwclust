## ----setup, include=FALSE------------------------------------------------
library("dtwclust")
library("ggplot2")
data("dtwclustTimings")

dist_single_results <- dtwclustTimings$dist$single
dist_multiple_results <- dtwclustTimings$dist$multiple
cent_results <- dtwclustTimings$cent
clus_tadpole_results <- dtwclustTimings$tadpole
partitional_results <- dtwclustTimings$partitional

factor_columns <- c("series_length", "window_size", "k", "num_repetitions")
adjust_factors <- function(df) {
    for (column in factor_columns) {
        if (column %in% colnames(df)) {
            df[[column]] <- factor(df[[column]])
        }
    }
    df
}

dist_single_results <- adjust_factors(dist_single_results)
dist_multiple_results <- adjust_factors(dist_multiple_results)
cent_results <- adjust_factors(cent_results)
clus_tadpole_results <- adjust_factors(clus_tadpole_results)
partitional_results$dtwlb_vs_dtwbasic$pam <- adjust_factors(partitional_results$dtwlb_vs_dtwbasic$pam)
partitional_results$dtwlb_vs_dtwbasic$pam_vs_reps <- adjust_factors(partitional_results$dtwlb_vs_dtwbasic$pam_vs_reps)
partitional_results$dtwlb_vs_dtwbasic$dba <- adjust_factors(partitional_results$dtwlb_vs_dtwbasic$dba)
partitional_results$sparse_pam_k$non_symmetric <- adjust_factors(partitional_results$sparse_pam_k$non_symmetric)
partitional_results$sparse_pam_k$symmetric <- adjust_factors(partitional_results$sparse_pam_k$symmetric)

# knitr defaults
knitr::opts_chunk$set(echo = FALSE, comment = "#>")

## ----dist-single-plot, fig.height=7, fig.cap="*Note the different vertical scales for each facet, even though they are all in microseconds.*"----
ggplot(dist_single_results,
       aes(x = series_length,
           y = median_time_us,
           group = window_size,
           colour = window_size)) +
    geom_line() +
    facet_wrap(~distance, scales = "free_y") +
    theme_bw() +
    theme(legend.position = "bottom")

## ----dist-multiple-lb-plot, fig.cap="*The facets' columns indicate the number of parallel workers, whereas the rows indicate the distance that was used. The vertical scales are different for each row, but they are all in milliseconds.*"----
id_gg <- grepl("lb", dist_multiple_results$distance)
ggplot(dist_multiple_results[id_gg,],
       aes(x = num_total,
           y = median_time_ms,
           colour = num_y,
           shape = series_length)) +
    geom_point(size = 3) +
    facet_grid(distance ~ num_workers, scales = "free_y") +
    theme_bw() +
    theme(legend.position = "bottom")

## ----dist-multiple-sbd-plot, fig.height=5, fig.cap="*The facets' columns indicate the number of parallel workers.*"----
id_gg <- grepl("^sbd$", dist_multiple_results$distance)
ggplot(dist_multiple_results[id_gg,],
       aes(x = num_total,
           y = median_time_ms,
           colour = series_length)) +
    geom_point(size = 3) +
    facet_wrap(~num_workers) +
    theme_bw() +
    theme(legend.position = "bottom")

## ----dist-multiple-dtw-plot, fig.height=8, fig.cap="*The facets' columns indicate the number of parallel workers, whereas the rows indicate the length of the series being considered (all series had the same length for each case). All times are in milliseconds.*"----
id_gg <- grepl("^dtw_[um]", dist_multiple_results$distance)
ggplot(dist_multiple_results[id_gg,],
       aes(x = num_total,
           y = median_time_ms,
           colour = window_size,
           shape = distance)) +
    geom_point(size = 3) +
    scale_shape_manual(values = c(0, 3)) +
    facet_grid(series_length ~ num_workers) +
    theme_bw() +
    theme(legend.position = "bottom")

## ----dist-multiple-gak-plot, fig.height=8, fig.cap="*The facets' columns indicate the number of parallel workers, whereas the rows indicate the length of the series being considered (all series had the same length for each case). All times are in milliseconds.*"----
id_gg <- grepl("^gak_", dist_multiple_results$distance)
ggplot(dist_multiple_results[id_gg,],
       aes(x = num_total,
           y = median_time_ms,
           colour = window_size)) +
    geom_point(size = 3) +
    facet_grid(series_length ~ num_workers) +
    theme_bw() +
    theme(legend.position = "bottom")

## ----cent-shape-plot, fig.height=5---------------------------------------
id_gg <- grepl("^shape_", cent_results$cent)
ggplot(cent_results[id_gg,],
       aes(x = num_series,
           y = median_time_ms,
           colour = series_length)) +
    geom_line() +
    facet_wrap(~cent) +
    theme_bw() + 
    theme(legend.position = "bottom")

## ----cent-dba-plot, fig.height=8, fig.cap="*The facets' rows indicate the length of the considered series. In the columns, 'byS' means 'by series' and 'byV' means 'by variable' (see DBA documentation). Note the different vertical scales for each row of the facets.*"----
id_gg <- grepl("^dba_", cent_results$cent)
ggplot(cent_results[id_gg,],
       aes(x = num_series,
           y = median_time_ms,
           colour = window_size)) +
    geom_line() +
    facet_grid(series_length ~ cent, scales = "free") +
    theme_bw() + 
    theme(legend.position = "bottom")

## ----clust-tadpole-plot, fig.cap="*In the facets, lbk stands for lb_keogh and lbi for lb_improved.*"----
ggplot(clus_tadpole_results,
       aes(x = num_series,
           y = median_time_s,
           colour = window_size)) +
    geom_line() +
    facet_wrap(~lb) +
    theme_bw() + 
    theme(legend.position = "bottom")

## ----clust-part-dtw-pam-plot, fig.height=5, fig.cap="*All values are in seconds.*"----
ggdf <- reshape2::melt(partitional_results$dtwlb_vs_dtwbasic$pam, 
                       id.vars=c("num_series", "k", "window_size"))
ggplot(ggdf,
       aes(x = num_series,
           y = value,
           colour = window_size)) +
    geom_line() +
    facet_wrap(~variable) +
    theme_bw() +
    theme(legend.position = "bottom")

## ----clust-part-dtw-pam-reps-plot, fig.height=5, fig.cap="*The window size was fixed at 20 for these experiments.*"----
cols <- setdiff(colnames(partitional_results$dtwlb_vs_dtwbasic$pam_vs_reps),
                "sparse_distmat_filled_percent")
ggdf <- reshape2::melt(partitional_results$dtwlb_vs_dtwbasic$pam_vs_reps[,cols], 
                       id.vars=c("num_series", "k", "num_repetitions"))
ggplot(ggdf,
       aes(x = num_series,
           y = value,
           colour = num_repetitions)) +
    geom_line() +
    facet_wrap(~variable) +
    theme_bw() +
    theme(legend.position = "bottom")

## ----clust-part-dtw-pam-reps-distmat-plot, fig.height=5------------------
ggplot(partitional_results$dtwlb_vs_dtwbasic$pam_vs_reps,
       aes(x = num_series,
           y = sparse_distmat_filled_percent,
           colour = num_repetitions)) +
    geom_line() +
    theme_bw() +
    theme(legend.position = "bottom")

## ----clust-part-dtw-dba-plot, fig.height=5, fig.cap="*All values are in seconds.*"----
ggdf <- reshape2::melt(partitional_results$dtwlb_vs_dtwbasic$dba, 
                       id.vars=c("num_series", "k", "window_size"))
ggplot(ggdf,
       aes(x = num_series,
           y = value,
           colour = window_size)) +
    geom_line() +
    facet_wrap(~variable) +
    theme_bw() +
    theme(legend.position = "bottom")

## ----clust-part-sparse-pam-k-nonsymmetric-plot, fig.height=5, fig.cap="*Note the different vertical scales.*"----
ggdf <- reshape2::melt(partitional_results$sparse_pam_k$non_symmetric, 
                       id.vars=c("num_series", "k"))
ggplot(ggdf,
       aes(x = num_series,
           y = value,
           colour = k)) +
    geom_line() +
    facet_wrap(~variable, scales = "free_y") +
    theme_bw() +
    theme(legend.position = "bottom")

## ----clust-part-sparse-pam-k-symmetric-plot, fig.height=5, fig.cap="*Note the different vertical scales.*"----
ggdf <- reshape2::melt(partitional_results$sparse_pam_k$symmetric, 
                       id.vars=c("num_series", "k"))
ggplot(ggdf,
       aes(x = num_series,
           y = value,
           colour = k)) +
    geom_line() +
    facet_wrap(~variable, scales = "free_y") +
    theme_bw() +
    theme(legend.position = "bottom")

