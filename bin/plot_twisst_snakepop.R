#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  args[i + 1]
}

weights_file <- get_arg("--weights")
window_file <- get_arg("--metadata")
topos_file <- get_arg("--topologies", NA)
out_prefix <- get_arg("--out-prefix")
formats <- strsplit(get_arg("--formats", "pdf,svg,png"), ",")[[1]]

width <- as.numeric(get_arg("--width", "14"))
height <- as.numeric(get_arg("--height", "8"))
top_n <- as.integer(get_arg("--top-n", "6"))
smooth_k <- as.integer(get_arg("--smooth-k", "15"))

if (is.null(weights_file) || is.null(window_file) || is.null(out_prefix)) {
  stop("Required: --weights FILE --metadata FILE --out-prefix PREFIX")
}

source("bin/twisst/plot_twisst.R")

open_device <- function(outfile, fmt, w, h) {
  if (fmt == "pdf") pdf(outfile, width = w, height = h)
  else if (fmt == "svg") svg(outfile, width = w, height = h)
  else if (fmt == "png") png(outfile, width = w, height = h, units = "in", res = 300)
  else stop(paste("Unsupported format:", fmt))
}

save_plot <- function(name, fmt, w, h, expr) {
  outfile <- paste0(out_prefix, ".", name, ".", fmt)
  open_device(outfile, fmt, w, h)
  tryCatch(
    force(expr),
    finally = dev.off()
  )
}

rolling_mean <- function(x, k = 15) {
  if (k <= 1 || length(x) < k) return(x)
  stats::filter(x, rep(1 / k, k), sides = 2)
}

message("Importing TWISST data...")

if (!is.na(topos_file) && file.exists(topos_file)) {
  twisst_data <- import.twisst(
    weights_files = weights_file,
    window_data_files = window_file,
    topos_file = topos_file,
    split_by_chrom = TRUE,
    recalculate_mid = TRUE
  )
} else {
  twisst_data <- import.twisst(
    weights_files = weights_file,
    window_data_files = window_file,
    split_by_chrom = TRUE,
    recalculate_mid = TRUE
  )
}

top_order <- order(twisst_data$weights_overall_mean, decreasing = TRUE)
top_n <- min(top_n, length(top_order))
top_idx <- top_order[seq_len(top_n)]

subset_topos <- function(obj, idx) {
  out <- obj
  out$weights <- lapply(obj$weights, function(x) x[, idx, drop = FALSE])
  out$weights_raw <- lapply(obj$weights_raw, function(x) x[, idx, drop = FALSE])
  out$weights_mean <- lapply(obj$weights_mean, function(x) x[idx])
  out$weights_overall_mean <- obj$weights_overall_mean[idx]
  out$topos <- obj$topos[idx]
  out
}

twisst_top <- subset_topos(twisst_data, top_idx)

make_genome_df <- function(obj, idx = NULL, smooth = FALSE, k = 15) {
  rows <- list()
  offset <- 0

  regions <- names(obj$weights)

  for (region in regions) {
    w <- obj$weights[[region]]
    wd <- obj$window_data[[region]]

    if (!is.null(idx)) w <- w[, idx, drop = FALSE]

    if (!("mid" %in% names(wd))) wd$mid <- (wd$start + wd$end) / 2

    pos <- wd$mid
    genome_pos <- pos + offset

    if (smooth) {
      for (j in seq_len(ncol(w))) {
        sm <- rolling_mean(w[, j], k = k)
        sm[is.na(sm)] <- w[is.na(sm), j]
        w[, j] <- sm
      }

      row_sums <- rowSums(w, na.rm = TRUE)
      row_sums[row_sums == 0] <- 1
      w <- w / row_sums
    }

    df <- data.frame(
      region = region,
      start = wd$start,
      end = wd$end,
      mid = wd$mid,
      genome_pos = genome_pos,
      w,
      check.names = FALSE
    )

    rows[[region]] <- df
    offset <- max(genome_pos, na.rm = TRUE) + 1
  }

  do.call(rbind, rows)
}

plot_genomewide <- function(obj, idx = NULL, smooth = FALSE, stacked = FALSE,
                            k = 15, main = "") {
  df <- make_genome_df(obj, idx = idx, smooth = smooth, k = k)

  topo_cols_use <- topo_cols
  topo_names <- colnames(df)[!(colnames(df) %in% c("region", "start", "end", "mid", "genome_pos"))]

  if (length(topo_names) > length(topo_cols_use)) {
    topo_cols_use <- rainbow(length(topo_names))
  } else {
    topo_cols_use <- topo_cols_use[seq_along(topo_names)]
  }

  par(mar = c(4, 4, 2, 1))
  plot(
    NA,
    xlim = range(df$genome_pos, na.rm = TRUE),
    ylim = c(0, 1),
    xlab = "Genome position",
    ylab = "Topology weighting",
    main = main,
    bty = "n",
    xaxt = "n"
  )

  regions <- unique(df$region)
  centers <- tapply(df$genome_pos, df$region, mean, na.rm = TRUE)

  abline(v = tapply(df$genome_pos, df$region, min, na.rm = TRUE), col = "grey85", lty = 3)

  if (stacked) {
    mat <- as.matrix(df[, topo_names, drop = FALSE])
    upper <- t(apply(mat, 1, cumsum))
    lower <- upper - mat

    for (j in rev(seq_along(topo_names))) {
      polygon(
        c(df$genome_pos, rev(df$genome_pos)),
        c(upper[, j], rev(lower[, j])),
        col = topo_cols_use[j],
        border = NA
      )
    }
  } else {
    for (j in seq_along(topo_names)) {
      lines(
        df$genome_pos,
        df[[topo_names[j]]],
        col = topo_cols_use[j],
        lwd = 1.2
      )
    }
  }

  axis(1, at = centers, labels = names(centers), las = 2, cex.axis = 0.45)
  legend(
    "topright",
    legend = topo_names,
    col = topo_cols_use,
    lwd = 2,
    bty = "n",
    cex = 0.7
  )
}

for (fmt in formats) {
  fmt <- trimws(fmt)
  message("Writing ", fmt, " plots...")

  save_plot("summary_barplot", fmt, 10, 6, {
    plot.twisst.summary(twisst_data, lwd = 3, cex = 0.7)
  })

  save_plot("summary_boxplot", fmt, 10, 6, {
    plot.twisst.summary.boxplot(twisst_data, lwd = 3, cex = 0.7, outline = FALSE)
  })

  save_plot("summary_topN_barplot", fmt, 10, 6, {
    plot.twisst.summary(twisst_top, lwd = 3, cex = 0.7)
  })

  save_plot("summary_topN_boxplot", fmt, 10, 6, {
    plot.twisst.summary.boxplot(twisst_top, lwd = 3, cex = 0.7, outline = FALSE)
  })

  save_plot("topologies_only", fmt, width, height, {
    par(mfrow = c(1, length(twisst_top$topos)), mar = c(1, 1, 2, 1), xpd = NA)
    for (i in seq_along(twisst_top$topos)) {
      plot.phylo(
        twisst_top$topos[[i]],
        type = "clad",
        edge.color = topo_cols[i],
        edge.width = 5,
        label.offset = 0.2,
        cex = 0.8
      )
      mtext(side = 3, text = names(twisst_top$topos)[i], col = topo_cols[i])
    }
  })

  save_plot("genomewide_raw_all_overlay", fmt, width, height, {
    plot_genomewide(twisst_data, smooth = FALSE, stacked = FALSE, main = "TWISST raw topology weights")
  })

  save_plot("genomewide_raw_all_stacked", fmt, width, height, {
    plot_genomewide(twisst_data, smooth = FALSE, stacked = TRUE, main = "TWISST raw topology weights, stacked")
  })

  save_plot("genomewide_smooth_all_overlay", fmt, width, height, {
    plot_genomewide(twisst_data, smooth = TRUE, stacked = FALSE, k = smooth_k, main = "TWISST smoothed topology weights")
  })

  save_plot("genomewide_smooth_all_stacked", fmt, width, height, {
    plot_genomewide(twisst_data, smooth = TRUE, stacked = TRUE, k = smooth_k, main = "TWISST smoothed topology weights, stacked")
  })

  save_plot("genomewide_raw_topN_overlay", fmt, width, height, {
    plot_genomewide(twisst_data, idx = top_idx, smooth = FALSE, stacked = FALSE, main = "TWISST raw Top-N topology weights")
  })

  save_plot("genomewide_raw_topN_stacked", fmt, width, height, {
    plot_genomewide(twisst_data, idx = top_idx, smooth = FALSE, stacked = TRUE, main = "TWISST raw Top-N topology weights, stacked")
  })

  save_plot("genomewide_smooth_topN_overlay", fmt, width, height, {
    plot_genomewide(twisst_data, idx = top_idx, smooth = TRUE, stacked = FALSE, k = smooth_k, main = "TWISST smoothed Top-N topology weights")
  })

  save_plot("genomewide_smooth_topN_stacked", fmt, width, height, {
    plot_genomewide(twisst_data, idx = top_idx, smooth = TRUE, stacked = TRUE, k = smooth_k, main = "TWISST smoothed Top-N topology weights, stacked")
  })
}

message("Done.")
