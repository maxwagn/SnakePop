#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default=NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  args[i + 1]
}

weights_file   <- get_arg("--weights")
metadata_file  <- get_arg("--metadata")
topologies_file <- get_arg("--topologies")
out_prefix     <- get_arg("--out-prefix")

top_n          <- get_arg("--top-n", "3")
smooth_n       <- as.integer(get_arg("--smooth", "1"))
combine_other  <- tolower(get_arg("--combine-other", "true")) %in% c("true","1","yes")
width          <- as.numeric(get_arg("--width", "14"))
height         <- as.numeric(get_arg("--height", "7"))
formats        <- strsplit(get_arg("--formats", "pdf,svg,png"), ",")[[1]]

# Optional: source Simon Martin's original plotting functions if present
# Put original plot_twisst.R here:
# bin/twisst/plot_twisst.R
plot_fun <- "bin/twisst/plot_twisst.R"
if (file.exists(plot_fun)) {
  source(plot_fun)
}

# -----------------------------
# Read SnakePop / Twisst output
# -----------------------------

weights <- read.table(
  weights_file,
  header = TRUE,
  comment.char = "#",
  stringsAsFactors = FALSE
)

metadata <- read.table(
  metadata_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

topologies <- readLines(topologies_file)
topologies <- topologies[nchar(topologies) > 0]

if (nrow(weights) != nrow(metadata)) {
  stop(
    "weights and metadata have different number of rows: ",
    nrow(weights), " vs ", nrow(metadata)
  )
}

topo_cols <- colnames(weights)

# -----------------------------
# Select topologies
# -----------------------------

mean_weights <- colMeans(weights, na.rm = TRUE)
mean_weights <- sort(mean_weights, decreasing = TRUE)

if (tolower(top_n) == "all") {
  selected <- names(mean_weights)
} else {
  selected <- names(mean_weights)[seq_len(min(as.integer(top_n), length(mean_weights)))]
}

plot_weights <- weights[, selected, drop = FALSE]

if (combine_other && tolower(top_n) != "all") {
  other <- setdiff(colnames(weights), selected)
  if (length(other) > 0) {
    plot_weights$Other <- rowSums(weights[, other, drop = FALSE], na.rm = TRUE)
  }
}

if (smooth_n > 1) {
  smooth_vec <- function(x, n) {
    stats::filter(x, rep(1 / n, n), sides = 2)
  }
  plot_weights <- as.data.frame(lapply(plot_weights, smooth_vec, n = smooth_n))
  plot_weights[is.na(plot_weights)] <- weights[is.na(plot_weights)]
}

# -----------------------------
# Coordinates
# -----------------------------

metadata$scaffold <- as.character(metadata$scaffold)
metadata$start <- as.numeric(metadata$start)
metadata$end <- as.numeric(metadata$end)
metadata$mid <- (metadata$start + metadata$end) / 2

chroms <- unique(metadata$scaffold)
offsets <- numeric(length(chroms))
names(offsets) <- chroms

running <- 0
centers <- numeric(length(chroms))
bounds <- numeric(length(chroms))

for (i in seq_along(chroms)) {
  chr <- chroms[i]
  chr_max <- max(metadata$end[metadata$scaffold == chr], na.rm = TRUE)
  offsets[chr] <- running
  centers[i] <- running + chr_max / 2
  running <- running + chr_max
  bounds[i] <- running
}

x <- metadata$mid + offsets[metadata$scaffold]

# -----------------------------
# Colours similar to Twisst examples
# -----------------------------

cols <- c(
  "#1f78b4", "#ff7f00", "#33a02c", "#e31a1c",
  "#6a3d9a", "#b15928", "#a6cee3", "#fdbf6f",
  "#b2df8a", "#fb9a99", "#cab2d6", "#ffff99"
)

while (length(cols) < ncol(plot_weights)) cols <- c(cols, cols)
cols <- cols[seq_len(ncol(plot_weights))]

# -----------------------------
# Labels
# -----------------------------

topo_label <- function(col) {
  if (col == "Other") return("Other")
  idx <- as.integer(gsub("[^0-9]", "", col))
  if (is.na(idx) || idx < 1 || idx > length(topologies)) return(col)
  topologies[idx]
}

legend_labels <- sapply(colnames(plot_weights), function(z) {
  if (z == "Other") return("Other")
  sprintf("%s (%.1f%%)", z, mean(weights[[z]], na.rm = TRUE) * 100)
})

# -----------------------------
# Plot function
# -----------------------------

make_plot <- function(outfile, device_fun) {
  device_fun(outfile, width = width, height = height)

  layout(
    matrix(c(1, 2, 3), ncol = 1),
    heights = c(1.1, 3.2, 2.2)
  )

  par(mar = c(0.5, 4, 2, 1), xpd = NA)

  # Topology text panel
  plot.new()
  title("Twisst topology weights", adj = 0, font.main = 2)

  n <- ncol(plot_weights)
  xs <- seq(0.1, 0.9, length.out = n)

  for (i in seq_len(n)) {
    rect(xs[i] - 0.035, 0.62, xs[i] + 0.035, 0.72, col = cols[i], border = NA)
    text(xs[i], 0.50, legend_labels[i], cex = 0.8)
    text(xs[i], 0.25, topo_label(colnames(plot_weights)[i]), cex = 0.55, family = "mono")
  }

  # Stacked area
  par(mar = c(0.5, 4, 0.5, 1))
  plot(
    range(x), c(0, 1),
    type = "n",
    xlab = "",
    ylab = "Topology weight",
    xaxt = "n",
    las = 1
  )

  y0 <- rep(0, nrow(plot_weights))
  for (i in seq_len(ncol(plot_weights))) {
    y1 <- y0 + plot_weights[[i]]
    polygon(
      c(x, rev(x)),
      c(y0, rev(y1)),
      col = cols[i],
      border = NA
    )
    y0 <- y1
  }

  abline(v = bounds[-length(bounds)], col = "white", lwd = 1)

  legend(
    "top",
    legend = legend_labels,
    fill = cols,
    horiz = TRUE,
    bty = "n",
    cex = 0.8,
    inset = -0.12
  )

  # Line plot
  par(mar = c(5, 4, 0.5, 1))
  plot(
    range(x), c(0, 1),
    type = "n",
    xlab = "Genomic position",
    ylab = "Topology weight",
    xaxt = "n",
    las = 1
  )

  for (i in seq_len(ncol(plot_weights))) {
    lines(x, plot_weights[[i]], col = cols[i], lwd = 1)
  }

  abline(v = bounds[-length(bounds)], col = "grey85", lwd = 0.8)
  axis(1, at = centers, labels = chroms, las = 2, cex.axis = 0.7)

  dev.off()
}

for (fmt in formats) {
  fmt <- tolower(trimws(fmt))
  outfile <- paste0(out_prefix, ".", fmt)

  if (fmt == "pdf") {
    make_plot(outfile, pdf)
  } else if (fmt == "svg") {
    make_plot(outfile, svg)
  } else if (fmt == "png") {
    make_plot(outfile, function(file, width, height) {
      png(file, width = width, height = height, units = "in", res = 300)
    })
  } else {
    warning("Unknown format skipped: ", fmt)
  }
}
