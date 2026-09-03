# 单细胞富集分析条形图（通路+基因）
# Reviewed Gallery plotting template; synthetic preview verified in this repository.
#
# Input: one row per cluster-pathway result. See ../data_schema.yml.
# Output: a ggplot object. No file is written unless save_path is supplied.

plot_clustered_pathway_bars <- function(
    enrichment,
    clusters = NULL,
    top_n = 5L,
    rank_by = "pvalue",
    cluster_col = "Cluster",
    description_col = "Description",
    pvalue_col = "pvalue",
    gene_col = "geneID",
    cluster_colors = NULL,
    title = "单细胞富集分析条形图（通路+基因）",
    x_label = expression(-Log[10](P)),
    bar_width = 0.62,
    group_gap = 0.85,
    label_x_fraction = 0.015,
    base_size = 11,
    save_path = NULL,
    width = 5.2,
    height = 5.1,
    dpi = 300) {
  required_packages <- c("dplyr", "ggplot2", "grid")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages) > 0L) {
    stop(
      "Missing required R packages: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }

  required_columns <- c(
    cluster_col,
    description_col,
    pvalue_col,
    gene_col,
    rank_by
  )
  missing_columns <- setdiff(required_columns, names(enrichment))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n < 1) {
    stop("top_n must be one positive integer.", call. = FALSE)
  }
  if (!is.numeric(enrichment[[pvalue_col]])) {
    stop(pvalue_col, " must be numeric.", call. = FALSE)
  }
  if (!is.numeric(enrichment[[rank_by]])) {
    stop(rank_by, " must be numeric.", call. = FALSE)
  }

  plot_data <- enrichment
  plot_data[[cluster_col]] <- as.character(plot_data[[cluster_col]])
  plot_data[[description_col]] <- as.character(plot_data[[description_col]])
  plot_data[[gene_col]] <- as.character(plot_data[[gene_col]])

  if (is.null(clusters)) {
    clusters <- unique(plot_data[[cluster_col]])
  }
  clusters <- as.character(clusters)
  if (length(clusters) == 0L) {
    stop("At least one cluster must be selected.", call. = FALSE)
  }

  absent_clusters <- setdiff(clusters, unique(plot_data[[cluster_col]]))
  if (length(absent_clusters) > 0L) {
    stop(
      "Selected clusters absent from the data: ",
      paste(absent_clusters, collapse = ", "),
      call. = FALSE
    )
  }

  plot_data <- plot_data[
    plot_data[[cluster_col]] %in% clusters &
      !is.na(plot_data[[description_col]]) &
      nzchar(plot_data[[description_col]]) &
      !is.na(plot_data[[pvalue_col]]) &
      is.finite(plot_data[[pvalue_col]]) &
      plot_data[[pvalue_col]] <= 1 &
      !is.na(plot_data[[rank_by]]) &
      is.finite(plot_data[[rank_by]]),
    ,
    drop = FALSE
  ]
  if (nrow(plot_data) == 0L) {
    stop("No valid enrichment rows remain after filtering.", call. = FALSE)
  }

  non_positive <- plot_data[[pvalue_col]] <= 0
  if (any(non_positive)) {
    positive_values <- plot_data[[pvalue_col]][plot_data[[pvalue_col]] > 0]
    replacement <- if (length(positive_values) > 0L) {
      min(positive_values) * 0.1
    } else {
      .Machine$double.xmin
    }
    warning(
      "Replaced ", sum(non_positive), " non-positive p-value(s) with ",
      format(replacement, scientific = TRUE), " before -log10 transformation."
    )
    plot_data[[pvalue_col]][non_positive] <- replacement
  }

  plot_data <- dplyr::as_tibble(plot_data) |>
    dplyr::mutate(
      .cluster = factor(.data[[cluster_col]], levels = clusters),
      .description = .data[[description_col]],
      .genes = gsub("/", ", ", .data[[gene_col]], fixed = TRUE),
      .score = -log10(.data[[pvalue_col]]),
      .rank = .data[[rank_by]]
    ) |>
    dplyr::group_by(.cluster) |>
    dplyr::slice_min(order_by = .rank, n = as.integer(top_n), with_ties = FALSE) |>
    dplyr::arrange(.cluster, dplyr::desc(.score), .description) |>
    dplyr::ungroup()

  rows_per_cluster <- table(plot_data$.cluster, useNA = "no")
  cluster_sizes <- as.integer(rows_per_cluster[clusters])
  if (anyNA(cluster_sizes) || any(cluster_sizes == 0L)) {
    stop("Every selected cluster must retain at least one pathway.", call. = FALSE)
  }

  plot_data$.row_in_cluster <- ave(
    seq_len(nrow(plot_data)),
    plot_data$.cluster,
    FUN = seq_along
  )
  block_starts <- c(0, head(cumsum(cluster_sizes), -1L)) +
    seq(0, by = group_gap, length.out = length(clusters))
  names(block_starts) <- clusters
  plot_data$.y <- unname(block_starts[as.character(plot_data$.cluster)]) +
    plot_data$.row_in_cluster

  cluster_labels <- data.frame(
    .cluster = clusters,
    .y = unname(block_starts[clusters]) + (cluster_sizes + 1) / 2,
    stringsAsFactors = FALSE
  )

  if (is.null(cluster_colors)) {
    default_colors <- c("#6BB9D2", "#D55640", "#7A9E4E", "#8E6BBE")
    if (length(clusters) > length(default_colors)) {
      cluster_colors <- grDevices::hcl.colors(length(clusters), "Dark 3")
    } else {
      cluster_colors <- default_colors[seq_along(clusters)]
    }
    names(cluster_colors) <- clusters
  } else {
    if (is.null(names(cluster_colors))) {
      if (length(cluster_colors) != length(clusters)) {
        stop(
          "An unnamed cluster_colors vector must match the number of clusters.",
          call. = FALSE
        )
      }
      names(cluster_colors) <- clusters
    }
    missing_colors <- setdiff(clusters, names(cluster_colors))
    if (length(missing_colors) > 0L) {
      stop(
        "Missing colours for clusters: ",
        paste(missing_colors, collapse = ", "),
        call. = FALSE
      )
    }
    cluster_colors <- cluster_colors[clusters]
  }

  max_score <- max(plot_data$.score, na.rm = TRUE)
  if (!is.finite(max_score) || max_score <= 0) {
    stop("The transformed p-values do not define a positive x range.", call. = FALSE)
  }
  label_x <- max_score * label_x_fraction
  cluster_label_x <- -0.055 * max_score
  x_upper <- max_score * 1.06
  gene_colors <- grDevices::adjustcolor(
    cluster_colors,
    red.f = 0.72,
    green.f = 0.72,
    blue.f = 0.72
  )

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .score, y = .y, fill = .cluster)
  ) +
    ggplot2::geom_col(
      width = bar_width,
      alpha = 0.88,
      orientation = "y"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = label_x, label = .description),
      hjust = 0,
      vjust = 0.1,
      size = base_size / ggplot2::.pt,
      colour = "black",
      lineheight = 0.95
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = label_x,
        label = .genes,
        colour = .cluster
      ),
      hjust = 0,
      vjust = 2.2,
      size = (base_size - 1) / ggplot2::.pt,
      fontface = "italic",
      lineheight = 0.95
    ) +
    ggplot2::geom_text(
      data = cluster_labels,
      ggplot2::aes(x = cluster_label_x, y = .y, label = .cluster),
      inherit.aes = FALSE,
      angle = 90,
      size = base_size / ggplot2::.pt,
      hjust = 0.5,
      vjust = 0.5,
      colour = "black"
    ) +
    ggplot2::scale_fill_manual(values = cluster_colors, guide = "none") +
    ggplot2::scale_colour_manual(values = gene_colors, guide = "none") +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_y_reverse(
      limits = c(max(plot_data$.y) + 0.9, min(plot_data$.y) - 0.65),
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, x_upper),
      clip = "off"
    ) +
    ggplot2::labs(x = x_label, y = NULL, title = title) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = base_size + 1),
      axis.text.x = ggplot2::element_text(size = base_size),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.ticks.length = grid::unit(1.8, "mm"),
      plot.title = ggplot2::element_text(
        size = base_size + 1,
        face = "bold",
        hjust = 0
      ),
      plot.margin = ggplot2::margin(t = 5.5, r = 18, b = 5.5, l = 42)
    )

  if (!is.null(save_path)) {
    ggplot2::ggsave(
      filename = save_path,
      plot = p,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
    if (!file.exists(save_path) || file.info(save_path)$size <= 0L) {
      stop("Plot export failed or produced an empty file: ", save_path)
    }
  }

  p
}

# Example using the bundled synthetic table:
# enrichment <- utils::read.csv("data/example.csv", check.names = FALSE)
# p <- plot_clustered_pathway_bars(
#   enrichment = enrichment,
#   clusters = c("CD14+ monocyte", "CD16+ monocyte"),
#   cluster_colors = c(
#     "CD14+ monocyte" = "#6BB9D2",
#     "CD16+ monocyte" = "#D55640"
#   ),
#   title = "KEGG pathway enrichment",
#   save_path = "preview.png"
# )
# p


# ---- Bundled synthetic example ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[[1]], winslash = "/", mustWork = TRUE)) else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
root <- if (basename(script_dir) == "code") dirname(script_dir) else script_dir
input <- utils::read.csv(file.path(root, "data", "example.csv"), check.names = FALSE, stringsAsFactors = FALSE)
out_arg <- sub("^--output-dir=", "", args[grepl("^--output-dir=", args)])
out_dir <- if (length(out_arg) && nzchar(out_arg[[1]])) out_arg[[1]] else file.path(root, "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
output <- file.path(out_dir, "render.png")
p <- plot_clustered_pathway_bars(
  enrichment = input,
  clusters = c("CD14+ monocyte", "CD16+ monocyte"),
  top_n = 5L, rank_by = "pvalue",
  cluster_colors = c("CD14+ monocyte" = "#6BB9D2", "CD16+ monocyte" = "#D55640"),
  title = "Single-cell pathway enrichment", save_path = output,
  width = 5.2, height = 5.1, dpi = 300
)
if (!inherits(p, "ggplot")) stop("The plotting function did not return a ggplot object.", call. = FALSE)
message("wrote ", output)
