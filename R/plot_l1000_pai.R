.require_plot_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Package `", package, "` is required for plotting. Please install it.")
  }
}

.plot_if_null <- function(x, y) {
  if (is.null(x)) {
    y
  } else {
    x
  }
}

.default_plot_p_col <- function(df) {
  candidates <- c(
    "padj_coherent", "padj_combined", "p.adjust", "p_coherent",
    "p_combined", "pvalue"
  )
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    stop("No p-value column found. Please provide `p_col` explicitly.")
  }
  hit[[1]]
}

.result_library <- function(x, library = NULL) {
  if (!is.null(library)) {
    return(library)
  }

  if (methods::is(x, "L1000PAIResult")) {
    return(x@metadata$library)
  }

  NULL
}

.plot_label_col <- function(df) {
  if ("compound_name" %in% colnames(df)) {
    return("compound_name")
  }
  if ("base_ligand" %in% colnames(df)) {
    return("base_ligand")
  }
  if ("target" %in% colnames(df)) {
    return("target")
  }
  if ("perturbagen" %in% colnames(df)) {
    return("perturbagen")
  }
  if ("term" %in% colnames(df)) {
    return("term")
  }
  stop("Could not find a label column in the result table.")
}

.plot_activity <- function(
    x,
    library = NULL,
    p_col = NULL,
    score_col = "activity_score",
    significance_cutoff = 0.05,
    activity_threshold = 1,
    top_n = 200,
    label_top = 10,
    title = NULL
) {
  .require_plot_package("ggplot2")

  lib <- .result_library(x, library = library)
  df <- if (!is.null(lib) && lib %in% available_l1000_annotation_catalogs()) {
    annotate_l1000_result(x, library = lib)
  } else {
    .as_result_df(x)
  }

  if (!is.character(score_col) || length(score_col) != 1L ||
      !score_col %in% colnames(df)) {
    stop("`score_col` must identify an activity-score column in the result table.")
  }
  if (!"activity_class" %in% colnames(df)) {
    df$activity_class <- "Unclassified"
  }

  p_col <- .plot_if_null(p_col, .default_plot_p_col(df))
  label_col <- .plot_label_col(df)
  df$plot_p <- suppressWarnings(as.numeric(df[[p_col]]))
  df$neg_log10_p <- -log10(pmax(df$plot_p, 1e-300))
  df$plot_score <- suppressWarnings(as.numeric(df[[score_col]]))
  df$plot_label <- as.character(df[[label_col]])

  df <- df[!is.na(df$plot_score) & !is.na(df$neg_log10_p), , drop = FALSE]
  df <- df[order(df$plot_p, -abs(df$plot_score), df$plot_label), , drop = FALSE]
  if (!is.null(top_n) && is.finite(top_n) && nrow(df) > top_n) {
    plot_df <- utils::head(df, top_n)
  } else {
    plot_df <- df
  }

  label_n <- if (is.null(label_top) || length(label_top) == 0 || is.na(label_top) || label_top <= 0) {
    0L
  } else {
    min(floor(label_top), nrow(plot_df))
  }
  label_df <- plot_df[seq_len(label_n), , drop = FALSE]

  class_values <- c(
    Activated = "#D55E00",
    Inhibited = "#009E73",
    Discordant = "#7570B3",
    Not_significant = "grey78",
    Positive_enrichment = "#D55E00",
    Negative_enrichment = "#009E73",
    Unclassified = "grey60"
  )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = plot_score,
      y = neg_log10_p,
      color = activity_class
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, color = "grey82") +
    ggplot2::geom_vline(xintercept = c(-activity_threshold, activity_threshold), linewidth = 0.35, color = "grey88", linetype = "dotted") +
    ggplot2::geom_hline(yintercept = -log10(significance_cutoff), linewidth = 0.35, color = "grey76", linetype = "dashed") +
    ggplot2::geom_point(alpha = 0.82, size = 2.2) +
    ggplot2::scale_color_manual(values = class_values, drop = FALSE) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.18))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.24))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = .plot_if_null(title, "L1000PAI activity plot"),
      subtitle = paste0("x = ", score_col, "; y = -log10(", p_col, ")"),
      x = score_col,
      y = paste0("-log10(", p_col, ")"),
      color = "Activity class"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )

  if (label_top > 0 && nrow(label_df) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p +
        ggrepel::geom_text_repel(
          data = label_df,
          ggplot2::aes(label = plot_label),
          size = 3,
          box.padding = 0.55,
          point.padding = 0.35,
          min.segment.length = 0,
          force = 4,
          force_pull = 0.05,
          max.time = 2,
          max.iter = 10000,
          max.overlaps = Inf,
          seed = 42,
          segment.alpha = 0.75,
          segment.size = 0.25,
          show.legend = FALSE
        )
    } else {
      p <- p +
        ggplot2::geom_text(
          data = label_df,
          ggplot2::aes(label = plot_label),
          size = 3,
          vjust = -0.6,
          show.legend = FALSE
        )
    }
  }

  p
}

.default_category_table_name <- function(library) {
  switch(
    library,
    chem = "big_class_summary",
    ligand = "ligand_group_summary",
    kinase_gpcr = "signal_family_summary",
    NULL
  )
}

.blank_category_plot <- function(title, message) {
  ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x, y)) +
    ggplot2::geom_text(ggplot2::aes(label = message), size = 4, color = "grey35") +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 11)
}

.plot_category <- function(
    x,
    library = NULL,
    category_table = NULL,
    top_n = 12,
    title = NULL
) {
  .require_plot_package("ggplot2")

  lib <- .result_library(x, library = library)
  if (is.null(category_table)) {
    if (is.null(lib) || !lib %in% available_l1000_annotation_catalogs()) {
      stop("`category_table` is required when no supported result library can be inferred.")
    }
    table_name <- .default_category_table_name(lib)
    tables <- if (methods::is(x, "L1000PAIResult")) {
      annotation_tables(x)
    } else {
      build_l1000_annotation_tables(x, library = lib)
    }
    category_table <- tables[[table_name]]
    if (is.null(category_table)) {
      stop("Could not find category table `", table_name, "`.")
    }
  }

  if (!is.data.frame(category_table) || nrow(category_table) == 0) {
    return(.blank_category_plot(.plot_if_null(title, "L1000PAI category summary"), "No category rows available."))
  }

  if ("n_significant" %in% colnames(category_table)) {
    df <- category_table[order(-category_table$n_significant, -abs(category_table$mean_activity_score), category_table$category), , drop = FALSE]
    df <- utils::head(df, top_n)
    df$category <- stats::reorder(df$category, df$n_significant)

    return(
      ggplot2::ggplot(df, ggplot2::aes(x = n_significant, y = category, fill = mean_activity_score)) +
        ggplot2::geom_col(width = 0.72, color = "white", linewidth = 0.2) +
        ggplot2::geom_text(ggplot2::aes(label = n_significant), hjust = -0.15, size = 3) +
        ggplot2::scale_fill_gradient2(low = "#009E73", mid = "white", high = "#D55E00", midpoint = 0, na.value = "grey85") +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
        ggplot2::labs(
          title = .plot_if_null(title, "L1000PAI category summary"),
          subtitle = "Bars show significant result counts; fill shows mean activity score.",
          x = "Significant results",
          y = NULL,
          fill = "Mean activity"
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          legend.position = "bottom"
        )
    )
  }

  if (all(c("group_meta_z", "group_direction", "n_meta_sig", "category") %in% colnames(category_table))) {
    df <- category_table[!is.na(category_table$group_meta_z), , drop = FALSE]
    if (nrow(df) == 0 || all(df$n_meta_sig == 0, na.rm = TRUE)) {
      return(.blank_category_plot(
        .plot_if_null(title, "L1000PAI weighted category summary"),
        "No category-level support at the current activity threshold."
      ))
    }

    df <- df[order(df$group_meta_fdr, -abs(df$group_meta_z), -df$n_meta_sig, df$category), , drop = FALSE]
    df <- utils::head(df, top_n)
    df$category <- stats::reorder(df$category, df$group_meta_z)

    return(
      ggplot2::ggplot(
        df,
        ggplot2::aes(
          x = group_meta_z,
          y = category,
          size = pmax(n_meta_sig, 1),
          fill = group_direction
        )
      ) +
        ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, color = "grey82") +
        ggplot2::geom_point(shape = 21, color = "grey25", stroke = 0.8) +
        ggplot2::geom_text(ggplot2::aes(label = n_meta_sig), size = 3, fontface = "bold") +
        ggplot2::scale_fill_manual(
          values = c(Activated = "#D55E00", Inhibited = "#009E73", Unresolved = "grey75"),
          drop = FALSE
        ) +
        ggplot2::scale_size_continuous(range = c(3, 10)) +
        ggplot2::labs(
          title = .plot_if_null(title, "L1000PAI weighted category summary"),
          subtitle = "x = group meta-Z; point label = number of meta-significant entities.",
          x = "Group meta-Z",
          y = NULL,
          fill = "Direction",
          size = "Meta-sig entities"
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          legend.position = "bottom"
        )
    )
  }

  stop("`category_table` must contain either `n_significant` or weighted `group_meta_z` columns.")
}

#' Plot L1000PAI activity or annotation summaries
#'
#' Create a quick diagnostic plot from an `L1000PAIResult` object or compatible
#' result data.frame. `type = "activity"` shows perturbagen activity scores
#' against statistical support. `type = "category"` uses packaged annotation
#' summaries for `chem`, `ligand`, or `kinase_gpcr` results.
#'
#' @param x An `L1000PAIResult` object or result data.frame.
#' @param type One of `"activity"` or `"category"`.
#' @param library Optional library name when `x` is not an `L1000PAIResult`.
#' @param category_table Optional precomputed category table for
#'   `type = "category"`.
#' @param p_col Optional p-value column for `type = "activity"`.
#' @param score_col Activity-score column for `type = "activity"`. The
#'   default is `"activity_score"`; `"activity_gls_z"` can be used for joint
#'   GLS results.
#' @param significance_cutoff Horizontal reference line for activity plots.
#' @param activity_threshold Vertical reference lines for activity plots.
#' @param top_n Maximum rows/categories to draw.
#' @param label_top Number of top activity rows to label.
#' @param title Optional plot title.
#'
#' @return A `ggplot` object.
#' @export
plot_l1000_pai <- function(
    x,
    type = c("activity", "category"),
    library = NULL,
    category_table = NULL,
    p_col = NULL,
    significance_cutoff = 0.05,
    activity_threshold = 1,
    top_n = 200,
    label_top = 10,
    title = NULL,
    score_col = "activity_score"
) {
  type <- match.arg(type)

  if (type == "activity") {
    return(.plot_activity(
      x = x,
      library = library,
      p_col = p_col,
      score_col = score_col,
      significance_cutoff = significance_cutoff,
      activity_threshold = activity_threshold,
      top_n = top_n,
      label_top = label_top,
      title = title
    ))
  }

  .plot_category(
    x = x,
    library = library,
    category_table = category_table,
    top_n = top_n,
    title = title
  )
}
