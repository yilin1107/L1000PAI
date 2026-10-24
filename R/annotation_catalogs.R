# Internal cache for packaged annotation catalogs.
.l1000_annotation_cache <- new.env(parent = emptyenv())

.supported_annotation_catalogs <- function() {
  c("chem", "ligand", "kinase_gpcr")
}

.annotation_file_name <- function(library) {
  switch(
    library,
    chem = "chem_annotation_catalog.csv.gz",
    ligand = "ligand_annotation_catalog.csv.gz",
    kinase_gpcr = "kinase_gpcr_annotation_catalog.csv.gz"
  )
}

.annotation_catalog_path <- function(library) {
  file.path(
    system.file("extdata", "annotation_catalogs", package = "L1000PAI"),
    .annotation_file_name(library)
  )
}

.read_annotation_catalog <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

.resolve_annotation_library <- function(x, library = NULL) {
  if (!is.null(library)) {
    return(match.arg(library, .supported_annotation_catalogs()))
  }

  if (methods::is(x, "L1000PAIResult")) {
    inferred <- x@metadata$library
    if (is.null(inferred) || !inferred %in% .supported_annotation_catalogs()) {
      stop(
        "Annotation catalogs are currently available only for: ",
        paste(.supported_annotation_catalogs(), collapse = ", "),
        "."
      )
    }
    return(inferred)
  }

  stop("`library` must be provided when `x` is not an `L1000PAIResult` object.")
}

.as_result_df <- function(x) {
  if (methods::is(x, "L1000PAIResult")) {
    return(result_table(x))
  }
  if (is.data.frame(x)) {
    return(x)
  }
  stop("`x` must be an `L1000PAIResult` object or a data.frame.")
}

.merge_catalog <- function(df, catalog, key_col, catalog_key = key_col) {
  idx <- match(as.character(df[[key_col]]), as.character(catalog[[catalog_key]]))
  ann_cols <- setdiff(colnames(catalog), catalog_key)
  out <- cbind(df, catalog[idx, ann_cols, drop = FALSE])
  rownames(out) <- NULL
  out
}

.prepare_term_source <- function(df) {
  if ("perturbagen" %in% colnames(df)) {
    return(as.character(df$perturbagen))
  }
  if ("term" %in% colnames(df)) {
    return(strip_direction_suffix(df$term))
  }
  stop("The result table must contain either `perturbagen` or `term`.")
}

.extract_cell_line <- function(x) {
  has_cell <- grepl("\\([^)]*\\)$", x)
  out <- rep(NA_character_, length(x))
  out[has_cell] <- sub("^.*\\(([^)]*)\\)$", "\\1", x[has_cell])
  out
}

.prepare_chem_result <- function(df) {
  term_source <- .prepare_term_source(df)
  df$compound_name <- trimws(strip_direction_suffix(term_source))
  df
}

.prepare_ligand_result <- function(df) {
  term_source <- trimws(strip_direction_suffix(.prepare_term_source(df)))
  df$base_ligand <- trimws(sub("\\([^)]*\\)$", "", term_source))
  df$cell_line <- .extract_cell_line(term_source)
  df
}

.prepare_kinase_gpcr_result <- function(df) {
  term_source <- trimws(strip_direction_suffix(.prepare_term_source(df)))
  df$target <- trimws(sub("\\([^)]*\\)$", "", term_source))
  df$cell_line <- .extract_cell_line(term_source)
  df$FineCategory <- if ("TargetSubclass" %in% colnames(df)) df$TargetSubclass else NA_character_
  df
}

.safe_two_sided_z <- function(p_value, direction_sign) {
  if (is.na(p_value) || is.na(direction_sign) || direction_sign == 0) {
    return(NA_real_)
  }

  p_value <- min(max(as.numeric(p_value), 1e-300), 1 - 1e-16)
  direction_sign * stats::qnorm(1 - p_value / 2)
}

.safe_one_sided_z <- function(p_value, direction_sign) {
  if (is.na(p_value) || is.na(direction_sign) || direction_sign == 0) {
    return(NA_real_)
  }

  p_value <- min(max(as.numeric(p_value), 1e-300), 1 - 1e-16)
  direction_sign * stats::qnorm(1 - p_value)
}

.preferred_activity_score_col <- function(df) {
  if ("activity_gls_z" %in% colnames(df)) {
    return("activity_gls_z")
  }
  if ("activity_score" %in% colnames(df)) {
    return("activity_score")
  }
  stop("An `activity_gls_z` or `activity_score` column is required.")
}

.empty_category_summary <- function(category_level = character()) {
  data.frame(
    category_level = category_level,
    category = character(),
    n_results = integer(),
    n_significant = integer(),
    n_activated = integer(),
    n_inhibited = integer(),
    n_discordant = integer(),
    n_not_significant = integer(),
    activity_score_source = character(),
    mean_activity_score = numeric(),
    mean_abs_activity_score = numeric(),
    top_activated = character(),
    top_inhibited = character(),
    stringsAsFactors = FALSE
  )
}

.summarize_category_table <- function(df, category_col, id_col) {
  if (!category_col %in% colnames(df)) {
    stop("Category column `", category_col, "` was not found.")
  }

  x <- df[!is.na(df[[category_col]]) & nzchar(df[[category_col]]), , drop = FALSE]
  if (nrow(x) == 0) {
    return(.empty_category_summary(category_col))
  }

  categories <- sort(unique(as.character(x[[category_col]])))
  score_col <- .preferred_activity_score_col(x)
  out <- lapply(categories, function(cat_name) {
    sub_df <- x[as.character(x[[category_col]]) == cat_name, , drop = FALSE]
    sig <- if ("is_significant" %in% colnames(sub_df)) as.logical(sub_df$is_significant) else !is.na(sub_df$activity_class) & sub_df$activity_class != "Not_significant"
    act_class <- as.character(sub_df$activity_class)
    act_score <- suppressWarnings(as.numeric(sub_df[[score_col]]))

    data.frame(
      category_level = category_col,
      category = cat_name,
      n_results = nrow(sub_df),
      n_significant = sum(sig, na.rm = TRUE),
      n_activated = sum(act_class == "Activated", na.rm = TRUE),
      n_inhibited = sum(act_class == "Inhibited", na.rm = TRUE),
      n_discordant = sum(act_class == "Discordant", na.rm = TRUE),
      n_not_significant = sum(act_class == "Not_significant", na.rm = TRUE),
      activity_score_source = score_col,
      mean_activity_score = if (all(is.na(act_score))) NA_real_ else mean(act_score, na.rm = TRUE),
      mean_abs_activity_score = if (all(is.na(act_score))) NA_real_ else mean(abs(act_score), na.rm = TRUE),
      top_activated = paste(utils::head(
        sub_df[[id_col]][order(act_score, decreasing = TRUE)], 5
      ), collapse = "; "),
      top_inhibited = paste(utils::head(
        sub_df[[id_col]][order(act_score, decreasing = FALSE)], 5
      ), collapse = "; "),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(-out$n_significant, -abs(out$mean_activity_score), out$category), , drop = FALSE]
}

.weighted_entity_summary <- function(df, entity_col) {
  p_raw_col <- if ("p_coherent" %in% colnames(df)) {
    "p_coherent"
  } else if ("p_combined" %in% colnames(df)) {
    "p_combined"
  } else if ("pvalue" %in% colnames(df)) {
    "pvalue"
  } else {
    stop("A raw p-value column is required.")
  }
  p_adj_col <- if ("padj_coherent" %in% colnames(df)) {
    "padj_coherent"
  } else if ("padj_combined" %in% colnames(df)) {
    "padj_combined"
  } else if ("p.adjust" %in% colnames(df)) {
    "p.adjust"
  } else {
    p_raw_col
  }
  score_col <- .preferred_activity_score_col(df)

  x <- df
  act_class <- as.character(x$activity_class)
  x$term_direction_sign <- ifelse(
    act_class == "Activated", 1,
    ifelse(act_class == "Inhibited", -1, 0)
  )
  if (score_col == "activity_gls_z") {
    x$signed_term_z <- ifelse(
      x$term_direction_sign == 0,
      NA_real_,
      suppressWarnings(as.numeric(x[[score_col]]))
    )
  } else if (p_raw_col == "p_coherent") {
    x$signed_term_z <- mapply(
      .safe_one_sided_z,
      x[[p_raw_col]],
      x$term_direction_sign
    )
  } else {
    x$signed_term_z <- mapply(
      .safe_two_sided_z,
      x[[p_raw_col]],
      x$term_direction_sign
    )
  }
  x$is_significant_internal <- if ("is_significant" %in% colnames(x)) {
    as.logical(x$is_significant)
  } else {
    !is.na(x[[p_adj_col]]) & as.numeric(x[[p_adj_col]]) < 0.05
  }

  entities <- sort(unique(as.character(x[[entity_col]])))
  out <- lapply(entities, function(entity_name) {
    sub_df <- x[as.character(x[[entity_col]]) == entity_name, , drop = FALSE]
    act_class_sub <- as.character(sub_df$activity_class)
    n_activated_sig <- sum(sub_df$is_significant_internal & act_class_sub == "Activated", na.rm = TRUE)
    n_inhibited_sig <- sum(sub_df$is_significant_internal & act_class_sub == "Inhibited", na.rm = TRUE)
    n_discordant_sig <- sum(sub_df$is_significant_internal & act_class_sub == "Discordant", na.rm = TRUE)

    supporting_direction <- if (n_activated_sig > n_inhibited_sig) {
      "Activated"
    } else if (n_inhibited_sig > n_activated_sig) {
      "Inhibited"
    } else if (n_activated_sig == 0 && n_inhibited_sig == 0 && n_discordant_sig > 0) {
      "Discordant"
    } else if (n_activated_sig == n_inhibited_sig && n_activated_sig > 0) {
      "Discordant"
    } else {
      "No clear direction"
    }

    n_supporting <- if (supporting_direction == "Activated") {
      n_activated_sig
    } else if (supporting_direction == "Inhibited") {
      n_inhibited_sig
    } else {
      0L
    }

    n_opposing <- if (supporting_direction == "Activated") {
      n_inhibited_sig
    } else if (supporting_direction == "Inhibited") {
      n_activated_sig
    } else {
      0L
    }

    meta_z <- if (sum(!is.na(sub_df$signed_term_z)) > 0) {
      sum(sub_df$signed_term_z, na.rm = TRUE) / sqrt(sum(!is.na(sub_df$signed_term_z)))
    } else {
      NA_real_
    }

    data.frame(
      entity = entity_name,
      n_tested = nrow(sub_df),
      tested_contexts = paste(sort(unique(sub_df$cell_line)), collapse = "; "),
      n_non_discordant = sum(act_class_sub %in% c("Activated", "Inhibited"), na.rm = TRUE),
      n_sig = sum(sub_df$is_significant_internal, na.rm = TRUE),
      n_activated_sig = n_activated_sig,
      n_inhibited_sig = n_inhibited_sig,
      n_discordant_sig = n_discordant_sig,
      supporting_direction = supporting_direction,
      n_supporting = n_supporting,
      n_opposing = n_opposing,
      support_ratio = n_supporting / nrow(sub_df),
      direction_consistency = if ((n_activated_sig + n_inhibited_sig) > 0) {
        n_supporting / (n_activated_sig + n_inhibited_sig)
      } else {
        NA_real_
      },
      activity_score_source = score_col,
      weighted_mean_activity = if (all(is.na(sub_df[[score_col]]))) {
        NA_real_
      } else {
        mean(sub_df[[score_col]], na.rm = TRUE)
      },
      meta_z = meta_z,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  out$meta_p <- ifelse(is.na(out$meta_z), NA_real_, 2 * stats::pnorm(-abs(out$meta_z)))
  out$meta_fdr <- stats::p.adjust(ifelse(is.na(out$meta_p), 1, out$meta_p), method = "BH")
  out$meta_direction <- ifelse(
    is.na(out$meta_z), "Unresolved",
    ifelse(out$meta_z > 0, "Activated", ifelse(out$meta_z < 0, "Inhibited", "Unresolved"))
  )
  out$neg_log10_meta_fdr <- -log10(pmax(out$meta_fdr, 1e-300))
  rownames(out) <- NULL
  out[order(out$meta_fdr, -abs(out$meta_z), -out$n_supporting, out$entity), , drop = FALSE]
}

.empty_weighted_category_summary <- function(category_level = character()) {
  data.frame(
    category_level = category_level,
    category = character(),
    n_entities = integer(),
    n_meta_sig = integer(),
    n_meta_activated = integer(),
    n_meta_inhibited = integer(),
    mean_support_ratio = numeric(),
    mean_direction_consistency = numeric(),
    median_meta_z = numeric(),
    mean_abs_meta_z = numeric(),
    group_meta_z = numeric(),
    group_meta_p = numeric(),
    group_meta_fdr = numeric(),
    group_direction = character(),
    top_entities = character(),
    stringsAsFactors = FALSE
  )
}

.summarize_weighted_category_table <- function(entity_df, category_col, entity_col = "entity") {
  if (!category_col %in% colnames(entity_df)) {
    stop("Category column `", category_col, "` was not found.")
  }

  x <- entity_df[!is.na(entity_df[[category_col]]) & nzchar(entity_df[[category_col]]), , drop = FALSE]
  if (nrow(x) == 0) {
    return(.empty_weighted_category_summary(category_col))
  }

  categories <- sort(unique(as.character(x[[category_col]])))
  out <- lapply(categories, function(cat_name) {
    sub_df <- x[as.character(x[[category_col]]) == cat_name, , drop = FALSE]
    group_meta_z <- if (sum(!is.na(sub_df$meta_z)) > 0) {
      sum(sub_df$meta_z, na.rm = TRUE) / sqrt(sum(!is.na(sub_df$meta_z)))
    } else {
      NA_real_
    }

    data.frame(
      category_level = category_col,
      category = cat_name,
      n_entities = nrow(sub_df),
      n_meta_sig = sum(sub_df$meta_fdr < 0.05, na.rm = TRUE),
      n_meta_activated = sum(sub_df$meta_fdr < 0.05 & sub_df$meta_direction == "Activated", na.rm = TRUE),
      n_meta_inhibited = sum(sub_df$meta_fdr < 0.05 & sub_df$meta_direction == "Inhibited", na.rm = TRUE),
      mean_support_ratio = if (all(is.na(sub_df$support_ratio))) NA_real_ else mean(sub_df$support_ratio, na.rm = TRUE),
      mean_direction_consistency = if (all(is.na(sub_df$direction_consistency))) NA_real_ else mean(sub_df$direction_consistency, na.rm = TRUE),
      median_meta_z = if (all(is.na(sub_df$meta_z))) NA_real_ else stats::median(sub_df$meta_z, na.rm = TRUE),
      mean_abs_meta_z = if (all(is.na(sub_df$meta_z))) NA_real_ else mean(abs(sub_df$meta_z), na.rm = TRUE),
      group_meta_z = group_meta_z,
      top_entities = paste(
        utils::head(
          sub_df[[entity_col]][order(sub_df$meta_fdr, -abs(sub_df$meta_z), -sub_df$n_supporting)],
          5
        ),
        collapse = "; "
      ),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  out$group_meta_p <- ifelse(is.na(out$group_meta_z), NA_real_, 2 * stats::pnorm(-abs(out$group_meta_z)))
  out$group_meta_fdr <- stats::p.adjust(ifelse(is.na(out$group_meta_p), 1, out$group_meta_p), method = "BH")
  out$group_direction <- ifelse(
    is.na(out$group_meta_z), "Unresolved",
    ifelse(out$group_meta_z > 0, "Activated", ifelse(out$group_meta_z < 0, "Inhibited", "Unresolved"))
  )
  rownames(out) <- NULL
  out[order(out$group_meta_fdr, -abs(out$group_meta_z), -out$n_meta_sig, out$category), , drop = FALSE]
}

.build_chem_annotation_tables <- function(annotated_result) {
  list(
    annotated_result = annotated_result,
    big_class_summary = .summarize_category_table(annotated_result, "BigClass", "compound_name"),
    mechanism_class_summary = .summarize_category_table(annotated_result, "MechanismClass", "compound_name"),
    target_family_summary = .summarize_category_table(annotated_result, "TargetFamily", "compound_name")
  )
}

.build_ligand_annotation_tables <- function(annotated_result) {
  base_summary <- .weighted_entity_summary(annotated_result, "base_ligand")
  ligand_catalog <- get_l1000_annotation_catalog("ligand")
  ligand_catalog_sub <- ligand_catalog[, setdiff(colnames(ligand_catalog), c("n_terms", "n_cell_lines", "n_exact_ligand_entries")), drop = FALSE]
  names(base_summary)[names(base_summary) == "entity"] <- "base_ligand"
  base_summary <- .merge_catalog(base_summary, ligand_catalog_sub, "base_ligand")

  list(
    annotated_result = annotated_result,
    base_ligand_summary = base_summary,
    ligand_group_summary = .summarize_weighted_category_table(base_summary, "LigandGroup", "base_ligand"),
    functional_theme_summary = .summarize_weighted_category_table(base_summary, "FunctionalTheme", "base_ligand")
  )
}

.build_kinase_gpcr_annotation_tables <- function(annotated_result) {
  target_summary <- .weighted_entity_summary(annotated_result, "target")
  kg_catalog <- get_l1000_annotation_catalog("kinase_gpcr")
  kg_catalog_sub <- kg_catalog[, setdiff(colnames(kg_catalog), c("n_terms", "n_cell_lines")), drop = FALSE]
  names(target_summary)[names(target_summary) == "entity"] <- "target"
  target_summary <- .merge_catalog(target_summary, kg_catalog_sub, "target")
  target_summary$FineCategory <- ifelse(
    target_summary$TargetClass == "Kinase / kinase-like" & !is.na(target_summary$KinaseFamily) & nzchar(target_summary$KinaseFamily),
    target_summary$KinaseFamily,
    ifelse(
      target_summary$TargetClass == "GPCR" & !is.na(target_summary$GPCRFamily) & nzchar(target_summary$GPCRFamily),
      target_summary$GPCRFamily,
      ifelse(!is.na(target_summary$TargetSubclass) & nzchar(target_summary$TargetSubclass), target_summary$TargetSubclass, target_summary$TargetClass)
    )
  )

  list(
    annotated_result = annotated_result,
    target_summary = target_summary,
    signal_family_summary = .summarize_weighted_category_table(target_summary, "SignalFamily", "target"),
    target_class_summary = .summarize_weighted_category_table(target_summary, "TargetClass", "target"),
    fine_category_summary = .summarize_weighted_category_table(target_summary, "FineCategory", "target")
  )
}

#' Available packaged annotation catalogs
#'
#' @return A character vector of annotation catalogs bundled in `L1000PAI`.
#' @export
available_l1000_annotation_catalogs <- function() {
  .supported_annotation_catalogs()
}

#' Get a packaged annotation catalog for an L1000 library
#'
#' Load the slim packaged annotation catalog used to annotate bundled `chem`,
#' `ligand`, or `kinase_gpcr` results.
#'
#' @param library One of `available_l1000_annotation_catalogs()`.
#'
#' @return A data.frame annotation catalog.
#' @export
get_l1000_annotation_catalog <- function(
    library = c("chem", "ligand", "kinase_gpcr")
) {
  library <- match.arg(library)

  if (exists(library, envir = .l1000_annotation_cache, inherits = FALSE)) {
    return(get(library, envir = .l1000_annotation_cache, inherits = FALSE))
  }

  path <- .annotation_catalog_path(library)
  if (!file.exists(path)) {
    stop(
      "Packaged annotation catalog for `", library, "` was not found. ",
      "Expected file: ", path
    )
  }

  catalog <- .read_annotation_catalog(path)
  assign(library, catalog, envir = .l1000_annotation_cache)
  catalog
}

#' Annotate an L1000 result table with packaged library classifications
#'
#' Add packaged `chem`, `ligand`, or `kinase_gpcr` annotations to an
#' `L1000PAIResult` object or result data.frame.
#'
#' @param x An `L1000PAIResult` object or a result data.frame.
#' @param library Optional library name when `x` is not an `L1000PAIResult`
#'   object.
#'
#' @return A data.frame with annotation columns appended.
#' @export
annotate_l1000_result <- function(x, library = NULL) {
  library <- .resolve_annotation_library(x, library = library)
  df <- .as_result_df(x)
  catalog <- get_l1000_annotation_catalog(library)

  if (library == "chem") {
    df <- .prepare_chem_result(df)
    return(.merge_catalog(df, catalog, "compound_name"))
  }

  if (library == "ligand") {
    df <- .prepare_ligand_result(df)
    return(.merge_catalog(df, catalog, "base_ligand"))
  }

  df <- .prepare_kinase_gpcr_result(df)
  .merge_catalog(df, catalog, "target")
}

#' Build packaged annotation tables from an L1000 result
#'
#' Create library-specific annotation tables directly from a result object.
#'
#' @param x An `L1000PAIResult` object or a result data.frame.
#' @param library Optional library name when `x` is not an `L1000PAIResult`
#'   object.
#'
#' @return A named list of data.frames.
#' @export
build_l1000_annotation_tables <- function(x, library = NULL) {
  library <- .resolve_annotation_library(x, library = library)
  annotated_result <- annotate_l1000_result(x, library = library)

  switch(
    library,
    chem = .build_chem_annotation_tables(annotated_result),
    ligand = .build_ligand_annotation_tables(annotated_result),
    kinase_gpcr = .build_kinase_gpcr_annotation_tables(annotated_result)
  )
}

#' Extract annotation tables from an `L1000PAIResult`
#'
#' If annotation tables were not precomputed during `run_l1000_pai()`, they are
#' generated on demand for supported libraries.
#'
#' @param x An `L1000PAIResult` object.
#'
#' @return A named list of data.frames, or `NULL` for unsupported libraries.
#' @export
annotation_tables <- function(x) {
  if (!methods::is(x, "L1000PAIResult")) {
    stop("`x` must be an L1000PAIResult object.")
  }

  cached <- x@metadata$annotation_tables
  if (!is.null(cached)) {
    return(cached)
  }

  library <- x@metadata$library
  if (is.null(library) || !library %in% .supported_annotation_catalogs()) {
    return(NULL)
  }

  build_l1000_annotation_tables(x, library = library)
}

#' Extract the annotated result table
#'
#' @param x An `L1000PAIResult` object.
#'
#' @return A data.frame with packaged annotations appended, or `NULL` when no
#'   packaged annotation catalog is available for the result library.
#' @export
annotated_result_table <- function(x) {
  tables <- annotation_tables(x)
  if (is.null(tables)) {
    return(NULL)
  }
  tables$annotated_result
}

#' Write packaged annotation tables to disk
#'
#' Export the library-specific annotation tables returned by
#' [build_l1000_annotation_tables()] as CSV files.
#'
#' @param x An `L1000PAIResult` object or a result data.frame.
#' @param out_dir Output directory for CSV files.
#' @param library Optional library name when `x` is not an `L1000PAIResult`
#'   object.
#' @param prefix Optional filename prefix. Defaults to the library name.
#'
#' @return Invisibly returns the named list of tables that were written.
#' @export
write_l1000_annotation_tables <- function(
    x,
    out_dir,
    library = NULL,
    prefix = NULL
) {
  library <- .resolve_annotation_library(x, library = library)
  if (!is.character(out_dir) || length(out_dir) != 1 || !nzchar(out_dir)) {
    stop("`out_dir` must be a single non-empty character string.")
  }

  tables <- build_l1000_annotation_tables(x, library = library)
  if (is.null(prefix)) {
    prefix <- library
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(tables)) {
    utils::write.csv(
      tables[[name]],
      file = file.path(out_dir, paste0(prefix, "_", name, ".csv")),
      row.names = FALSE,
      na = ""
    )
  }

  invisible(tables)
}
