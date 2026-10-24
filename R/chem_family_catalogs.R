# Internal cache for the versioned Chem-family taxonomy.
.chem_family_catalog_cache <- new.env(parent = emptyenv())

.supported_chem_family_versions <- function() "2.0.0"

.chem_family_component_files <- function() {
  c(
    term_registry = "chem_term_registry_v2.csv",
    family_registry = "chem_family_registry_v2.csv",
    membership = "chem_family_membership_v2.csv",
    source_registry = "chem_source_registry_v2.csv",
    manifest = "chem_taxonomy_v2_manifest.csv"
  )
}

.chem_family_catalog_path <- function(component, version) {
  if (version != "2.0.0") {
    stop("Unsupported Chem-family taxonomy version: ", version)
  }
  file.path(
    system.file(
      "extdata", "annotation_catalogs", "v2",
      package = "L1000PAI"
    ),
    unname(.chem_family_component_files()[[component]])
  )
}

#' Available Chem-family taxonomy versions
#'
#' @return A character vector of packaged Chem-family taxonomy versions.
#' @export
available_chem_family_versions <- function() {
  .supported_chem_family_versions()
}

#' Get a versioned Chem-family catalog component
#'
#' The v2 Chem-family taxonomy is separate from the legacy one-row-per-compound
#' annotation catalog returned by [get_l1000_annotation_catalog()]. It stores a
#' complete exact-term registry plus curated many-to-many family memberships.
#'
#' @param component One of `term_registry`, `family_registry`, `membership`,
#'   `source_registry`, or `manifest`.
#' @param version Taxonomy version. Currently `2.0.0`.
#'
#' @return A data.frame containing the requested frozen taxonomy component.
#' @export
get_chem_family_catalog <- function(
    component = c(
      "term_registry", "family_registry", "membership",
      "source_registry", "manifest"
    ),
    version = "2.0.0"
) {
  component <- match.arg(component)
  version <- match.arg(version, .supported_chem_family_versions())
  cache_key <- paste(version, component, sep = "::")

  if (exists(cache_key, envir = .chem_family_catalog_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .chem_family_catalog_cache, inherits = FALSE))
  }

  path <- .chem_family_catalog_path(component, version)
  if (!file.exists(path)) {
    stop(
      "Packaged Chem-family component `", component,
      "` was not found for version ", version, ". Expected file: ", path
    )
  }

  out <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = ""
  )
  assign(cache_key, out, envir = .chem_family_catalog_cache)
  out
}

.empty_chem_membership_row <- function(membership) {
  out <- membership[0, , drop = FALSE]
  out[1, ] <- NA
  out
}

#' Annotate Chem results with versioned many-to-many families
#'
#' Expands each Chem result into long-format family memberships without
#' silently choosing the first family. No activity aggregation or enrichment is
#' performed. Use `analysis_unit_id` to enforce one frozen compound vote during
#' downstream family analysis.
#'
#' @param x An `L1000PAIResult` object or result data.frame containing a
#'   `perturbagen` or `term` column.
#' @param family_level One of `direct_moa`, `pathway_context`, or `all`.
#' @param evidence_tier Membership evidence tiers to retain. Defaults to A and B.
#' @param include_unmapped If `TRUE`, retain input rows without a membership at
#'   the requested level.
#' @param version Chem-family taxonomy version.
#'
#' @return A long-format data.frame. A result can occur more than once when its
#'   analysis unit belongs to multiple families.
#' @export
annotate_chem_families <- function(
    x,
    family_level = c("direct_moa", "pathway_context", "all"),
    evidence_tier = c("A", "B"),
    include_unmapped = TRUE,
    version = "2.0.0"
) {
  family_level <- match.arg(family_level)
  version <- match.arg(version, .supported_chem_family_versions())
  allowed_tiers <- c("A", "B", "C", "U")
  if (!is.character(evidence_tier) || length(evidence_tier) == 0L ||
      !all(evidence_tier %in% allowed_tiers)) {
    stop("`evidence_tier` must contain one or more of: A, B, C, U.")
  }
  if (!is.logical(include_unmapped) || length(include_unmapped) != 1L ||
      is.na(include_unmapped)) {
    stop("`include_unmapped` must be TRUE or FALSE.")
  }

  df <- .prepare_chem_result(.as_result_df(x))
  df$result_row_id <- seq_len(nrow(df))

  term_registry <- get_chem_family_catalog("term_registry", version)
  membership <- get_chem_family_catalog("membership", version)
  family_registry <- get_chem_family_catalog("family_registry", version)
  source_registry <- get_chem_family_catalog("source_registry", version)

  term_idx <- match(df$compound_name, term_registry$library_term)
  term_ann <- term_registry[term_idx, setdiff(names(term_registry), "library_term"), drop = FALSE]
  names(term_ann)[names(term_ann) == "taxonomy_version"] <- "term_taxonomy_version"
  base <- cbind(df, library_term = df$compound_name, term_ann)
  rownames(base) <- NULL

  family_level_lookup <- stats::setNames(
    family_registry$family_level,
    family_registry$family_id
  )
  membership$family_level_internal <- unname(family_level_lookup[membership$family_id])
  membership <- membership[
    membership$evidence_tier %in% evidence_tier &
      (family_level == "all" | membership$family_level_internal == family_level),
    ,
    drop = FALSE
  ]
  names(membership)[names(membership) == "taxonomy_version"] <- "membership_taxonomy_version"

  membership_index <- split(
    seq_len(nrow(membership)),
    membership$analysis_unit_id,
    drop = TRUE
  )
  membership_output <- membership[
    ,
    setdiff(names(membership), c("analysis_unit_id", "family_level_internal")),
    drop = FALSE
  ]
  names(membership_output)[names(membership_output) == "curation_status"] <-
    "membership_curation_status"
  empty_membership <- .empty_chem_membership_row(membership_output)

  expanded <- lapply(seq_len(nrow(base)), function(i) {
    unit <- base$analysis_unit_id[[i]]
    idx <- if (is.na(unit)) integer() else membership_index[[unit]]
    if (is.null(idx) || length(idx) == 0L) {
      if (!include_unmapped) return(NULL)
      mem <- empty_membership
    } else {
      mem <- membership_output[idx, , drop = FALSE]
    }
    cbind(base[rep(i, nrow(mem)), , drop = FALSE], mem)
  })
  expanded <- Filter(Negate(is.null), expanded)
  if (length(expanded) == 0L) {
    return(data.frame())
  }
  out <- do.call(rbind, expanded)
  rownames(out) <- NULL

  fam_idx <- match(out$family_id, family_registry$family_id)
  fam_ann <- family_registry[
    fam_idx,
    setdiff(names(family_registry), c("family_id", "source_id", "taxonomy_version")),
    drop = FALSE
  ]
  names(fam_ann)[names(fam_ann) == "family_level"] <- "family_level"
  out <- cbind(out, fam_ann)

  source_idx <- match(out$source_id, source_registry$source_id)
  out$membership_source_citation <- source_registry$citation[source_idx]
  out$membership_source_url <- source_registry$url[source_idx]
  out$mapping_status <- ifelse(
    is.na(out$library_term_id),
    "unmapped_library_term",
    ifelse(is.na(out$family_id), "no_membership_at_requested_level", "matched")
  )

  out <- out[order(out$result_row_id, out$family_level, out$family_id, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.chem_family_forbidden_columns <- function(x) {
  pattern <- paste0(
    "(^|_)(activity_score|activity_gls_z|activity_class|nes(_up|_down)?|",
    "pvalue|p_value|p_combined|padj.*|fdr|qvalue|neg_log10.*|evaluated|",
    "recovered|recovery.*|is_significant|n_recovered.*|best_recovered.*)(_|$)"
  )
  grep(pattern, names(x), ignore.case = TRUE, perl = TRUE, value = TRUE)
}

#' Verify the packaged Chem-family taxonomy
#'
#' Runs schema, key, foreign-key, pair-coverage, alias-collapse and outcome-
#' firewall checks. SHA-256 verification is optional because it requires the
#' suggested `digest` package.
#'
#' @param verify_hash Recompute packaged-table SHA-256 hashes.
#' @param stop_on_error Stop when any check fails.
#' @param version Chem-family taxonomy version.
#'
#' @return A data.frame containing one row per integrity check.
#' @export
verify_chem_family_catalogs <- function(
    verify_hash = FALSE,
    stop_on_error = TRUE,
    version = "2.0.0"
) {
  version <- match.arg(version, .supported_chem_family_versions())
  term <- get_chem_family_catalog("term_registry", version)
  family <- get_chem_family_catalog("family_registry", version)
  membership <- get_chem_family_catalog("membership", version)
  source <- get_chem_family_catalog("source_registry", version)
  manifest <- get_chem_family_catalog("manifest", version)

  checks <- list()
  add_check <- function(id, pass, detail) {
    checks[[length(checks) + 1L]] <<- data.frame(
      check_id = id,
      pass = isTRUE(pass),
      detail = as.character(detail),
      stringsAsFactors = FALSE
    )
  }

  add_check("term_count", nrow(term) == 5425L, paste("rows=", nrow(term)))
  add_check("term_id_unique", anyDuplicated(term$library_term_id) == 0L, "library_term_id")
  add_check("library_term_unique", anyDuplicated(term$library_term) == 0L, "exact library_term")
  add_check("family_id_unique", anyDuplicated(family$family_id) == 0L, "family_id")
  membership_key <- paste(membership$analysis_unit_id, membership$family_id, sep = "\r")
  add_check("membership_key_unique", anyDuplicated(membership_key) == 0L, "analysis_unit_id + family_id")
  add_check("pair_complete", all(term$pair_complete & term$up_present & term$dn_present), paste("complete=", sum(term$pair_complete)))
  add_check("membership_unit_fk", all(membership$analysis_unit_id %in% term$analysis_unit_id), "analysis_unit_id")
  add_check("membership_family_fk", all(membership$family_id %in% family$family_id), "family_id")
  add_check("membership_source_fk", all(membership$source_id %in% source$source_id), "source_id")
  preferred_n <- tapply(as.integer(term$preferred_for_analysis), term$analysis_unit_id, sum)
  add_check("one_preferred_term_per_unit", all(preferred_n == 1L), paste("units=", length(preferred_n)))
  forbidden <- unique(unlist(lapply(list(term, family, membership, source), .chem_family_forbidden_columns)))
  add_check("outcome_firewall", length(forbidden) == 0L, if (length(forbidden)) paste(forbidden, collapse = ";") else "no outcome-derived columns")
  add_check("curated_term_count", sum(term$curation_status == "mechanism_curated_v2") == 80L, paste("curated=", sum(term$curation_status == "mechanism_curated_v2")))
  add_check("plx4720_present", "PLX-4720" %in% term$library_term, "PLX-4720 exact paired term")
  alias_units <- stats::setNames(term$analysis_unit_id, term$library_term)
  add_check("tanespimycin_alias", identical(alias_units[["17-AAG"]], alias_units[["Tanespimycin"]]), "17-AAG/Tanespimycin")
  add_check("ixazomib_program", identical(alias_units[["Ixazomib"]], alias_units[["MLN-9708"]]), "Ixazomib/MLN-9708")

  if (isTRUE(verify_hash)) {
    if (!requireNamespace("digest", quietly = TRUE)) {
      add_check("sha256", FALSE, "suggested package `digest` is unavailable")
    } else {
      frozen <- manifest[manifest$artifact_role == "frozen_taxonomy_table", , drop = FALSE]
      observed <- vapply(frozen$artifact_id, function(id) {
        component <- switch(id,
          chem_term_registry_v2 = "term_registry",
          chem_family_registry_v2 = "family_registry",
          chem_family_membership_v2 = "membership",
          chem_source_registry_v2 = "source_registry"
        )
        tolower(digest::digest(
          file = .chem_family_catalog_path(component, version),
          algo = "sha256",
          serialize = FALSE
        ))
      }, character(1))
      add_check("sha256", identical(unname(observed), tolower(frozen$sha256)), paste("files=", length(observed)))
    }
  }

  out <- do.call(rbind, checks)
  rownames(out) <- NULL
  if (isTRUE(stop_on_error) && any(!out$pass)) {
    failed <- out$check_id[!out$pass]
    stop("Chem-family taxonomy verification failed: ", paste(failed, collapse = ", "))
  }
  out
}
