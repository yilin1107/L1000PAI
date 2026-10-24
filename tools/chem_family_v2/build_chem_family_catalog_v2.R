taxonomy_version <- "2.0.0"
schema_version <- "1.0.0"
lock_date <- "2026-08-26"
lock_utc <- "2026-08-26T00:00:00Z"

`%||%` <- function(x, y) if (is.null(x) || !nzchar(x)) y else x

script_path <- function() {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg) == 0L) return(normalizePath(".", winslash = "/", mustWork = TRUE))
  normalizePath(sub("^--file=", "", arg[[1]]), winslash = "/", mustWork = TRUE)
}

arg_value <- function(flag, default = "") {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) default else args[[hit + 1L]]
}

repo_dir <- normalizePath(
  file.path(dirname(script_path()), "..", ".."),
  winslash = "/",
  mustWork = TRUE
)

default_library <- file.path(
  dirname(repo_dir), "L1000PAIdata", "inst", "extdata",
  "LINCS_L1000_Chem_Pert_Consensus_Sigs.txt"
)
chem_library <- arg_value(
  "--chem-library",
  Sys.getenv("L1000PAI_CHEM_LIBRARY", unset = default_library)
)
chem_library <- normalizePath(chem_library, winslash = "/", mustWork = TRUE)

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("The build requires the suggested package `digest` for SHA-256 manifests.")
}

sha256_file <- function(path) {
  tolower(digest::digest(file = path, algo = "sha256", serialize = FALSE))
}

slug <- function(x) {
  out <- toupper(gsub("[^A-Za-z0-9]+", "_", x))
  gsub("(^_+|_+$)", "", out)
}

norm_term <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))

read_chem_pairs <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  fields <- strsplit(lines, "\t", fixed = TRUE)
  signature_term <- vapply(fields, function(x) x[[1]], character(1))
  direction <- ifelse(
    grepl(" up$", signature_term), "up",
    ifelse(grepl(" dn$", signature_term), "dn", NA_character_)
  )
  if (anyNA(direction)) stop("Every Chem signature must end in ` up` or ` dn`.")
  library_term <- sub(" (up|dn)$", "", signature_term)
  n_genes <- vapply(
    fields,
    function(x) sum(nzchar(x[-1L])),
    integer(1)
  )
  key <- paste(library_term, direction, sep = "\r")
  if (anyDuplicated(key)) stop("Duplicate directional Chem signatures were found.")

  terms <- sort(unique(library_term), method = "radix")
  up_key <- setNames(n_genes[direction == "up"], library_term[direction == "up"])
  dn_key <- setNames(n_genes[direction == "dn"], library_term[direction == "dn"])

  data.frame(
    library_term_id = sprintf("CHEM_V2_%05d", seq_along(terms)),
    library_term = terms,
    library_term_norm = norm_term(terms),
    canonical_compound_id = sprintf("CHEM_V2_%05d", seq_along(terms)),
    canonical_compound_name = terms,
    analysis_unit_id = sprintf("CHEM_V2_%05d", seq_along(terms)),
    analysis_unit_name = terms,
    term_relation = "unreviewed_library_term",
    preferred_for_analysis = TRUE,
    normalization_collision_resolution = "none",
    identity_evidence_tier = "U",
    identity_source_id = "SRC_L1000_LIBRARY_TERM",
    curation_status = "unreviewed",
    multi_target_flag = NA,
    tool_compound_flag = NA,
    up_present = terms %in% names(up_key),
    dn_present = terms %in% names(dn_key),
    n_up_genes = unname(up_key[terms]),
    n_dn_genes = unname(dn_key[terms]),
    pair_complete = terms %in% names(up_key) & terms %in% names(dn_key),
    taxonomy_version = taxonomy_version,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

source_registry <- data.frame(
  source_id = c(
    "SRC_L1000_LIBRARY_TERM", "SRC_CURATION_POLICY",
    "SRC_NCI_IDASANUTLIN", "SRC_NCI_ERLOTINIB", "SRC_NCI_SELUMETINIB",
    "SRC_NCI_VEMURAFENIB", "SRC_NCI_SIROLIMUS", "SRC_PUBMED_AZD8055",
    "SRC_PUBMED_NVP_BEZ235", "SRC_NCI_ZSTK474", "SRC_NCI_TANESPIMYCIN",
    "SRC_NCI_IXAZOMIB", "SRC_PMC_ONX0914", "SRC_NCI_ALISERTIB",
    "SRC_NCI_DOCETAXEL", "SRC_PUBMED_NOCODAZOLE", "SRC_PUBCHEM_EQUILIN",
    "SRC_FDA_STATIN", "SRC_NCI_VANDETANIB", "SRC_NCI_PONATINIB"
  ),
  citation = c(
    "LINCS L1000 Chem Perturbation Consensus Signatures (local provider artifact)",
    "L1000PAI Chem family v2 curation policy",
    "NCI Drug Dictionary: idasanutlin", "NCI Drug Dictionary: erlotinib hydrochloride",
    "NCI Drug Dictionary: selumetinib", "NCI Drug Dictionary: vemurafenib",
    "NCI Drug Dictionary: sirolimus", "Chresta et al. Cancer Research 2010: AZD8055",
    "Maira et al. Molecular Cancer Therapeutics 2008: NVP-BEZ235",
    "NCI Drug Dictionary: PI3K inhibitor ZSTK474", "NCI Drug Dictionary: tanespimycin",
    "NCI Drug Dictionary: ixazomib citrate", "ONX-0914 immunoproteasome review",
    "NCI Drug Dictionary: alisertib", "NCI Drug Dictionary: docetaxel",
    "De Brabander et al. Cell Biology International Reports 1976: nocodazole",
    "PubChem Compound Summary: Equilin", "U.S. FDA cholesterol medicines guide",
    "NCI Drug Dictionary: vandetanib", "NCI Drug Dictionary: ponatinib hydrochloride"
  ),
  url = c(
    "", "",
    "https://www.cancer.gov/publications/dictionaries/cancer-drug/def/idasanutlin",
    "https://www.cancer.gov/publications/dictionaries/cancer-drug/def/erlotinib-hydrochloride",
    "https://www.cancer.gov/publications/dictionaries/cancer-drug/def/selumetinib",
    "https://www.cancer.gov/publications/dictionaries/cancer-drug/def/vemurafenib",
    "https://www.cancer.gov/publications/dictionaries/cancer-drug/def/sirolimus",
    "https://pubmed.ncbi.nlm.nih.gov/20028854/",
    "https://pubmed.ncbi.nlm.nih.gov/18829560/",
    "https://www.cancer.gov/publications/dictionaries/cancer-drug/def/pi3k-inhibitor-zstk474",
    "https://www.cancer.gov/publications/dictionaries/cancer-drug/def/tanespimycin",
    "https://www.cancer.gov/publications/dictionaries/cancer-drug/def/ixazomib-citrate",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC9688407/",
    "https://www.cancer.gov/publications/dictionaries/cancer-drug/def/alisertib",
    "https://www.cancer.gov/publications/dictionaries/cancer-drug/def/docetaxel",
    "https://pubmed.ncbi.nlm.nih.gov/85702/",
    "https://pubchem.ncbi.nlm.nih.gov/compound/Equilin",
    "https://www.fda.gov/consumers/womens-health-topics/cholesterol-medicines-guide",
    "https://www.cancer.gov/publications/dictionaries/cancer-drug/def/vandetanib",
    "https://www.cancer.gov/publications/dictionaries/cancer-drug/def/ponatinib-hydrochloride"
  ),
  source_type = c(
    "local_signature_library", "curation_policy",
    rep("government_drug_dictionary", 5), "primary_literature", "primary_literature",
    "government_drug_dictionary", "government_drug_dictionary",
    "government_drug_dictionary", "peer_reviewed_review",
    "government_drug_dictionary", "government_drug_dictionary", "primary_literature",
    "government_database", "regulatory_guidance", "government_drug_dictionary",
    "government_drug_dictionary"
  ),
  accessed_date = lock_date,
  redistribution_status = c(
    "provider_artifact_not_redistributed", "project_authored",
    rep("link_only", 18)
  ),
  taxonomy_version = taxonomy_version,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

family_row <- function(id, label, level, parent, definition, inclusion, exclusion,
                       role, tier, source_id) {
  data.frame(
    family_id = id,
    family_label = label,
    family_level = level,
    parent_family_id = parent,
    definition = definition,
    inclusion_rule = inclusion,
    exclusion_rule = exclusion,
    polarity_model = "positive Chem score denotes concordance with compound-induced transcription",
    analysis_role = role,
    minimum_evidence_tier = tier,
    source_id = source_id,
    taxonomy_version = taxonomy_version,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

family_registry <- do.call(rbind, list(
  family_row("MDM2_ANTAGONIST", "MDM2 antagonists", "direct_moa", "P53_CHECKPOINT_RESPONSE", "Compounds that disrupt the MDM2-p53 interaction.", "Direct MDM2 antagonist with documented p53-reactivating pharmacology.", "Downstream p53 modulators without direct MDM2 antagonism.", "strict_primary", "A", "SRC_NCI_IDASANUTLIN"),
  family_row("ERBB_TKI_CORE", "ERBB-family tyrosine kinase inhibitors", "direct_moa", "ERBB_RTK_AXIS", "Compounds whose primary pharmacology directly inhibits EGFR/ERBB-family kinases.", "Direct EGFR, dual ERBB, or pan-ERBB kinase inhibitor.", "Downstream MEK inhibitors and broad multikinase agents without core ERBB selectivity.", "strict_primary", "A", "SRC_NCI_ERLOTINIB"),
  family_row("MEK1_2_INHIBITOR", "MEK1/2 inhibitors", "direct_moa", "RAS_RAF_MEK_ERK_AXIS", "Compounds that directly inhibit MEK1 and/or MEK2.", "Direct MEK1/2 pharmacology.", "Upstream EGFR/ERBB or RAF inhibition.", "strict_primary", "A", "SRC_NCI_SELUMETINIB"),
  family_row("SELECTIVE_BRAF_INHIBITOR", "Selective BRAF/RAF inhibitors", "direct_moa", "RAS_RAF_MEK_ERK_AXIS", "Selective or primary RAF inhibitors used to suppress oncogenic RAF signaling.", "Direct RAF inhibitor with RAF as the primary intended target.", "RAF-containing multikinase drugs, downstream MEK inhibitors, and BCR::ABL1 inhibitors.", "strict_primary", "A", "SRC_NCI_VEMURAFENIB"),
  family_row("RAPALOG_FKBP12_MTORC1", "Rapalog FKBP12-mTORC1 inhibitors", "direct_moa", "PI3K_AKT_MTOR_AXIS", "Rapamycin and analogues that inhibit mTORC1 through FKBP12.", "Documented rapalog mechanism.", "Catalytic mTOR, PI3K, or dual PI3K-mTOR inhibitors.", "strict_primary", "A", "SRC_NCI_SIROLIMUS"),
  family_row("MTOR_KINASE_INHIBITOR", "Catalytic mTOR kinase inhibitors", "direct_moa", "PI3K_AKT_MTOR_AXIS", "ATP-competitive inhibitors of mTOR kinase complexes.", "Direct catalytic mTORC1/2 inhibition.", "Allosteric rapalogs and PI3K-only inhibitors.", "strict_primary", "A", "SRC_PUBMED_AZD8055"),
  family_row("DUAL_PI3K_MTOR_INHIBITOR", "Dual PI3K-mTOR inhibitors", "direct_moa", "PI3K_AKT_MTOR_AXIS", "Compounds that directly inhibit both PI3K and mTOR catalytic activity.", "Documented dual PI3K and mTOR pharmacology.", "PI3K-only, mTOR-only, or allosteric rapalog pharmacology.", "strict_primary", "A", "SRC_PUBMED_NVP_BEZ235"),
  family_row("CLASS_I_PI3K_INHIBITOR", "Class I PI3K inhibitors", "direct_moa", "PI3K_AKT_MTOR_AXIS", "Compounds whose primary pharmacology inhibits class I PI3K.", "Direct class I PI3K inhibition.", "Dual PI3K-mTOR and mTOR-only inhibition.", "strict_primary", "A", "SRC_NCI_ZSTK474"),
  family_row("HSP90_INHIBITOR", "HSP90 inhibitors", "direct_moa", "HSP90_CHAPERONE_PROTEOSTASIS", "Compounds that directly inhibit HSP90 chaperone function.", "Direct HSP90 inhibitor.", "General proteostasis perturbagens without HSP90 pharmacology.", "strict_primary", "A", "SRC_NCI_TANESPIMYCIN"),
  family_row("CONSTITUTIVE_PROTEASOME_INHIBITOR", "Constitutive proteasome inhibitors", "direct_moa", "UBIQUITIN_PROTEASOME_AXIS", "Compounds targeting constitutive proteasome catalytic activity, including PSMB5-directed agents.", "Clinically established or selective constitutive proteasome inhibitor.", "Immunoproteasome-selective inhibitors and broad protease tools.", "strict_primary", "A", "SRC_NCI_IXAZOMIB"),
  family_row("IMMUNOPROTEASOME_PSMB8_INHIBITOR", "Immunoproteasome PSMB8 inhibitors", "direct_moa", "UBIQUITIN_PROTEASOME_AXIS", "Compounds selective for immunoproteasome PSMB8/LMP7.", "Documented PSMB8-selective inhibition.", "Constitutive proteasome inhibitors and broad protease tools.", "strict_primary", "B", "SRC_PMC_ONX0914"),
  family_row("AURORA_KINASE_INHIBITOR", "Aurora kinase inhibitors", "direct_moa", "MITOTIC_SPINDLE_PERTURBATION", "Compounds directly inhibiting one or more Aurora kinases.", "Direct Aurora-family kinase inhibition.", "Microtubule-directed agents without Aurora inhibition.", "strict_primary", "A", "SRC_NCI_ALISERTIB"),
  family_row("MICROTUBULE_STABILIZER", "Microtubule stabilizers", "direct_moa", "MITOTIC_SPINDLE_PERTURBATION", "Compounds that bind tubulin and stabilize microtubules.", "Direct microtubule stabilization.", "Aurora inhibitors and microtubule destabilizers.", "strict_primary", "A", "SRC_NCI_DOCETAXEL"),
  family_row("MICROTUBULE_DESTABILIZER", "Microtubule destabilizers", "direct_moa", "MITOTIC_SPINDLE_PERTURBATION", "Compounds that disrupt microtubule polymerization or stability.", "Direct microtubule destabilization.", "Aurora inhibitors and microtubule stabilizers.", "strict_primary", "A", "SRC_PUBMED_NOCODAZOLE"),
  family_row("ESTROGENIC_ER_LIGAND", "Estrogenic ER agonist ligands", "direct_moa", "ESTROGEN_RECEPTOR_RESPONSE", "Estrogenic ligands that elicit an estrogen-receptor agonist transcriptional response.", "Documented estrogenic ER ligand activity.", "ER antagonists, degraders, and unrelated steroid hormones.", "strict_primary", "B", "SRC_PUBCHEM_EQUILIN"),
  family_row("HMGCR_INHIBITOR_STATIN", "HMGCR inhibitors (statins)", "direct_moa", "CHOLESTEROL_BIOSYNTHESIS", "Statins that directly inhibit HMG-CoA reductase.", "Recognized statin-class HMGCR inhibitor.", "Non-statin cholesterol-pathway modulators.", "strict_primary", "A", "SRC_FDA_STATIN"),
  family_row("P53_CHECKPOINT_RESPONSE", "p53 checkpoint response", "pathway_context", "", "Transcriptional context associated with p53 checkpoint activation.", "Perturbagens expected to engage the p53 checkpoint axis.", "None; context membership is not a direct-MoA claim.", "exploratory_context", "B", "SRC_CURATION_POLICY"),
  family_row("ERBB_RTK_AXIS", "ERBB/RTK signaling context", "pathway_context", "", "Transcriptional context of ERBB and related RTK signaling perturbation.", "Direct ERBB inhibitors and explicitly retained ERBB-containing multikinase context drugs.", "None; context membership is not a strict same-mechanism claim.", "exploratory_context", "B", "SRC_CURATION_POLICY"),
  family_row("RAS_RAF_MEK_ERK_AXIS", "RAS-RAF-MEK-ERK signaling context", "pathway_context", "", "Transcriptional context of MAPK-axis perturbation.", "ERBB, RAF, MEK, and explicitly retained multikinase context compounds.", "None; context membership is not a strict same-mechanism claim.", "exploratory_context", "B", "SRC_CURATION_POLICY"),
  family_row("PI3K_AKT_MTOR_AXIS", "PI3K-AKT-mTOR signaling context", "pathway_context", "", "Transcriptional context of PI3K-AKT-mTOR pathway perturbation.", "PI3K, mTOR, dual PI3K-mTOR, and rapalog compounds.", "None; context membership is not a strict same-mechanism claim.", "exploratory_context", "B", "SRC_CURATION_POLICY"),
  family_row("HSP90_CHAPERONE_PROTEOSTASIS", "HSP90 chaperone/proteostasis context", "pathway_context", "", "Transcriptional context linked to HSP90 chaperone inhibition.", "HSP90 inhibitors.", "None; context membership is not a strict same-mechanism claim.", "exploratory_context", "B", "SRC_CURATION_POLICY"),
  family_row("UBIQUITIN_PROTEASOME_AXIS", "Ubiquitin-proteasome context", "pathway_context", "", "Transcriptional context linked to proteasome and protease-tool perturbation.", "Constitutive proteasome, immunoproteasome, and explicitly flagged protease tools.", "None; context membership is not a strict same-mechanism claim.", "exploratory_context", "B", "SRC_CURATION_POLICY"),
  family_row("MITOTIC_SPINDLE_PERTURBATION", "Mitotic-spindle perturbation context", "pathway_context", "", "Transcriptional context shared by Aurora and microtubule perturbation.", "Aurora inhibitors and microtubule stabilizing/destabilizing agents.", "None; context membership is not a strict same-mechanism claim.", "exploratory_context", "B", "SRC_CURATION_POLICY"),
  family_row("ESTROGEN_RECEPTOR_RESPONSE", "Estrogen-receptor response context", "pathway_context", "", "Transcriptional context associated with estrogen-receptor agonism.", "Estrogenic ER ligands.", "None; context membership is not a strict same-mechanism claim.", "exploratory_context", "B", "SRC_CURATION_POLICY"),
  family_row("CHOLESTEROL_BIOSYNTHESIS", "Cholesterol-biosynthesis context", "pathway_context", "", "Transcriptional context associated with HMGCR inhibition and sterol homeostasis.", "Statins/HMGCR inhibitors.", "None; context membership is not a strict same-mechanism claim.", "exploratory_context", "B", "SRC_CURATION_POLICY")
))
rownames(family_registry) <- NULL

groups <- list(
  MDM2_ANTAGONIST = c("AMG-232", "CGM-097", "Nutlin-3", "RG7112", "RG-7388", "SAR405838"),
  ERBB_TKI_CORE = c("Afatinib", "Canertinib", "Dacomitinib", "Erlotinib", "Gefitinib", "Lapatinib", "Neratinib", "PD-158780"),
  MEK1_2_INHIBITOR = c("Cobimetinib", "MEK-162", "PD-198306", "Refametinib", "Ro-4987655", "Selumetinib", "TAK-733"),
  SELECTIVE_BRAF_INHIBITOR = c("PLX-4720", "Dabrafenib", "Vemurafenib", "SB-590885"),
  RAPALOG_FKBP12_MTORC1 = c("Sirolimus", "Everolimus", "Temsirolimus"),
  MTOR_KINASE_INHIBITOR = c("AZD-8055", "OSI-027"),
  DUAL_PI3K_MTOR_INHIBITOR = c("NVP-BEZ235", "PI-103", "GSK-1059615"),
  CLASS_I_PI3K_INHIBITOR = c("ZSTK-474", "Taselisib"),
  HSP90_INHIBITOR = c("17-AAG", "Tanespimycin", "Geldanamycin", "Alvespimycin", "Ganetespib", "NVP-AUY922", "BIIB-021", "PU-H71", "Radicicol"),
  CONSTITUTIVE_PROTEASOME_INHIBITOR = c("Bortezomib", "Carfilzomib", "Delanzomib", "Oprozomib", "Ixazomib", "MLN-9708"),
  IMMUNOPROTEASOME_PSMB8_INHIBITOR = "ONX-0914",
  AURORA_KINASE_INHIBITOR = c("Alisertib", "MLN-8054", "AMG-900", "AT-9283", "Danusertib", "Tozasertib"),
  MICROTUBULE_STABILIZER = c("Docetaxel", "Cabazitaxel"),
  MICROTUBULE_DESTABILIZER = c("Nocodazole", "Parbendazole"),
  ESTROGENIC_ER_LIGAND = c("Estradiol", "Estrone", "Equilin", "Alpha-Estradiol"),
  HMGCR_INHIBITOR_STATIN = c("Atorvastatin", "Fluvastatin", "Lovastatin", "Pitavastatin", "Pravastatin", "Rosuvastatin", "Simvastatin")
)

context_only <- c(
  "Vandetanib", "GW-5074", "RAF-265", "Sorafenib", "Regorafenib",
  "Ponatinib", "MG-132", "z-leu3-VS"
)
curated_terms <- sort(unique(c(unlist(groups, use.names = FALSE), context_only)), method = "radix")

term_registry <- read_chem_pairs(chem_library)
if (!all(curated_terms %in% term_registry$library_term)) {
  stop("Curated terms missing from the Chem library: ", paste(setdiff(curated_terms, term_registry$library_term), collapse = ", "))
}
if (length(curated_terms) != 80L) stop("Expected exactly 80 curated library terms.")

curated_idx <- match(curated_terms, term_registry$library_term)
term_registry$canonical_compound_id[curated_idx] <- paste0("CMP_", slug(curated_terms))
term_registry$canonical_compound_name[curated_idx] <- curated_terms
term_registry$analysis_unit_id[curated_idx] <- paste0("UNIT_", slug(curated_terms))
term_registry$analysis_unit_name[curated_idx] <- curated_terms
term_registry$term_relation[curated_idx] <- "exact_library_name"
term_registry$identity_evidence_tier[curated_idx] <- "C"
term_registry$curation_status[curated_idx] <- "mechanism_curated_v2"
term_registry$multi_target_flag[curated_idx] <- FALSE
term_registry$tool_compound_flag[curated_idx] <- FALSE

set_override <- function(term, canonical_id, canonical_name, unit_id, unit_name,
                         relation, preferred, tier, source_id) {
  i <- match(term, term_registry$library_term)
  if (is.na(i)) stop("Unknown override term: ", term)
  term_registry$canonical_compound_id[i] <<- canonical_id
  term_registry$canonical_compound_name[i] <<- canonical_name
  term_registry$analysis_unit_id[i] <<- unit_id
  term_registry$analysis_unit_name[i] <<- unit_name
  term_registry$term_relation[i] <<- relation
  term_registry$preferred_for_analysis[i] <<- preferred
  term_registry$identity_evidence_tier[i] <<- tier
  term_registry$identity_source_id[i] <<- source_id
}

set_override("RG-7388", "CMP_IDASANUTLIN", "Idasanutlin", "UNIT_IDASANUTLIN", "Idasanutlin", "code_name", TRUE, "A", "SRC_NCI_IDASANUTLIN")
set_override("AMG-232", "CMP_NAVTEMADLIN", "Navtemadlin", "UNIT_NAVTEMADLIN", "Navtemadlin", "code_name", TRUE, "B", "SRC_CURATION_POLICY")
set_override("MEK-162", "CMP_BINIMETINIB", "Binimetinib", "UNIT_BINIMETINIB", "Binimetinib", "code_name", TRUE, "B", "SRC_CURATION_POLICY")
set_override("NVP-AUY922", "CMP_LUMINESPIB", "Luminespib", "UNIT_LUMINESPIB", "Luminespib", "code_name", TRUE, "B", "SRC_CURATION_POLICY")
set_override("NVP-BEZ235", "CMP_DACTOLISIB", "Dactolisib", "UNIT_DACTOLISIB", "Dactolisib", "code_name", TRUE, "B", "SRC_CURATION_POLICY")
set_override("Alpha-Estradiol", "CMP_17ALPHA_ESTRADIOL", "17alpha-Estradiol", "UNIT_17ALPHA_ESTRADIOL", "17alpha-Estradiol", "name_variant", TRUE, "B", "SRC_CURATION_POLICY")
set_override("Sirolimus", "CMP_SIROLIMUS_RAPAMYCIN", "Sirolimus (rapamycin)", "UNIT_SIROLIMUS_RAPAMYCIN", "Sirolimus (rapamycin)", "exact_synonym", TRUE, "A", "SRC_NCI_SIROLIMUS")
set_override("17-AAG", "CMP_TANESPIMYCIN", "Tanespimycin", "UNIT_TANESPIMYCIN", "Tanespimycin/17-AAG", "code_name", TRUE, "A", "SRC_NCI_TANESPIMYCIN")
set_override("Tanespimycin", "CMP_TANESPIMYCIN", "Tanespimycin", "UNIT_TANESPIMYCIN", "Tanespimycin/17-AAG", "exact_synonym", FALSE, "A", "SRC_NCI_TANESPIMYCIN")
set_override("Ixazomib", "CMP_IXAZOMIB", "Ixazomib", "UNIT_IXAZOMIB_PROGRAM", "Ixazomib drug program", "active_drug_name", TRUE, "A", "SRC_NCI_IXAZOMIB")
set_override("MLN-9708", "CMP_IXAZOMIB_CITRATE", "Ixazomib citrate (MLN-9708)", "UNIT_IXAZOMIB_PROGRAM", "Ixazomib drug program", "prodrug_program", FALSE, "A", "SRC_NCI_IXAZOMIB")

multitarget_terms <- c("Vandetanib", "RAF-265", "Sorafenib", "Regorafenib", "Ponatinib", "AT-9283", "Danusertib", "Tozasertib")
tool_terms <- c("PD-158780", "PLX-4720", "SB-590885", "GW-5074", "PD-198306", "PI-103", "MG-132", "z-leu3-VS", "MLN-8054", "Radicicol")
term_registry$multi_target_flag[match(multitarget_terms, term_registry$library_term)] <- TRUE
term_registry$tool_compound_flag[match(tool_terms, term_registry$library_term)] <- TRUE

collision_norms <- names(which(table(term_registry$library_term_norm) > 1L))
term_registry$normalization_collision_resolution[
  term_registry$library_term_norm %in% collision_norms
] <- "distinct_analysis_units_pending_identity_review"

source_for_family <- setNames(family_registry$source_id, family_registry$family_id)

membership_rows <- list()
add_membership <- function(terms, family_id, role, action, selectivity,
                           evidence_tier = "B", source_id = source_for_family[[family_id]],
                           curation_status = "curated_v2") {
  units <- unique(term_registry$analysis_unit_id[match(terms, term_registry$library_term)])
  if (anyNA(units)) stop("Membership includes an unknown library term.")
  membership_rows[[length(membership_rows) + 1L]] <<- data.frame(
    analysis_unit_id = units,
    family_id = family_id,
    membership_role = role,
    pharmacologic_action = action,
    selectivity_class = selectivity,
    evidence_tier = evidence_tier,
    source_id = source_id,
    curation_status = curation_status,
    taxonomy_version = taxonomy_version,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

add_membership(groups$MDM2_ANTAGONIST, "MDM2_ANTAGONIST", "primary", "MDM2-p53 interaction antagonism", "direct")
add_membership(groups$ERBB_TKI_CORE, "ERBB_TKI_CORE", "primary", "ERBB-family kinase inhibition", "selective_to_pan_ERBB")
add_membership(groups$MEK1_2_INHIBITOR, "MEK1_2_INHIBITOR", "primary", "MEK1/2 inhibition", "direct")
add_membership(groups$SELECTIVE_BRAF_INHIBITOR, "SELECTIVE_BRAF_INHIBITOR", "primary", "RAF kinase inhibition", "selective_or_primary_RAF")
add_membership(groups$RAPALOG_FKBP12_MTORC1, "RAPALOG_FKBP12_MTORC1", "primary", "FKBP12-dependent allosteric mTORC1 inhibition", "rapalog")
add_membership(groups$MTOR_KINASE_INHIBITOR, "MTOR_KINASE_INHIBITOR", "primary", "ATP-competitive mTOR kinase inhibition", "mTORC1_2")
add_membership(groups$DUAL_PI3K_MTOR_INHIBITOR, "DUAL_PI3K_MTOR_INHIBITOR", "primary", "dual PI3K and mTOR inhibition", "dual_target")
add_membership(groups$CLASS_I_PI3K_INHIBITOR, "CLASS_I_PI3K_INHIBITOR", "primary", "class I PI3K inhibition", "PI3K")
add_membership(groups$HSP90_INHIBITOR, "HSP90_INHIBITOR", "primary", "HSP90 chaperone inhibition", "direct")
add_membership(groups$CONSTITUTIVE_PROTEASOME_INHIBITOR, "CONSTITUTIVE_PROTEASOME_INHIBITOR", "primary", "constitutive proteasome inhibition", "constitutive_proteasome")
add_membership(groups$IMMUNOPROTEASOME_PSMB8_INHIBITOR, "IMMUNOPROTEASOME_PSMB8_INHIBITOR", "primary", "PSMB8/LMP7 inhibition", "immunoproteasome", "B")
add_membership(groups$AURORA_KINASE_INHIBITOR, "AURORA_KINASE_INHIBITOR", "primary", "Aurora kinase inhibition", "AURKA_to_pan_Aurora")
add_membership(groups$MICROTUBULE_STABILIZER, "MICROTUBULE_STABILIZER", "primary", "microtubule stabilization", "taxane")
add_membership(groups$MICROTUBULE_DESTABILIZER, "MICROTUBULE_DESTABILIZER", "primary", "microtubule destabilization", "direct")
add_membership(groups$ESTROGENIC_ER_LIGAND, "ESTROGENIC_ER_LIGAND", "primary", "estrogen-receptor agonist signaling", "estrogenic_ligand", "B")
add_membership(groups$HMGCR_INHIBITOR_STATIN, "HMGCR_INHIBITOR_STATIN", "primary", "HMG-CoA reductase inhibition", "statin")

add_membership(groups$MDM2_ANTAGONIST, "P53_CHECKPOINT_RESPONSE", "context", "p53 checkpoint reactivation", "pathway_context", "B", "SRC_CURATION_POLICY")
add_membership(c(groups$ERBB_TKI_CORE, "Vandetanib"), "ERBB_RTK_AXIS", "context", "ERBB/RTK signaling perturbation", "pathway_context", "B", "SRC_CURATION_POLICY")
add_membership(c(groups$ERBB_TKI_CORE, groups$MEK1_2_INHIBITOR, groups$SELECTIVE_BRAF_INHIBITOR, "GW-5074", "RAF-265", "Sorafenib", "Regorafenib", "Ponatinib", "Vandetanib"), "RAS_RAF_MEK_ERK_AXIS", "context", "MAPK-axis perturbation", "pathway_context", "B", "SRC_CURATION_POLICY")
add_membership(c(groups$RAPALOG_FKBP12_MTORC1, groups$MTOR_KINASE_INHIBITOR, groups$DUAL_PI3K_MTOR_INHIBITOR, groups$CLASS_I_PI3K_INHIBITOR), "PI3K_AKT_MTOR_AXIS", "context", "PI3K-AKT-mTOR-axis perturbation", "pathway_context", "B", "SRC_CURATION_POLICY")
add_membership(groups$HSP90_INHIBITOR, "HSP90_CHAPERONE_PROTEOSTASIS", "context", "HSP90/proteostasis perturbation", "pathway_context", "B", "SRC_CURATION_POLICY")
add_membership(c(groups$CONSTITUTIVE_PROTEASOME_INHIBITOR, groups$IMMUNOPROTEASOME_PSMB8_INHIBITOR, "MG-132", "z-leu3-VS"), "UBIQUITIN_PROTEASOME_AXIS", "context", "ubiquitin-proteasome perturbation", "pathway_context", "B", "SRC_CURATION_POLICY")
add_membership(c(groups$AURORA_KINASE_INHIBITOR, groups$MICROTUBULE_STABILIZER, groups$MICROTUBULE_DESTABILIZER), "MITOTIC_SPINDLE_PERTURBATION", "context", "mitotic-spindle perturbation", "phenotypic_context", "B", "SRC_CURATION_POLICY")
add_membership(groups$ESTROGENIC_ER_LIGAND, "ESTROGEN_RECEPTOR_RESPONSE", "context", "estrogen-receptor response", "pathway_context", "B", "SRC_CURATION_POLICY")
add_membership(groups$HMGCR_INHIBITOR_STATIN, "CHOLESTEROL_BIOSYNTHESIS", "context", "cholesterol-biosynthesis perturbation", "pathway_context", "B", "SRC_CURATION_POLICY")

membership <- do.call(rbind, membership_rows)
membership <- membership[order(membership$analysis_unit_id, membership$family_id, method = "radix"), , drop = FALSE]
rownames(membership) <- NULL

forbidden_pattern <- paste0(
  "(^|_)(activity_score|activity_gls_z|activity_class|nes(_up|_down)?|",
  "pvalue|p_value|p_combined|padj.*|fdr|qvalue|neg_log10.*|evaluated|",
  "recovered|recovery.*|is_significant|n_recovered.*|best_recovered.*)(_|$)"
)

validate_tables <- function() {
  if (nrow(term_registry) != 5425L) stop("Expected 5,425 Chem library terms.")
  if (anyDuplicated(term_registry$library_term_id)) stop("Duplicate library_term_id.")
  if (anyDuplicated(term_registry$library_term)) stop("Duplicate exact library_term.")
  if (anyDuplicated(family_registry$family_id)) stop("Duplicate family_id.")
  key <- paste(membership$analysis_unit_id, membership$family_id, sep = "\r")
  if (anyDuplicated(key)) stop("Duplicate analysis_unit_id/family_id membership.")
  if (!all(term_registry$pair_complete == (term_registry$up_present & term_registry$dn_present))) stop("Pair completeness mismatch.")
  if (!all(term_registry$pair_complete)) stop("The Chem v2 source universe requires complete Up/Down pairs.")
  if (!all(membership$analysis_unit_id %in% term_registry$analysis_unit_id)) stop("Dangling membership analysis_unit_id.")
  if (!all(membership$family_id %in% family_registry$family_id)) stop("Dangling membership family_id.")
  if (!all(membership$source_id %in% source_registry$source_id)) stop("Dangling membership source_id.")
  if (!all(family_registry$source_id %in% source_registry$source_id)) stop("Dangling family source_id.")
  preferred_n <- tapply(as.integer(term_registry$preferred_for_analysis), term_registry$analysis_unit_id, sum)
  if (!all(preferred_n == 1L)) stop("Every analysis unit must have exactly one preferred term.")
  canonical_name_n <- tapply(term_registry$canonical_compound_name, term_registry$canonical_compound_id, function(x) length(unique(x)))
  if (!all(canonical_name_n == 1L)) stop("A canonical compound ID maps to multiple names.")
  if (length(collision_norms) != 3L) stop("Expected three normalized-name collision groups in the source library.")
  if (any(term_registry$normalization_collision_resolution[term_registry$library_term_norm %in% collision_norms] == "none")) stop("Normalized-name collision is unresolved.")
  for (tbl in list(term_registry, family_registry, membership, source_registry)) {
    bad <- grep(forbidden_pattern, names(tbl), ignore.case = TRUE, perl = TRUE, value = TRUE)
    if (length(bad)) stop("Outcome-derived columns are forbidden: ", paste(bad, collapse = ", "))
  }
  anchors <- list(
    "Erlotinib" = "ERBB_TKI_CORE", "Cobimetinib" = "MEK1_2_INHIBITOR",
    "PLX-4720" = "SELECTIVE_BRAF_INHIBITOR", "Alisertib" = "AURORA_KINASE_INHIBITOR",
    "Docetaxel" = "MICROTUBULE_STABILIZER", "Nocodazole" = "MICROTUBULE_DESTABILIZER",
    "Bortezomib" = "CONSTITUTIVE_PROTEASOME_INHIBITOR", "RG-7388" = "MDM2_ANTAGONIST"
  )
  for (term in names(anchors)) {
    unit <- term_registry$analysis_unit_id[match(term, term_registry$library_term)]
    if (!paste(unit, anchors[[term]], sep = "\r") %in% key) stop("Scientific anchor failed: ", term)
  }
  invisible(TRUE)
}

validate_tables()

data_raw_dir <- file.path(repo_dir, "data-raw", "chem_family_v2")
package_dir <- file.path(repo_dir, "inst", "extdata", "annotation_catalogs", "v2")
stage <- tempfile("chem-family-v2-stage-")
dir.create(stage, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)

tables <- list(
  chem_term_registry_v2 = term_registry,
  chem_family_registry_v2 = family_registry,
  chem_family_membership_v2 = membership,
  chem_source_registry_v2 = source_registry
)

for (nm in names(tables)) {
  utils::write.csv(
    tables[[nm]], file.path(stage, paste0(nm, ".csv")),
    row.names = FALSE, na = "", fileEncoding = "UTF-8"
  )
}

manifest_rows <- lapply(names(tables), function(nm) {
  path <- file.path(stage, paste0(nm, ".csv"))
  data.frame(
    artifact_id = nm,
    artifact_role = "frozen_taxonomy_table",
    relative_path_or_provider_key = file.path("inst", "extdata", "annotation_catalogs", "v2", paste0(nm, ".csv")),
    schema_version = schema_version,
    taxonomy_version = taxonomy_version,
    sha256 = sha256_file(path),
    n_rows = nrow(tables[[nm]]),
    primary_key = switch(nm,
      chem_term_registry_v2 = "library_term_id; unique(library_term)",
      chem_family_registry_v2 = "family_id",
      chem_family_membership_v2 = "analysis_unit_id + family_id",
      chem_source_registry_v2 = "source_id"
    ),
    created_utc = lock_utc,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})

manifest <- do.call(rbind, c(manifest_rows, list(
  data.frame(
    artifact_id = c("chem_signature_provider", "chem_family_v2_builder", "legacy_roster_provenance"),
    artifact_role = c("provider_source", "builder", "legacy_outcome_seen_source"),
    relative_path_or_provider_key = c(
      "L1000PAIdata::LINCS_L1000_Chem_Pert_Consensus_Sigs.txt",
      "tools/chem_family_v2/build_chem_family_catalog_v2.R",
      "legacy:Targeted_drug_10GEO_figure/targeted_drug_10geo_same_mechanism_drug_roster.csv"
    ),
    schema_version = schema_version,
    taxonomy_version = taxonomy_version,
    sha256 = c(
      sha256_file(chem_library),
      sha256_file(script_path()),
      "8d823fca5def3cc48eedb7a94e53c48d3bdbc3c1725432a899810c12db07ca19"
    ),
    n_rows = c(10850L, NA_integer_, 86L),
    primary_key = c("directional signature term", "file", "legacy case-term row"),
    created_utc = lock_utc,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
)))
rownames(manifest) <- NULL
utils::write.csv(manifest, file.path(stage, "chem_taxonomy_v2_manifest.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")

for (out_dir in c(data_raw_dir, package_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(stage, full.names = TRUE)
  ok <- file.copy(files, file.path(out_dir, basename(files)), overwrite = TRUE)
  if (!all(ok)) stop("Failed to publish Chem family v2 artifacts to ", out_dir)
}

message("Chem family v2 build complete")
message("  terms: ", nrow(term_registry), " (curated: ", sum(term_registry$curation_status == "mechanism_curated_v2"), ")")
message("  families: ", nrow(family_registry))
message("  memberships: ", nrow(membership))
message("  library SHA-256: ", sha256_file(chem_library))
