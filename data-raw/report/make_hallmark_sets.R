# Hallmark E2F targets and G2M checkpoint (MSigDB, CC BY 4.0) for the proliferation check of
# pai_deg_report(). Run from the package root; written with msigdbr 26.1.1 (MSigDB 2026.1.Hs).
h <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
h <- h[h$gs_name %in% c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT"), ]
d <- unique(data.frame(gene_set = h$gs_name, gene_symbol = h$gene_symbol, msigdb_version = h$db_version, stringsAsFactors = FALSE))
d <- d[order(d$gene_set, d$gene_symbol), ]
write.csv(d, "inst/extdata/report/hallmark_proliferation_sets.csv", row.names = FALSE)
