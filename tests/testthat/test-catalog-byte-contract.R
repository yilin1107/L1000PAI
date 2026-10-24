test_that("frozen catalog hashes cover canonical shipped bytes on every platform", {
  for (version in available_chem_family_versions()) {
    manifest <- get_chem_family_catalog("manifest", version = version)
    frozen <- manifest[manifest$artifact_role %in%
      c("frozen_taxonomy_table", "shared_identity_table"), , drop = FALSE]
    expect_gt(nrow(frozen), 0L)
    for (i in seq_len(nrow(frozen))) {
      relative <- sub("^inst/", "", frozen$relative_path_or_provider_key[i])
      path <- system.file(relative, package = "L1000PAI")
      expect_true(nzchar(path), info = relative)
      bytes <- readBin(path, what = "raw", n = file.info(path)$size)
      if (grepl("\\.csv$", path)) {
        expect_false(any(bytes == as.raw(13)), info = paste(relative, "must use LF"))
      }
      expect_identical(
        digest::digest(bytes, algo = "sha256", serialize = FALSE),
        frozen$sha256[i],
        info = paste(version, relative)
      )
    }
    expect_true(all(verify_chem_family_catalogs(
      version = version, verify_hash = TRUE, stop_on_error = FALSE
    )$pass))
  }
})
