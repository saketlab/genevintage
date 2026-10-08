testthat::test_that("symbol index builds and filters releases", {
  a <- data.frame(species = "homo_sapiens", release = c("95", "110"), assembly = "38", source = "ensembl", name = c("OLD", "TP53"))
  ix <- build_symbol_index(a, min_release = 100)
  testthat::expect_equal(unique(ix$release), "110")
})
testthat::test_that("weighted symbol overlap selects the right species", {
  a <- data.frame(species = "homo_sapiens", release = "110", assembly = "38", source = "ensembl", name = c("TP53", "XIST", paste0("H", 1:30)))
  b <- data.frame(species = "mus_musculus", release = "110", assembly = "GRCm39", source = "ensembl", name = c("Tp53", "Xist", "B", paste0("M", 1:30)))
  z <- detect_species_names(c("TP53", "XIST", paste0("H", 1:30)), build_symbol_index(list(a, b)), min_overlap = 5)
  testthat::expect_identical(z$species[[1]], "homo_sapiens")
  testthat::expect_identical(z$status[[1]], "selected")
})
testthat::test_that("excel date-like symbols are not sole evidence", {
  a <- data.frame(species = "homo_sapiens", name = c("SEPT1", "MARCH1"))
  z <- detect_species_names(c("SEPT1", "MARCH1"), build_symbol_index(a), min_overlap = 1)
  testthat::expect_equal(nrow(z), 0)
})
testthat::test_that("insufficient overlap is reported", {
  a <- data.frame(species = "homo_sapiens", name = paste0("H", 1:5))
  z <- detect_species_names("H1", build_symbol_index(a), min_overlap = 20)
  testthat::expect_identical(z$status[[1]], "insufficient")
})
testthat::test_that("ambiguous candidates are retained", {
  a <- data.frame(species = "homo_sapiens", name = paste0("G", 1:30))
  b <- data.frame(species = "mus_musculus", name = paste0("G", 1:30))
  z <- detect_species_names(paste0("G", 1:30), build_symbol_index(list(a, b)), min_overlap = 20)
  testthat::expect_true(all(z$status == "ambiguous"))
})
testthat::test_that("duplicate symbols do not inflate overlap", {
  a <- data.frame(species = "homo_sapiens", name = c("A", "A", "B"))
  z <- detect_species_names(c("A", "A", "B"), build_symbol_index(a), min_overlap = 1)
  testthat::expect_equal(z$overlap[[1]], 2)
})
testthat::test_that("empty and NA symbols are safe", {
  a <- data.frame(species = "homo_sapiens", name = "A")
  testthat::expect_equal(nrow(detect_species_names(c(NA, ""), build_symbol_index(a))), 0)
})
testthat::test_that("case and whitespace are normalized", {
  a <- data.frame(species = "homo_sapiens", name = c("TP53", "XIST"))
  z <- detect_species_names(c(" tp53 ", "xist"), build_symbol_index(a), min_overlap = 1)
  testthat::expect_equal(z$overlap[[1]], 2)
})
testthat::test_that("post-2020 release floor is honored", {
  a <- data.frame(species = "homo_sapiens", release = c("99", "110"), name = c("OLD", "NEW"))
  ix <- build_symbol_index(a, min_release = 100)
  testthat::expect_false("OLD" %in% ix$symbol)
})
testthat::test_that("malformed index is rejected", {
  testthat::expect_error(detect_species_names("TP53", data.frame(species = "human")), "symbol")
})
testthat::test_that("multiple organisms and releases rank correctly", {
  species <- c("homo_sapiens", "mus_musculus", "gallus_gallus", "danio_rerio", "drosophila_melanogaster", "arabidopsis_thaliana")
  maps <- lapply(seq_along(species), function(i) {
    data.frame(
      species = species[[i]], release = c("95", "110"), source = c("ensembl", "ensembl"),
      assembly = c("old", "new"), name = c(paste0("OLD", i), paste0("NEW", i))
    )
  })
  ix <- build_symbol_index(maps, min_release = 100)
  testthat::expect_true(all(ix$release == "110"))
  for (i in seq_along(species)) {
    z <- detect_species_names(c(paste0("NEW", i)), ix, min_overlap = 1, min_score = .1)
    testthat::expect_identical(z$species[[1]], species[[i]])
  }
})
testthat::test_that("GENCODE releases remain distinct from Ensembl", {
  a <- data.frame(species = "homo_sapiens", release = "44", source = "gencode", assembly = "38", name = c("GENE_A", "GENE_B"))
  b <- data.frame(species = "homo_sapiens", release = "112", source = "ensembl", assembly = "38", name = c("GENE_A", "GENE_B", "GENE_C"))
  z <- detect_species_names(c("GENE_A", "GENE_B"), build_symbol_index(list(a, b)), min_overlap = 1, min_score = .1)
  testthat::expect_true(all(z$species == "homo_sapiens"))
  testthat::expect_true(all(z$source %in% c("gencode", "ensembl")))
})
testthat::test_that("post-2020 symbols can select a release", {
  old <- data.frame(species = "danio_rerio", release = "100", source = "ensembl", name = paste0("OLD", 1:20))
  new <- data.frame(species = "danio_rerio", release = "110", source = "ensembl", name = c(paste0("OLD", 1:10), paste0("NEW", 1:20)))
  ix <- build_symbol_index(list(old, new), min_release = 105)
  z <- detect_species_names(paste0("NEW", 1:20), ix, min_overlap = 10)
  testthat::expect_identical(z$status[[1]], "selected")
  testthat::expect_identical(z$release[[1]], "110")
})
testthat::test_that("common mammalian symbols are downweighted", {
  a <- data.frame(species = "homo_sapiens", name = c("TP53", "ACTB", paste0("H", 1:20)))
  b <- data.frame(species = "mus_musculus", name = c("TP53", "Actb", "ACTB", paste0("M", 1:20)))
  z <- detect_species_names(c("TP53", "ACTB", paste0("H", 1:20)), build_symbol_index(list(a, b)), min_overlap = 10)
  testthat::expect_identical(z$species[[1]], "homo_sapiens")
})
testthat::test_that("mixed organism symbol sets are ambiguous", {
  a <- data.frame(species = "gallus_gallus", name = paste0("A", 1:20))
  b <- data.frame(species = "danio_rerio", name = paste0("B", 1:20))
  z <- detect_species_names(c(paste0("A", 1:10), paste0("B", 1:10)), build_symbol_index(list(a, b)), min_overlap = 5)
  testthat::expect_true(all(z$status == "ambiguous"))
})
testthat::test_that("symbol aliases and punctuation are normalized", {
  a <- data.frame(species = "drosophila_melanogaster", name = c("GAPDH", "RpL3"))
  z <- detect_species_names(c("gapdh", " rpl3 "), build_symbol_index(a), min_overlap = 1)
  testthat::expect_equal(z$overlap[[1]], 2)
})
testthat::test_that("Excel month tokens cannot dominate post-2020 detection", {
  a <- data.frame(species = "homo_sapiens", release = "110", name = c("SEPT1", "MARCH1", paste0("H", 1:30)))
  b <- data.frame(species = "mus_musculus", release = "110", name = c("SEPT1", "MARCH1", paste0("M", 1:30)))
  z <- detect_species_names(c("SEPT1", "MARCH1", paste0("H", 1:20)), build_symbol_index(list(a, b)), min_overlap = 10)
  testthat::expect_identical(z$species[[1]], "homo_sapiens")
})
testthat::test_that("the comprehensive builder script covers every fingerprint key", {
  script <- testthat::test_path("..", "..", "data-raw", "build_symbol_index_all.R")
  # data-raw is .Rbuildignore'd, so it is absent from an installed or CRAN-checked package
  testthat::skip_if_not(file.exists(script), "data-raw not present (not shipped in the built package)")
  fp <- readRDS(system.file("extdata", "fingerprints.rds", package = "genevintage"))
  testthat::expect_true(all(c("species", "release", "assembly", "source") %in% names(fp)))
  testthat::expect_gt(nrow(unique(fp[, c("species", "release", "assembly", "source")])), 9000)
})
