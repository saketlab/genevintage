# Orthologs. The network paths are skipped off-CRAN; the rest is pure logic.

test_that("a species cannot be its own ortholog target", {
  expect_error(fetch_orthologs("homo_sapiens", "homo_sapiens", 116), "must differ")
})

test_that("a release Compara predates fails with a usable message", {
  skip_net()
  # r85 has a tsv/ tree but no ensembl-compara/ in it. The floor differs by
  # division -- Ensembl Genomes starts far earlier -- so it is read off the
  # server rather than hardcoded.
  expect_error(
    suppressWarnings(fetch_orthologs("homo_sapiens", "mus_musculus", 85)),
    "no homologies"
  )
})

test_that("species in different divisions are refused before any download", {
  skip_no_index()
  expect_error(
    fetch_orthologs("arabidopsis_thaliana", "homo_sapiens", 63),
    "different Ensembl divisions"
  )
})

test_that("an unknown target species is refused", {
  skip_no_index()
  expect_error(
    fetch_orthologs("homo_sapiens", "not_a_species", 116),
    "not a species in the index"
  )
})

test_that("stream_homologies keeps only cross-species orthologs", {
  # paralogs share the file with orthologs, and only the requested partner
  # species is wanted -- neither should survive
  wide <- c(
    "gene_stable_id", "protein_stable_id", "species", "identity",
    "homology_type", "homology_gene_stable_id", "homology_protein_stable_id",
    "homology_species", "homology_identity", "is_high_confidence"
  )
  row <- function(id, type, oid, sp, conf = "1") {
    c(id, "P1", "homo_sapiens", "90.0", type, oid, "P2", sp, "88.0", conf)
  }
  u <- hom_file(row("ENSG1", "ortholog_one2one", "ENSMUSG1", "mus_musculus"),
    row("ENSG2", "ortholog_one2many", "ENSMUSG2", "mus_musculus", "0"),
    row("ENSG3", "other_paralog", "ENSG9", "homo_sapiens"),
    row("ENSG4", "ortholog_one2one", "ENSRNOG1", "rattus_norvegicus"),
    header = wide
  )

  out <- stream_homologies(u, "mus_musculus")
  expect_equal(nrow(out), 2L)
  expect_setequal(out$id, c("ENSG1", "ENSG2"))
  expect_true(all(startsWith(out$homology_type, "ortholog")))
  expect_equal(out$high_confidence, c(TRUE, FALSE))
  expect_equal(out$identity, c(90, 90))
})

test_that("a homology file missing a needed column is rejected", {
  u <- hom_file(c("ENSG1", "mus_musculus"),
    header = c("gene_stable_id", "homology_species")
  )
  expect_error(stream_homologies(u, "mus_musculus"), "missing column")
})

test_that("an older homology file without is_high_confidence still reads", {
  old <- c(
    "gene_stable_id", "homology_type", "homology_gene_stable_id",
    "homology_species", "identity", "homology_identity"
  )
  u <- hom_file(c("ENSG1", "ortholog_one2one", "ENSMUSG1", "mus_musculus", "90", "88"),
    header = old
  )
  out <- stream_homologies(u, "mus_musculus")
  expect_equal(out$ortholog_id, "ENSMUSG1")
  expect_true(is.na(out$high_confidence))
})

test_that("no ortholog in the target species gives an empty frame, not an error", {
  u <- hom_file(c(
    "ENSG1", "ortholog_one2one", "ENSRNOG1", "rattus_norvegicus",
    "90", "88", "1"
  ))
  out <- stream_homologies(u, "mus_musculus")
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0L)
})

test_that("orthologs() finds mouse orthologs of human genes", {
  skip_net()
  o <- suppressMessages(orthologs(c("ENSG00000141510", "ENSG00000012048"),
    to = "mus_musculus"
  ))
  expect_true(all(c("id", "ortholog_id", "ortholog_name") %in% names(o)))
  expect_true("Trp53" %in% o$ortholog_name) # human TP53
  expect_true(all(startsWith(o$ortholog_id, "ENSMUSG")))
})

test_that("a versioned identifier is matched and reported back as given", {
  skip_net()
  o <- suppressMessages(orthologs("ENSG00000141510.16", to = "mus_musculus"))
  expect_equal(o$id, "ENSG00000141510.16")
})

test_that("any two indexed species can be paired, not just human", {
  skip_net()
  # mouse -> rat exercises a pair where neither side is human, and the mouse
  # file carries the pair directly (no swap)
  o <- suppressMessages(orthologs(c("ENSMUSG00000020886", "ENSMUSG00000021466"),
    to = "rattus_norvegicus", species = "mus_musculus"
  ))
  expect_gt(nrow(o), 0)
  expect_true(all(startsWith(o$ortholog_id, "ENSRNOG")))
  expect_true("Ptch1" %in% o$ortholog_name)
})

test_that("ortholog_species lists what Compara actually covers", {
  skip_net()
  sp <- ortholog_species(116)
  expect_gt(length(sp), 200)
  expect_true(all(c("homo_sapiens", "mus_musculus", "danio_rerio") %in% sp))
})

test_that("both sides of the pair carry an identifier and a name", {
  skip_net()
  o <- suppressMessages(orthologs(c("ENSG00000141510", "ENSG00000012048"),
    to = "mouse"
  ))
  expect_true(all(c("id", "name", "ortholog_id", "ortholog_name") %in% names(o)))
  expect_equal(o$name[o$id == "ENSG00000141510"], "TP53")
  expect_equal(o$ortholog_name[o$id == "ENSG00000141510"], "Trp53")
})

test_that("gene names work as input and give what the identifiers give", {
  skip_net()
  by_id <- suppressMessages(orthologs(c("ENSG00000141510", "ENSG00000012048"),
    to = "mouse", species = "human"
  ))
  by_name <- suppressMessages(orthologs(c("TP53", "BRCA1"),
    to = "mouse", species = "human"
  ))
  # names resolve to identifiers first, so the two agree row for row
  expect_equal(by_name[order(by_name$id), ], by_id[order(by_id$id), ],
    ignore_attr = "row.names"
  )
  # a name is not case-sensitive here: no two human symbols differ only by case
  lower <- suppressMessages(orthologs(c("tp53"), to = "mouse", species = "human"))
  expect_equal(lower$ortholog_name, "Trp53")
})

test_that("names need a species, since they do not carry one", {
  skip_no_index()
  expect_error(
    orthologs(c("TP53", "BRCA1"), to = "mus_musculus"),
    "do not say which species"
  )
})

test_that("an unmatched or ambiguous name is reported, not swallowed", {
  skip_net()
  expect_warning(
    suppressMessages(orthologs(c("TP53", "NOT_A_GENE"), to = "mouse", species = "human")),
    "matched no gene"
  )
  # Y_RNA names 758 human genes; all are kept rather than the first taken
  expect_warning(
    suppressMessages(orthologs(c("TP53", "Y_RNA"),
      to = "mouse", species = "human",
      names = FALSE, one2one = FALSE
    )),
    "match more than one gene"
  )
})

test_that("the swap path turns the whole row around, not just the ids", {
  # Compara records human->mouse only in the mouse file, so the pair is read
  # from there and reversed: ids, identities and direction all have to follow.
  u <- hom_file(c(
    "ENSMUSG1", "ortholog_one2many", "ENSG1", "homo_sapiens",
    "70", "80", "1"
  ))
  out <- stream_homologies(u, "homo_sapiens", swap = TRUE)
  expect_equal(out$id, "ENSG1")
  expect_equal(out$ortholog_id, "ENSMUSG1")
  expect_equal(out$homology_type, "ortholog_many2one") # one2many, other way up
  expect_equal(out$identity, 80)
  expect_equal(out$ortholog_identity, 70)
})

test_that("an empty result has the same columns as a populated one", {
  skip_no_index()
  # code reading $ortholog_name would otherwise break only when nothing matched
  e <- .empty_homology(TRUE)
  expect_true(all(c("id", "name", "ortholog_id", "ortholog_name") %in% names(e)))
  expect_equal(nrow(e), 0L)
  expect_false("name" %in% names(.empty_homology(FALSE)))
})

test_that("a completely empty homology file is reported as such", {
  f <- tempfile(fileext = ".tsv.gz")
  on.exit(unlink(f))
  writeLines(character(), gzfile(f)) # no header at all, so hom_file cannot make it
  expect_error(
    stream_homologies(paste0("file://", f), "mus_musculus"),
    "empty homology file"
  )
})

test_that("a source species with no homology file still finds the pair", {
  # Compara records a pair in one direction only. If the source species has no
  # homology directory of its own, the target's file -- read the other way
  # round -- is the whole answer, so the first missing directory must not be
  # fatal.
  hom <- data.frame(
    id = "ENSG1", ortholog_id = "ENSMUSG1",
    homology_type = "ortholog_one2one", identity = 90,
    ortholog_identity = 88, high_confidence = TRUE,
    stringsAsFactors = FALSE
  )
  testthat::local_mocked_bindings(
    .collections = function(root, release, kind, species) {
      if (species == "homo_sapiens") stop("Compara has no homologies for homo_sapiens")
      "file:///dev/null"
    },
    stream_homologies = function(url, to, swap = FALSE, chunk = 2e5) {
      if (swap) hom else hom[0, ]
    }
  )
  out <- genevintage:::.homologies_for_pair(
    "root", "116", "protein",
    "homo_sapiens", "mus_musculus"
  )
  expect_equal(nrow(out), 1L)
  expect_equal(out$ortholog_id, "ENSMUSG1")
})

test_that("neither species having a homology file is still an error", {
  testthat::local_mocked_bindings(
    .collections = function(root, release, kind, species) stop("no homologies for ", species)
  )
  expect_error(
    genevintage:::.homologies_for_pair(
      "root", "116", "protein",
      "homo_sapiens", "mus_musculus"
    ),
    "no homologies"
  )
})

test_that("a pair genuinely absent from both files is an empty answer", {
  testthat::local_mocked_bindings(
    .collections = function(root, release, kind, species) "file:///dev/null",
    stream_homologies = function(url, to, swap = FALSE, chunk = 2e5) .empty_homology()
  )
  out <- genevintage:::.homologies_for_pair(
    "root", "116", "protein",
    "homo_sapiens", "mus_musculus"
  )
  expect_equal(nrow(out), 0L)
  expect_true(all(c("id", "ortholog_id") %in% names(out)))
})
