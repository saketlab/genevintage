# gene_ids(), annotate_genes() and compare_releases(), offline against a
# hand-built annotation.

ann <- data.frame(
  id = c("ENSG1", "ENSG2", "ENSG3", "ENSG4"),
  name = c("TP53", "BRCA1", "Y_RNA", "Y_RNA"),
  chr = c("17", "17", "1", "2"),
  biotype = c("protein_coding", "protein_coding", "misc_RNA", "misc_RNA"),
  start = c(1L, 10L, 20L, 30L), end = c(5L, 15L, 25L, 35L),
  strand = c("-", "+", "+", "-"), span = c(5L, 6L, 6L, 6L),
  stringsAsFactors = FALSE
)

test_that("names map back to identifiers, keeping every match", {
  out <- gene_ids(c("TP53", "Y_RNA"), mapping = ann, quiet = TRUE)
  # a symbol is not a key: Y_RNA names two genes here and both are kept
  expect_equal(nrow(out), 3L)
  expect_equal(out$input, c("TP53", "Y_RNA", "Y_RNA"))
  expect_setequal(out$id[out$input == "Y_RNA"], c("ENSG3", "ENSG4"))
})

test_that("an input matching nothing is kept with NA rather than dropped", {
  out <- gene_ids(c("TP53", "NOT_A_GENE"), mapping = ann, quiet = TRUE)
  # a lookup returning fewer rows than it was given cannot be joined back
  expect_equal(out$input, c("TP53", "NOT_A_GENE"))
  expect_true(is.na(out$id[2]))
  expect_named(out, c("input", "id", "name"))
})

test_that("all-digit input is read as Entrez without being told", {
  local_mocked_bindings(
    .locate_annotation = function(...) {
      list(
        species = "homo_sapiens",
        source = "ensembl", release = "116",
        assembly = "38", map = ann
      )
    },
    fetch_xrefs = function(...) {
      data.frame(
        id = "ENSG1", xref = "7157",
        stringsAsFactors = FALSE
      )
    }
  )
  out <- gene_ids("7157", species = "human", quiet = TRUE)
  expect_equal(out$id, "ENSG1")
  expect_equal(out$name, "TP53")
})

test_that("a supplied mapping fetches nothing, even for Entrez input", {
  local_mocked_bindings(fetch_xrefs = function(...) stop("must not fetch"))
  # the mapping carries names, not cross-references
  expect_error(
    gene_ids("7157", mapping = ann, quiet = TRUE),
    "carries no cross-references"
  )
  # unless it carries them itself
  withx <- ann
  withx$entrez <- c("7157", "672", NA, NA)
  out <- gene_ids(c("7157", "999"), mapping = withx, quiet = TRUE)
  expect_equal(out$input, c("7157", "999"))
  expect_equal(out$id, c("ENSG1", NA))
})

test_that("from = 'name' is honoured even when the input looks numeric", {
  local_mocked_bindings(fetch_xrefs = function(...) stop("must not fetch xrefs"))
  out <- gene_ids("7157", from = "name", mapping = ann, quiet = TRUE)
  expect_true(is.na(out$id))
})

test_that("annotate_genes returns the whole record, input order kept", {
  out <- annotate_genes(c("BRCA1", "ENSG1", "NOPE"), mapping = ann, quiet = TRUE)
  expect_equal(out$input, c("BRCA1", "ENSG1", "NOPE"))
  expect_equal(out$id, c("ENSG2", "ENSG1", NA))
  expect_true(all(c("chr", "biotype", "start", "end", "strand", "span") %in% names(out)))
  expect_equal(names(out)[1:3], c("input", "id", "name"))
})

test_that("annotate_genes omits coordinates the annotation does not carry", {
  bare <- ann[, c("id", "name", "chr", "biotype")]
  out <- annotate_genes("TP53", mapping = bare, quiet = TRUE)
  expect_named(out, c("input", "id", "name", "chr", "biotype"))
})

test_that("compare_releases names what changed and stays silent about what did not", {
  a <- data.frame(
    id = c("G1", "G2", "G3"), name = c("AAA", "BBB", NA),
    chr = "1", biotype = c("protein_coding", "lncRNA", "lncRNA"),
    stringsAsFactors = FALSE
  )
  b <- data.frame(
    id = c("G1", "G3", "G4"), name = c("AAA", "CCC", "DDD"),
    chr = "1", biotype = c("protein_coding", "lncRNA", "lncRNA"),
    stringsAsFactors = FALSE
  )
  local_mocked_bindings(
    fetch_mapping = function(species, release, ...) if (release == 1) a else b
  )

  d <- compare_releases("human", from = 1, to = 2, quiet = TRUE)
  expect_setequal(d$change, c("removed", "renamed", "added"))
  expect_equal(d$change[d$id == "G2"], "removed") # in `from`, not in `to`
  expect_equal(d$change[d$id == "G4"], "added")
  expect_equal(d$change[d$id == "G3"], "renamed") # NA -> CCC is a change
  expect_false("G1" %in% d$id) # identical rows are omitted
})

test_that("compare_releases separates a biotype change from a symbol change", {
  a <- data.frame(
    id = "G1", name = "AAA", chr = "1", biotype = "lncRNA",
    stringsAsFactors = FALSE
  )
  b <- data.frame(
    id = "G1", name = "AAA", chr = "1", biotype = "protein_coding",
    stringsAsFactors = FALSE
  )
  local_mocked_bindings(
    fetch_mapping = function(species, release, ...) if (release == 1) a else b
  )
  d <- compare_releases("human", from = 1, to = 2, quiet = TRUE)
  expect_equal(d$change, "retyped")
  expect_equal(d$biotype_from, "lncRNA")
  expect_equal(d$biotype_to, "protein_coding")
})

test_that("compare_releases reports a summary unless told not to", {
  a <- data.frame(
    id = "G1", name = "AAA", chr = "1", biotype = "lncRNA",
    stringsAsFactors = FALSE
  )
  local_mocked_bindings(
    fetch_mapping = function(species, release, ...) if (release == 1) a else a[0, ]
  )
  expect_message(compare_releases("human", from = 1, to = 2), "1 removed")
})

test_that("a repeated input yields a row each time", {
  # the result is joined back positionally, so collapsing duplicates breaks it
  out <- gene_ids(c("TP53", "NOPE", "TP53", "NOPE"),
    from = "name",
    mapping = ann, quiet = TRUE
  )
  expect_equal(out$input, c("TP53", "NOPE", "TP53", "NOPE"))
  expect_equal(out$id, c("ENSG1", NA, "ENSG1", NA))
})

test_that("an unmatched row does not change a column's type", {
  # padding with NA_character_ turned start, end and span into strings
  a <- annotate_genes(c("TP53", "NOPE"), mapping = ann, quiet = TRUE)
  expect_type(a$start, "integer")
  expect_type(a$end, "integer")
  expect_type(a$span, "integer")
  # and when nothing matches at all
  b <- annotate_genes("NOPE", mapping = ann, quiet = TRUE)
  expect_type(b$start, "integer")
  expect_equal(nrow(b), 1L)
})

test_that("a name matching several genes keeps input order across the expansion", {
  out <- gene_ids(c("Y_RNA", "TP53"), from = "name", mapping = ann, quiet = TRUE)
  expect_equal(out$input, c("Y_RNA", "Y_RNA", "TP53"))
  expect_equal(out$id[3], "ENSG1")
})

test_that("a non-character input is refused before anything else", {
  expect_error(gene_ids(c(NA, NA), mapping = ann), "must be a character vector")
  expect_error(gene_ids(1:3, mapping = ann), "must be a character vector")
})

test_that("empty and all-NA input come back the same shape", {
  for (x in list(character(), NA_character_, c(NA_character_, NA_character_))) {
    out <- gene_ids(x, from = "name", mapping = ann, quiet = TRUE)
    expect_equal(nrow(out), length(x))
    expect_named(out, c("input", "id", "name"))
  }
})

test_that("compare_releases does not depend on the order of duplicate ids", {
  # Ensembl gene ids are unique within a release, so this is a guard on
  # malformed input rather than a case the fetch path produces
  a1 <- data.frame(
    id = c("G1", "G1"), name = c("A", "B"), chr = "1",
    biotype = "protein_coding", stringsAsFactors = FALSE
  )
  b1 <- data.frame(
    id = "G1", name = "A", chr = "1", biotype = "protein_coding",
    stringsAsFactors = FALSE
  )
  local_mocked_bindings(
    fetch_mapping = function(species, release, ...) if (release == 1) a1 else b1
  )
  fwd <- compare_releases("human", from = 1, to = 2, quiet = TRUE)
  local_mocked_bindings(
    fetch_mapping = function(species, release, ...) if (release == 1) a1[2:1, ] else b1
  )
  rev <- compare_releases("human", from = 1, to = 2, quiet = TRUE)
  expect_equal(nrow(fwd) + nrow(rev), 1L) # exactly one order reports a rename
})

test_that("a rename and a retype at once are reported as a rename", {
  a <- data.frame(
    id = "G1", name = "A", chr = "1", biotype = "lncRNA",
    stringsAsFactors = FALSE
  )
  b <- data.frame(
    id = "G1", name = "B", chr = "1", biotype = "protein_coding",
    stringsAsFactors = FALSE
  )
  local_mocked_bindings(
    fetch_mapping = function(species, release, ...) if (release == 1) a else b
  )
  d <- compare_releases("human", from = 1, to = 2, quiet = TRUE)
  expect_equal(d$change, "renamed")
  # the biotype columns still carry both sides, so the retype is not lost
  expect_equal(c(d$biotype_from, d$biotype_to), c("lncRNA", "protein_coding"))
})
