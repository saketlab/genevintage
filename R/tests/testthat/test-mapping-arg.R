# A caller's own annotation: a local GTF, a custom build, or a species the
# index does not cover. Nothing is detected and nothing is downloaded, which
# makes all of this offline.

fake <- data.frame(
  id = c("SAL0001", "SAL0002", "SAL0003", "SAL0004", "SAL0005", "SAL0006"),
  name = c("COX1", "RPL10", "MRPL1", "XIST-like", "tRNA-Ala", "5S_rRNA"),
  chr = c("MT", "3", "3", "X", "7", "7"),
  biotype = c(
    "protein_coding", "protein_coding", "protein_coding",
    "lncRNA", "tRNA", "rRNA"
  ),
  stringsAsFactors = FALSE
)

test_that("a supplied mapping is used instead of fetching one", {
  # a salamander-shaped identifier no scheme knows, on a species not indexed
  expect_equal(mito_genes(mapping = fake, quiet = TRUE)$name, "COX1")
  expect_equal(sex_genes(mapping = fake, quiet = TRUE)$name, "XIST-like")
  expect_equal(trna_genes(mapping = fake, quiet = TRUE)$name, "tRNA-Ala")
  expect_equal(rrna_genes(mapping = fake, quiet = TRUE)$name, "5S_rRNA")
  expect_equal(lncrna_genes(mapping = fake, quiet = TRUE)$name, "XIST-like")
  expect_equal(nrow(coding_genes(mapping = fake, quiet = TRUE)), 3L)
  expect_equal(nrow(noncoding_genes(mapping = fake, quiet = TRUE)), 3L)
})

test_that("the ribosomal rules apply to a supplied mapping too", {
  expect_equal(ribosomal_genes(mapping = fake, quiet = TRUE)$name, "RPL10")
  expect_equal(ribosomal_genes(mapping = fake, which = "mito", quiet = TRUE)$name, "MRPL1")
  expect_setequal(
    ribosomal_genes(mapping = fake, which = "all", quiet = TRUE)$name,
    c("RPL10", "MRPL1", "5S_rRNA")
  )
})

test_that("ids still subset a supplied mapping, by identifier or by name", {
  expect_equal(
    mito_genes(c("SAL0001", "SAL0002"), mapping = fake, quiet = TRUE)$input,
    "SAL0001"
  )
  expect_equal(mito_genes("COX1", mapping = fake, quiet = TRUE)$input, "COX1")
  expect_equal(mito_genes("COX1", mapping = fake, quiet = TRUE)$id, "SAL0001")
})

test_that("gene_biotypes and genes_by_biotype read a supplied mapping", {
  tb <- gene_biotypes(mapping = fake, quiet = TRUE)
  expect_named(tb, c("biotype", "class", "n"))
  expect_equal(sum(tb$n), nrow(fake))
  expect_equal(nrow(genes_by_biotype("lncRNA", mapping = fake, quiet = TRUE)), 1L)
  expect_error(
    genes_by_biotype("no_such_biotype", mapping = fake, quiet = TRUE),
    "no biotype no_such_biotype"
  )
})

test_that("a mapping missing what these functions read is refused", {
  expect_error(mito_genes(mapping = data.frame(x = 1)), "needs column")
  expect_error(
    mito_genes(mapping = data.frame(id = "a", name = "b")),
    "chr, biotype"
  )
  expect_error(mito_genes(mapping = "not a frame"), "must be a data frame")
})

test_that("no network or index is touched when a mapping is supplied", {
  # if either were consulted these would error, since the species is unknown
  testthat::local_mocked_bindings(
    fetch_mapping = function(...) stop("fetch_mapping must not be called"),
    detect_species = function(...) stop("detect_species must not be called")
  )
  expect_equal(nrow(mito_genes(mapping = fake, quiet = TRUE)), 1L)
})

test_that("coordinates come through when the annotation carries them", {
  m <- cbind(fake,
    start = seq(100, 600, by = 100), end = seq(200, 700, by = 100),
    strand = "+", span = 101L
  )
  g <- coding_genes(mapping = m, quiet = TRUE)
  expect_true(all(c("start", "end", "strand", "span") %in% names(g)))
  expect_equal(unique(g$span), 101L)
  # and are simply absent when it does not
  expect_false("span" %in% names(coding_genes(mapping = fake, quiet = TRUE)))
})

# The gaps below are all reachable with mapping =, so they test the branch
# logic offline rather than needing a fetched annotation.

test_that("which= is read case-insensitively and can name both chromosomes", {
  m <- data.frame(
    id = c("g1", "g2", "g3"), name = c("XIST", "SRY", "TP53"),
    chr = c("X", "Y", "17"), biotype = "protein_coding", stringsAsFactors = FALSE
  )
  expect_equal(sex_genes(mapping = m, which = "y", quiet = TRUE)$name, "SRY")
  expect_equal(sex_genes(mapping = m, which = "Y", quiet = TRUE)$name, "SRY")
  expect_setequal(
    sex_genes(mapping = m, which = c("X", "Y"), quiet = TRUE)$name,
    c("XIST", "SRY")
  )
  expect_error(sex_genes(mapping = m, which = "W", quiet = TRUE), "no chromosome W")
})

test_that("a missing or empty biotype is classified, not dropped", {
  m <- data.frame(
    id = c("g1", "g2", "g3", "g4"), name = c("A", "B", "C", "D"),
    chr = "1", biotype = c("protein_coding", NA, "", "TEC"),
    stringsAsFactors = FALSE
  )
  # NA and "" are unclassifiable, not non-coding -- "other" is the honest bucket
  expect_setequal(
    genes_by_biotype("other", mapping = m, quiet = TRUE)$name,
    c("B", "C", "D")
  )
  expect_equal(genes_by_biotype("coding", mapping = m, quiet = TRUE)$name, "A")
  expect_equal(nrow(noncoding_genes(mapping = m, quiet = TRUE)), 0L)
})

test_that("an annotation that does not draw the lncRNA class can be told to stay quiet", {
  # fly and yeast file everything non-coding under a generic ncRNA biotype
  m <- data.frame(
    id = c("g1", "g2"), name = c("A", "B"), chr = "1",
    biotype = c("ncRNA", "protein_coding"), stringsAsFactors = FALSE
  )
  expect_message(lncrna_genes(mapping = m), "ncRNA")
  expect_silent(g <- lncrna_genes(mapping = m, quiet = TRUE))
  expect_equal(nrow(g), 0L)
})

test_that("a name that only looks ribosomal is excluded", {
  # the S6 kinase is protein-coding and matches ^RP[SL]; the protein/mito/all
  # branches are covered above, so this adds only the name rule and which="rrna"
  m <- rbind(fake, data.frame(
    id = "SAL0007", name = "RPS6KA1", chr = "3",
    biotype = "protein_coding", stringsAsFactors = FALSE
  ))
  expect_equal(ribosomal_genes(mapping = m, quiet = TRUE)$name, "RPL10")
  expect_equal(ribosomal_genes(mapping = m, which = "rrna", quiet = TRUE)$name, "5S_rRNA")
})

test_that("hemoglobin excludes the lookalikes and the pseudogenes", {
  m <- data.frame(
    id = paste0("g", 1:6),
    name = c("HBB", "HBA1", "HBEGF", "HBS1L", "HBAP1", "HBP1"),
    chr = "11",
    biotype = c(
      "protein_coding", "protein_coding", "protein_coding",
      "protein_coding", "unprocessed_pseudogene", "protein_coding"
    ),
    stringsAsFactors = FALSE
  )
  # HBEGF, HBS1L and HBP1 fall outside the subunit letters; HBAP1 matches the
  # name and is removed by the biotype
  expect_setequal(hemoglobin_genes(mapping = m, quiet = TRUE)$name, c("HBB", "HBA1"))
})

test_that("hemoglobin is found under each clade's spelling", {
  for (nm in c("HBA1", "Hba-a1", "hbaa1", "Hbb-bs", "hbbe2")) {
    m <- data.frame(
      id = "g1", name = nm, chr = "11", biotype = "protein_coding",
      stringsAsFactors = FALSE
    )
    expect_equal(nrow(hemoglobin_genes(mapping = m, quiet = TRUE)), 1L, label = nm)
  }
})

test_that("the MHC is found under each clade's name and nowhere else", {
  for (nm in c("HLA-A", "HLA-DRB1", "H2-K1", "H2-Ab1", "mhc1uba", "mhc2b")) {
    m <- data.frame(
      id = "g1", name = nm, chr = "6", biotype = "protein_coding",
      stringsAsFactors = FALSE
    )
    expect_equal(nrow(mhc_genes(mapping = m, quiet = TRUE)), 1L, label = nm)
  }
  # the hyphen keeps H2- clear of the H2A histones
  for (nm in c("H2AC1", "H2BC1", "MHCX")) {
    m <- data.frame(
      id = "g1", name = nm, chr = "6", biotype = "protein_coding",
      stringsAsFactors = FALSE
    )
    expect_equal(nrow(mhc_genes(mapping = m, quiet = TRUE)), 0L, label = nm)
  }
})

test_that("MHC pseudogenes and lncRNAs carrying an HLA name are excluded", {
  m <- data.frame(
    id = paste0("g", 1:3), name = c("HLA-A", "HLA-H", "HLA-F-AS1"), chr = "6",
    biotype = c("protein_coding", "unprocessed_pseudogene", "lncRNA"),
    stringsAsFactors = FALSE
  )
  expect_equal(mhc_genes(mapping = m, quiet = TRUE)$name, "HLA-A")
})

test_that("a species with no MHC says so rather than implying an empty one", {
  m <- data.frame(
    id = "g1", name = "YDL246C", chr = "IV",
    biotype = "protein_coding", stringsAsFactors = FALSE
  )
  expect_message(mhc_genes(mapping = m, species = "yeast"), "no MHC genes")
  expect_equal(nrow(mhc_genes(mapping = m, quiet = TRUE)), 0L)
})
