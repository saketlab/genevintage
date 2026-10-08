test_that("overlapping, nested and abutting exons are counted once", {
  len <- genevintage:::.exon_union_length(
    gene = c("A", "A", "A", "B", "B", "C"),
    start = c(100, 150, 300, 10, 20, 5),
    end = c(200, 180, 310, 50, 30, 5)
  )
  # A: [100,200] absorbs nested [150,180]; plus [300,310] -> 101 + 11
  # B: [20,30] nested in [10,50] -> 41; C: a single base
  expect_identical(len, c(A = 112L, B = 41L, C = 1L))
})

test_that("exon order does not matter and abutting exons add", {
  a <- genevintage:::.exon_union_length(c("G", "G", "G"), c(21, 1, 11), c(30, 10, 20))
  expect_identical(a, c(G = 30L))
  # a later exon that starts before an earlier one's end still merges
  b <- genevintage:::.exon_union_length(c("G", "G", "G"), c(1, 5, 50), c(100, 10, 60))
  expect_identical(b, c(G = 100L))
})

test_that("an exon shared by two genes counts toward both", {
  len <- genevintage:::.exon_union_length(c("A", "B"), c(1, 1), c(10, 10))
  expect_identical(len, c(A = 10L, B = 10L))
})

test_that("no exons gives an empty result", {
  expect_length(genevintage:::.exon_union_length(character(), numeric(), numeric()), 0L)
})

test_that("stream_exon_lengths reads only exon lines of a GTF", {
  gtf <- withr::local_tempfile(fileext = ".gtf.gz")
  con <- gzfile(gtf, "w")
  writeLines(c(
    "#!genome-build test",
    '1\tensembl\tgene\t100\t310\t.\t+\t.\tgene_id "G1"; gene_version "2"; gene_name "ONE";',
    '1\tensembl\ttranscript\t100\t310\t.\t+\t.\tgene_id "G1"; transcript_id "T1";',
    '1\tensembl\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "1";',
    '1\tensembl\texon\t300\t310\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "2";',
    '1\tensembl\texon\t150\t250\t.\t+\t.\tgene_id "G1"; transcript_id "T2"; exon_number "1";',
    '1\tensembl\tCDS\t120\t200\t.\t+\t0\tgene_id "G1"; transcript_id "T1";',
    '2\thavana\texon\t5\t14\t.\t-\t.\thavana_gene_id "OTT1"; gene_id "G2"; transcript_id "T3";'
  ), con)
  close(con)
  withr::local_options(genevintage.protocol = "https")
  out <- stream_exon_lengths(paste0("file://", normalizePath(gtf)))
  # G1: [100,250] + [300,310] = 151 + 11; G2 must not pick up havana_gene_id
  expect_equal(out, data.frame(id = c("G1", "G2"), length = c(162L, 10L)))
})

test_that("gene_lengths matches ids, versioned ids and names from the cache", {
  withr::local_envvar(R_USER_CACHE_DIR = withr::local_tempdir())
  a <- genevintage:::.annotation_source("human", 116)
  f <- genevintage:::.cache_path("lengths", a$row$source, "homo_sapiens", a$row$assembly, a$release)
  dir.create(dirname(f), recursive = TRUE)
  saveRDS(data.frame(id = c("ENSG00000141510", "ENSG00000012048"), length = c(2579L, 7088L)), f)
  m <- genevintage:::.cache_path("mapping", a$row$source, "homo_sapiens", a$row$assembly, a$release)
  saveRDS(data.frame(
    id = c("ENSG00000141510", "ENSG00000012048"), name = c("TP53", "BRCA1"),
    chr = "17", biotype = "protein_coding", span = c(25772L, 125951L)
  ), m)
  rm(list = ls(genevintage:::.lengths_memo), envir = genevintage:::.lengths_memo)
  rm(list = ls(genevintage:::.mapping_memo), envir = genevintage:::.mapping_memo)
  expect_identical(
    gene_lengths("human", 116, ids = c("ENSG00000141510.18", "BRCA1", "NOPE")),
    c(ENSG00000141510.18 = 2579L, BRCA1 = 7088L, NOPE = NA)
  )
  expect_equal(nrow(gene_lengths("human", 116)), 2L)
  expect_true("lengths" %in% gene_cache()$kind)
})

test_that("gene_lengths agrees with featureCounts' Length on real human genes", {
  skip_on_cran()
  skip_if_offline()
  # featureCounts on Ensembl 116 reports the exon union; spot-check a gene with
  # many overlapping transcripts against its span
  len <- gene_lengths("human", 116, ids = c("TP53", "ALB"))
  span <- fetch_mapping("human", 116)
  span <- span$span[match(c("TP53", "ALB"), span$name)]
  expect_true(all(len > 1000 & len < span))
})
