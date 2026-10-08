test_that("stream_gtf pulls gene records without touching disk", {
  skip_on_cran()
  skip_if_offline()
  g <- stream_gtf(paste0(
    "https://ftp.ensembl.org/pub/release-115/gtf/",
    "saccharomyces_cerevisiae/",
    "Saccharomyces_cerevisiae.R64-1-1.115.gtf.gz"
  ))
  expect_gt(nrow(g), 6000)
  expect_true(all(grepl("^Y|^Q0|^R00", g$id[1:50])))
  expect_true(any(!is.na(g$name)))
  expect_named(g, c(
    "id", "name", "chr", "start", "end", "strand", "span",
    "id_version", "biotype"
  ))
  # chr is column 1 of the GTF -- yeast numbers its chromosomes in Roman and
  # calls its mitochondrion Mito, which is what mito_genes() reads
  expect_true(all(c("I", "II", "Mito") %in% g$chr))
})

test_that("a streamed mapping feeds gene_names directly", {
  skip_on_cran()
  skip_if_offline()
  skip_no_index()
  # release 63 predates gene feature rows (exercises the attribute-derived
  # path); yeast uses main-site vertebrate numbering, hence 63.
  m <- fetch_mapping("saccharomyces_cerevisiae", 63)
  expect_gt(nrow(m), 6000)
  out <- gene_names(c("YDL246C", "YDR387C", "NOT_A_GENE"), m, warn = FALSE)
  expect_length(out, 3)
  expect_equal(as.vector(out)[1], "SOR2") # has a symbol
  # YDR387C has no symbol at this release: the GTF repeats the identifier, which
  # is how Ensembl says "no name". Both it and an unknown id come back unchanged.
  expect_equal(as.vector(out)[2:3], c("YDR387C", "NOT_A_GENE"))
  expect_equal(attr(out, "mapped"), c(TRUE, FALSE, FALSE))
})

test_that("a release number from the wrong division is refused", {
  skip_no_index()
  # 115 is a vertebrates release; Arabidopsis is numbered by Ensembl Genomes,
  # whose releases stop at 63. (Yeast and fly are not the example to use here:
  # Ensembl's main site carries them too, so they are indexed under vertebrate
  # release numbering and 115 is a perfectly good release for them.)
  expect_error(fetch_mapping("arabidopsis_thaliana", 115), "no index row")
})

test_that("a mapping carries coordinates and strand", {
  skip_net()
  m <- fetch_mapping("saccharomyces_cerevisiae", 115)
  expect_true(all(c("start", "end", "strand", "span") %in% names(m)))
  expect_setequal(unique(m$strand), c("+", "-"))
  expect_false(any(is.na(m$start)))
  expect_true(all(m$end >= m$start))
  expect_equal(m$span, m$end - m$start + 1L)
  # a genomic span, so kilobases -- not the handful of bases a bad parse gives
  expect_gt(median(m$span), 500)
})

test_that("span is the genomic extent, not the exonic length", {
  skip_net()
  # this distinction matters: TPM and FPKM divide by union-exon length, which
  # is shorter than the span wherever a gene has introns. Human TP53 spans
  # ~25kb of genome and has ~2.6kb of exons.
  m <- fetch_mapping("homo_sapiens", 116, assembly = "38")
  tp53 <- m[m$name == "TP53" & !is.na(m$name), ]
  expect_equal(nrow(tp53), 1L)
  expect_gt(tp53$span, 19000) # the span, introns included
  expect_equal(tp53$chr, "17")
  expect_equal(tp53$strand, "-")
})
