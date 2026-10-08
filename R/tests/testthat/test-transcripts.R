# An explicit table keeps these offline; the fetch path is exercised at the end.
tx_map <- data.frame(
  transcript_id = c(
    "ENSMUST00000000001", "ENSMUST00000000002", "ENSMUST00000000003",
    "ENSMUST00000000004", "ENSMUST00000000005"
  ),
  id = c(
    "ENSMUSG00000000001", "ENSMUSG00000000001", "ENSMUSG00000000002",
    "ENSMUSG00000000003", "ENSMUSG00000000004"
  ),
  name = c("Gnai3", "Gnai3", "ENSMUSG00000000002", "Dup", "Dup") # unnamed; shared name
)

gtf_url <- function(lines, env = parent.frame()) {
  f <- withr::local_tempfile(fileext = ".gtf.gz", .local_envir = env)
  con <- gzfile(f, "w")
  writeLines(lines, con)
  close(con)
  paste0("file://", normalizePath(f))
}

test_that("stream_transcripts reads transcript lines", {
  withr::local_options(genevintage.protocol = "https")
  url <- gtf_url(c(
    "#!genome-build test",
    '1\tensembl\tgene\t1\t50\t.\t+\t.\tgene_id "G1"; gene_name "ONE";',
    '1\tensembl\ttranscript\t1\t50\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "ONE";',
    '1\tensembl\texon\t1\t50\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "ONE";',
    '1\tensembl\ttranscript\t1\t40\t.\t+\t.\tgene_id "G1"; transcript_id "T2"; gene_name "ONE";',
    '2\thavana\ttranscript\t5\t9\t.\t-\t.\thavana_gene_id "OTT1"; gene_id "G2"; transcript_id "T3";'
  ))
  expect_equal(stream_transcripts(url), data.frame(
    transcript_id = c("T1", "T2", "T3"), id = c("G1", "G1", "G2"),
    name = c("ONE", "ONE", NA)
  ))
})

test_that("a GTF without transcript lines is read from its exons, once each", {
  withr::local_options(genevintage.protocol = "https")
  # Ensembl before release 75 wrote only exon and CDS lines
  url <- gtf_url(c(
    '1\tprotein_coding\texon\t1\t20\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "ONE";',
    '1\tprotein_coding\tCDS\t5\t20\t.\t+\t0\tgene_id "G1"; transcript_id "T1"; gene_name "ONE";',
    '1\tprotein_coding\texon\t30\t40\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "ONE";',
    '1\tprotein_coding\texon\t1\t20\t.\t+\t.\tgene_id "G1"; transcript_id "T2"; gene_name "ONE";'
  ))
  out <- stream_transcripts(url, chunk = 2) # duplicates span chunk boundaries
  expect_equal(out$transcript_id, c("T1", "T2"))
  expect_equal(out$id, c("G1", "G1"))
})

test_that("transcript2gene ignores versions and keeps input order", {
  out <- transcript2gene(
    c("ENSMUST00000000003.1", "ENSMUST00000000001.4", "NOPE", "ENSMUST00000000002"),
    mapping = tx_map
  )
  expect_equal(out$input, c("ENSMUST00000000003.1", "ENSMUST00000000001.4", "NOPE", "ENSMUST00000000002"))
  expect_equal(out$gene_id, c("ENSMUSG00000000002", "ENSMUSG00000000001", NA, "ENSMUSG00000000001"))
  # an unnamed gene is named by its id, as gene_names() keeps it
  expect_equal(out$gene_name, c("ENSMUSG00000000002", "Gnai3", NA, "Gnai3"))
})

test_that("a GENCODE PAR_Y transcript stays with its own gene", {
  gc <- data.frame(
    transcript_id = c("ENST00000381192.10", "ENST00000381192.10_PAR_Y"),
    gene_id = c("ENSG00000002586.20", "ENSG00000002586.20_PAR_Y"), gene_name = "CD99"
  )
  out <- transcript2gene(c("ENST00000381192.10_PAR_Y", "ENST00000381192.9", "ENST00000381192_PAR_Y"), mapping = gc)
  expect_equal(out$gene_id, c("ENSG00000002586.20_PAR_Y", "ENSG00000002586.20", "ENSG00000002586.20_PAR_Y"))
  # Ensembl has no PAR_Y copy; falls back to the bare id
  ens <- transcript2gene("ENST00000381192.10_PAR_Y", mapping = gc[1, ])
  expect_equal(ens$gene_id, "ENSG00000002586.20")
})

test_that("a handful of transcript ids takes the newest release", {
  skip_no_index()
  # two ids fit no release; dating them picks the smallest
  ids <- c("ENSGALT00000000003", "ENSGALT00000000004")
  seen <- character()
  local_mocked_bindings(.transcript_table = function(a, refresh = FALSE) {
    seen <<- c(seen, a$release)
    data.frame(transcript_id = ids, id = c("G3", "G4"), name = c("PANX2", "RFKL"))
  })
  out <- transcript2gene(ids, species = "chicken", quiet = TRUE)
  newest <- .transcript_candidates("gallus_gallus")$release[1]
  expect_equal(attr(out, "release"), newest)
  expect_equal(seen, newest)
  expect_error(transcript2gene("ENSMUST00000000001", species = "mouse", source = "gencode"), "no gencode release")
})

test_that("transcript ids pick the transcript index", {
  expect_true(.is_transcript(c("ENSMUST00000178537.2", "FBtr0070129", "ENSG00000141510")))
  expect_false(.is_transcript(c("ENSG00000141510", "ENSMUSG00000051951", "ENST00000269305")))
  expect_false(.is_transcript(character()))
})

test_that("transcript2gene accepts its own output as `mapping`", {
  first <- transcript2gene(tx_map$transcript_id, mapping = tx_map)
  again <- transcript2gene("ENSMUST00000000001.9", mapping = first)
  expect_equal(again$gene_name, "Gnai3")
  expect_error(transcript2gene("x", mapping = data.frame(a = 1)), "transcript_id")
})

test_that("a dense matrix sums to genes in first-seen order", {
  m <- mk(c("ENSMUST00000000001.1", "ENSMUST00000000003.2", "ENSMUST00000000002.1"))
  out <- aggregate_transcripts(m, mapping = tx_map, quiet = TRUE)
  expect_equal(rownames(out), c("ENSMUSG00000000001", "ENSMUSG00000000002"))
  expect_equal(colnames(out), c("s1", "s2"))
  expect_equal(unname(out[1, ]), c(1L + 3L, 4L + 6L))
  expect_equal(unname(out[2, ]), c(2L, 5L))
})

test_that("unmapped rows are kept under their id, or dropped", {
  m <- mk(c("ENSMUST00000000001", "ENSMUST00000999999", "ENSMUST00000000002"))
  kept <- aggregate_transcripts(m, mapping = tx_map, quiet = TRUE)
  expect_equal(rownames(kept), c("ENSMUSG00000000001", "ENSMUST00000999999"))
  dropped <- aggregate_transcripts(m, mapping = tx_map, unmapped = "drop", quiet = TRUE)
  expect_equal(rownames(dropped), "ENSMUSG00000000001")
  expect_equal(sum(dropped), sum(m[c(1, 3), ]))
  expect_message(aggregate_transcripts(m, mapping = tx_map), "1 of 3 transcripts have no gene")
})

test_that("labelling by name disambiguates genes that share one", {
  m <- mk(c("ENSMUST00000000004", "ENSMUST00000000005", "ENSMUST00000000001"))
  out <- aggregate_transcripts(m, mapping = tx_map, label = "name", quiet = TRUE)
  expect_equal(rownames(out), c("Dup_ENSMUSG00000000003", "Dup_ENSMUSG00000000004", "Gnai3"))
  expect_equal(nrow(out), 3L) # two genes named Dup stay two rows
})

test_that("mostly unmapped input warns with a catchable class", {
  m <- mk(c("ENSMUST00000999998", "ENSMUST00000999999", "ENSMUST00000000001"))
  expect_warning(
    aggregate_transcripts(m, mapping = tx_map, quiet = TRUE),
    class = "genevintage_unmapped_transcripts"
  )
})

test_that("a sparse matrix comes back sparse and summed", {
  skip_if_not_installed("Matrix")
  sm <- Matrix::Matrix(c(1, 0, 2, 0, 3, 0),
    nrow = 3, sparse = TRUE,
    dimnames = list(c("ENSMUST00000000001", "ENSMUST00000000003", "ENSMUST00000000002"), c("a", "b"))
  )
  out <- aggregate_transcripts(sm, mapping = tx_map, quiet = TRUE)
  expect_s4_class(out, "dgCMatrix")
  expect_true(methods::validObject(out))
  expect_equal(rownames(out), c("ENSMUSG00000000001", "ENSMUSG00000000002"))
  expect_equal(as.matrix(out), matrix(c(3, 0, 0, 3), 2,
    dimnames = list(rownames(out), c("a", "b"))
  ))
})

test_that("data frames sum too, and non-numeric ones are refused", {
  df <- data.frame(s1 = 1:2, s2 = 3:4, row.names = c("ENSMUST00000000001", "ENSMUST00000000002"))
  out <- aggregate_transcripts(df, mapping = tx_map, quiet = TRUE)
  expect_s3_class(out, "data.frame")
  expect_equal(out$s1, 3L)
  df$s2 <- c("a", "b")
  expect_error(aggregate_transcripts(df, mapping = tx_map, quiet = TRUE), "numeric")
})

test_that("bad input is refused, not guessed", {
  expect_error(aggregate_transcripts(1:3, mapping = tx_map), "transcript2gene")
  m <- mk(c("ENSMUST00000000001", "ENSMUST00000000002"))
  expect_error(aggregate_transcripts(m, row_ids = "x", mapping = tx_map), "2 rows")
  skip_no_index()
  expect_error(transcript2gene("ENSGALT00000000003", species = "chicken", source = "gencode"), "no gencode release")
})

test_that("candidate releases are the newest per assembly, newest first", {
  skip_no_index()
  cand <- genevintage:::.transcript_candidates("gallus_gallus")
  expect_equal(cand$assembly[1:2], c("7", "6"))
  expect_true(as.numeric(cand$release[1]) > as.numeric(cand$release[2]))
})

test_that("transcript2gene maps real mouse and chicken transcripts", {
  skip_net()
  mm <- transcript2gene(c("ENSMUST00000178537.2", "ENSMUST00000178862.2"), quiet = TRUE)
  expect_equal(attr(mm, "species"), "mus_musculus")
  expect_true(all(startsWith(mm$gene_id, "ENSMUSG")))
  # pre-GRCg7b ids, renumbered at release 107, so an older assembly is chosen
  gg <- transcript2gene(c("ENSGALT00000000003", "ENSGALT00000000004"), quiet = TRUE)
  expect_equal(attr(gg, "species"), "gallus_gallus")
  expect_false(attr(gg, "assembly") == "7")
  expect_equal(gg$gene_name, c("PANX2", "RFKL"))
})
