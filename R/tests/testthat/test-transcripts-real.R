# Real transcript-level quantifications; provenance in transcripts/README.md.
tx_dir <- test_path("transcripts")
ttruth <- if (file.exists(file.path(tx_dir, "truth.csv"))) {
  utils::read.csv(file.path(tx_dir, "truth.csv"), colClasses = "character")
} else {
  NULL
}
skip_no_tx <- function() skip_if_not(!is.null(ttruth), "no transcript fixtures")
skip_no_tx_index <- function(species) {
  skip_if_not(any(.tx_fp()$species == species), paste("no transcript index for", species))
}

read_tx <- function(f) {
  read_fx(f, tx_dir, function(k) {
    utils::read.delim(gzfile(k), check.names = FALSE, stringsAsFactors = FALSE)
  })
}

counts <- function(x, cols) {
  m <- as.matrix(x[cols])
  rownames(m) <- x[[1]]
  m
}

# the declared release is feasible and scores with the best
expect_declared_release <- function(ids, row) {
  rel <- suppressWarnings(detect_release(ids, species = row$species))
  k <- which(rel$source == row$source & rel$release == row$release)
  expect_true(length(k) == 1L, label = paste(row$file, "indexed"))
  expect_true(isTRUE(rel$feasible[k]), label = paste(row$file, "feasible"))
  expect_lte(rel$dist[k], rel$dist[1] + 1e-9)
}

test_that("every fixture's species is detected from its transcript ids", {
  skip_no_tx()
  skip_no_index()
  for (i in seq_len(nrow(ttruth))) {
    ids <- read_tx(ttruth$file[i])[[1]]
    expect_equal(detect_species(ids)$species[1], ttruth$species[i], label = ttruth$file[i])
  }
})

test_that("the declared annotation is detected from the transcript ids alone", {
  skip_no_tx()
  dated <- ttruth[nzchar(ttruth$release), ]
  for (i in seq_len(nrow(dated))) {
    skip_no_tx_index(dated$species[i])
    expect_declared_release(read_tx(dated$file[i])[[1]], dated[i, ])
  }
  for (f in c("gencode_v27_salmon.tsv.gz", "geuvadis_gencode_v12.tsv.gz")) {
    row <- ttruth[ttruth$file == f, ]
    best <- suppressWarnings(detect_release(read_tx(f)[[1]]))[1, ]
    expect_equal(unlist(best[c("source", "release", "assembly")]),
      unlist(row[c("source", "release", "assembly")]),
      label = f
    )
  }
})

test_that("a filtered transcriptome is placed by its newest id", {
  skip_no_tx()
  skip_no_tx_index("mus_musculus")
  # a filtered set whose top id is M23's qmax
  ids <- read_tx("tasic_mouse_salmon.tsv.gz")$Name
  best <- suppressWarnings(detect_release(ids))[1, ]
  expect_equal(
    unlist(best[c("source", "release", "assembly")]),
    c(source = "gencode", release = "M23", assembly = "38")
  )
})

test_that("summing RSEM's isoforms reproduces RSEM's own gene counts", {
  skip_no_tx()
  iso <- read_tx("gencode_v27_rsem_ERR188021.tsv.gz")
  gen <- read_tx("gencode_v27_rsem_ERR188021_genes.tsv.gz")
  m <- counts(iso, "expected_count")
  out <- aggregate_transcripts(m, mapping = iso[c("transcript_id", "gene_id")], quiet = TRUE)
  expect_setequal(rownames(out), gen$gene_id)
  expect_equal(sum(out), sum(m))
  # both files round to two decimals
  n_iso <- as.vector(table(iso$gene_id)[gen$gene_id])
  drift <- abs(out[gen$gene_id, 1] - gen$expected_count)
  expect_true(all(drift <= 0.005 * (n_iso + 1) + 1e-9))
  skip_if_not_installed("Matrix")
  sp <- aggregate_transcripts(Matrix::Matrix(m, sparse = TRUE),
    mapping = iso[c("transcript_id", "gene_id")], quiet = TRUE
  )
  expect_equal(as.matrix(sp), out)
})

test_that("GENCODE 27 transcripts are dated and summed to RSEM's genes", {
  skip_no_tx()
  skip_net()
  skip_no_tx_index("homo_sapiens")
  row <- ttruth[ttruth$file == "gencode_v27_salmon.tsv.gz", ]
  s <- read_tx(row$file)
  iso <- read_tx("gencode_v27_rsem_ERR188021.tsv.gz")

  tab <- transcript2gene(s$Name, quiet = TRUE)
  expect_equal(
    unlist(attributes(tab)[c("species", "source", "release", "assembly")]),
    c(species = "homo_sapiens", source = "gencode", release = "27", assembly = "38")
  )
  # RSEM on the same transcriptome gives the truth
  truth <- iso$gene_id[match(s$Name, iso$transcript_id)]
  expect_false(anyNA(truth))
  expect_identical(tab$gene_id, truth)

  m <- counts(s, c("ERR188021", "ERR188088"))
  by_id <- aggregate_transcripts(m, quiet = TRUE)
  expect_equal(by_id, rowsum(m, truth, reorder = FALSE))
  expect_equal(colSums(by_id), colSums(m))

  by_name <- aggregate_transcripts(m, label = "name", quiet = TRUE)
  expect_equal(unname(by_name), unname(by_id))
  gapdh <- unique(truth[tab$gene_name == "GAPDH"])
  expect_length(gapdh, 1L)
  expect_equal(by_name["GAPDH", ], colSums(m[truth == gapdh, ]))
  par <- grepl("_PAR_Y$", s$Name)
  expect_true(any(par))
  expect_true(all(endsWith(tab$gene_id[par], "_PAR_Y")))
})

test_that("Ensembl 98 fly transcripts map to the release's genes", {
  skip_no_tx()
  skip_net()
  skip_no_tx_index("drosophila_melanogaster")
  row <- ttruth[ttruth$file == "ensembl98_fly_salmon.tsv.gz", ]
  q <- read_tx(row$file)

  tab <- transcript2gene(q$Name, quiet = TRUE)
  expect_equal(attr(tab, "species"), "drosophila_melanogaster")
  expect_equal(attr(tab, "source"), "ensembl")
  known <- !is.na(q$gene_id)
  expect_identical(tab$gene_id[known], q$gene_id[known])

  m <- counts(q, "NumReads")
  out <- aggregate_transcripts(m, quiet = TRUE)
  gid <- .coalesce(tab$gene_id, q$Name)
  expect_equal(out, rowsum(m, gid, reorder = FALSE))
  expect_equal(sum(out), sum(m))
})

test_that("Tasic mouse transcripts map to the authors' genes", {
  skip_no_tx()
  skip_net()
  q <- read_tx("tasic_mouse_salmon.tsv.gz")
  tab <- transcript2gene(q$Name, quiet = TRUE)
  # CHR drops scaffold transcripts; ALL keeps them
  expect_equal(
    unlist(attributes(tab)[c("species", "source", "release", "assembly")]),
    c(species = "mus_musculus", source = "gencode_all", release = "M23", assembly = "38")
  )
  # in the salmon index but in no M23 or Ensembl 98 GTF
  miss <- c("ENSMUST00000229179.1", "ENSMUST00000230581.1")
  ok <- !is.na(tab$gene_id)
  expect_setequal(q$Name[!ok], miss)
  expect_identical(.bare(tab$gene_id[ok]), q$gene_id[ok])

  m <- counts(q, c("SRR7311351", "SRR7311386"))
  by_name <- aggregate_transcripts(m, label = "name", quiet = TRUE)
  expect_equal(colSums(by_name), colSums(m))
  expect_equal(nrow(by_name), length(unique(q$gene_id[ok])) + length(miss))
  expect_true(all(miss %in% rownames(by_name)))
  expect_equal(by_name["Snap25", ], colSums(m[tab$gene_name %in% "Snap25", ]))
})

test_that("GEUVADIS GENCODE 12 transcripts land on GRCh37 and sum to the provider's genes", {
  skip_no_tx()
  skip_net()
  skip_no_tx_index("homo_sapiens")
  row <- ttruth[ttruth$file == "geuvadis_gencode_v12.tsv.gz", ]
  g <- read_tx(row$file)

  tab <- transcript2gene(g$TargetID, quiet = TRUE)
  expect_equal(attr(tab, "assembly"), "37")
  expect_equal(attr(tab, "source"), "gencode")
  expect_identical(tab$gene_id, g$gene_id)

  m <- counts(g, names(g)[3:4])
  expect_equal(aggregate_transcripts(m, quiet = TRUE), rowsum(m, g$gene_id, reorder = FALSE))
})
