# Real pipeline outputs with a provider-declared annotation vintage.
pipe_dir <- test_path("pipelines")
ptruth <- if (file.exists(file.path(pipe_dir, "truth.csv"))) {
  utils::read.csv(file.path(pipe_dir, "truth.csv"), colClasses = "character")
} else {
  NULL
}

pipe_ids <- function(f) read_fx(f, dir = pipe_dir)

skip_no_pipelines <- function() {
  skip_if_not(!is.null(ptruth) && nrow(ptruth) > 0, "no pipeline fixtures")
}

test_that("pipeline matrices resolve to the right species", {
  skip_no_pipelines()
  for (i in seq_len(nrow(ptruth))) {
    expect_equal(detect_species(pipe_ids(ptruth$file[i]))$species[1],
      ptruth$species[i],
      label = ptruth$file[i]
    )
  }
})

test_that("the declared annotation is never ruled out and ranks near the top", {
  skip_no_pipelines()
  for (i in seq_len(nrow(ptruth))) {
    out <- suppressMessages(suppressWarnings(
      detect_release(pipe_ids(ptruth$file[i]), species = ptruth$species[i])
    ))
    k <- which(out$source == ptruth$source[i] & out$release == ptruth$release[i])
    skip_if(!length(k), paste("not indexed:", ptruth$source[i], ptruth$release[i]))
    expect_true(out$feasible[k], label = paste(ptruth$file[i], "feasible"))
    # Not rank 1: adjacent vintages are often indistinguishable by gene count
    # alone -- GENCODE M23 and M25 differ by 20 genes out of 55000. Top 5 is
    # what the fingerprint can honestly promise.
    expect_lte(k, 5L, label = paste(ptruth$file[i], "rank"))
  }
})

test_that("the annotation family is identified, not just the vintage", {
  skip_no_pipelines()
  for (i in seq_len(nrow(ptruth))) {
    out <- suppressMessages(suppressWarnings(
      detect_release(pipe_ids(ptruth$file[i]), species = ptruth$species[i])
    ))
    expect_equal(out$source[1], ptruth$source[i],
      label = paste(ptruth$file[i], "source")
    )
  }
})

test_that("versioned pipeline identifiers parse", {
  skip_no_pipelines()
  v <- ptruth[as.numeric(ptruth$versioned) > 0.9, ]
  skip_if(nrow(v) == 0)
  for (i in seq_len(nrow(v))) {
    expect_gt(mean(!is.na(.id_version(pipe_ids(v$file[i])))), 0.9, label = v$file[i])
  }
})
