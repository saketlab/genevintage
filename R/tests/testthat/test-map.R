test_that("names are returned, ids are retained, order never moves", {
  m <- data.frame(
    id = c("ENSG00000141510", "ENSG00000284733"),
    name = c("TP53", "ENSG00000284733")
  ) # 2nd has no symbol
  x <- c("ENSG00000284733", "ENSG00000141510.16", "ENSG00000999999", "TP53")
  out <- gene_names(x, m, warn = FALSE)
  expect_equal(as.vector(out), c("ENSG00000284733", "TP53", "ENSG00000999999", "TP53"))
  expect_equal(attr(out, "mapped"), c(FALSE, TRUE, FALSE, FALSE))
  expect_length(out, length(x))
})

test_that("a name equal to its id counts as no name", {
  m <- data.frame(id = "ENSG00000284733", name = "ENSG00000284733")
  expect_false(attr(gene_names("ENSG00000284733", m, warn = FALSE), "mapped"))
  m2 <- data.frame(id = c("A", "B"), name = c("", NA))
  expect_equal(as.vector(gene_names(c("A", "B"), m2, warn = FALSE)), c("A", "B"))
})

test_that("versioned input matches an unversioned mapping and vice versa", {
  expect_equal(as.vector(gene_names("ENSG00000141510.16",
    data.frame(id = "ENSG00000141510", name = "TP53"),
    warn = FALSE
  )), "TP53")
  expect_equal(as.vector(gene_names("ENSG00000141510",
    data.frame(id = "ENSG00000141510.16", name = "TP53"),
    warn = FALSE
  )), "TP53")
})

test_that("collisions warn, and `unique` resolves them without touching ids", {
  m <- data.frame(id = c("ENSG1", "ENSG2"), name = c("DUP", "DUP"))
  expect_warning(gene_names(c("ENSG1", "ENSG2"), m), "duplicate")
  u <- gene_names(c("ENSG1", "ENSG2"), m, unique = TRUE, warn = FALSE)
  expect_equal(anyDuplicated(u), 0L)
  # unmapped ids are already unique -- they must not be suffixed
  m2 <- data.frame(id = "ENSG1", name = "DUP")
  expect_equal(
    as.vector(gene_names(c("ENSG1", "ENSGX"), m2, unique = TRUE, warn = FALSE)),
    c("DUP", "ENSGX")
  )
})

test_that("the collision warning carries a class a caller can catch alone", {
  m <- data.frame(id = c("ENSG1", "ENSG2"), name = c("DUP", "DUP"))
  caught <- character()
  withCallingHandlers(
    gene_names(c("ENSG1", "ENSG2"), m),
    genevintage_duplicate_names = function(w) {
      caught <<- c(caught, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(caught, 1L)
  expect_match(caught, "duplicate")
})

test_that("a wrong-species mapping warns instead of silently returning ids", {
  m <- data.frame(id = "ENSMUSG00000051951", name = "Xkr4")
  expect_warning(gene_names(c("ENSG00000141510", "ENSG00000284733"), m), "wrong species")
})

# Real references: half the pairs become the mapping, all ids go in.
test_that("round-trips against real 10x references", {
  skip_if_not(dir.exists(truth_dir))
  fs <- list.files(truth_dir, pattern = "\\.tsv\\.gz$")
  skip_if(length(fs) == 0)
  for (f in fs) {
    tr <- read_truth(f)
    half <- seq_len(nrow(tr)) %% 2 == 1
    out <- gene_names(tr$id, tr[half, ], warn = FALSE)
    real <- .is_name(tr$name, tr$id)

    expect_length(out, nrow(tr))
    # in the half we supplied, every real name is recovered
    expect_equal(as.vector(out[half & real]), tr$name[half & real], label = f)
    # everything else keeps its id, exactly as given
    expect_equal(as.vector(out[!half | !real]), tr$id[!half | !real], label = f)
  }
})

test_that("references really do repeat the id when there is no symbol", {
  skip_if_not(dir.exists(truth_dir))
  fs <- list.files(truth_dir, pattern = "\\.tsv\\.gz$")
  skip_if(length(fs) == 0)
  frac <- vapply(fs, function(f) {
    tr <- read_truth(f)
    mean(tr$id == tr$name)
  }, numeric(1))
  expect_true(any(frac > 0.3)) # human/chicken/macaque/arabidopsis do
  expect_true(all(vapply(fs, function(f) {
    tr <- read_truth(f)
    mean(!nzchar(tr$name))
  }, numeric(1)) < 0.01)) # blanks are not the idiom
})

test_that("factor mapping columns are coerced, not silently turned into codes", {
  m <- data.frame(id = factor("ENSG1"), name = factor("TP53"))
  expect_equal(as.vector(gene_names("ENSG1", m, warn = FALSE)), "TP53")
})

test_that("conflicting duplicate mapping rows warn instead of first-winning silently", {
  m <- data.frame(id = c("ENSG1", "ENSG1"), name = c("OLD", "NEW"))
  expect_warning(gene_names("ENSG1", m), "conflicting")
  # agreeing duplicates are normal and must stay quiet
  ok <- data.frame(id = c("ENSG1", "ENSG1"), name = c("TP53", "TP53"))
  expect_silent(gene_names("ENSG1", ok, warn = TRUE))
})

test_that("unique = TRUE really is unique, even with repeated input ids", {
  m <- data.frame(id = "ENSG1", name = "TP53")
  out <- gene_names(c("ENSG1", "ENSG1", "ENSGX", "ENSGX"), m, unique = TRUE, warn = FALSE)
  expect_equal(anyDuplicated(out), 0L)
})

test_that("a PAR_Y tagged id resolves against an untagged mapping", {
  m <- data.frame(id = "ENSG00000182378", name = "ZBED1")
  expect_equal(as.vector(gene_names("ENSG00000182378.14_PAR_Y", m, warn = FALSE)), "ZBED1")
})

test_that("no input still returns an empty character vector", {
  # ifelse() takes its type from the test, so an empty input used to come back
  # as logical(0) -- a silent change to the documented return type
  for (m in list(
    data.frame(id = "a", name = "b"),
    data.frame(id = character(), name = character())
  )) {
    out <- gene_names(character(), m, warn = FALSE)
    expect_type(out, "character")
    expect_length(out, 0L)
    expect_type(attr(out, "mapped"), "logical")
  }
  expect_type(gene_names(character(), data.frame(id = character(), name = character()),
    unique = TRUE
  ), "character")
})
