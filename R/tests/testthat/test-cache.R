# The cache lives in the user's filespace, so its behaviour there is a contract:
# it must not appear until something is written, and it must be removable.

# Point the cache somewhere disposable for the whole file. tools::R_user_dir()
# reads the environment variable each call, so this is enough.
local_cache <- function(env = parent.frame()) {
  d <- withr::local_tempdir(.local_envir = env)
  withr::local_envvar(c(R_USER_CACHE_DIR = d), .local_envir = env)
  tools::R_user_dir("genevintage", "cache")
}

# a cache directory holding the named files, ready to inspect
seed_cache <- function(names, env = parent.frame()) {
  d <- local_cache(env)
  dir.create(d, recursive = TRUE)
  for (n in names) saveRDS(data.frame(id = "x"), file.path(d, n))
  d
}

test_that("looking at an empty cache does not create one", {
  d <- local_cache()
  expect_false(dir.exists(d))
  out <- gene_cache()
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0L)
  # merely reporting must leave nothing behind in the user's filespace
  expect_false(dir.exists(d))
})

test_that("cached files are described, not just listed", {
  seed_cache(c(
    "mapping-ensembl-homo_sapiens-38-116.rds",
    "mapping-gencode_all-homo_sapiens-38-44.rds",
    "mapping-ensembl-anas_platyrhynchos-1-116.rds",
    "orthologs-homo_sapiens-mus_musculus-116-protein.rds"
  ))
  x <- gene_cache()
  x <- x[order(x$file), ]
  # species and sources both carry underscores, so the separator must be a
  # character the fields cannot contain -- otherwise the pair of an ortholog
  # file collapses into the species column and a source has to be enumerated
  expect_equal(x$species, c(
    "anas_platyrhynchos", "homo_sapiens",
    "homo_sapiens", "homo_sapiens"
  ))
  expect_equal(x$to, c(NA, NA, NA, "mus_musculus"))
  expect_equal(x$release, c("116", "116", "44", "116"))
  expect_equal(x$kind, c("mapping", "mapping", "mapping", "orthologs"))
  expect_true(all(x$bytes > 0))
})

test_that("clear removes everything, or only what is stale", {
  d <- seed_cache(c(
    "mapping-ensembl-homo_sapiens-38-116.rds",
    "mapping-ensembl-mus_musculus-39-116.rds"
  ))
  old <- file.path(d, "mapping-ensembl-homo_sapiens-38-116.rds")
  new <- file.path(d, "mapping-ensembl-mus_musculus-39-116.rds")
  Sys.setFileTime(old, Sys.time() - 100 * 86400)

  expect_message(gene_cache(clear = 90), "removed 1 of 2")
  expect_equal(gene_cache()$file, basename(new))

  expect_message(gene_cache(clear = TRUE), "removed 1 of 1")
  expect_equal(nrow(gene_cache()), 0L)
})

test_that("clear rejects anything that is not TRUE, FALSE or days", {
  local_cache()
  expect_error(gene_cache(clear = "everything"), "TRUE, FALSE or a number")
})

test_that("clearing drops the in-session copy too", {
  d <- local_cache()
  f <- file.path(d, "mapping-ensembl-homo_sapiens-38-116.rds")
  x <- data.frame(
    id = "ENSG1", name = "A", chr = "1", biotype = "protein_coding",
    span = 10L
  )
  genevintage:::.cache_put(f, genevintage:::.mapping_memo, x)
  # .cache_put creates the directory; nothing else does
  expect_true(dir.exists(d))
  expect_equal(nrow(gene_cache()), 1L)

  suppressMessages(gene_cache(clear = TRUE))
  # a stale memo would keep serving the table that was just deleted
  expect_null(genevintage:::.cache_get(f, genevintage:::.mapping_memo, names(x)))
})

test_that("an interrupted download is visible and removable", {
  d <- seed_cache("mapping-ensembl-homo_sapiens-38-116.rds")
  # .cache_put writes <f>.tmp then renames; a crash between the two leaves this
  file.create(file.path(d, "mapping-ensembl-homo_sapiens-38-115.rds.tmp"))

  x <- gene_cache()
  # matching only "\\.rds$" hid the orphan here and left clear unable to touch it
  expect_equal(nrow(x), 2L)
  expect_true("partial" %in% x$kind)
  expect_equal(x$species[x$kind == "partial"], "homo_sapiens")
  expect_equal(x$release[x$kind == "partial"], "115")

  suppressMessages(gene_cache(clear = TRUE))
  expect_equal(length(list.files(d, all.files = TRUE, no.. = TRUE)), 0L)
})


test_that("a full clear leaves no directory behind", {
  d <- local_cache()
  genevintage:::.cache_put(
    file.path(d, "mapping-ensembl-homo_sapiens-38-116.rds"),
    genevintage:::.mapping_memo, data.frame(id = "x")
  )
  expect_true(dir.exists(d))
  suppressMessages(gene_cache(clear = TRUE))
  expect_false(dir.exists(d))
})

test_that("a file the package did not write is not described as a species", {
  d <- seed_cache(c("mapping-ensembl-homo_sapiens-38-116.rds", "notes.rds"))
  # sub() no-ops on a non-match, which used to report "notes" as the species
  x <- gene_cache()
  expect_equal(nrow(x), 1L)
  expect_equal(x$species, "homo_sapiens")
  suppressMessages(gene_cache(clear = TRUE))
  expect_true(file.exists(file.path(d, "notes.rds"))) # and is not deleted
})
