# Protocol preference: FTP first, HTTPS to catch it.

test_that("a URL converts between the two protocols", {
  expect_equal(
    genevintage:::.as_ftp("https://ftp.ensembl.org/pub/x"),
    "ftp://ftp.ensembl.org/pub/x"
  )
  expect_equal(
    genevintage:::.as_ftp("http://ftp.ensembl.org/pub/x"),
    "ftp://ftp.ensembl.org/pub/x"
  )
  expect_equal(
    genevintage:::.as_https("ftp://ftp.ebi.ac.uk/pub/x"),
    "https://ftp.ebi.ac.uk/pub/x"
  )
  # already the target protocol: unchanged, not doubled
  expect_equal(genevintage:::.as_ftp("ftp://h/x"), "ftp://h/x")
  expect_equal(genevintage:::.as_https("https://h/x"), "https://h/x")
})

test_that("the protocol option pins the choice without probing", {
  u <- "https://ftp.ensembl.org/pub/x"
  withr::with_options(
    list(genevintage.protocol = "https"),
    expect_equal(genevintage:::.url_candidates(u), u)
  )
  withr::with_options(
    list(genevintage.protocol = "ftp"),
    expect_equal(
      genevintage:::.url_candidates(u),
      "ftp://ftp.ensembl.org/pub/x"
    )
  )
})

test_that("an HTML index yields bare entry names", {
  n <- genevintage:::.listing_names(c(
    '<a href="?C=N;O=D">Name</a>',
    '<a href="/pub/release-116/">Parent Directory</a>',
    '<td><a href="gene.txt.gz">gene.txt.gz</a></td>',
    '<td><a href="homo_sapiens/">homo_sapiens/</a></td>'
  ))
  # the sort links, the parent link and the trailing slash are all gone
  expect_setequal(n, c("gene.txt.gz", "homo_sapiens"))
})

test_that("an FTP ls -l listing yields the same names", {
  n <- genevintage:::.listing_names(c(
    "drwxr-xr-x 2 ftp ftp    4096 Jan 10  2024 homo_sapiens",
    "-rw-r--r-- 1 ftp ftp 2842270 May 15  2017 gene.txt.gz"
  ))
  expect_setequal(n, c("gene.txt.gz", "homo_sapiens"))
  # a name containing spaces survives: the split is on field count, not on space
  expect_equal(genevintage:::.listing_names(
    "-rw-r--r-- 1 ftp ftp 55 May 15  2017 my file.txt"
  ), "my file.txt")
})

test_that("a bare NLST listing is understood too", {
  expect_setequal(
    genevintage:::.listing_names(c("gene.txt.gz", "homo_sapiens/")),
    c("gene.txt.gz", "homo_sapiens")
  )
})

test_that("both listing shapes give a caller the same answer", {
  # this is the point of normalising: neither caller learns which protocol ran
  html <- genevintage:::.listing_names(c(
    '<a href="Compara.116.protein_default.homologies.tsv.gz">x</a>',
    '<a href="Compara.116.ncrna_default.homologies.tsv.gz">x</a>'
  ))
  ftp <- genevintage:::.listing_names(c(
    "-rw-r--r-- 1 ftp ftp 109478724 Jan 1 2025 Compara.116.protein_default.homologies.tsv.gz",
    "-rw-r--r-- 1 ftp ftp   1234567 Jan 1 2025 Compara.116.ncrna_default.homologies.tsv.gz"
  ))
  expect_setequal(html, ftp)
  rx <- "^Compara\\.[0-9]+\\.protein_[a-z]+\\.homologies\\.tsv\\.gz$"
  expect_length(grep(rx, html, value = TRUE), 1L)
  expect_equal(grep(rx, html, value = TRUE), grep(rx, ftp, value = TRUE))
})

test_that("only the last protocol tried is retried", {
  # a blocked protocol must fail fast; the one that works earns the retries
  tried <- character()
  f <- function(u) {
    tried <<- c(tried, u)
    stop("nope")
  }
  withr::with_options(list(genevintage.protocol = "auto"), {
    # force two candidates without probing the network
    cand <- c("ftp://h/x", "https://h/x")
    for (i in seq_along(cand)) {
      v <- genevintage:::.attempt(function() f(cand[i]),
        tries = if (i == length(cand)) 3L else 1L,
        wait = 0
      )
      expect_s3_class(v, "try-error")
    }
  })
  expect_equal(sum(tried == "ftp://h/x"), 1L)
  expect_equal(sum(tried == "https://h/x"), 3L)
})

test_that("a working call is not retried", {
  n <- 0L
  v <- genevintage:::.attempt(function() {
    n <<- n + 1L
    "ok"
  })
  expect_equal(v, "ok")
  expect_equal(n, 1L)
})

test_that("a malformed listing yields no entries rather than junk ones", {
  n <- genevintage:::.listing_names
  # markup with no links is an error page or an empty directory; returning the
  # raw HTML as a filename is worse than returning nothing
  expect_length(n("<html><body>No files found</body></html>"), 0L)
  # a truncated ls -l line is not a filename
  expect_length(n("-rw-r--r-- gene.txt.gz"), 0L)
  # a symlink entry is the name, not "name -> target"
  expect_equal(
    n("lrwxrwxrwx 1 ftp ftp 7 Jan 1  2025 current -> release-116"),
    "current"
  )
  expect_length(n(character()), 0L)
  expect_length(n(c("", "   ")), 0L)
})

test_that("a real listing still parses after the hardening", {
  n <- genevintage:::.listing_names
  expect_equal(n("-rw-r--r-- 1 ftp ftp 2842270 May 15  2017 gene.txt.gz"), "gene.txt.gz")
  expect_equal(n("-rw-r--r--. 1 ftp ftp 55 May 15  2017 selinux.txt"), "selinux.txt")
  expect_setequal(n(c("gene.txt.gz", "homo_sapiens/")), c("gene.txt.gz", "homo_sapiens"))
})
