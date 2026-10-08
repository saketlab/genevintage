# Shared fixture access. testthat sources helper*.R before any test file.
fx_dir <- test_path("fixtures")
truth_dir <- test_path("truth")
manifest <- if (file.exists(file.path(fx_dir, "manifest.csv"))) {
  utils::read.csv(file.path(fx_dir, "manifest.csv"))
} else {
  NULL
}

# The corpus is walked by several tests; inflating 9.5MB once is enough. The
# golden and pipeline directories are read the same way, so one reader serves
# all three rather than each test file rolling its own.
read_fx <- local({
  cache <- new.env(parent = emptyenv())
  function(f, dir = fx_dir, read = function(k) readLines(gzfile(k), warn = FALSE)) {
    key <- file.path(dir, f)
    if (is.null(cache[[key]])) cache[[key]] <- read(key)
    cache[[key]]
  }
})

read_truth <- function(f) {
  utils::read.delim(gzfile(file.path(truth_dir, f)), colClasses = "character")
}

skip_no_index <- function() skip_if_not(nrow(.fp()) > 0, "no fingerprint index")

# Every network test opens with the same three guards.
skip_net <- function() {
  skip_on_cran()
  skip_if_offline()
  skip_no_index()
}

# A minimal annotation holding the given gene names.
nom <- function(...) {
  name <- c(...)
  data.frame(id = paste0("g", seq_along(name)), name = name, chr = "1", biotype = "protein_coding")
}

# A two-column numeric matrix with the given row names, for rename/correct tests.
mk <- function(ids) {
  matrix(seq_len(2 * length(ids)),
    nrow = length(ids),
    dimnames = list(ids, c("s1", "s2"))
  )
}

# A Compara homology file, for the offline stream tests. Rows are character
# vectors of fields; returns the file:// URL and cleans up with the test.
HOM_HDR <- c(
  "gene_stable_id", "homology_type", "homology_gene_stable_id",
  "homology_species", "identity", "homology_identity", "is_high_confidence"
)

hom_file <- function(..., header = HOM_HDR, env = parent.frame()) {
  f <- tempfile(fileext = ".tsv.gz")
  withr::defer(unlink(f), envir = env)
  rows <- vapply(list(...), paste, character(1), collapse = "\t")
  writeLines(c(paste(header, collapse = "\t"), rows), gzfile(f))
  paste0("file://", f)
}
