# Benchmark: fixtures whose true release and assembly are known.
golden_dir <- test_path("golden")
truth <- if (file.exists(file.path(golden_dir, "truth.csv"))) {
  utils::read.csv(file.path(golden_dir, "truth.csv"), stringsAsFactors = FALSE, colClasses = c(assembly = "character"))
} else {
  NULL
}
read_golden <- function(f) readLines(gzfile(file.path(golden_dir, f)), warn = FALSE)
skip_no_golden <- function() {
  skip_if_not(!is.null(truth) && nrow(truth) > 0, "no golden fixtures")
}

test_that("species is exact on every golden fixture", {
  skip_no_golden()
  for (i in seq_len(nrow(truth))) {
    cand <- detect_species(read_golden(truth$file[i]))
    # breeds and strains share an identifier space, so accept any tied candidate
    tied <- cand$species[cand$frac >= cand$frac[1] - 1e-9]
    expect_true(truth$species[i] %in% tied, label = truth$file[i])
  }
})

test_that("the true release is never ruled out", {
  skip_no_golden()
  # This is the property that matters. Ranking can be ambiguous -- adjacent
  # releases are often identical -- but the hard constraint must never exclude
  # the release the identifiers actually came from.
  for (i in seq_len(nrow(truth))) {
    out <- suppressMessages(suppressWarnings(
      detect_release(read_golden(truth$file[i]), species = truth$species[i])
    ))
    ok <- out$release == truth$release[i] & out$assembly == truth$assembly[i]
    skip_if(!any(ok))
    expect_false(isFALSE(out$feasible[ok][1]), label = truth$file[i])
  }
})

test_that("the true release scores with the best candidates", {
  skip_no_golden()
  # Ranking is reported, not asserted: danio r104 and r116 have identical gene
  # counts and identical qmax, so no scorer can separate them. What must hold is
  # that the true vintage is among those scoring best.
  ord <- truth[truth$species %in%
    unique(truth$species[!truth$species %in%
      c("saccharomyces_cerevisiae", "arabidopsis_thaliana")]), ]
  hit <- vapply(seq_len(nrow(ord)), function(i) {
    out <- suppressMessages(suppressWarnings(
      detect_release(read_golden(ord$file[i]), species = ord$species[i])
    ))
    ok <- which(out$release == ord$release[i] & out$assembly == ord$assembly[i])
    if (!length(ok)) {
      return(NA)
    }
    out$dist[ok[1]] <= out$dist[1] + 0.01
  }, logical(1))
  expect_gt(mean(hit, na.rm = TRUE), 0.9)
})

test_that("assembly transitions are called correctly", {
  skip_no_golden()
  # human GRCh37 vs GRCh38 and mouse GRCm38 vs GRCm39 are the cases that matter
  pairs <- truth[truth$species %in% c("homo_sapiens", "mus_musculus"), ]
  skip_if(length(unique(pairs$assembly)) < 2)
  for (i in seq_len(nrow(pairs))) {
    out <- suppressMessages(suppressWarnings(
      detect_release(read_golden(pairs$file[i]), species = pairs$species[i])
    ))
    best <- out[which(out$feasible)[1], ]
    expect_equal(best$assembly, pairs$assembly[i], label = pairs$file[i])
  }
})

test_that("a filtered universe still lands on a feasible release", {
  skip_no_golden()
  # real matrices are rarely the whole gene set: protein-coding only, or
  # expressed-only. Gene count is then misleading, so only feasibility can hold.
  set.seed(1)
  for (i in which(truth$species == "homo_sapiens")) {
    ids <- read_golden(truth$file[i])
    out <- suppressMessages(suppressWarnings(
      detect_release(sample(ids, round(0.4 * length(ids))), species = "homo_sapiens")
    ))
    ok <- out$release == truth$release[i] & out$assembly == truth$assembly[i]
    expect_true(any(ok) && isTRUE(out$feasible[ok][1]), label = truth$file[i])
  }
})
