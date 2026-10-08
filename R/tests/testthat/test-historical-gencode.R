test_that("the 2014 Monocle single-cell IDs select their documented GENCODE 17", {
  skip_if_not(
    file.exists(test_path("pipelines", "gse52529_gencode_v17.txt.gz")),
    "no pipeline fixture"
  )
  ids <- read_fx("gse52529_gencode_v17.txt.gz", dir = test_path("pipelines"))
  expect_length(ids, 47192L)
  expect_true(all(!is.na(.id_version(ids))))
  out <- suppressWarnings(detect_release(ids))
  expect_identical(out$source[1], "gencode")
  expect_identical(out$release[1], "17")
  expect_identical(out$assembly[1], "37")
  expect_true(out$feasible[1])
})

test_that("historical GENCODE references are indexed on the correct assembly", {
  fp <- .fp()
  human <- fp[fp$species == "homo_sapiens" & fp$source == "gencode", ]
  old <- human[match(as.character(17:19), human$release), ]
  expect_false(anyNA(old$release))
  expect_true(all(old$assembly == "37"))
  expect_true(all(is.na(old$frozen_release)))
  expect_true(all(human$assembly[.rel_num(human$release) >= 20] == "38"))
})

test_that("fetch_mapping routes historical GENCODE to its own GTF", {
  seen <- character()
  withr::local_envvar(R_USER_CACHE_DIR = tempfile())
  local_mocked_bindings(
    stream_gtf = function(url, ...) {
      seen <<- c(seen, url)
      data.frame(
        id = "ENSG00000000003.10", name = "TSPAN6",
        chr = "chrX", biotype = "protein_coding", span = 1L
      )
    }
  )
  for (rel in 17:19) {
    expect_identical(fetch_mapping("human", rel,
      source = "gencode",
      assembly = "37"
    )$name, "TSPAN6")
  }
  expect_identical(seen, sprintf(paste0(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/",
    "release_%s/gencode.v%s.annotation.gtf.gz"
  ), 17:19, 17:19))
})
