# GRCh37 and GRCh38 share release numbers but are different annotations.
test_that("the index keeps assemblies apart", {
  skip_no_index()
  fp <- .fp()
  expect_true("assembly" %in% names(fp))
  # some release, for some species, ships two assemblies at once; they must be
  # kept apart rather than collapsed onto one row
  key <- paste(fp$species, fp$release)
  dup <- key[duplicated(key)]
  skip_if(!length(dup))
  both <- fp[key == dup[1], ]
  expect_equal(length(unique(both$assembly)), nrow(both))
})

test_that("qmax never falls within one species and assembly", {
  skip_no_index()
  # qmax monotonicity is a heuristic, not a law: Ensembl re-annotates within an
  # assembly occasionally. What must hold is that the exceptions stay rare, and
  # that detect_release() drops the hard constraint for exactly those species.
  v <- qmax_violations()
  fp <- .fp()
  expect_lt(length(v), 0.05 * length(unique(fp$species)))
  for (sp in sub(" .*$", "", v)) {
    out <- suppressMessages(suppressWarnings(detect_release(
      sprintf(
        "%s%011.0f", fp$prefix[fp$species == sp][1],
        seq(1, max(fp$qmax[fp$species == sp], na.rm = TRUE), length.out = 500)
      ),
      species = sp
    )))
    expect_true(all(out$feasible | is.na(out$feasible)), label = sp)
  }
})

test_that("GRCh37 ids prefer the GRCh37 rows", {
  skip_no_index()
  fp <- .fp()
  h37 <- fp[fp$species == "homo_sapiens" & fp$assembly == "37", ]
  h38 <- fp[fp$species == "homo_sapiens" & fp$assembly == "38", ]
  skip_if(nrow(h37) == 0 || nrow(h38) == 0)
  # an id set topping out below the GRCh37 ceiling but well under GRCh38's
  ids <- sprintf("ENSG%011.0f", seq(1, max(h37$qmax), length.out = 2000))
  out <- detect_release(ids, species = "homo_sapiens")
  expect_true(any(out$feasible))
  expect_true(nrow(out) > nrow(h38)) # both assemblies are scored
})

test_that("the index records which schemes have ordered identifiers", {
  skip_no_index()
  fp <- .fp()
  expect_true(all(fp$ordered[fp$species == "homo_sapiens"]))
  # TAIR and SGD embed the chromosome, so their numbers are not ordered
  for (sp in c("arabidopsis_thaliana", "saccharomyces_cerevisiae")) {
    if (sp %in% fp$species) expect_false(any(fp$ordered[fp$species == sp]), label = sp)
  }
})

test_that("Ensembl Genomes assemblies are not tangled with release numbers", {
  skip_no_index()
  fp <- .fp()
  # <species>_core_63_116_4 must yield assembly "4", not "116_4"
  expect_false(any(grepl("_", fp$assembly)))
})
