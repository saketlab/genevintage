# Shared setup for the correct_genenames() regression tests.

# padding; keeps a short list below the numeric-list limit
plain <- paste0("GENE", 1:60)

# Current human names of the genes that Excel damages.
hs_now <- function() {
  nom(
    paste0("SEPTIN", c(1:12, 14)), "SEPTIN7P2", paste0("MARCHF", 1:11),
    "MTARC1", "MTARC2", "DELEC1", "SELENOF", "SCN1A", "MARS2", "AGO1", "TP53"
  )
}

# Current mouse names of the same genes. MGI lists Sep2 as an old name of Apoa1.
mm_now <- function() {
  nom(
    paste0("Septin", c(1:12, 14)), paste0("Marchf", 1:11), "Mtarc1", "Mtarc2",
    "Selenof", "Apoa1", "Mars2", "Mei1", "Mei4", "Trp53"
  )
}

# repair x among the padding; the result keeps the padding
run <- function(x, mapping = hs_now(), ...) {
  correct_genenames(c(x, plain), mapping = mapping, quiet = TRUE, ...)
}

fixed <- function(x, mapping = hs_now(), ...) {
  as.vector(run(x, mapping, ...))[seq_along(x)]
}
