# The naming rules, exercised against the real vocabulary of many organisms.
#
# mito_genes(), sex_genes() and trna_genes() all turn on strings that differ by
# clade. These fixtures hold each species' distinct seqnames and biotypes, taken
# from the actual annotation, so the rules can be checked over the whole range
# without the network.

vocab_dir <- test_path("vocab")
vocab_files <- if (dir.exists(vocab_dir)) list.files(vocab_dir, full.names = TRUE) else character()
vocab <- lapply(vocab_files, readRDS)
names(vocab) <- vapply(vocab, function(v) paste(v$source, v$species, v$release), "")

skip_no_vocab <- function() skip_if(length(vocab) == 0, "no vocabulary fixtures")

# Which species have a mitochondrial sequence or a sex chromosome is not
# guessed here -- Ensembl ships no MT for the dog, cat, duck, zebra finch,
# lamprey or fugu gene sets, and no sex chromosome for zebrafish or yeast. The
# fixtures are the ground truth, so these tests assert the rule holds wherever
# the feature exists, and that it is found often enough to be doing real work.
mito_of <- function(v) {
  v$chr$chr[grepl(genevintage:::MITO_SEQNAMES, v$chr$chr, ignore.case = TRUE)]
}

test_that("the vocabulary fixtures cover more than one clade", {
  skip_no_vocab()
  sp <- unique(vapply(vocab, function(v) v$species, ""))
  expect_gt(length(sp), 15)
  expect_true(any(vapply(vocab, function(v) v$source != "ensembl", logical(1))))
})

test_that("at most one mitochondrial sequence is identified per annotation", {
  skip_no_vocab()
  # the danger is a rule loose enough to also match a nuclear chromosome
  for (v in vocab) {
    hit <- v$chr$chr[grepl(genevintage:::MITO_SEQNAMES, v$chr$chr, ignore.case = TRUE)]
    expect_lte(length(hit), 1L, label = paste(
      v$species, "mito seqnames:",
      paste(hit, collapse = ",")
    ))
  }
})

test_that("a found mitochondrial sequence always looks like one", {
  skip_no_vocab()
  found <- 0L
  for (v in vocab) {
    hit <- mito_of(v)
    if (!length(hit)) next
    found <- found + 1L
    # a mitochondrial genome is small: tens of genes, never thousands. A rule
    # that matched a nuclear chromosome would show up here as a huge count.
    n <- v$chr$Freq[v$chr$chr == hit]
    expect_gt(n, 10)
    expect_lt(n, 500)
  }
  # and it is found in most of them, so the rule is not quietly matching nothing
  expect_gt(found, length(vocab) * 0.5)
})

test_that("sex chromosomes are only ever X, Y, Z or W", {
  skip_no_vocab()
  for (v in vocab) {
    got <- sub("^chr", "", genevintage:::.sex_seqnames(
      rep(v$chr$chr, pmin(v$chr$Freq, 1))
    ))
    expect_true(all(got %in% c("X", "Y", "Z", "W")),
      label = paste(v$species, paste(got, collapse = ","))
    )
  }
  # zebrafish sex is polygenic and yeast mating type is a locus: both correctly
  # yield nothing, so an empty answer is a real one
  for (sp in c("danio_rerio", "saccharomyces_cerevisiae")) {
    v <- Filter(function(x) x$species == sp, vocab)
    for (x in v) expect_length(genevintage:::.sex_seqnames(x$chr$chr), 0L)
  }
  found <- sum(vapply(vocab, function(x) {
    length(genevintage:::.sex_seqnames(x$chr$chr)) > 0
  }, logical(1)))
  expect_gt(found, length(vocab) * 0.5)
})

test_that("birds get Z and W, mammals X and Y, and never both systems", {
  skip_no_vocab()
  sexof <- function(sp, src = "ensembl") {
    v <- vocab[[paste(src, sp, if (src == "ensembl") "116" else "44")]]
    if (is.null(v)) {
      return(NULL)
    }
    sub("^chr", "", genevintage:::.sex_seqnames(v$chr$chr))
  }
  # the assertion is the *system*, not its completeness: a W is small and
  # heterochromatic and is simply absent from the duck and zebra finch
  # assemblies, as a Y is from the horse, cat and opossum ones
  for (bird in c("gallus_gallus", "taeniopygia_guttata", "anas_platyrhynchos")) {
    s <- sexof(bird)
    if (!is.null(s)) {
      expect_true("Z" %in% s, label = bird)
      expect_true(all(s %in% c("Z", "W")), label = bird)
    }
  }
  for (mam in c("homo_sapiens", "mus_musculus", "bos_taurus")) {
    s <- sexof(mam)
    if (!is.null(s)) expect_true(all(s %in% c("X", "Y")) && "X" %in% s, label = mam)
  }
  # no annotation should ever report a mammalian and an avian system together
  for (v in vocab) {
    got <- sub("^chr", "", genevintage:::.sex_seqnames(v$chr$chr))
    expect_false(any(c("X", "Y") %in% got) && any(c("Z", "W") %in% got),
      label = v$species
    )
  }
})

test_that("every biotype an organism uses lands in exactly one class", {
  skip_no_vocab()
  for (v in vocab) {
    cls <- genevintage:::.biotype_class(v$biotype$biotype)
    expect_true(all(cls %in% genevintage:::BIOTYPE_CLASSES),
      label = paste(v$species, paste(unique(cls), collapse = ","))
    )
    expect_equal(length(cls), nrow(v$biotype))
    # every organism has protein-coding genes, and they are the majority class
    # or close to it -- an annotation with none means the biotype was not read
    n_coding <- sum(v$biotype$Freq[cls == "coding"])
    expect_gt(n_coding, 1000)
  }
})

test_that("pseudogenes never leak into the non-coding class", {
  skip_no_vocab()
  for (v in vocab) {
    cls <- genevintage:::.biotype_class(v$biotype$biotype)
    nc <- v$biotype$biotype[cls == "noncoding"]
    expect_false(any(grepl("pseudogene", nc)), label = v$species)
    expect_false(any(grepl("^(IG|TR)_", nc)), label = v$species)
  }
})

test_that("tRNA biotypes are recognised in every organism that has them", {
  skip_no_vocab()
  seen <- 0L
  for (v in vocab) {
    t <- v$biotype$biotype[grepl(genevintage:::TRNA_BIOTYPES, v$biotype$biotype)]
    expect_true(all(t %in% c("tRNA", "Mt_tRNA")), label = v$species)
    # any biotype whose name mentions tRNA must be matched, or the rule is too
    # narrow for some clade
    mentions <- v$biotype$biotype[grepl("tRNA", v$biotype$biotype)]
    expect_setequal(t, mentions)
    if (length(t)) seen <- seen + 1L
  }
  expect_gt(seen, 10L)
})

test_that("the ribosomal name rules hold across every naming convention", {
  skip_no_vocab()
  for (v in vocab) {
    rp <- v$rp_names
    if (!length(rp)) next
    keep <- grepl(genevintage:::RIBO_PROTEIN_RX, rp, ignore.case = TRUE) &
      !grepl(genevintage:::RIBO_PROTEIN_NOT_RX, rp, ignore.case = TRUE)
    mito <- grepl(genevintage:::MITORIBO_PROTEIN_RX, rp, ignore.case = TRUE)
    # RP and MRP are disjoint in every organism: nothing is both
    expect_equal(length(intersect(rp[keep], rp[mito])), 0L, label = v$species)
    # and no S6 kinase survives, whatever the species capitalises it as
    expect_false(any(grepl("^RPS6K", rp[keep], ignore.case = TRUE)), label = v$species)
  }
})

test_that("the RNA-class rules never claim the same biotype twice", {
  skip_no_vocab()
  for (v in vocab) {
    b <- v$biotype$biotype
    rr <- grepl(genevintage:::RRNA_BIOTYPES, b)
    tr <- grepl(genevintage:::TRNA_BIOTYPES, b)
    ln <- grepl(genevintage:::LNCRNA_BIOTYPES, b)
    expect_true(all(rr + tr + ln <= 1), label = v$species)
    # every one of them is non-coding, in every organism
    cls <- genevintage:::.biotype_class(b)
    expect_true(all(cls[rr | tr | ln] == "noncoding"), label = v$species)
  }
})

test_that("every organism's rRNA and lncRNA vocabulary is recognised", {
  skip_no_vocab()
  seen_r <- 0L
  seen_l <- 0L
  for (v in vocab) {
    b <- v$biotype$biotype
    # anything named rRNA is either matched or is explicitly a pseudogene
    mentions_r <- b[grepl("rRNA", b)]
    unmatched <- mentions_r[!grepl(genevintage:::RRNA_BIOTYPES, mentions_r)]
    expect_true(all(grepl("pseudogene", unmatched)),
      label = paste(v$species, paste(unmatched, collapse = ","))
    )
    if (any(grepl(genevintage:::RRNA_BIOTYPES, b))) seen_r <- seen_r + 1L
    # lncRNA is matched wherever the annotation draws the class at all
    if (any(grepl("lncRNA|lincRNA|antisense", b))) {
      expect_true(any(grepl(genevintage:::LNCRNA_BIOTYPES, b)), label = v$species)
      seen_l <- seen_l + 1L
    }
  }
  expect_gt(seen_r, 10L)
  expect_gt(seen_l, 10L)
})
