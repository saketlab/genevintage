# Gene sets defined by position: mitochondrial and sex-linked.

test_that("the sex chromosomes of a clade are read off the annotation", {
  # mammals XY, birds ZW, C. elegans X with no Y
  expect_setequal(genevintage:::.sex_seqnames(c("1", "2", "X", "Y", "MT")), c("X", "Y"))
  expect_setequal(genevintage:::.sex_seqnames(c("1", "2", "Z", "W", "MT")), c("Z", "W"))
  expect_setequal(genevintage:::.sex_seqnames(c("I", "II", "III", "IV", "V", "X")), "X")
  expect_setequal(genevintage:::.sex_seqnames(c("chr1", "chrX", "chrY")), c("chrX", "chrY"))
})

test_that("a Roman-numbered chromosome X is not mistaken for a sex chromosome", {
  # yeast numbers I..XVI, so its "X" is chromosome 10 and it has no sex
  # chromosome at all. Roman numbering always passes IX and XI on the way.
  yeast <- c(
    "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
    "XI", "XII", "XIII", "XIV", "XV", "XVI", "Mito"
  )
  expect_length(genevintage:::.sex_seqnames(yeast), 0L)
  # C. elegans stops at V, so its X is real
  expect_equal(genevintage:::.sex_seqnames(c("I", "II", "III", "IV", "V", "X")), "X")
})

test_that("every spelling of the mitochondrion is recognised", {
  rx <- genevintage:::MITO_SEQNAMES
  # MT vertebrates, Mito yeast, MtDNA worm, Mt Arabidopsis, the long fly one,
  # and chrM for references written the UCSC way
  for (v in c("MT", "Mito", "MtDNA", "Mt", "mitochondrion_genome", "chrM", "M")) {
    expect_true(grepl(rx, v, ignore.case = TRUE), label = v)
  }
  # and nothing else is
  for (v in c("1", "X", "Y", "Z", "W", "MTND1", "KI270728.1", "Pt", "IV")) {
    expect_false(grepl(rx, v, ignore.case = TRUE), label = v)
  }
})

test_that("a species must be named when it cannot be detected", {
  skip_no_index()
  expect_error(mito_genes(), "name a `species`")
  expect_error(mito_genes(c("MT-CO1", "MT-ND1")), "do not say which species")
})

test_that("human mitochondrial genes are the canonical 37", {
  skip_net()
  mt <- mito_genes(species = "human", quiet = TRUE)
  # 13 protein-coding + 22 tRNA + 2 rRNA is the vertebrate mitochondrial genome
  expect_equal(nrow(mt), 37L)
  expect_true(all(mt$chr == "MT"))
  expect_true(all(c("MT-CO1", "MT-ND1", "MT-RNR1") %in% mt$name))
  expect_equal(sum(mt$biotype == "protein_coding"), 13L)
})

test_that("the mitochondrion is found whatever the annotation calls it", {
  skip_net()
  for (case in list(
    c("chicken", "MT"), c("fly", "mitochondrion_genome"),
    c("yeast", "Mito"), c("worm", "MtDNA")
  )) {
    mt <- mito_genes(species = case[1], quiet = TRUE)
    expect_gt(nrow(mt), 20)
    expect_true(all(mt$chr == case[2]), label = case[1])
  }
})

test_that("sex chromosomes follow the clade, not a constant", {
  skip_net()
  expect_setequal(unique(sex_genes(species = "human", quiet = TRUE)$chr), c("X", "Y"))
  expect_setequal(unique(sex_genes(species = "chicken", quiet = TRUE)$chr), c("Z", "W"))
  expect_equal(unique(sex_genes(species = "worm", quiet = TRUE)$chr), "X") # X but no Y
  # zebrafish sex is polygenic; yeast mating type is a locus, not a chromosome
  expect_equal(nrow(sex_genes(species = "zebrafish", quiet = TRUE)), 0L)
  expect_equal(nrow(sex_genes(species = "yeast", quiet = TRUE)), 0L)
})

test_that("asking for a chromosome a species lacks says which it has", {
  skip_net()
  expect_error(
    sex_genes(species = "chicken", which = "Y", quiet = TRUE),
    "has no chromosome Y; it has Z, W"
  )
  expect_gt(nrow(sex_genes(species = "human", which = "Y", quiet = TRUE)), 100)
})

test_that("a matrix's own row names come back, versions and all", {
  skip_net()
  ids <- c("ENSG00000210049.1", "ENSG00000198804.2", "ENSG00000141510.16")
  mt <- mito_genes(ids, species = "human", quiet = TRUE)
  expect_setequal(mt$input, ids[1:2]) # TP53 is nuclear, and excluded
  expect_true(all(mt$input %in% ids)) # so `ids %in% mt$input` subsets directly
  # the version is echoed back, while id stays the bare stable identifier
  expect_true(all(grepl("\\.[0-9]+$", mt$input)))
  expect_false(any(grepl("\\.[0-9]+$", mt$id)))
})

test_that("gene names work as input, resolving to their genes", {
  skip_net()
  mt <- mito_genes(c("MT-CO1", "TP53", "MT-ND1"), species = "human", quiet = TRUE)
  expect_setequal(mt$name, c("MT-CO1", "MT-ND1"))
})

test_that("biotypes group into the classes people ask for", {
  cls <- genevintage:::.biotype_class
  expect_equal(cls("protein_coding"), "coding")
  expect_equal(
    cls(c("lncRNA", "miRNA", "Mt_tRNA", "rRNA")),
    rep("noncoding", 4)
  )
  expect_equal(
    cls(c("processed_pseudogene", "unitary_pseudogene")),
    rep("pseudogene", 2)
  )
  expect_equal(cls(c("IG_V_gene", "TR_J_gene")), rep("immune", 2))
  # tested before the IG prefix, so an IG pseudogene is a pseudogene
  expect_equal(cls("IG_V_pseudogene"), "pseudogene")
  expect_equal(cls(c("TEC", NA, "")), rep("other", 3))
})

test_that("non-coding does not quietly mean not-protein-coding", {
  skip_net()
  # for human the difference is 15,205 pseudogenes plus TEC entries
  nc <- noncoding_genes(species = "human", quiet = TRUE)
  expect_false(any(grepl("pseudogene", nc$biotype)))
  expect_false(any(nc$biotype == "TEC"))
  ps <- genes_by_biotype("pseudogene", species = "human", quiet = TRUE)
  expect_gt(nrow(ps), 10000)
  expect_equal(length(intersect(nc$id, ps$id)), 0L)
})

test_that("protein-coding counts match the published figures", {
  skip_net()
  # these are the numbers quoted for each genome, so a drift here is a real bug
  for (case in list(c("human", 20131), c("fly", 13986), c("yeast", 6600))) {
    n <- nrow(coding_genes(species = case[1], quiet = TRUE))
    expect_equal(n, as.integer(case[2]), label = case[1])
  }
})

test_that("a biotype can be asked for exactly or as a class", {
  skip_net()
  expect_true(all(genes_by_biotype("lncRNA",
    species = "human",
    quiet = TRUE
  )$biotype == "lncRNA"))
  two <- genes_by_biotype(c("miRNA", "snoRNA"), species = "human", quiet = TRUE)
  expect_setequal(unique(two$biotype), c("miRNA", "snoRNA"))
  expect_error(
    genes_by_biotype("not_a_biotype", species = "human", quiet = TRUE),
    "no biotype not_a_biotype"
  )
  expect_error(genes_by_biotype(character()), "at least one biotype")
})

test_that("gene_biotypes reports the vocabulary with its classes", {
  skip_net()
  tb <- gene_biotypes(species = "human", quiet = TRUE)
  expect_named(tb, c("biotype", "class", "n"))
  expect_gt(nrow(tb), 20)
  expect_equal(tb$biotype[1], "lncRNA") # commonest first
  expect_true(all(tb$class %in% genevintage:::BIOTYPE_CLASSES))
  expect_equal(sum(tb$n), 78941L) # every gene accounted for
})

test_that("ribosomal proteins are separated from the pseudogenes that share their name", {
  skip_net()
  # ^RP[SL] alone hits 1775 human genes; only 102 are protein-coding, the rest
  # are ribosomal protein pseudogenes and lncRNAs carrying the name
  m <- fetch_mapping("homo_sapiens", 116)
  naive <- sum(grepl("^RP[SL]", ifelse(is.na(m$name), "", m$name)))
  rb <- ribosomal_genes(species = "human", quiet = TRUE)
  expect_gt(naive, 1000)
  expect_equal(nrow(rb), 91L)
  expect_true(all(rb$biotype == "protein_coding"))
  expect_false(any(grepl("pseudogene", rb$biotype)))
})

test_that("RP-named genes that are not ribosomal proteins are excluded", {
  skip_net()
  rb <- ribosomal_genes(species = "human", quiet = TRUE)
  # the S6 kinases are signalling enzymes and RPS19BP1 binds a ribosomal
  # protein rather than being one -- all protein-coding, so the biotype filter
  # alone leaves them in
  expect_false(any(grepl("^RPS6K", rb$name)))
  expect_false("RPS19BP1" %in% rb$name)
  # the paralogues are real ribosomal proteins and must survive
  expect_true(all(c("RPL10L", "RPL3L", "RPS27L", "RPSA", "RPLP0", "RPS4Y1")
  %in% rb$name))
})

test_that("the three ribosomal sets are distinct and add up", {
  skip_net()
  p <- ribosomal_genes(species = "human", which = "protein", quiet = TRUE)
  r <- ribosomal_genes(species = "human", which = "rrna", quiet = TRUE)
  m <- ribosomal_genes(species = "human", which = "mito", quiet = TRUE)
  a <- ribosomal_genes(species = "human", which = "all", quiet = TRUE)
  expect_equal(nrow(a), nrow(p) + nrow(r) + nrow(m))
  expect_length(intersect(p$id, m$id), 0L) # RP and MRP do not overlap
  expect_true(all(r$biotype %in% c("rRNA", "Mt_rRNA")))
})

test_that("mitoribosomal proteins are nuclear, not mitochondrial", {
  skip_net()
  # a real confusion: MRPL/MRPS encode the mitochondrion's ribosome but are
  # transcribed from the nucleus, so mito_genes() must not return them
  mr <- ribosomal_genes(species = "human", which = "mito", quiet = TRUE)
  expect_gt(nrow(mr), 50)
  expect_false(any(mr$chr == "MT"))
  expect_length(intersect(mr$id, mito_genes(species = "human", quiet = TRUE)$id), 0L)
})

test_that("ribosomal proteins are found whatever the species capitalises", {
  skip_net()
  # RPL human/yeast, Rpl mouse, RpL fly, rpl- worm
  for (case in list(
    c("mouse", "Rpl"), c("fly", "RpL"), c("worm", "rpl-"),
    c("yeast", "RPL")
  )) {
    rb <- ribosomal_genes(species = case[1], quiet = TRUE)
    expect_gt(nrow(rb), 50)
    expect_true(any(startsWith(rb$name, case[2])), label = case[1])
  }
})


test_that("tRNA genes are found under both names Ensembl gives them", {
  skip_net()
  # vertebrate GTFs carry only the mitochondrial set: 22 tRNAs is the complete
  # mitochondrial genome complement, and nuclear tRNAs are predicted separately
  h <- trna_genes(species = "human", quiet = TRUE)
  expect_equal(nrow(h), 22L)
  expect_true(all(h$biotype == "Mt_tRNA"))
  expect_true(all(h$chr == "MT"))
  # yeast, worm and Arabidopsis annotate the nuclear ones
  for (case in list(c("yeast", 299), c("worm", 634), c("arabidopsis", 689))) {
    t <- trna_genes(species = case[1], quiet = TRUE)
    expect_equal(nrow(t), as.integer(case[2]), label = case[1])
    expect_true(all(t$biotype == "tRNA"), label = case[1])
  }
})

test_that("tRNA and rRNA are different sets", {
  skip_net()
  t <- trna_genes(species = "human", quiet = TRUE)
  r <- ribosomal_genes(species = "human", which = "rrna", quiet = TRUE)
  expect_length(intersect(t$id, r$id), 0L)
  # both are mitochondrial here, but they are not the same genes
  expect_true(all(c("MT-TF", "MT-TV") %in% t$name))
  expect_true(all(c("MT-RNR1", "MT-RNR2") %in% r$name))
})

test_that("rRNA genes are the RNA of the ribosome, not its proteins", {
  skip_net()
  r <- rrna_genes(species = "human", quiet = TRUE)
  expect_equal(nrow(r), 55L) # 53 nuclear + 2 mitochondrial
  expect_setequal(unique(r$biotype), c("rRNA", "Mt_rRNA"))
  # the 497 rRNA pseudogenes are a separate question and are not included
  expect_false(any(grepl("pseudogene", r$biotype)))
  expect_gt(nrow(genes_by_biotype("rRNA_pseudogene", species = "human", quiet = TRUE)), 400)
  # and they are disjoint from the ribosomal proteins
  expect_length(intersect(r$id, ribosomal_genes(species = "human", quiet = TRUE)$id), 0L)
})

test_that("lncRNA is found under every spelling the annotations use", {
  skip_net()
  # modern Ensembl says lncRNA
  expect_gt(nrow(lncrna_genes(species = "human", quiet = TRUE)), 30000)
  # zebrafish at the same release still uses the older split vocabulary
  z <- lncrna_genes(species = "zebrafish", quiet = TRUE)
  expect_gt(nrow(z), 1000)
  expect_true(any(z$biotype %in% c("antisense", "lincRNA")))
  # chimpanzee says lincRNA
  expect_gt(nrow(lncrna_genes(species = "pan_troglodytes", quiet = TRUE)), 1000)
})

test_that("an annotation with no lncRNA class says so rather than returning nothing", {
  skip_net()
  # Drosophila files everything non-coding under a generic ncRNA biotype
  expect_message(g <- lncrna_genes(species = "fly"), "generic 'ncRNA' biotype")
  expect_equal(nrow(g), 0L)
  expect_gt(nrow(genes_by_biotype("ncRNA", species = "fly", quiet = TRUE)), 1000)
})

test_that("the RNA classes do not overlap one another", {
  skip_net()
  sets <- list(
    rrna = rrna_genes(species = "human", quiet = TRUE),
    trna = trna_genes(species = "human", quiet = TRUE),
    lnc = lncrna_genes(species = "human", quiet = TRUE),
    rp = ribosomal_genes(species = "human", quiet = TRUE)
  )
  for (i in seq_along(sets)) {
    for (j in seq_along(sets)) {
      if (i < j) {
        expect_equal(length(intersect(sets[[i]]$id, sets[[j]]$id)), 0L,
          label = paste(names(sets)[i], names(sets)[j])
        )
      }
    }
  }
  # all three RNA classes are non-coding, none is protein-coding
  for (nm in c("rrna", "trna", "lnc")) {
    expect_false(any(sets[[nm]]$biotype == "protein_coding"), label = nm)
  }
})

test_that("every gene-set function works on GENCODE as well as Ensembl", {
  skip_net()
  # GENCODE names chromosomes the UCSC way (chr1, chrM, chrX) and writes the
  # biotype as gene_type rather than gene_biotype. Both were silent failures:
  # the biotype-driven functions returned zero rows on GENCODE annotations.
  fns <- c(
    "mito_genes", "sex_genes", "coding_genes", "noncoding_genes",
    "ribosomal_genes", "trna_genes", "rrna_genes", "lncrna_genes"
  )
  for (f in fns) {
    g <- do.call(f, list(
      species = "human", release = "44", source = "gencode",
      quiet = TRUE
    ))
    expect_gt(nrow(g), 0, label = paste(f, "on GENCODE"))
    expect_true(all(c("id", "name", "chr", "biotype") %in% names(g)), label = f)
  }
  # the counts agree with the Ensembl release GENCODE v44 corresponds to
  for (f in c("mito_genes", "trna_genes", "ribosomal_genes")) {
    e <- nrow(do.call(f, list(
      species = "human", release = "112",
      assembly = "38", quiet = TRUE
    )))
    g <- nrow(do.call(f, list(
      species = "human", release = "44",
      source = "gencode", quiet = TRUE
    )))
    expect_equal(g, e, label = f)
  }
})

test_that("GENCODE's UCSC-style chromosome names are understood", {
  skip_net()
  m <- fetch_mapping("homo_sapiens", "44", source = "gencode")
  expect_true(all(startsWith(m$chr[seq_len(100)], "chr")))
  expect_true("chrM" %in% m$chr)
  expect_setequal(genevintage:::.sex_seqnames(m$chr), c("chrX", "chrY"))
  expect_false(any(is.na(m$biotype))) # gene_type, not gene_biotype
})

test_that("both GENCODE flavours and both species are indexed", {
  skip_no_index()
  fp <- genevintage:::.fp()
  g <- fp[fp$source %in% c("gencode", "gencode_all"), ]
  expect_setequal(unique(g$species), c("homo_sapiens", "mus_musculus"))
  for (src in c("gencode", "gencode_all")) {
    # Test coverage of the original range, allowing newly archived releases.
    expect_true(all(as.character(20:49) %in%
      g$release[g$source == src & g$species == "homo_sapiens"]))
    expect_true(all(paste0("M", 10:37) %in%
      g$release[g$source == src & g$species == "mus_musculus"]))
  }
  # human is numbered plainly, mouse with an M prefix
  expect_true(all(grepl("^M", g$release[g$species == "mus_musculus"])))
  expect_false(any(grepl("^M", g$release[g$species == "homo_sapiens"])))
})

test_that("a caller's own token comes back, whether id or gene name", {
  m <- data.frame(
    id = c("ENSG00000141510", "ENSG00000198804"),
    name = c("TP53", "MT-CO1"), chr = c("17", "MT"),
    biotype = "protein_coding", stringsAsFactors = FALSE
  )
  # a gene name used to come back as the Ensembl id, so subsetting a matrix
  # with symbol row names silently matched nothing
  r <- genevintage:::.restrict(m, "MT-CO1")
  expect_equal(r$input, "MT-CO1")
  expect_equal(r$id, "ENSG00000198804") # the stable id is still there
  expect_equal(nrow(genevintage:::.restrict(m, "MT-CO")), 0L) # no partial match
})

test_that("a vector mixing identifiers and gene names keeps both", {
  m <- data.frame(
    id = c("ENSG00000141510", "ENSG00000198804"),
    name = c("TP53", "MT-CO1"), chr = c("17", "MT"),
    biotype = "protein_coding", stringsAsFactors = FALSE
  )
  # deciding for the vector as a whole dropped whichever kind lost the vote
  o <- genevintage:::.restrict(m, c("ENSG00000141510", "MT-CO1"))
  expect_equal(nrow(o), 2L)
  expect_setequal(o$input, c("ENSG00000141510", "MT-CO1"))
  expect_setequal(o$name, c("TP53", "MT-CO1"))
})

test_that("the input column appears only when there was an input", {
  skip_net()
  # input is present exactly when the caller supplied something
  expect_false("input" %in% names(mito_genes(species = "human", quiet = TRUE)))
  got <- names(mito_genes(c("MT-CO1"), species = "human", quiet = TRUE))
  expect_equal(got[1], "input")
  expect_true(all(c("id", "name", "chr", "biotype") %in% got))
})

test_that("an unknown biotype is caught even when a class is named too", {
  skip_net()
  # the class match used to satisfy the check and the typo was dropped silently
  expect_error(
    genes_by_biotype(c("coding", "protien_coding"),
      species = "human",
      quiet = TRUE
    ),
    "no biotype protien_coding"
  )
})
