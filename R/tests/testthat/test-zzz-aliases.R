test_that("every PascalCase alias is identical to its snake_case target", {
  pairs <- list(
    c("AggregateTranscripts", "aggregate_transcripts"),
    c("AnnotateGenes", "annotate_genes"),
    c("BuildSymbolIndex", "build_symbol_index"),
    c("BuildSymbolIndexFromGTF", "build_symbol_index_from_gtf"),
    c("CodingGenes", "coding_genes"),
    c("CompareReleases", "compare_releases"),
    c("CorrectGeneNames", "correct_genenames"),
    c("DetectRelease", "detect_release"),
    c("DetectSpecies", "detect_species"),
    c("DetectSpeciesNames", "detect_species_names"),
    c("EntrezIDs", "entrez_ids"),
    c("FetchMapping", "fetch_mapping"),
    c("FetchOrthologs", "fetch_orthologs"),
    c("FetchTranscripts", "fetch_transcripts"),
    c("FetchXrefs", "fetch_xrefs"),
    c("GeneBiotypes", "gene_biotypes"),
    c("GeneCache", "gene_cache"),
    c("GeneIDs", "gene_ids"),
    c("GeneNames", "gene_names"),
    c("GeneID2Name", "geneid2name"),
    c("GenesByBiotype", "genes_by_biotype"),
    c("GuessScheme", "guess_scheme"),
    c("HemoglobinGenes", "hemoglobin_genes"),
    c("lncRNAGenes", "lncrna_genes"),
    c("MHCGenes", "mhc_genes"),
    c("MitoGenes", "mito_genes"),
    c("NoncodingGenes", "noncoding_genes"),
    c("OrthologSpecies", "ortholog_species"),
    c("Orthologs", "orthologs"),
    c("ParseEns", "parse_ens"),
    c("QmaxViolations", "qmax_violations"),
    c("ResolveSpecies", "resolve_species"),
    c("RibosomalGenes", "ribosomal_genes"),
    c("rRNAGenes", "rrna_genes"),
    c("SexGenes", "sex_genes"),
    c("SpeciesTable", "species_table"),
    c("StreamGTF", "stream_gtf"),
    c("StreamHomologies", "stream_homologies"),
    c("StreamTranscripts", "stream_transcripts"),
    c("StreamXrefs", "stream_xrefs"),
    c("Transcript2Gene", "transcript2gene"),
    c("tRNAGenes", "trna_genes")
  )
  for (p in pairs) {
    expect_identical(get(p[1]), get(p[2]), label = p[1])
  }
})

test_that("every exported function has a PascalCase alias", {
  ns <- getNamespaceExports("genevintage")
  snake <- grep("_|^[a-z0-9]+$", ns, value = TRUE)
  pascal <- setdiff(ns, snake)
  # each snake_case export's PascalCase form differs only in case/underscores
  norm <- function(x) tolower(gsub("_", "", x))
  expect_setequal(norm(snake), norm(pascal))
})
