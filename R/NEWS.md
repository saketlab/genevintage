# genevintage 0.1.3

- Added `transcript2gene()` to map Ensembl/GENCODE transcript identifiers to
  gene ids or names 
- Added `aggregate_transcripts()` to sum a transcript-level matrix 
  to gene level
- `detect_species()` and `detect_release()` accept transcript identifiers
- Add support for `gene_lengths()`/`stream_exon_lengths()` to fetch exon-union gene lengths
- Added `stream_homologies()` to read older Compara dumps that have no
  `is_high_confidence` column; `high_confidence` is then `NA`.
- `entrez_ids()` no longer fails on identifiers dated before release 85, which
  has no cross-reference files; it reads them from the newest release.

# genevintage 0.1.2

- `geneid2name()` and `correct_genenames()` (renamed from `correct_names()`)
  also accept a matrix, sparse matrix or data frame, replacing row names in
  place instead of returning a vector.
- Every exported function also has a PascalCase alias (`GeneID2Name()`,
  `MitoGenes()`, `CorrectGeneNames()`, and so on).
- Add functions to infer build from gene names

# genevintage 0.1.1

- `detect_species_names()` infers species and annotation build from gene
  symbols alone, no identifiers needed, using an index from
  `build_symbol_index()`.
- Index human GENCODE 17-19 on GRCh37.

# genevintage 0.1.0

First release.

