# 0.1.2

- `geneid2name()` and `correct_genenames()` (renamed from `correct_names()`)
  also accept a matrix, sparse matrix or data frame, replacing row names in
  place instead of returning a vector.
- Every exported function also has a PascalCase alias (`GeneID2Name()`,
  `MitoGenes()`, `CorrectGeneNames()`, and so on).
- Add functions to infer build from gene names

# 0.1.1

- `detect_species_names()` infers species and annotation build from gene
  symbols alone, no identifiers needed, using an index from
  `build_symbol_index()`.
- Index human GENCODE 17-19 on GRCh37.

# 0.1.0

First release.

