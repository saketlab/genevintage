# genevintage <img src="man/figures/logo.png" align="right" height="139" alt="" />

[![R-CMD-check](https://github.com/saketlab/genevintage/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/saketlab/genevintage/actions/workflows/R-CMD-check.yaml)

Genevintage supports easy conversion between gene ids and gene names for any genome
and GTF build.

```r
ids <- c("ENSG00000141510.16", "ENSG00000284733", "ENSG00000012048")
```

Using a reverse index, `genevintage` can detect the actual release:

```r
head(DetectRelease(ids), 3)
#>        species  source release assembly feasible   dist n_index
#> 1 homo_sapiens gencode      26       38     TRUE 9.6486   58219
#> 2 homo_sapiens gencode      27       38     TRUE 9.6505   58288
#> 3 homo_sapiens gencode      28       38     TRUE 9.6534   58381
```

Gene names can be retrieved with a single command:

```r
GeneID2Name(ids)
#> Best fingerprint match: homo_sapiens, gencode release 26 (assembly 38) -- 79 candidates score within 10%
#> [1] "TP53"            "ENSG00000284733" "BRCA1"
```

## Also does

```r
DetectSpeciesNames(c("TP53", "BRCA1", "MYC"), idx)    # idx <- BuildSymbolIndex(...); detect from symbols
m <- GeneID2Name(m)                                   # rownames(m), id -> name
CorrectGeneNames(c("9-Sep", "1-Mar"), species = "human") # undo Excel's date/float autocorrect
GeneIDs(c("TP53", "BRCA1"), species = "human")        # name -> id, and back
AnnotateGenes(rownames(counts))                       # full record: symbol, biotype, coords
MitoGenes(rownames(counts)); SexGenes(species = "human"); RibosomalGenes(species = "human")
SexGenes(species = "chicken", which = "W")            # per-clade: X/Y, Z/W, or none (yeast)
GenesByBiotype("lncRNA", species = "mouse")           # any Ensembl/GENCODE biotype
HemoglobinGenes(species = "human"); MHCGenes(species = "mouse")
CompareReleases("human", from = 110, to = 116)        # added / removed / renamed genes
Orthologs(ids, to = "mouse")                          # common names resolve too
```

Every function above is also available under its snake_case name
(`GeneID2Name()` is also `geneid2name()`, `MitoGenes()` is also `mito_genes()`,
and so on) for callers who prefer that style.

## Install

```r
# install.packages("remotes")
remotes::install_github("saketlab/genevintage")
```
