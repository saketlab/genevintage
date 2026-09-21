<img src="R/man/figures/logo.png" align="right" height="160" alt="genevintage logo" />

# genevintage

Given a vector of gene IDs, find out which species and Ensembl/GENCODE release
they came from, then convert them to gene names. 

| Client | Status | Path |
|---|---|---|
| R | available | [`R/`](R/) |

Docs: <https://genevintage.saketlab.org/R/>

## R

```r
# install.packages("remotes")
remotes::install_github("saketlab/genevintage", subdir = "R")
```

```r
library(genevintage)
ids <- c("ENSG00000141510.16", "ENSG00000284733", "ENSG00000012048")

detect_release(ids)          # ranked releases, Ensembl and GENCODE
#>        species  source release assembly feasible   dist n_index
#> 1 homo_sapiens gencode      26       38     TRUE 9.6486   58219

geneid2name(ids)             # detect + map in one step
#> [1] "TP53"            "ENSG00000284733" "BRCA1"
```

Also: gene classes (`mito_genes()`, `sex_genes()`, `coding_genes()`,
`ribosomal_genes()`, `lncrna_genes()`), `orthologs()`, `compare_releases()`,
`annotate_genes()`, `correct_genenames()` (undo Excel's `9-Sep` autocorrect). See the
[vignette](genevintage.saketlab.org) for details.

## License

MIT
