# genevintage 0.1.4

- `correct_genenames()` has these changes:
  - no longer changes a list that is mostly numbers, such as Entrez IDs or array
    spot coordinates. `42248` does not become `SEPTIN1`.
  - no longer change a five-digit number if the list has a number that is not
    a gene or a date. In `c("42248", "7157", "TP53", "MYC", "GAPDH")`, `42248`
    stays unchanged.
  - does not read a date as damaged unless the time is midnight.
  - also converts month names from other Excel locales (`mars-02`, `ago-01`). When the
    month comes first, the day must have two digits.
  - can read dates in the form `2019/9/9`. A four-digit year must start with 19
    or 20.
  - It reads `2003-09-01` as `SEPT3`, not `SEPT1`, if the list also has a date
    such as `2005-09-01`. Some Excel files store `SEPT3` as September 2003, and
    `SEPT5` as September 2005. The year gives the gene number. If the year gives
    no gene for any of these dates, it reads the day.
  - It keeps the day reading if the same year and month has other days, such as
    `2012-09-01` with `2012-09-02`. A merge of two files can give this.
  - It gives two new attributes for each input. `category` is the form of the
    damage, such as `day_month` or `date_year`. `locale` is the languages that
    use the month word, such as `nl` for `mei-01`.
  - In the year reading, it also tries other languages. `2002-03-01` can be
    `MARS2`, and `2001-05-01` can be `MEI1`.
- An input with `species` and no `release` does not stop when gene names include a
  digit-only, `LOC` or RefSeq symbol.
- `orthologs()` reads a digit-only, `LOC` or RefSeq symbol as a gene name. It does
  not remove these symbols from the result without a warning.

# genevintage 0.1.3

- New function `transcript2gene()` maps Ensembl and GENCODE transcript
  identifiers to gene IDs or gene names.
- New function `aggregate_transcripts()` sums a transcript-level matrix to gene
  level.
- `detect_species()` and `detect_release()` accept transcript identifiers.
- New functions `gene_lengths()` and `stream_exon_lengths()` get exon-union gene
  lengths.
- New function `stream_homologies()` reads older Compara dumps that have no
  `is_high_confidence` column. For these dumps, `high_confidence` is `NA`.
- `entrez_ids()` does not fail on identifiers from before release 85. These
  releases have no cross-reference files.

# genevintage 0.1.2

- `geneid2name()` and `correct_genenames()` accept a matrix, sparse matrix or
  data frame. They replace the row names and return the same object.
- The new name of `correct_names()` is `correct_genenames()`.
- Each exported function has a PascalCase alias, such as `GeneID2Name()`,
  `MitoGenes()` and `CorrectGeneNames()`.
- New functions find the annotation build from gene names.

# genevintage 0.1.1

- `detect_species_names()` finds the species and the annotation build from gene
  symbols. It does not need identifiers. It uses an index from
  `build_symbol_index()`.
- The index includes human GENCODE 17-19 on GRCh37.

# genevintage 0.1.0

First release.
