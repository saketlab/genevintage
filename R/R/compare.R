# What changed between two annotation vintages.

#' Compare two releases of one annotation
#'
#' Joins two vintages on the stable identifier and reports what differs: genes
#' the newer release added, genes it dropped, and genes whose symbol or biotype
#' changed. Rows that are identical in both are not returned.
#'
#' Added and removed genes are identified by stable identifier presence.
#' Finding a successor to a retired identifier requires a separate lookup.
#'
#' @param species Species in any form [resolve_species()] accepts.
#' @param from,to The two releases, in either order.
#' @param assembly,source Passed to [fetch_mapping()]. Both releases are read
#'   from the same assembly and annotation family to isolate release changes.
#' @param quiet Suppress the summary message.
#' @return A data frame of `id`, `change` (`"added"`, `"removed"`, `"renamed"`
#'   or `"retyped"`), `name_from`, `name_to`, `biotype_from` and `biotype_to`.
#' @seealso [detect_release()] to find out which vintages you have.
#' @export
#' @examples
#' \dontrun{
#' d <- compare_releases("human", from = 110, to = 116)
#' table(d$change)
#' d[d$change == "renamed", ]
#' }
compare_releases <- function(species, from, to, assembly = NULL, source = NULL,
                             quiet = FALSE) {
  species <- resolve_species(species)
  a <- fetch_mapping(species, from, assembly = assembly, source = source)
  b <- fetch_mapping(species, to, assembly = assembly, source = source)
  a <- a[!duplicated(a$id), ]
  b <- b[!duplicated(b$id), ]

  ids <- union(a$id, b$id)
  ia <- match(ids, a$id)
  ib <- match(ids, b$id)
  out <- data.frame(
    id = ids,
    name_from = a$name[ia], name_to = b$name[ib],
    biotype_from = a$biotype[ia], biotype_to = b$biotype[ib],
    stringsAsFactors = FALSE
  )

  # paired missing values compare equal
  same <- function(x, y) (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)
  out$change <- ifelse(is.na(ia), "added",
    ifelse(is.na(ib), "removed",
      ifelse(!same(out$name_from, out$name_to), "renamed",
        ifelse(!same(out$biotype_from, out$biotype_to), "retyped", NA_character_)
      )
    )
  )
  out <- out[
    !is.na(out$change),
    c("id", "change", "name_from", "name_to", "biotype_from", "biotype_to")
  ]
  row.names(out) <- NULL

  if (!quiet) {
    n <- table(factor(out$change, c("added", "removed", "renamed", "retyped")))
    message(sprintf(
      "%s %s -> %s: %d added, %d removed, %d renamed, %d retyped",
      species, from, to, n[["added"]], n[["removed"]],
      n[["renamed"]], n[["retyped"]]
    ))
  }
  out
}
