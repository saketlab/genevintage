# What changed between two annotation vintages.

#' Compare two releases of one annotation
#'
#' This function joins two releases on the stable identifier and reports the
#' differences. It reports genes that are only in `to` (`"added"`), genes that
#' are only in `from` (`"removed"`), and genes with a changed name or biotype. It
#' does not return rows that are identical in the two releases.
#'
#' The function finds added and removed genes from the presence of the stable
#' identifier. To find the successor of a retired identifier, you must do a
#' separate lookup.
#'
#' @param species Species in any form that [resolve_species()] accepts.
#' @param from,to The two releases. If `to` is the older release, `"added"` and
#'   `"removed"` change places.
#' @param assembly,source The function passes these to [fetch_mapping()]. It
#'   reads the two releases from the same assembly and annotation family, to
#'   isolate the release changes.
#' @param quiet Do not show the summary message.
#' @return A data frame of `id`, `change` (`"added"`, `"removed"`, `"renamed"`
#'   or `"retyped"`), `name_from`, `name_to`, `biotype_from` and `biotype_to`.
#' @seealso [detect_release()] to find which releases you have.
#' @export
#' @examplesIf interactive()
#' d <- compare_releases("human", from = 110, to = 116)
#' table(d$change)
#' d[d$change == "renamed", ]
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
