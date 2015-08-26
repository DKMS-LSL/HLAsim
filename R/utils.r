is.empty <- function(x) {
  is.null(x) || length(x) == 0L || (length(x) == 1L && !nzchar(x))
}

nunique <- function(x, ...) {
  if (is.factor(x)) {
    length(levels(x))
  } else {
    length(unique(x, ...))
  }
}

chin <- data.table::`%chin%`

`%ni%` <- Negate(`%in%`)

`%.%` <- purrr::compose

`%||%` <- function (a, b) if (is.empty(a)) b else a

hla_burst <- function(a) {
  strsplit(a, split = ":", fixed = TRUE)
}

hla_collapse <- function(a) {
  paste0(a, collapse = ":")
}

strip_field_four <- function(x, unique = TRUE) {
  if (is.null(x)) {
    return(NULL)
  }
  if (all(is.na(x))) {
    return(NA_character_)
  }
  fun <- function(x) {
    x <- x[1:3]
    hla_collapse(x[!is.na(x)])
  }
  ans <- vapply(hla_burst(x), fun, FUN.VALUE = character(1))
  if (unique) {
    unique(ans)
  } else ans
}

allele1 <- function(x) {
  vapply(strsplit(x, split = "/", fixed = TRUE), `[`, 1, FUN.VALUE = "")
}

allele2 <- function(x) {
  vapply(strsplit(x, split = "/", fixed = TRUE), `[`, 2, FUN.VALUE = "")
}

is_empty <- function(x) length(x) == 0

none_empty <- function(...) {
  all(vapply(list(...), function(x) length(x) != 0L, FUN.VALUE = FALSE))
}

maybe_duplicated <- function(x) {
  length(x) > 1 && any(duplicated(x))
}

colon <- function(...) paste0(..., collapse = ":")

slash <- function(...) paste0(..., collapse = "/")

is.hla <- function(x) {
  p <- "((\\d\\d\\d?):([[:alnum:]]{2,})(:[[:alnum:]]{2,})?(:[[:alnum:]]{2,})?)|(NEW)"
  grepl(p, x)
}

is_nmdp <- function(a) {
  p <- "^\\d{2, }:[A-Z]+$"
  grepl(p, a)
}

is_g_code <- function(a) {
  p <- "^\\d\\d+:\\d\\d+:\\d\\d+G$"
  grepl(p, a)
}

is_ambiguous <- function(a) {
  p <- "(^\\d{2, }:[A-Z]+$)|(^\\d\\d+:\\d\\d+:\\d\\d+G$)"
  grepl(p, a)
}

is.genotype <- function(x) {
  ptn <- "^(\\d\\d\\d?:)+[[:alnum:]]+/(\\d\\d\\d?:)+[[:alnum:]]+$"
  grepl(ptn, x)
}

is.homozygous <- function(a1, a2) {
  ifelse(a1 == a2, TRUE, FALSE)
}

is_homozygous <- function(gtp, sep = "/") {
  sgtp <- strsplit(gtp, sep, fixed = TRUE)
  vapply(sgtp, function(x) x[1] == x[2], FUN.VALUE = logical(1))
}

is.heterozygous <- function(a1, a2) {
  ifelse(a1 != a2, TRUE, FALSE)
}

is_heterozygous <- function(gtp, sep = "/") {
  sgtp <- strsplit(gtp, sep, fixed = TRUE)
  vapply(sgtp, function(x) x[1] != x[2], FUN.VALUE = logical(1))
}

maybe_exon_shuffling <- function(a) {
  is_ambiguous(allele(a, 1)) && is_ambiguous(allele(a, 2))
}

match_hla_gene <- function(gene) {
  gene <- match.arg(gene, c("A", "B", "C", "DPB1", "DRB1", "DQB1"))
  paste0('HLA-', gene)
}

valid_date <- function(x) {
  if (!grepl("\\d\\d/\\d\\d/\\d{4}", x)) {
    stop(sQuote(x), " is no valid date in format 'DD/MM/YYYY'",
         call. = FALSE)
  }
  TRUE
}

strip_date <- function(x, fmt = "%Y-%m-%d") {
  sx <- strftime(x, format = fmt)
  as.POSIXct(sx, format = fmt)
}

at_least_two <- function(x) {
  if (length(x) == 1L) rep(x, 2) else x[order(x)]
}

max_table <- function(x) {
  rs <- tabulate(x)
  levels(x)[rs == max(rs)]
}

allele2string <- function(x, split = "/", ...) {
  purrr::map_v(x, ~ paste0(., collapse = split))
}

string2allele <- function(x, split = "/", ...) {
    purrr::map(strsplit(x, split = split, fixed = TRUE), ~ rep_allele(.[1], .[2]))
}

strsplitN <- function(x, split, n, from = "start", collapse = split, ...) {
  from <- match.arg(from, c("start", "end"))
  xs <- strsplit(x, split, ...)
  end <- vapply(xs, length, 0L)
  if (from == "end") {
    end <- end + 1L
    n <- lapply(end, `-`, n)
    n <- .mapply(`[<-`, list(x = n, i = lapply(n, `<`, 0), value = 0L), NULL)
  } else {
    n <- lapply(rep.int(0, length(xs)), `+`, n)
    n <- .mapply(`[<-`, list(x = n, i = Map(`>`, n, end), value = end), NULL)
  }
  n <- lapply(n, sort %.% unique)
  unlist(.mapply(function(x, n) paste0(x[n], collapse = collapse), list(x = xs, n = n), NULL))
}

starts_with <- function(p, s, ignore.case = FALSE) {
  p <- paste0("^", p, ".*$")
  grepl(p, s, fixed = FALSE, perl = FALSE, ignore.case = ignore.case)
}

map_nextype_basis_to_dna_version <- function() {
  stmt <- '
  select distinct nextype_basis_id, dna_version_id
  from ngsrep.nextype_alleles_per_eag'
  rs <- tbl_dt(orcl::ora_query(stmt))
  setkeyv(rs, "nextype_basis_id")
  rs
}

map_dna_version_to_nextype_basis <- function() {
  rs1 <- map_nextype_basis_to_dna_version()
  stmt <- '
  select nextype_basis_id,
    count(lims_donor_id) as count
  from ngsrep.res_nmdp_codes
  group by nextype_basis_id
  order by count(lims_donor_id)'
  rs2 <- tbl_dt(orcl::ora_query(stmt))
  setkeyv(rs2, "nextype_basis_id")
  left_join(rs1, rs2, by = "nextype_basis_id") %>%
    group_by(dna_version_id) %>%
    dplyr::filter(min_rank(desc(count)) == 1) %>%
    ungroup() %>%
    dplyr::select(dna_version_id, nextype_basis_id)
}

#' @export
dna2basis <- function(id) {
  dna2basis_[dna_version_id == id, nextype_basis_id]
}

#' @export
basis2dna <- function(id) {
  basis2dna_[nextype_basis_id == id, dna_version_id]
}
