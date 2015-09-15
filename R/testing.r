## Testing code ####

make_exon_getter <- function(n) {
  n <- as.character(n)
  function(x) {
    eag_num <- x[exon == n, unique(eag_num)]
    if (length(eag_num) == 2) {
      eag_num
    } else if (length(eag_num) == 1) {
      rep(eag_num, 2)
    } else stop("More than three EAGs for ", x[, unique(lims_donor_id)])
  }
}

exon2 <- make_exon_getter(2)

exon3 <- make_exon_getter(3)

get_nums <- function(id, testdata = HLAdata::testdata_dpb1_1267) {
  rs <- testdata[lims_donor_id == id]
  eag_numbers(
    e2 = rs[exon == 2 & trimws(main_eag) %in% 1:2, unique(eag_num)] %||% NA_character_,
    e3 = rs[exon == 3 & trimws(main_eag) %in% 1:2, unique(eag_num)] %||% NA_character_
  )
}

LDID_checker <- function(ldids, lookup, tdata) {
  iID <- iterators::iter(ldids)
  map <- make_mapper(lookup)
  remap <- make_remapper(lookup)
  jkrs <- lookup$jkr_tbl$eag_allele
  id <- NULL
  td <- NULL
  function(rerun = FALSE, remove_row = NULL, debug_mapper = FALSE, debug_remapper = FALSE, show_jokers = TRUE) {
    if (!rerun) {
      id <<- iterators::nextElem(iID)
      td <<- tdata[lims_donor_id == id, ]
    }
    if (!is.null(remove_row)) {
      td <- td[-remove_row]
    }
    if (debug_mapper) {
      base::debug(map)
      on.exit(base::undebug(map))
    }
    if (debug_remapper) {
      base::debug(remap)
      on.exit(base::undebug(remap), add = TRUE)
    }
    ## Testing
    num <- get_nums(id, td)
    rep <- map(num)
    new_num <- remap(rep)
    new_rep <- map(new_num)

    ## Printing
    print(td)
    cat("\nEAG nums \t\t\tRemapped EAG nums", sep = "\n")
    fmt <- " a:%-9sc:%-9s \t a:%-9sc:%-9s \n b:%-9sd:%-9s \t b:%-9sd:%-9s"
    cat(sprintf(fmt, num[1], num[3], new_num[1], new_num[3],
                num[2], num[4], new_num[2], new_num[4]), sep = "\n")

    cat("\nMapped allele \t\t\tBackmapped allele", sep = "\n")
    fmt <- " a: %s \t\t\t a: %s \n b: %s \t\t\t b: %s"
    cat(sprintf(fmt, allele(rep, 1), allele(new_rep, 1),
                allele(rep, 2), allele(new_rep, 2)), sep = "\n")

    if (show_jokers) {
      cat("\n")
      cat(sQuote(allele(rep, 1)), " is ", if (allele(rep, 1) %in% jkrs) "" else "not ", "in the joker group\n", sep = "")
      cat(sQuote(allele(rep, 2)), " is ", if (allele(rep, 2) %in% jkrs) "" else "not ", "in the joker group\n", sep = "")
    }
    invisible(list(num = num, rep = rep, new_num = new_num, new_rep = new_rep))
  }
}

#' Retrieve testing data from ngsrep for a gene
#'
#' @param gene A HLA gene.
#' @param nextype_basis_id A neXtype Basis ID.
#'
#' @return
#' A table with the fields:
#' \describe{
#'  \item{lims_donor_id}{}
#'  \item{nextype_run_id}{}
#'  \item{orig_nmdp1}{}
#'  \item{orig_nmdp2}{}
#'  \item{exon}{}
#'  \item{eag_num}{}
#'  \item{main_eag}{}
#'  \item{total_reads}{}
#'  \item{best_reads}{}
#' }
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' testdata <- get_testdata("DPB1", "1267")
#' }
get_testdata <- function(gene, nextype_basis_id) {
  stopifnot(requireNamespace("orcl"))
  gene <- match_hla_gene(gene)
  fmt <- "
  SELECT b.LIMS_DONOR_ID, b.NEXTYPE_RUN_ID, a.ORIG_NMDP1,
    a.ORIG_NMDP2, TRIM(b.EXON) AS EXON, b.EAG_NUM, b.MAIN_EAG,
    b.TOTAL_READS, b.BEST_READS, b.BEST_READS / b.TOTAL_READS AS BEST_TO_TOTAL
  FROM RES_NMDP_CODES a
  INNER JOIN RES_MATRIX_MASTER b
  ON a.GENE               = b.GENE
  AND a.NEXTYPE_RUN_ID    = b.NEXTYPE_RUN_ID
  AND a.NEXTYPE_SAMPLE_ID = b.NEXTYPE_SAMPLE_ID
  AND a.LIMS_DONOR_ID     = b.LIMS_DONOR_ID
  WHERE (b.RESULT = 'E' OR b.RESULT = 'N')
  AND a.NEXTYPE_BASIS_ID  = '%s'
  AND a.GENE              = '%s'"
  db <- orcl::Ngsrep()
  testdata_all <- db$query(sprintf(fmt, nextype_basis_id, gene))
  db$disconnect()
  clean_testing_data(td = testdata_all)
}

clean_testing_data <- function(td) {
  pat <- "^\\[2\\]\\d+/\\d+:\\[3\\]\\d+/\\d+$"
  ## Remove NAs, NotAvailable, and unresolved ambiguities.
  td <- td[!is.na(orig_nmdp1) & !is.na(orig_nmdp2)]
  td <- td[orig_nmdp1 != "NotAvailable" & orig_nmdp2 != "NotAvailable"]
  td <- td[!grepl("/", orig_nmdp1, fixed = TRUE) & !grepl("/", orig_nmdp2, fixed = TRUE)]
  ## If a LIMS_DONOR_ID was used across multiple runs, take only data from the last
  ## run
  setkeyv(td, c("lims_donor_id", "nextype_run_id"))
  non_dups2 <- td[, .N, by = .(lims_donor_id, nextype_run_id)][
    !duplicated(lims_donor_id, fromLast = TRUE)][
      , .(lims_donor_id, nextype_run_id)]
  td2 <- td[non_dups2]
  rs <- td2[, paste0(sort(at_least_two(unique(eag_num))), collapse = "/"),
            by = .(lims_donor_id, exon)][
              , paste0(paste0("[", exon, "]", V1), collapse = ":"), by = lims_donor_id]
  rs <- rs[grepl(pattern = pat, x = V1)]
  non_dups <- rs[!duplicated(rs[, V1]), .(lims_donor_id)]
  setkeyv(non_dups, "lims_donor_id")
  td3 <- td2[non_dups]
  #   dt <- td3 %>%
  #     dplyr::group_by(lims_donor_id) %>%
  #     dplyr::summarise(lower_limit = median(total_reads) - 3*mad(total_reads))
  #   td3 %>%
  #     right_join(dt, by = "lims_donor_id") %>%
  #     dplyr::filter(total_reads > lower_limit & best_reads > 30) %>%
  #     dplyr::select(-lower_limit) %>%
  #     data.table()
  td3
}

#' Make a tester
#'
#' Generate a list of functions \code{$test()}, \code{$get_successes(all = FALSE)},
#' \code{$get_failures(all = FALSE)}, \code{$reset()}, \code{$get_run_ldids()},
#' that allow testing the performance of the mapper against NeXtype results fro
#' a test dataset.
#'
#' @param testdata A dataset generated by a call to \code{\link{get_testdata}}.
#' @param lookup A <\code{\link{lookup_table}}>.
#' @return A list of functions.
#'
#' @keywords internal
#' @export
#' @examples
#' lookup <- lookup_list(rep_dpb1, eag)
#' tester <- make_test_instance(testdata, lookup)
#' tester$test(limit = 20)
#' tester$get_successes()
#' tester$get_failures()
#' tester$test(limit = 20)
#' tester$get_failures()
#'
#' run_tests_on("ID11732060", testdata, lookup)
#'
make_test_instance <- function(testdata, lookup) {
  ## global funcs
  `%do%` <- foreach::`%do%`
  map    <- make_mapper(lookup)
  remap  <- make_remapper(lookup)
  ## global variables
  lims_donor_ids <- testdata[, unique(lims_donor_id)]
  tested         <- character(0)
  successes      <- character(0)
  ambiguities    <- character(0)
  failures       <- character(0)
  ## return funcs
  list(
    test = function(limit = -1, ignore_remapping = FALSE, exclude_missing_exon3 = FALSE, suppress_msgs = FALSE) {
      if (exclude_missing_exon3) {
        ldid <- testdata[, .N, by = .(lims_donor_id, exon)][duplicated(lims_donor_id), lims_donor_id]
        testdata <- testdata[lims_donor_id == ldid]
        lims_donor_ids <- testdata[, unique(lims_donor_id)]
      }
      ## reset the global variables
      successes      <<- character(0)
      ambiguities    <<- character(0)
      failures       <<- character(0)
      ldids_to_test <- setdiff(lims_donor_ids, tested)
      n <- if (limit < 0) length(ldids_to_test) else limit
      all_run_ldids <- ldids_to_test[1:n]
      ## set up progress bar
      pbar <- utils::txtProgressBar(min = 1, max = length(all_run_ldids), style = 3)
      cat("\nRunning tests on ", length(all_run_ldids), " allelic combinations\n\n", sep = "")
      on.exit(close(pbar))
      iCnt <- iterators::icount()
      iLDID <- iterators::iter(all_run_ldids)

      null <- foreach::foreach(ldid = iLDID, cnt = iCnt) %do% {
        tested <<- unique(c(tested, ldid))
        rs <- testdata[lims_donor_id == ldid]
        if (length(unique(rs[, nextype_run_id])) > 1) {
          warning("Fix duplicate runs in ", sQuote(ldid), call. = FALSE, immediate. = TRUE)
          failures <<- unique(c(failures, ldid))
        }
        else {
          nums <- eag_numbers(
            e2 = rs[exon == 2 & trimws(main_eag) %in% 1:2, unique(eag_num)] %||% NA_character_,
            e3 = rs[exon == 3 & trimws(main_eag) %in% 1:2, unique(eag_num)] %||% NA_character_
          )
          al1 <- rs[, unique(orig_nmdp1)]
          al2 <- rs[, unique(orig_nmdp2)]
          if (is.na(al1) || is.na(al2)) {
            warning("Fix missing allele in ", sQuote(ldid), call. = FALSE, immediate. = TRUE)
            failures <<- unique(c(failures, ldid))
          }
          else {
            ## 1) Map original EAG nums to allele.
            als <- at_least_two(union(al1, al2))
            rep <- rep_allele(a1 = als[1], a2 = als[2])
            rep_new <- tryCatch(map(nums), error = function(e) {
              warning("Mapper failure in ", sQuote(ldid), ": ", e$message, call. = FALSE, immediate. = TRUE)
              failures <<- unique(c(failures, ldid))
            })
            if (is.na(rep_new) || rep != rep_new) {
              msg <- sprintf("Mapping error in %s: was:  %s/%s; should have been:  %s/%s",
                             sQuote(ldid), allele(rep_new, 1), allele(rep_new, 2), allele(rep, 1), allele(rep, 2))
              warning(msg, immediate. = TRUE, call. = FALSE)
              failures <<- unique(c(failures, ldid))
            }
            else if (!ignore_remapping) {
              ## 2) Remap the EAG nums based on the new allele.
              nums_new <- tryCatch(remap(rep_new), error = function(e) {
                warning("Remapper failure in ", sQuote(ldid), ": ", e$message, call. = FALSE, immediate. = TRUE)
                failures <<- unique(c(failures, ldid))
              })
              if (is.na(nums) || is.na(nums_new) || nums != nums_new) {
                ambiguities <<- unique(c(ambiguities, ldid))
                if (!suppress_msgs) {
                  msg <- sprintf("Remapping ambiguity in %s: was: %s/%s; should have been: %s/%s",
                                 sQuote(ldid), colon(exon(nums_new, 2)), colon(exon(nums_new, 3)),
                                 colon(exon(nums, 2)), colon(exon(nums, 3)))
                  message(msg, appendLF = TRUE)
                }
              }
              ## 3) Map the remapped EAG nums to allele.
              rep_back <- tryCatch(map(nums_new), error = function(e) {
                warning("Backmapper failure in ", sQuote(ldid), ": ", e$message, call. = FALSE, immediate. = TRUE)
                failures <<- unique(c(failures, ldid))
              })
              if (is.na(rep_back) || rep != rep_back) {
                msg <- sprintf("Backmapping error in %s: was: %s/%s; should have been: %s/%s",
                               sQuote(ldid), allele(rep_back, 1), allele(rep_back, 2), allele(rep, 1), allele(rep, 2))
                warning(msg, immediate. = TRUE, call. = FALSE)
                failures <<- unique(c(failures, ldid))
              }
              else {
                utils::setTxtProgressBar(pb = pbar, value = cnt)
                successes <<- unique(c(successes, ldid))
              }
            }
          }
        }
      }
    },
    get_successes = function() {
      successes
    },
    get_failures = function() {
      failures
    },
    get_ambigs = function() {
      ambiguities
    },
    get_tested = function() {
      tested
    },
    reset = function() {
      successes <<- character(0)
      failures <<- character(0)
      ambiguities <<- character(0)
      tested <<- character(0)
    },
    get_all_ids = function() {
      lims_donor_ids
    }
  )
}

#' @param alleles A vector of alleles
#' @param mapper A <mapper> function.
#' @param remapper A <remapper> fucntion.
#' @rdname make_test_instance
#' @keywords internal
make_test_instance2 <- function(alleles, mapper, remapper) {
  ## global funcs
  `%do%` <- foreach::`%do%`
  map    <- make_mapper(lookup)
  ## global objects
  aa <- combn(alleles, 2, function(x) rep_allele(x[1], x[2]), simplify = FALSE)
  all_successes  <- list()
  all_failures   <- list()
  run_successes  <- list()
  run_failures   <- list()
  ## return funcs
  list(
    test = function(limit = -1) {
      ## reset the run successes
      run_successes <<- list()
      run_failures  <<- list()
      alleles_done <- length(all_successes) + length(all_failures) + 1
      n <- if (limit < 0) length(aa) else alleles_done + limit
      iAA <- iterators::iter(aa[alleles_done:n])
      null <- foreach::foreach(a = iAA) %do% {
        nums <- tryCatch(remapper(a), error = function(e) {
          warning("Remapper error in ", allele(a, 1), "/", allele(a, 2), ": ", e$message,
                  call. = FALSE, immediate. = TRUE)
          run_failures <<- c(run_failures, list(a))
          all_failures <<- c(all_failures, list(a))
          eag_numbers(NA_character_, NA_character_)
        })

        if (!all(is.na(nums))) {
          anew <- tryCatch(mapper(nums), error = function(e) {
            warning("Mapper error in ", allele(a, 1), "/", allele(a, 2), ": ", e$message,
                    call. = FALSE, immediate. = TRUE)
            run_failures <<- c(run_failures, list(a))
            all_failures <<- c(all_failures, list(a))
            eag_numbers(NA_character_, NA_character_)
          })

          if (a == anew) {
            run_successes <<- c(run_successes, list(a))
            all_successes <<- c(all_successes, list(a))
            cat(".", sep = "")
          }
          else {
            msg <- sprintf("Mapping error in: %s/%s\tgave: %s/%s",
                           allele(a, 1), allele(a, 2), allele(anew, 1), allele(anew, 2))
            warning(msg, immediate. = TRUE, call. = FALSE)
            run_failures <<- c(run_failures, list(a))
            all_failures <<- c(all_failures, list(a))
          }
        }
      }
    },
    get_successes = function(all = FALSE) {
      if (all) all_successes else run_successes
    },
    get_failures = function(all = FALSE) {
      if (all) all_failures else run_failures
    },
    reset = function() {
      all_successes <<- character(0)
      all_failures <<- character(0)
    }
  )
}

update_testdata <- function(ldid, rnid, remove_row, td, overwrite = FALSE) {
  td1 <- td[lims_donor_id == ldid & nextype_run_id == rnid][-remove_row]
  td2 <- td[lims_donor_id != ldid]
  td <- rbind(td1, td2)
  setkeyv(td, c("lims_donor_id", "nextype_run_id"))
  if (overwrite) {
    testdata <- td
    devtools::use_data(testdata, overwrite = TRUE)
  }
  td
}

