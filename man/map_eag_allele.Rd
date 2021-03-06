% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapping.r
\name{map_eag_allele}
\alias{map_eag_allele}
\title{Map eag numbers to an eag allele}
\usage{
map_eag_allele(nums, lookup)
}
\arguments{
\item{nums}{EAG numbers for exons 2 and 3 as
an <\code{\link{eag_numbers}}> object.}

\item{lookup}{A <\code{\link{lookup_list}}> object.}
}
\value{
An <\code{\link{eag_allele}}>.
}
\description{
Map eag numbers to an eag allele
}
\examples{
nums <- eag_numbers(e2 = c(6673401, 6673404), e3 = c(6672138, 6670791))
lookup <- lookup_list(rep_dpb1, eag)

map_eag_allele(nums, lookup)

}
\keyword{internal}

