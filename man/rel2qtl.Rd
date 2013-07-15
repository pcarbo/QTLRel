\name{rel2qtl}
\alias{rel2qtl}

\title{Convert data from QTLRel to R/qtl format}

\description{
  Convert the data for a QTL mapping experiment from the QTLRel format
  to that used by R/qtl (\url{http://www.rqtl.org}).
}

\usage{
rel2qtl(gdat, pdat, gmap)
}

\arguments{
 \item{gdat}{Genotype data}
 \item{pdat}{Phenotype data}
 \item{gmap}{Genetic map}
}

\value{
  A cross object for the R/qtl package (\url{http://www.rqtl.org}).
}

\details{
  Pedigree information is ignored, and X chromosome data is omitted.

  The data are treated as an intercross.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(miscEx)
f8 <- rel2qtl(gdatF8, pdatF8, gmapF8)
summary(f8)
}

\seealso{ \code{\link{qtl2rel}} }

\keyword{manip}
