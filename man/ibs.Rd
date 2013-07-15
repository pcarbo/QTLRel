\name{ibs}
\alias{ibs}
\title{
Estimate Jacquard condensed identity coefficients
}
\description{
Estimate Jacquard condensed identity coefficients by identity-by-state (IBS) from genotypic data.
}
\usage{
ibs(x)
}
\arguments{
  \item{x}{
genotype data with genotypes ("AA", "AB", "BB", or, 1, 2, 3) and without missing data, or probabilities for these genotypes (e.g., obtained by using \code{\link{genoProb}}). In case of genotype data, rows represent individuals and columns represent SNPs.
}
}
\value{
A matrix G with G[,j] being the j-th Jacquard identity coefficients.
}
\note{
Currently only support the two-allele data.
}
\seealso{
\code{\link{genMatrix}}
}

