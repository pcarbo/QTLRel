\name{scanTwo}
\alias{scanTwo}
\title{
Genome Scan for Epistasis
}
\description{
Evaluate log-likelihood ratio test statistic for epistasis (QTL by QTL interaction).
}
\usage{
scanTwo(y, x, gdat, prdat=NULL, vc=NULL, numGeno=FALSE,
   minorGenoFreq=0, rmv=TRUE)
}
\arguments{
  \item{y}{
a numeric vector or a numeric matrix of one column (representing a phenotype for instance).
}
  \item{x}{
a data frame or matrix, representing covariates if not missing.
}
  \item{gdat}{
genotype data. Should be a matrix or a data frame, with each row representing an observation and each column a marker locus. The column names should be marker names. Optional if an object \code{prdat} from \code{\link{genoProb}} is used as an argument.
}
  \item{prdat}{
an object from \code{\link{genoProb}}.
}
  \item{vc}{
an object from \code{\link{estVC}} or \code{\link{aicVC}}, or an estimated variance-covariance matrix induced by relatedness and environment.
}
  \item{numGeno}{
whether to treat numeric coding of genotypes as numeric.
}
  \item{minorGenoFreq}{
specify the minimum tolerable minor genotype frequency at a scanning locus if \code{gdat} is used.
}
  \item{rmv}{
a logical variable. If true, then the scanning locus will be skipped if the minor genotype frequency at the locus is smaller than \code{minorGenoFreq}. Otherwise, the scanning process will stop and return with NULL.
}
}
\value{
A matrix whose entry in the upper triangle is the log-likelihood test statistic for epistatic effect.
}

\seealso{
\code{\link{scanOne}}.
}

