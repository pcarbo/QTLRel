\name{blup}
\alias{blup}
\title{
Best Linear Unbiased Prediction
}
\description{
Estimate the best linear unbiased prediction (BLUP) for various effects in the model.
}
\usage{
blup(object)
}
\arguments{
  \item{object}{
an object from \code{\link{estVC}} or \code{\link{aicVC}}.
}
}
\value{
\item{fixed}{BLUP for fixed effects.}
\item{AA,DD,\dots}{BLUP for random (genetic) variance components "AA", "DD", \dots}
\item{EE}{BLUP for residual effect.}
}

\seealso{
\code{\link{estVC}} and \code{\link{aicVC}}.
}
\examples{
data(miscEx)

\dontrun{
# only consider additive genetic variance component
o<- estVC(y=pdatF8$bwt, x=pdatF8$sex, v=list(AA=gmF8$AA,DD=gmF8$DD,
	HH=NULL, AD=NULL, MH=NULL, EE=diag(length(pdatF8$bwt))))
b<- blup(o)
}
}

