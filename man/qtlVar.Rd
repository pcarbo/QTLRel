\name{qtlVar}
\alias{qtlVar}
\title{
QTL Variance
}
\description{
Estimate variance in a quantitative trait induced by QTL.
}
\usage{
qtlVar(lrt,prdat,simulation=FALSE,nsim=25)
}
\arguments{
  \item{lrt}{
a data frame (a, d, ...), where 'a' and 'd' are respectively additive and dominance effects.
}
  \item{prdat}{
a 3-D array that provides probabilities of genotypes "AA", "AB" and "BB". If \code{prDat} is an object of \code{\link{genoProb}}, then \code{prdat} can be \code{prDat$pr}. 
}
  \item{simulation}{
whether to use simulations to estimate the variance explained by QTL.
}
  \item{nsim}{
number of simulations to perform if \code{simulation} is TRUE.
}
}
\value{
A vector displaying the estimated variance at each loci.
}
\note{
Correlations among observations are ignored, and this function should be used with caution.
}
\seealso{
\code{\link{scanOne}} and \code{\link{genoProb}}
}
\examples{
data(miscEx)

\dontrun{
# impute missing genotypes
gdtmp<- gdatF8
   gdtmp<- replace(gdtmp,is.na(gdtmp),0)
# rung 'genoProb'
prDat<- genoProb(gdat=gdtmp, gmap=gmapF8, step=Inf,
   gr=8, method="Haldane", verbose=TRUE)
# estimate variance components
o<- estVC(y=pdatF8$bwt, x=pdatF8$sex, v=list(AA=gmF8$AA,DD=gmF8$DD,
	HH=NULL, AD=NULL, MH=NULL, EE=diag(length(pdatF8$bwt))))

# genome scan
pv.hk<- scanOne(y=pdatF8$bwt, x=pdatF8$sex, prdat=prDat, vc=o)

# run 'qtlVar'
qef<- NULL
for(n in 1:length(llk.hk$par))
   qef<- rbind(qef,llk.hk$par[[n]][c("a","d")])
   qef<- as.data.frame(qef)
qtlVar(qef,prDat$pr)[1:3]
}
}

