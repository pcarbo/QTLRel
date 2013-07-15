\name{lodci}
\alias{lodci}
\title{
Estimate LOD Support Intervals
}
\description{
Estimate LOD support intervals.
}
\usage{
lodci(llk,cv=0,lod=1.5,drop=3)
}
\arguments{
  \item{llk}{
a data frame with components (chr, dist, y, ...), where "chr" is the chromosome on which the scanning locus is located, "dist" is the genetic or physical position of the scanning locus, and "y" is the test statistic.
}
  \item{cv}{
threshold. Reported support intervals cover at least one scanning locus where \code{llk$y > cv}.
}
  \item{lod}{
\code{lod} (1.5 by default) LOD support intervals are reported when \code{llk$y} is converted to LOD score.
}
  \item{drop}{
3 by default. See "details".
}
}
\details{
In case of multiple peaks on a chromosome, a peak has to satisfy: a) above the threshold \code{cv}; b) drops, e.g., 3 LOD on both sides except chromosome ends. So if two peaks close to each other but LOD between them doesn't drop, e.g., 3 LOD, only one of them is considered.
}
\value{
A data frame with the following components:
\item{chr}{the chromosome}
\item{lower}{the lower bound}
\item{upper}{the upper bound}
\item{index}{indicates which scanning loci}
}
\examples{
data(miscEx)

\dontrun{
# impute missing genotypes
gdtmp<- gdatF8
   gdtmp<- replace(gdtmp,is.na(gdtmp),0)
# run 'genoProb'
prDat<- genoProb(gdat=gdtmp, gmap=gmapF8, step=Inf,
   gr=8, method="Haldane", verbose=TRUE)
# estimate variance components
o<- estVC(y=pdatF8$bwt, x=pdatF8$sex, v=list(AA=gmF8$AA,DD=gmF8$DD,
	HH=NULL, AD=NULL, MH=NULL, EE=diag(length(pdatF8$bwt))))

# genome scan
llk.hk<- scanOne(y=pdatF8$bwt, x=pdatF8$sex, vc=o, prdat=prDat)

# extract LOD support intervals
tmp<- data.frame(y=llk.hk$p, chr=llk.hk$chr, dist=llk.hk$dist)
lodci(tmp, cv=10, lod=1.5, drop=3)
}
}

