\name{scanOne}
\alias{scanOne}
\title{
Genome Scan for QTL
}
\description{
Evaluate log-likelihood ratio test statistics or P-values at scanning loci along the genome.
}
\usage{
scanOne(y, x, gdat, prdat=NULL, vc=NULL, intcovar=NULL,
   numGeno=FALSE, test=c("None","F","Chisq"),
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
an object from \code{\link{genoProb}}, or in the same form.
}
  \item{vc}{
an object from \code{\link{estVC}} or \code{\link{aicVC}}, or an estimated variance-covariance matrix induced by relatedness and environment.
}
  \item{intcovar}{
Covariates that interact with QTL.
}
  \item{numGeno}{
whether to treat numeric coding of genotypes as numeric.
}
  \item{test}{
"None", "F" or "Chisq".
}
  \item{minorGenoFreq}{
specify the minimum tolerable minor genotype frequency at a scanning locus if \code{gdat} is used.
}
  \item{rmv}{
a logical variable. If true, then the scanning locus will be skipped if the minor genotype frequency at the locus is smaller than \code{minorGenoFreq}. Otherwise, the scanning process will stop and return with NULL.
}
}
\details{
The test at a scanning locus under the assumption of no QTL effect versus the assumption of QTL effect is performed by conditioning on the estimated polygenic genetic variance-covariance matrix. Normality is assumed for the random effects.

It is possible to extend the Haley-Knott approach to multiple-allelic cases under the assumption that allele effects are all additive. Then, \code{prdat} should be provided and be of class "addEff".
}
\value{
A list with at least the following components:
   \item{p}{P-value at the snp (marker) if \code{test} is "F" or "Chisq", or the log-likelihood ratio statistic at the SNP (marker) if "test" is "None"}
   \item{v}{Percentage of variation explained by QTL related effects at the snp (marker)}
   \item{parameters}{estimated parameters at all scanning loci, including additive effect \code{a} and dominance effect \code{d} if \code{prdat} is not NULL}
}
\references{
Haley, C. S., and S. A. Knott (1992). A simple regression method for mapping quantitative trait loci in line crosses using flanking markers. Heredity 69: 315-324.
}

\seealso{
\code{\link{genoImpute}} and  \code{\link{genoProb}}.
}
\examples{
data(miscEx)

\dontrun{
# impute missing genotypes
gdatTmp<- genoImpute(gdatF8, gmap=gmapF8, step=Inf,
   gr=8, na.str=NA)

# estimate variance components
o<- estVC(y=pdatF8$bwt, x=pdatF8$sex, v=list(AA=gmF8$AA,DD=gmF8$DD,
	HH=NULL, AD=NULL, MH=NULL, EE=diag(length(pdatF8$bwt))))

# genome scan and plotting
pv<- scanOne(y=pdatF8$bwt, x=pdatF8$sex, gdat=gdatTmp, vc=o)
pv
plot(pv,gmap=gmapF8)

# Haley-Knott method
gdtmp<- gdatF8; unique(unlist(gdtmp))
   gdtmp<- replace(gdtmp,is.na(gdtmp),0)
prDat<- genoProb(gdat=gdtmp, gmap=gmapF8, step=Inf,
   gr=8, method="Haldane", verbose=TRUE)
pv.hk<- scanOne(y=pdatF8$bwt, x=pdatF8$sex, prdat=prDat, vc=o)
pv.hk
plot(pv.hk)
}
}

