\name{mAIC}
\alias{mAIC}
\alias{mAIC.HK}
\title{
Multiple QTL AIC
}
\description{
Multiple QTL model selection by AIC criterion.
}
\usage{
mAIC(y,x,gdat,prdat=NULL,vc=NULL,chrIdx,xin,k=2,
   direction=c("both","backward","forward"),ext=FALSE,verbose=FALSE)
}
\arguments{
 \item{y}{
a numeric vector or a numeric matrix of one column (representing a phenotype for instance).
}
  \item{x}{
a data frame or matrix, representing covariates if not missing.
}
  \item{gdat}{
genotype data. Should be a matrix or a data frame, with each row representing an observation and each column a marker locus. The column names should be marker names. Numeric coding of genotype is treated as numeric. Ignored if \code{prdat} is an object from \code{\link{genoProb}}.
}
  \item{vc}{
an object from \code{\link{estVC}} or \code{\link{aicVC}}, or an estimated variance-covariance matrix induced by relatedness. The scan will assume no polygenic variation if \code{vc} is NULL.
}
  \item{prdat}{
an object from \code{\link{genoProb}}.
}
  \item{chrIdx}{
chromsome index of markers in columns of \code{gdat} if given. Ignored if \code{prdat} is an object from \code{\link{genoProb}}.
}
  \item{xin}{
vector indicating whether a locus is already in the model.
}
  \item{k}{
penalty on a parameter. The selection criterion is the known "AIC" if \code{k = 2} and is "BIC" if \code{k = log(n)} where "n" is the sample size.
}
  \item{direction}{
the mode of search: "both", "forward" or "backward" with default "both".
}
  \item{ext}{
a logical variable. True if ones wants more exhaustive search.
}
  \item{verbose}{
a logical variable. True if ones wants to track the process for monitoring purpose.
}
}
\details{
Makes use of "Haley-Knott" method (Haley and Knott 1992) if \code{prdat} is an object from \code{\link{genoProb}}.
}
\value{
a list with the following components:

   \code{model}: the resulting model;

   \code{aic}: AIC of the model;

   \code{snp}: selected SNPs.

   \code{xin}: vector indicating whether a SNP is selected.
}
\references{
Haley, C. S., and S. A. Knott (1992). A simple regression method for mapping quantitative trait loci in line crosses using flanking markers. Heredity 69: 315-324.
}
\note{
Currently only suitable for advanced intercross lines (or diallelic data).
}
\seealso{
\code{\link{optim}}, \code{\link{genoProb}} and \code{\link{aicVC}}.
}

\examples{
data(miscEx)

\dontrun{
# impute missing genotypes
gdat.imp<- genoImpute(gdatF8, gmap=gmapF8, step=Inf,
   gr=8, na.str=NA)
# estimate variance components
o<- estVC(y=pdatF8$bwt, x=pdatF8$sex, v=list(AA=gmF8$AA,DD=gmF8$DD,
	HH=NULL, AD=NULL, MH=NULL, EE=diag(length(pdatF8$bwt))))

# run 'genoProb'
gdtmp<- gdatF8
   gdtmp<- replace(gdtmp,is.na(gdtmp),0)
prDat<- genoProb(gdat=gdtmp, gmap=gmapF8, step=Inf,
   gr=8, method="Haldane", verbose=TRUE)

# genome scan
llk.hk<- scanOne(y=pdatF8$bwt, x=pdatF8$sex, prdat=prDat, vc=o)
xin<- llk.hk$p > 10

# run 'mAIC' based on genome scan results
mg<- mAIC(y=pdatF8$bwt, x=pdatF8$sex, prdat=prDat, vc=o, xin=xin,
   k=5, direction="back", verbose=TRUE)
mg$model$value # likelihood of the final model
}
}

