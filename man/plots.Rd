\name{plotit}
\alias{plot}
\alias{plotit}
\alias{plot.scanOne}
\alias{plot.scanTwo}
\title{
Plotting
}
\description{
Plot mapping results.
}
\usage{
\method{plot}{scanOne}(x,...)

plotit(lrt,cv,bychr=FALSE,chr.labels=TRUE,type="p",lty=NULL,
   col=NULL,pch=NULL,cex=NULL,...)
}
\arguments{
  \item{x}{
object from \code{\link{scanOne}} or \code{\link{scanTwo}}.
}
  \item{lrt}{
a data frame with (chr, dist, y,...) or (chr, dist, y, group,...), where "chr" represents chromosome, "dist" position on the chromosome, "y" the test statistic.
}
  \item{cv}{
threshold to be drawn on the plot.
}
  \item{cex}{
see \code{\link{par}}.
}
  \item{bychr}{
a logical variable. If true, the plot will be displayed per chromosomes.
}
  \item{chr.labels}{
a logical variable. If true, the chromosome names will be drawn.
}
  \item{type,lty,col,pch}{
See \code{\link{plot.default}}.
}
  \item{\dots}{
other options passed to R plot function. To call \code{\link{plot}} to plot results of \code{\link{scanOne}}, one may need to provide a genetic map \code{gmap} that should be data frame (snp, chr, dist, ...), where "snp" is the SNP (marker) name, "chr" is the chromosome where the "snp" is, and "dist" is the genetic distance in centi-Morgan (cM) from the leftmost SNP (marker) on the chromosome.

}
}
\note{
A genetic map 'gmap' may be needed to plot an object of \code{\link{scanOne}} or \code{\link{scanTwo}}. The color option may not give what is expected.
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

# genome scan
llk<- scanOne(y=pdatF8$bwt, x=pdatF8$sex, vc=o, gdat=gdat.imp)

# plotting
plot(llk, gmap=gmapF8) # gmap is needed

# plotting in another way
idx<- match(colnames(gdat.imp), gmapF8$snp)
tmp<- data.frame(chr=gmapF8$chr[idx],dist=gmapF8$dist[idx],y=llk$p)
plotit(tmp, main="Mapping Plot", xlab="Chromosome", ylab="LRT",
   col=as.integer(tmp$ch)\%\%2+2,type="p")
}
}

