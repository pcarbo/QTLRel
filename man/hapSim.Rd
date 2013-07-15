\name{hapSim}
\alias{hapSim}
\title{
Generate Genotypic Data
}
\description{
Simulate gametic data from a pedigree.
}
\usage{
hapSim(ped,gmap,ids,hap,method=c("Haldane","Kosambi"),
   recode.pedigree=FALSE)
}
\arguments{
  \item{ped}{
a pedigree, which is a data frame (id, sex, sire, dam, ...). In "sex", male should be "M", "Male" or 1, and female should be "F", "Female" or 2 (other than 0 and 1). If given, "generation" can be numeric 0, 1, 2, ... or non-numeric "F0", "F1", "F2", ... Note that 0 is reserved for missing values.
}
  \item{gmap}{
a genetic map. Should be data frame (snp, chr, dist, ...), where "snp" is the SNP (marker) name, "chr" is the chromosome where the "snp" is, and "dist" is the genetic distance in centi-Morgan (cM) from the leftmost SNP (marker) on the chromosome. If \code{gmap} is missing but \code{hap} not, all but the first two columns of \code{hap} are ignored.
}
  \item{ids}{
genotypic data are extracted only for individuals with IDs specified by \code{ids}. If missing, genotypic data are extracted for all individuals in the pedigree. If \code{ped} is an object of \code{\link{pedRecode}}, \code{ids} should be referred to "old" IDs.
}
  \item{hap}{
founders' haplotype data if not missing. Rows correspond to all founders, which should be in the first places in the pedigree \code{ped}, in the exact order and columns correspond to loci in the genetic map \code{gmap} in the exact order. For an individual, the haplotype should be (f1 m1 f2 m2 ...) where fi is the allele from father at the i-th locus and mi is the allele from mother at the i-th locus. Elements should be non-negative integers that are not larger than 16384. If missing, two founders with alleles 1 and 2 are assumed.
}
  \item{method}{
whether "Haldane" or "Kosambi" mapping function should be used. This will be ignored if the recombination rate \code{recRate} is a component of \code{gmap}.
}
  \item{recode.pedigree}{
a logical variable. True if the pedigree needs to be recoded.
}
}
\details{
The pedigree should be in the same format as an output of \code{\link{pedRecode}}.
}
\value{
a matrix giving haplotypes.
}
\seealso{
\code{\link{pedRecode}} for more information.
}
\examples{
data(miscEx)

\dontrun{
# prepare pedigree in desired format
pedR<- pedRecode(pedF8)
# fake founder haplotypes
hapDat<- rbind(rep(1:2,nrow(gmapF8)),rep(3:4,nrow(gmapF8)))
# simulate hyplotypes for F8 individuals
hd<- hapSim(pedR, gmapF8, ids=pedR$id[pedR$gen=="F8"],
   hap=hapDat, recode.pedigree=TRUE)
dim(hd)
hd[1:5,1:10]
}
}

