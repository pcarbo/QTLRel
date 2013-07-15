\name{kinship}
\alias{kinship}
\title{
Calculate kinship coefficients
}
\description{
Calculate kinship coefficients from a pedigree.
}
\usage{
kinship(ped, ids)
}
\arguments{
  \item{ped}{
a pedigree, which a data frame (id, sire, dam, ...). If given, "generation" can be numeric 0, 1, 2, ... or non-numeric "F0", "F1", "F2", ... If "sex" is included, male should be "M", "Male" or 1, and female should be "F", "Female" or 2 (other than 0 and 1). Note that 0 is reserved for missing values.
}
  \item{ids}{
IDs of the individuals. If given, kinship coefficients are extracted for individuals with ID \code{ids}; otherwise, kinship coefficients are provided for all individuals in the pedigree.
}
}
\value{
A matrix giving kinship coefficients.
}
\examples{
data(miscEx)

ids<- sample(pedF8$id,10)
\dontrun{
ksp<- kinship(pedF8,ids=ids)
}
}
