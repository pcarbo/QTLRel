\name{rem}
\alias{rem}
\title{
Random effect matrices
}
\description{
Construct matrices associated with random effects.
}
\usage{
rem(formula,data)
}
\arguments{
 \item{formula}{
a formula of the form: ~ Z | G1/.../Gk + ..., corresponding to random effects Z*G_i + Z*G_\{ij\} + ... in a mixed effect model. If Z=1 as in most cases, then it can be ~ G1/.../Gk + ...
}
  \item{data}{
a data frame that contains all the variables in the formula.
}
}
\value{
a list of matrices that are associated with random effects.
}

\examples{
\dontrun{
# make-up example
dat<- data.frame(
   group=c("A","A","A","A","A","A","B","B","B","B"),
   sex=c("F","F","F","M","M","M","F","F","M","M"),
   pass=c("Y","N","N","Y","Y","Y","Y","N","N","Y"),
   z=1:10)

# random effect pass, group and sex, where sex is nested
# within group:
# y_{ijk} = x_{ij}*b + group_i + sex_{ij} + z*pass_{ij}
#           + e_{ijk}
rem(~ group/sex + z|pass,data=dat)
}
}

