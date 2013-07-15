\name{aicVC}
\alias{aicVC}
\title{
AIC Model Selection
}
\description{
Select genetic variance components via Akaike's information criterion (AIC).
}
\usage{
aicVC(y,x,v=vector("list",6),initpar,k=2,init=6,keep=6,
   direction=c("forward","backward"),nit=25,verbose=FALSE,
   method=c("Nelder-Mead","BFGS","CG","SANN"),control=list(),
   hessian=FALSE)
}
\arguments{
  \item{y}{
a numeric vector or a numeric matrix of one column (representing a phenotype for instance).
}
  \item{x}{
a data frame or matrix, representing covariates if not missing.
}
  \item{v}{
a list of variance components (AA, DD, HH, AD, MH, EE, ...) where "AA" and "DD" are respectively additive and dominance genetic matrices, "HH", "AD" and "MH" are other genetic matrices that one may be interested in (see "details"), "EE" is the residual matrix that is usually assumed to be an identity matrix, and "..." are other random components of interest. If a genetic component is not considered, it should be set to NULL.
}
  \item{initpar}{
optional initial parameter values.
}
  \item{k}{
penalty on a parameter. The selection criterion is the known "AIC" if \code{k = 2} and is "BIC" if \code{k = log(n)} where "n" is the sample size.
}
  \item{init}{
indicator of the initial model pertaining to genetic variance components. For instance, c(1,2,6) indicates the initial model includes "AA", "DD" and "EE". By default, only "EE" is considered.
}
 \item{keep}{
indicator of which variance components should be forced into the final model. The default is "EE".
}
  \item{direction}{
the mode of search. Either "forward" or "backward" with default "forward".
}
  \item{nit}{
number of iterations to call \code{\link{optim}} for optimization.
}
  \item{verbose}{
a logical variable. True if ones wants to track the process for monitoring purpose.
}
  \item{method}{
the optimization method to be used. See \code{\link{optim}} for details.
}
  \item{control}{
a list of control parameters to be passed to \code{\link{optim}}.
}
  \item{hessian}{
logical. Should a numerically differentiated Hessian matrix be returned?
}
}
\details{
In genome-wide association studies (GWAS), random effects are usually added to a model to account for polygenic variation. Abney et al (2000) showed that five variance components including the most interesting additive and dominance variance components are potentially induced by polygenes. The above function is intended for selecting variance components that contribute "most" to a quantitative trait.

Function \code{\link{estVC}} is called by the above function to estimate the parameters and maximum likelihood in each model. Refer to \code{\link{estVC}} for more information.
}
\value{
\item{aic}{AIC of the final model.}
\item{model}{gives parameter estimates, log-likihood, and other information.}
\item{lik}{log-likelihood of the model selected at each intermediate step.}
\item{trace}{indicates which variance components were selected at each intermediate step.}
}
\references{
Abney, M., M. S. McPeek, and C. Ober (2000). Estimation of variance components of quantitative traits in inbred populations. Am. J. Hum. Genet. 141, 629-650.
}
\seealso{
\code{\link{estVC}} for more information.
}
\examples{
data(miscEx)

\dontrun{
# forward selection
# any variance component will be selected
# if AIC improve by 1e-5 or larger
v=list(AA=gmF8$AA, DD=gmF8$DD, HH=NULL, AD=NULL,
	MH=NULL, EE=diag(length(pdatF8$bwt)))
o<- aicVC(y=pdatF8$bwt, x=pdatF8$sex, k=0, v=v, verbose=TRUE)
o

# forward selection
of<- aicVC(y=pdatF8$bwt, x=pdatF8$sex, v=v, k=1/2, init=6,
	direction="for", verbose=TRUE)
of

# backward elimination
ob<- aicVC(y=pdatF8$bwt, x=pdatF8$sex, v=v, k=1/2, keep=6,
	direction="back", verbose=TRUE)
ob
}
}

