\name{gls}
\alias{gls}
\title{
Generalized Least Squares Estimates
}
\description{
Obtain estimates using generalized least squares (gls).
}
\usage{
gls(formula,data,vc=NULL)
}
\arguments{
  \item{formula}{
an object of class "formula": a symbolic description of the model to be fitted.
}
  \item{data}{
an data frame containing the variables in the model.
}
  \item{vc}{
an object from \code{\link{estVC}} or \code{\link{aicVC}} or an estimated variance-covariance matrix induced by relatedness and environment if not NULL.
}
}
\value{
A matrix with columns: "Estimate", "Std. Error", "t value" and "Pr(>|t|)".
}

\seealso{
\code{\link{lm}}.
}

