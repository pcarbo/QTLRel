\name{qqPlot}
\alias{qqPlot}
\title{
Quantile-Quantile Plots
}
\description{
Quantile-Quantile Plots With the Ability to Draw Confidence Bands.
}
\usage{
qqPlot(y,x="norm", ...,
   xlab=if(is.numeric(x)) deparse(substitute(x)) else x,
   ylab=deparse(substitute(y)),main="Q-Q Plot",
   type="p",col=1,lty=2,lwd=1,pch=1,cex=0.7,plot.it=TRUE,
   confidence=.95,qqline=c("observed","expected","none"))
}
\arguments{
  \item{y}{
a numeric vector of data values.
}
  \item{x}{
either a numeric vector of data values, or a character string naming a distribution function such as "norm".
}
  \item{...}{
parameters passed to the distribution specified by \code{x} (if non-numerical).
}
  \item{xlab}{
a label for the x axis.
}
  \item{ylab}{
a label for the y axis.
}
  \item{main}{
a main title for the plot.
}
  \item{type}{
1-character string giving the type of plot desired.
}
  \item{col}{
color for points and lines.
}
  \item{lty}{
line type.
}
  \item{lwd}{
line width.
}
  \item{pch}{
plotting character for points.
}
  \item{cex}{
factor for expanding the size of plotted symbols.
}
  \item{plot.it}{
whether or not to draw a plot. if plotting, points outside the confidence bands will be indicated by different a color.
}
  \item{confidence}{
confidence level for the confidence band, or \code{FALSE} for no band.
}
  \item{qqline}{
whether or not to draw a reference line. if "observed", the line passes through the first and third observed quartiles; if "expected", the point (x,y) is expected to fall on the line if x and y follow the same distribution; if "none", no reference line is drawn.
}
}
\details{
If \code{x} is numeric, a two-sample test of the null hypothesis that \code{x} and \code{y} were drawn from the same continuous distribution is performed. Alternatively, \code{x} can be a character string naming a continuous distribution function. In such a case, a one-sample test is carried out of the null that \code{y} was draw from distribution \code{x} with parameters specified by "...".
}
\value{
   \item{x}{quantiles of \code{x}}
   \item{y}{quantiles of \code{y}}
   \item{lower, upper}{lower and upper limits if \code{confidence} is specified}
}
\references{
George Marsaglia, Wai Wan Tsang and Jingbo Wang (2003), Evaluating Kolmogorov's distribution. Journal of Statistical Software 8 (18): 1-4.

Vijayan N. Nair (1982). Q-Q plots with confidence bands for comparing several populations.

William J. Conover (1971). Practical Nonparametric Statistics. New York: John Wiley & Sons.
}

\seealso{
\code{\link{ks.test}}.
}
\examples{
\dontrun{
par(mfrow=c(1,2))
x<- rnorm(200, mean=0.7,sd=2); y<- rnorm(200, sd=2)
qqPlot(y,x,qqline="exp")
qqPlot(y=y,x="norm",sd=2)
ks.test(x,y)
}
}

