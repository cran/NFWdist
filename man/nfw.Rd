\name{nfw}
\alias{NFW}
\alias{NFWdist}
\alias{nfw}
\alias{dnfw}
\alias{pnfw}
\alias{qnfw}
\alias{rnfw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Standard Distribution Functions for the 3D NFW Profile
}
\description{
Density, distribution function, quantile function and random generation for the 3D NFW profile
}
\usage{
dnfw(x, con = 5, log = FALSE)
pnfw(q, con = 5, log.p = FALSE)
qnfw(p, con = 5, log.p = FALSE)
rnfw(n, con = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,q}{
Vector of quantiles. This is scaled such that \option{x} and \option{q} are equal to R/Rvir for NFW. This means the PDF is only defined between 0 and 1.
}
  \item{p}{
Vector of probabilities.
}
  \item{n}{
Number of observations. If length(n) > 1, the length is taken to be the number required.
}
  \item{con}{
The NFW profile concentration parameter, where c=Rvir/Rs.
}
  \item{log, log.p}{
Logical; if TRUE, probabilities/densities p are returned as log(p).
}
}
\details{
The novel part of this package is the general solution for the CDF inversion (i.e. qnfw). As far as I can see this has not been published anywhere, and it is a useful function for populating halos in something like an HOD.

One of \code{lamW} (fastest) or \code{gsl} (easier to install) must be installed to use the \code{qnfw} and \code{rnfw} functions!. Try to install \code{lamW} first (since it is about four times faster), but if that is tricky due to \code{Rcpp} dependencies then use \code{gsl} instead.
}
\value{
dnfw gives the density, pnfw gives the distribution function, qnfw gives the quantile function, and rnfw generates random deviates.
}
\references{
Robotham & Howlett, 2018, arXiv 1805.09550
}
\author{
Aaron Robotham
}
\note{
This seems to work at least as efficiently as accept reject, but it is ultimately much more elegant code in any case.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{lambert_W0} (\code{gsl}) or \code{lambertW0} (\code{lamW}).
}
\examples{
#Both the PDF (dnfw) integrated up to x, and CDF at q (pnfw) should be the same:
#0.373, 0.562, 0.644, 0.712

for(con in c(1,5,10,20)){
  print(integrate(dnfw, lower=0, upper=0.5, con=con)$value)
  print(pnfw(0.5, con=con))
}

#The qnfw should invert the pnfw, returning the input vector (1:9)/10:
for(con in c(1,5,10,20)){
  print(qnfw(p=pnfw(q=(1:9)/10,con=con), con=con))
}

#The sampling from rnfw should recreate the expected PDF from dnfw:

for(con in c(1,5,10,20)){
  plot(density(rnfw(1e6,con=con), bw=0.01))
  lines(seq(0,1,len=1e3), dnfw(seq(0,1,len=1e3),con=con),col='red')
  legend('topright',legend=paste('con =',con))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ nfw }% use one of  RShowDoc("KEYWORDS")
