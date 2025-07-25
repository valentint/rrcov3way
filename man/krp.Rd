% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mtx-utils.R
\name{krp}
\alias{krp}
\title{The Khatri-Rao product of two matrices}
\usage{
krp(A, B)
}
\arguments{
\item{A}{Matrix of order I x K.}

\item{B}{Matrix of order J x K.}
}
\value{
The IJ x K matrix of columnwise Kronecker products.
}
\description{
The function \code{krp(A,B)} returns the Khatri-Rao product of two matrices \code{A} and \code{B}, of
 dimensions I x K and J x K respectively. The result is an IJ x K matrix formed by the matching
 column-wise Kronecker products, i.e. the k-th column of the Khatri-Rao product is
 defined as \code{kronecker(A[, k], B[, k])}.
}
\examples{
a <- matrix(1:12, 3, 4)
b <- diag(1:4)
krp(a, b)
krp(b, a)
}
\references{
Khatri, C. G., and Rao, C. Radhakrishna (1968).
     Solutions to Some Functional Equations and Their Applications to Characterization of Probability Distributions.
     Sankhya: Indian J. Statistics, Series A 30, 167-180.

 	Smilde, A., Bro R. and Gelardi, P. (2004). Multi-way Analysis: Applications in Chemical Sciences, Chichester:Wiley
}
\author{
Valentin Todorov, \email{<valentin@todorov.at>}
}
