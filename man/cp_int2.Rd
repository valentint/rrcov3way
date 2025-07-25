% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cp_int2.R
\name{cp_int2}
\alias{cp_int2}
\title{ATLD-ALS algorithm for Candecomp/Parafac (CP)}
\usage{
cp_int2(
  X,
  n,
  m,
  p,
  ncomp,
  initconv = 0.01,
  conv = 1e-06,
  const = "none",
  start = "random",
  maxit = 5000,
  trace = FALSE
)
}
\arguments{
\item{X}{A three-way array or a matrix. If \code{X} is a matrix
(matricised threeway array), \code{n}, \code{m} and \code{p} must be
given and are the number of A-, B- and C-mode entities respectively}

\item{n}{Number of A-mode entities}

\item{m}{Number of B-mode entities}

\item{p}{Number of C-mode entities}

\item{ncomp}{Number of components to extract}

\item{initconv}{Convergence criterion for the initialization phase (ATLD),
default is \code{conv=1e-2}.}

\item{conv}{Convergence criterion, default is \code{conv=1e-6}.}

\item{const}{Optional constraints for each mode. Can be a three element
character vector or a single character, one of \code{"none"} for no
constraints (default), \code{"orth"} for orthogonality constraints or
\code{"zerocor"} for zero correlation between the extracted factors.
For example, \code{const="orth"} means orthogonality constraints for
all modes, while \code{const=c("orth", "none", "none")} sets the
orthogonality constraint only for mode A.}

\item{start}{Initial values for the A, B and C components. Can be
\code{"svd"} for starting point of the algorithm from SVD's,
\code{"random"} for random starting point (orthonormalized
component matrices or nonnegative matrices in case of nonnegativity
constraint), or a list containing user specified components.}

\item{maxit}{Maximum number of iterations, default is \code{maxit=10000}.}

\item{trace}{Logical, provide trace output.}
}
\value{
The result of the decomposition as a list with the following
 elements:
   \itemize{
   \item \code{fit} Value of the loss function
   \item \code{fp} Fit value expressed as a percentage
   \item \code{ss} Sum of squares
   \item \code{A} Component matrix for the A-mode
   \item \code{B} Component matrix for the B-mode
   \item \code{C} Component matrix for the C-mode
   \item \code{iter} Number of iterations
   \item \code{tripcos} Minimal triple cosine between two components
       across the three component matrices, used to inspect degeneracy
   \item \code{mintripcos} Minimal triple cosine during the iterative
       algorithm observed at every 10 iterations, used to inspect degeneracy
   \item \code{ftiter} Matrix containing in each row the function value
       and the minimal triple cosine at every 10 iterations
   \item \code{const} Optional constraints (same as the input parameter
       \code{const})
   }
}
\description{
Integrated algorithm combining ATLD and ALS for the
 minimization of the Candecomp/Parafac (CP) loss function.
}
\note{
The argument \code{const} should be a three element character vector.
 Set \code{const[j]="none"} for unconstrained update in j-th mode weight
 matrix (the default),
 \code{const[j]="orth"} for orthogonal update in j-th mode weight matrix or
 \code{const[j]="zerocor"} for zero correlation between the extracted
     factors.
 The default is unconstrained update for all modes.

 The loss function to be minimized is \eqn{sum(k)|| X(k) - A D(k) B' ||^2},
 where \eqn{D(k)} is a diagonal matrix holding the \code{k}-th row of
 \code{C}.
}
\examples{

\dontrun{
## Example with the OECD data
 data(elind)
 dim(elind)

 res <- cp_int2(elind, ncomp=3)
 res$fp
 res$fp
 res$iter

 res <- cp_int2(elind, ncomp=3, const="orth")
 res$A
}
}
\references{
H.-L. Wu, M. Shibukawa, K. Oguma, An alternating trilinear decomposition
 algorithm with application to calibration of HPLC-DAD for
 simultaneous determination of overlapped chlorinated aromatic hydrocarbons,
 \emph{Journal of Chemometrics} \bold{12} (1998) 1--26.

Simonacci, V. and Gallo, M. (2020). An ATLD--ALS method for the trilinear decomposition
 of large third-order tensors, \emph{Soft Computing} 24 13535--13546.

Todorov, V. and Simonacci, V. and Gallo, M. and Trendafilov, N. (2023). A novel 
 estimation procedure for robust CANDECOMP/PARAFAC model fitting. 
 \emph{Econometrics and Statistics}. In press.
}
\author{
Valentin Todorov, \email{<valentin@todorov.at>}; Violetta Simonacci, \email{violetta.simonacci@unina.it}
}
