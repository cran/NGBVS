\name{Random value generation from the Generalized Inverse Gaussian Distribution}
\alias{rgig}
\title{Random value generation from the Generalized Inverse Gaussian Distribution}
\description{
Random value generation from the Generalized Inverse Gaussian (GIG) Distribution.
}
\usage{
rgig(n = 10, lambda = 1, chi = 1, psi = 1)
}
\arguments{
\item{n}{Number of observations.}
\item{lambda}{A shape and scale and parameter.}
\item{chi}{Shape parameter. Must be positive.}
\item{psi}{Scale parameter. Must be positive.}
}
 
\details{
\code{rgig} uses the code from the GIG-random number generator from
the \R package \emph{fBasics}. I copied the code from the "ghyp" package because it had not longer a maintainer.
}

\value{
A vector with random values from the GIG distrigution.
}

\references{

  The algorithm for simulating generalized inverse gaussian variates is
  copied from the \R package \emph{fBasics} from Diethelm Wuertz.

  Dagpunar, J.S. (1989). \emph{An easily implemented generalised inverse
  Gaussian generator.} Commun. Statist. -Simula., \bold{18}, 703--710.

  Raible S. (2000).
  \emph{Levy Processes in Finance: Theory, Numerics and Empirical Facts},
  PhD Thesis, University of Freiburg, Germany, 161 pages.

}

\author{David Luethi. Minor changes made by Abdulaziz Alenazi \email{a.alenazi@nbu.edu.sa}.}

\seealso{ \code{\link{integrate}},  \code{\link{uniroot}}, \code{\link{spline}} }

\examples{
x <- rgig(n = 10, lambda = 1, chi = 1, psi = 1)
}
