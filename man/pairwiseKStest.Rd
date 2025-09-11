\name{pairwiseKStest}
\alias{pairwiseKStest}
\title{Pairwise Kolmogorov--Smirnov tests for matrix columns}
\usage{
pairwiseKStest(X,
               alternative = c("greater","less","two.sided"),
               use = c("asis","pairwise.complete.obs"),
               exact = NULL)
}
\description{
Computes pairwise Kolmogorov--Smirnov (KS) tests between all columns of a
numeric matrix or data frame, returning the matrix of p-values. Typical inputs
include the \code{Bvals} matrix from \code{MTest}. For one-sided alternatives
(\code{"greater"} or \code{"less"}), the p-value matrix is directional:
rows correspond to \code{x} and columns to \code{y}.
}

\arguments{
  \item{X}{Numeric matrix or data frame. Columns are compared pairwise by KS tests.
  A common use is \code{Bvals} from \code{MTest}. If you want to compare only predictors,
  pass \code{Bvals[ , -1]} to exclude the \code{"global"} column.}
  \item{alternative}{Character string: \code{"greater"}, \code{"less"}, or \code{"two.sided"}.
  See \code{\link[stats]{ks.test}} for the meaning of these options.}
  \item{use}{Character string: \code{"asis"} (default; no NA filtering, replicates the original
  behavior), or \code{"pairwise.complete.obs"} (remove NAs pairwise for each test).}
  \item{exact}{Logical or \code{NULL}. Forwarded to \code{\link[stats]{ks.test}}.
  Leave as \code{NULL} to use the default decision of \code{ks.test} (as in the original function).}
}

\details{
The function performs a KS test for each ordered pair of columns \code{(i, j)}
using \code{ks.test(X[, i], X[, j], alternative = alternative, exact = exact)}.
For one-sided alternatives, the result is not symmetric, since rows play the
role of \code{x} and columns the role of \code{y}.

The returned \code{Suggestion} follows the same rule as the original function:
for \code{alternative = "greater"}, it sorts the row sums of the p-value matrix
(descending); for \code{"less"}, it sorts the column sums; for \code{"two.sided"},
no suggestion is returned.
}

\value{
A list of class \code{pairwiseKStest} with components:
  \item{KSpwMatrix}{Numeric matrix of p-values. Rows are \code{x}, columns are \code{y}.}
  \item{alternative}{Character string describing the alternative hypothesis used.}
  \item{Suggestion}{For \code{"greater"}: row sums of \code{KSpwMatrix} (sorted decreasing).
  For \code{"less"}: column sums (sorted decreasing). For \code{"two.sided"}: a message
  indicating no suggestions.}
}

\author{
Víctor Morales Oñate \email{vmorales.ppb@gmail.com} \cr
Bolívar Morales Oñate \email{bmoralesonate@gmail.com} \cr
\url{https://sites.google.com/site/moralesonatevictor/} \cr
\url{https://www.linkedin.com/in/vmoralesonate/}
}

\references{
Morales-Oñate, V., and Morales-Oñate, B. (2023).
\emph{MTest: a Bootstrap Test for Multicollinearity}. Revista Politécnica, 51(2), 53--62.
\doi{10.33333/rp.vol51n2.05}
}

\seealso{
\code{\link[stats]{ks.test}}, \code{\link{MTest}}
}

\examples{
## Typical workflow with MTest:
## (use small nboot for speed in examples)
set.seed(1)
data(simDataMTest, package = "MTest")
m1 <- stats::lm(y ~ ., data = simDataMTest)
boot.sol <- MTest(m1, nboot = 30, trace = FALSE, seed = 123)

## Compare only predictors (exclude "global"):
ks_res_greater <- pairwiseKStest(boot.sol$Bvals[, -1],
                                 alternative = "greater",
                                 use = "asis",    # same behavior as the original
                                 exact = NULL)    # let ks.test decide

ks_res_greater$KSpwMatrix
ks_res_greater$Suggestion

## Two-sided (no suggestion by design):
ks_res_twosided <- pairwiseKStest(boot.sol$Bvals[, -1],
                                  alternative = "two.sided")
ks_res_twosided$KSpwMatrix
}
\keyword{multicollinearity}
\keyword{nonparametric}
\keyword{ks-test}
