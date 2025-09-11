\name{MTest}
\alias{MTest}
\title{Bootstrap-based test for multicollinearity (Klein and VIF)}
\usage{
MTest(object, nboot = 100, nsam = NULL, trace = FALSE, seed = NULL,
      valor_vif = 0.9)
}
\description{
\code{MTest} implements a nonparametric (pairs) bootstrap to assess
multicollinearity by providing achieved significance levels (ASL)
for two widely used diagnostics: Klein's rule and the Variance Inflation
Factor (VIF). It returns bootstrap distributions of the global \eqn{R^2}
and the auxiliary \eqn{R_j^2} (from regressions of each predictor on the remaining
predictors), along with p-values for both rules.
}

\arguments{
  \item{object}{A fitted model, typically of class \code{"lm"}.
  The function uses its \code{model.frame} and \code{model.matrix} as the fixed design for resampling.}
  \item{nboot}{Integer. Number of bootstrap iterations (rows resampled with replacement).}
  \item{nsam}{Integer. Bootstrap sample size per iteration. Defaults to the original number of rows.}
  \item{trace}{Logical. If \code{TRUE}, shows a progress bar.}
  \item{seed}{Integer. Optional RNG seed for reproducibility.}
  \item{valor_vif}{Numeric in \eqn{[0,1)}. Threshold applied to \eqn{R_j^2} for the VIF rule:
  the ASL is \eqn{P(R_j^2 > c)}, where \eqn{c =} \code{valor_vif}. Since \eqn{\mathrm{VIF}_j = 1/(1-R_j^2)},
  \code{valor_vif = 0.9} corresponds roughly to \code{VIF >= 10}.}
}

\details{
\strong{Model.}
Consider the linear regression model
\deqn{Y_i = \beta_0 + \beta_1 X_{1i} + \cdots + \beta_p X_{pi} + u_i,\quad i=1,\ldots,n,}
and the \emph{auxiliary} regressions obtained by regressing each predictor \eqn{X_j}
on the remaining predictors \eqn{X_{-j}}. Let \eqn{R_g^2} be the global coefficient of
determination and \eqn{R_j^2} the coefficient of determination of the \eqn{j}-th auxiliary regression.

\strong{Diagnostics and achieved significance levels (ASL).}
\itemize{
  \item \emph{Klein's rule}: flag multicollinearity if \eqn{R_j^2 > R_g^2}.
        We estimate the ASL as \eqn{P(R_g^2 < R_j^2)} using the bootstrap distribution.
  \item \emph{VIF rule}: flag multicollinearity if \code{VIF} exceeds a threshold.
        Since \eqn{\mathrm{VIF}_j = 1/(1-R_j^2)}, this is equivalent to testing \eqn{R_j^2}
        against \code{valor_vif}. We estimate \eqn{P(R_j^2 >} \code{valor_vif}\eqn{)}.
}

\strong{Bootstrap scheme.}
The function resamples rows of the model frame (pairs bootstrap) and, for each bootstrap
sample, computes \eqn{R_g^2} and \eqn{R_j^2} (hence VIF) using the same expanded design
matrix as the original fit. This makes the procedure robust to transformed terms on either
side of the formula (e.g., \code{log(y)}, \code{I(X1^2)}, interactions, factors,
\code{poly()}, etc.).
}

\value{
An object of class \code{MTest}, which is a list containing:
\item{pval_vif}{Named numeric vector of ASL for the VIF rule,
\eqn{P(R_j^2 >} \code{valor_vif}\eqn{)}.}
\item{pval_klein}{Named numeric vector of ASL for Klein's rule,
\eqn{P(R_g^2 < R_j^2)}.}
\item{Bvals}{Numeric matrix of size \code{nboot x (p+1)} with columns
\code{"global"} for \eqn{R_g^2} and one column per predictor for \eqn{R_j^2}.}
\item{VIFvals}{Numeric matrix \code{nboot x p} with VIF values per predictor (one column per design column).}
\item{vif.tot}{Observed VIF per predictor from the original design.}
\item{R.tot}{Named numeric vector with observed \eqn{R_g^2} and all \eqn{R_j^2}.}
\item{nsam}{Bootstrap sample size actually used.}
\item{nboot}{Number of bootstrap iterations actually performed.}
}

\section{Interpretation}{
\itemize{
  \item Larger \code{pval_klein[j]} indicates stronger evidence that predictor \eqn{j}
  violates Klein's rule (\eqn{R_j^2} often exceeds \eqn{R_g^2}).
  \item Larger \code{pval_vif[j]} indicates that \eqn{R_j^2} frequently exceeds
  \code{valor_vif} (equivalently, \code{VIF} exceeds the implied threshold).
}
}

\section{Notes}{
For factor predictors, the underlying design includes multiple columns;
\code{VIFvals} and VIF-related summaries are returned \emph{per design column}.
In singular bootstrap samples some statistics may be \code{NA}.
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
\code{\link[car]{vif}} for classical VIF computation.
}

\examples{
## Minimal example (small nboot for speed)
set.seed(1)
data(simDataMTest, package = "MTest")
m1 <- stats::lm(y ~ ., data = simDataMTest)

boot.sol <- MTest(m1, nboot = 50, trace = FALSE, seed = 123, valor_vif = 0.90)

boot.sol$pval_vif
boot.sol$pval_klein
head(boot.sol$Bvals)
print(boot.sol)
}
\keyword{regression}
\keyword{multicollinearity}
\keyword{bootstrap}
