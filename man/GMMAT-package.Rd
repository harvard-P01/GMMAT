\name{GMMAT-package}
\alias{GMMAT-package}
\alias{GMMAT}
\docType{package}
\title{
Generalized linear Mixed Model Association Tests
}
\description{
GMMAT is an R package for performing association tests using generalized linear mixed models (GLMMs) in genome-wide association studies (GWAS). GMMAT first fits a GLMM with covariate adjustment and random effects to account for population structure and family or cryptic relatedness, and then performs score tests for each genetic variant in the GWAS. For candidate gene studies, GMMAT can also perform Wald tests to get the effect size estimate for each genetic variant.
}
\details{
\tabular{ll}{
Package: \tab GMMAT\cr
Type: \tab Package\cr
Version: \tab 0.7\cr
Date: \tab 2015-11-18\cr
License: \tab GPL-3\cr
}
}
\author{
Han Chen, Matthew P. Conomos

Maintainer: Han Chen <hanchen@hsph.harvard.edu>
}
\references{
Brent, R.P. (1973) "Chapter 4: An Algorithm with Guaranteed Convergence for Finding a Zero of a Function", Algorithms for Minimization without Derivatives, Englewood Cliffs, NJ: Prentice-Hall, ISBN 0-13-022335-2.

Breslow, N.E. and Clayton, D.G. (1993) Approximate Inference in Generalized Linear Mixed Models. Journal of the American Statistical Association 88, 9-25.

Chen, H., Wang, C., Conomos, M.P., Stilp, A.M., Li, Z., Sofer, T., Szpiro, A.A., Chen, W., Brehm, J.M., Celedón, J.C., Redline, S.S., Papanicolaou, G.J., Thornton, T.A., Laurie, C.C., Rice, K. and Lin, X. Control for population structure and relatedness for binary traits in genetic association studies using logistic mixed models. (Submitted)

Gilmour, A.R., Thompson, R. and Cullis, B.R. (1995) Average Information REML: An Efficient Algorithm for Variance Parameter Estimation in Linear Mixed Models. Biometrics 51, 1440-1450.

Nelder, J.A. and Mead, R. (1965) A simplex algorithm for function minimization. Computer Journal 7, 308-313.

Yang, J., Lee, S.H., Goddard, M.E. and Visscher, P.M. (2011) GCTA: A Tool for Genome-wide Complex Trait Analysis. The American Journal of Human Genetics 88, 76-82.

Zhou, X. and Stephens, M. (2012) Genome-wide efficient mixed-model analysis for association studies. Nature Genetics 44, 821-824.
}
\keyword{ package }
\keyword{ generalized linear mixed model }


