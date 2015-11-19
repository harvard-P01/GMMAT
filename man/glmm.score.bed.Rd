\name{glmm.score.bed}
\alias{glmm.score.bed}
\title{
Performing GLMM based score tests using plink binary genotypes
}
\description{
Use residuals and projection matrix from the null GLMM to perform score tests for association with genotypes in a plink bed file (binary genotypes).
}
\usage{
glmm.score.bed(res, P, plinkfiles, outfile, center = T, select = NULL, 
	missing.method = "impute2mean", nperbatch = 100)
}
\arguments{
  \item{res}{
rescaled residuals from the null GLMM. Note if you use \code{glmmkin( )} to fit the null model, rescaled residuals can be obtained by dividing raw residuals by the dispersion parameter.
}
  \item{P}{
projection matrix from the null GLMM. Direct output of \code{glmmkin( )}.
}
  \item{plinkfiles}{
the plink binary input file name. Note that only the prefix without .bed, .bim or .fam should be used. Only SNP major mode recognized in the binary file.
}
  \item{outfile}{
the output file name.
}
  \item{center}{
a logical switch for centering genotypes before tests. If TRUE, genotypes will be centered to have mean 0 before tests, otherwise raw values will be directly used in tests (default = TRUE).
}
  \item{select}{
an optional vector indicating the order of individuals in \code{infile}. If supplied, the length must match the number of individuals in \code{infile} (default = NULL). Individuals to be excluded should be coded 0. For example, \code{select = c(2, 3, 1, 0)} means the 1st individual in \code{infile} corresponds to the 2nd individual in \code{res} and \code{P}, the 2nd individual in \code{infile} corresponds to the 3rd individual in \code{res} and \code{P}, the 3rd individual in \code{infile} corresponds to the 1st individual in \code{res} and \code{P}, the 4th individual in \code{infile} is not included in \code{res} or \code{P}.
}
  \item{missing.method}{
method of handling missing genotypes. Either "impute2mean" or "omit" (default = "impute2mean").
}
  \item{nperbatch}{
an integer for how many SNPs should be tested in a batch (default = 100). The computational time can increase dramatically if this value is either small or large. The optimal value for best performance depends on the user's system.
}
}
\value{
computational time in seconds, excluding I/O time.
}
\references{
Chen, H., Wang, C., Conomos, M.P., Stilp, A.M., Li, Z., Sofer, T., Szpiro, A.A., Chen, W., Brehm, J.M., Celedón, J.C., Redline, S.S., Papanicolaou, G.J., Thornton, T.A., Laurie, C.C., Rice, K. and Lin, X. Control for population structure and relatedness for binary traits in genetic association studies using logistic mixed models. (Submitted)
}
\author{
Han Chen
}
\seealso{
\code{\link{glmmkin}}, \code{\link{glmm.wald.bed}}, \code{\link{glmm.score.text}}
}
\examples{
plinkfiles <- strsplit(system.file("extdata", "geno.bed", package = "GMMAT"), 
       ".bed", fixed = TRUE)[[1]]
outfile <- "glmm.score.bed.testoutfile.txt"
data(example)
attach(example)
model0 <- glmmkin(disease ~ age + sex, data = pheno, kins = GRM, 
       family = binomial(link = "logit"))
glmm.score.bed(model0$res/model0$theta[1], model0$P, plinkfiles = plinkfiles, 
	outfile = outfile)
}
\keyword{ generalized linear mixed model }
\keyword{ score test }
