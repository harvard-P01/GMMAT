\name{glmm.score.text}
\alias{glmm.score.text}
\title{
Performing GLMM based score tests
}
\description{
Use residuals and projection matrix from the null GLMM to perform score tests for association with genotypes in a plain text file (or compressed .gz or .bz2 file).
}
\usage{
glmm.score.text(res, P, infile, outfile, tol = 1e-05, center = T, 
	select = NULL, missing.method = "impute2mean", infile.nrow = NULL, 
	infile.nrow.skip = 0, infile.sep = "\t", infile.na = "NA", 
	infile.ncol.skip = 1, infile.ncol.print = 1, 
	infile.header.print = "SNP", nperbatch = 100)
}
\arguments{
  \item{res}{
rescaled residuals from the null GLMM. Note if you use \code{glmmkin( )} to fit the null model, rescaled residuals can be obtained by dividing raw residuals by the dispersion parameter.
}
  \item{P}{
projection matrix from the null GLMM. Direct output of \code{glmmkin( )}.
}
  \item{infile}{
the input file name. This is a plain text file with some delimiters (comma, space, tab or something else), with one row for each SNP and one column for each individual. SNPs should be coded as numeric values (0/1/2 or dosages allowed, A/C/G/T coding is not recognized). There can be additional rows and columns to skip at the beginning. The order of individuals can be different from \code{res} and \code{P} (see the argument \code{select}). Some compressed files (.gz and .bz2) also allowed.
}
  \item{outfile}{
the output file name.
}
  \item{tol}{
the threshold for determining monomorphism. If a SNP has value range less than the tolerance, it will be considered monomorphic and its association test p-value will be NA (default = 1e-5).
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
  \item{infile.nrow}{
number of rows to read in \code{infile}, including number of rows to skip at the beginning. If NULL, the program will determine how many rows there are in \code{infile} automatically and read all rows (default = NULL).
}
  \item{infile.nrow.skip}{
number of rows to skip at the beginning of \code{infile}. Must be nonnegative integers. Useful when header or comment lines are present (default = 0).
}
  \item{infile.sep}{
delimiter in \code{infile} (default = "\\t").
}
  \item{infile.na}{
symbol in \code{infile} to denote missing genotypes (default = "NA").
}
  \item{infile.ncol.skip}{
number of columns to skip before genotype data in \code{infile}. These columns can be SNP name, alleles and/or quality measures and should be placed at the beginning in each line. After skipping these columns, the program will read in genotype data and perform score tests. Must be nonnegative integers. It is recommended that SNP name should be included as the first column in \code{infile} and genotype data should start from the second column or later (default = 1).
}
  \item{infile.ncol.print}{
a vector indicating which column(s) in \code{infile} should be printed to the output directly. These columns can be SNP name, alleles and/or quality measures placed at the beginning in each line. Must be nonnegative integers, no greater than \code{infile.ncol.skip} and sorted numerically in ascending order. By default, it is assumed that the first column is SNP name and genotype data start from the second column, and SNP name should be carried over to the output (default = 1). Should be set to NULL if \code{infile.ncol.skip} is 0.
}
  \item{infile.header.print}{
a character vector indicating column name(s) of column(s) selected to print by \code{infile.ncol.print} (default = "SNP"). Should be set to NULL if \code{infile.ncol.skip} is 0.
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
\code{\link{glmmkin}}, \code{\link{glmm.wald.text}}, \code{\link{glmm.score.bed}}
}
\examples{
infile <- system.file("extdata", "geno.txt", package = "GMMAT")
outfile <- "glmm.score.text.testoutfile.txt"
data(example)
attach(example)
model0 <- glmmkin(disease ~ age + sex, data = pheno, kins = GRM, 
       family = binomial(link = "logit"))
glmm.score.text(model0$res/model0$theta[1], model0$P, infile = infile, 
	outfile = outfile, infile.nrow.skip = 5, infile.ncol.skip = 3, 
	infile.ncol.print = 1:3, infile.header.print = c("SNP", "Allele1", 
	"Allele2"))
}
\keyword{ generalized linear mixed model }
\keyword{ score test }
