\name{lme.EC}
\alias{lme.EC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{function for testing a single/pooled variant for continuous traits 
with family data using Linear Mixed Effects model}
\description{
   Fit linear mixed effects (LME) model to test a single/pooled variant 
for associations against a continuous phenotype with family data. The 
\code{lmekin} function from package \code{coxme} is used. 
}
\usage{
lme.EC(snp,phen,test.dat,covar,kmat,chr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snp}{ a numeric vector with genotype of a single/pooled variant }
  \item{phen}{ a character string for the phenotype name of a continuous 
trait of interest in \code{test.dat}  }
  \item{test.dat}{ the product of merging phenotype, genotype and pedigree data }
  \item{covar}{ a character vector for covariates in \code{test.dat} }
  \item{kmat}{ relationship coefficient (twice of kinship coefficient) 
matrix based on pedigree file }
  \item{chr}{ chromosome number }
}
\details{
   The \code{lme.EC} function fits a Linear Mixed Effects model (LME) that uses 
relationship coefficient matrix as within pedigree correlation matrix to test 
association between a continuous phenotype and a single/pooled genetic variant 
with additive model. The trait-SNP association test is carried out by the 
\code{lmekin} function from package \code{coxme}. This function is called in 
\code{lme.ped} function to test all single/pooled variants.
}
\value{
  \item{ntotal }{number of individuals with genotype, phenotype and covariates}
  \item{nmiss }{number of individuals with missing genotype among \code{ntotal}}
  \item{maf_ntotal }{minor allele frequency based on \code{ntotal}}
  \item{beta }{regression coefficient of single SNP test or burden test}
  \item{se }{standard error of \code{beta}}
  \item{Z }{Wald Z statistic}
  \item{remark }{additional information of the analysis }
  \item{p }{p-value of single SNP test or burden test}
  \item{MAC }{minor allele count}
  \item{n0 }{the number of individuals with 0 copy of coded alleles}
  \item{n1 }{the number of individuals with 1 copy of coded alleles}
  \item{n2 }{the number of individuals with 2 copies of coded alleles}
}
\references{ 
 coxme package: mixed-effects Cox models, sparse matrices, and modeling data from 
large pedigrees. Beth Atkinson (atkinson@mayo.edu) for pedigree functions.Terry 
Therneau (therneau@mayo.edu) for all other functions. 2007. Ref Type: Computer 
Program. http://cran.r-project.org/web/packages/coxme/. 

 Abecasis, G. R., Cardon, L. R., Cookson, W. O., Sham, P. C., & Cherny, S. S (2001). 
Association analysis in a variance components framework. \emph{Genet Epidemiol}, 
\bold{21} Suppl 1, S341-S346.
}

\author{ Ming-Huei Chen <mhchen@bu.edu> and Qiong Yang <qyang@bu.edu> }

\examples{
\dontrun{
lme.EC(snp=rsnps.dat$counts,phen="trait",test.dat=rsnps.dat,
covar=c("age","sex"),kmat=kmat,chr=1)}
}
