\name{glmm.EC}
\alias{glmm.EC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{function for testing a single/pooled variant for continuous traits 
with family data using generalized linear mixed effects model}
\description{
   Fit generalized linear mixed effects model (GLMM) with logistic link 
that treats each pedigree as a cluster to test a single/pooled variant for 
associations against a continuous phenotype with family data. The \code{glmer} 
function from package \code{lme4} is used. 
}
\usage{
glmm.EC(snp,phen,test.dat,covar,chr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snp}{ a numeric vector with genotype of a single/pooled variant  }
  \item{phen}{ a character string for the phenotype name of a binary trait of 
interest in \code{test.dat}  }
  \item{test.dat}{ the product of merging phenotype, genotype and pedigree data, 
should be ordered by "famid" }
  \item{covar}{ a character vector for covariates in \code{test.dat} }
  \item{chr}{ chromosome number }
}
\details{
   The \code{glmm.EC} function fits a generalized linear mixed effects model (GLMM) 
with logistic link that treats each pedigree as a cluster to test association between 
a binary trait and a single/pooled genetic variant with additive model. The 
trait-variant association test is carried out by the \code{glmer} function from package 
\code{lme4}. P-value from likelihood ratio test (LRT) is reported. This function is 
called in \code{glmm.ped} function to test all single/pooled variants.
}
\value{
  \item{ntotal }{number of individuals with genotype, phenotype and covariates}
  \item{nmiss }{number of individuals with missing genotype among \code{ntotal}}
  \item{maf_ntotal }{minor allele frequency based on \code{ntotal}}
  \item{beta }{regression coefficient of single SNP test or burden test}
  \item{se }{standard error of \code{beta}}
  \item{Z }{Z statistic based on signed LRT }
  \item{remark }{additional information of the analysis }
  \item{p }{LRT p-value of a single variant test or burden test}
  \item{MAC }{minor allele count}
  \item{n0 }{the number of individuals with 0 copy of coded alleles}
  \item{n1 }{the number of individuals with 1 copy of coded alleles}
  \item{n2 }{the number of individuals with 2 copies of coded alleles}
}
\references{ 
 Bates D, Maechler M, Bolker B and Walker S (2014). lme4: Linear mixed-effects 
models using Eigen and S4. R package version 1.1-7, 
http://CRAN.R-project.org/package=lme4. 
}

\author{ Ming-Huei Chen <mhchen@bu.edu> and Qiong Yang <qyang@bu.edu> }

\examples{
\dontrun{
glmm.EC(snp=rsnps.dat[,"snp1"],phen="disease",test.dat=rsnps.dat,
covar=c("sex"),chr=1)
}
}
