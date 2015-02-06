\name{coxph.EC}
\alias{coxph.EC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{function for testing a single/pooled variant for survival traits with family 
data using Cox proportional hazards regression model}
\description{
  Fit Cox proportional hazards regression model to test a single/pooled variant for 
associations against a survival phenotype with family data. The \code{coxph} function 
from package \code{survival} is used. 
}
\usage{
coxph.EC(snp,phen,test.dat,covar,chr,time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snp}{ a numeric vector with genotype of a single/pooled variant  }
  \item{phen}{ a character string for the phenotype name of a continuous trait of 
interest in \code{test.dat}  }
  \item{test.dat}{ the product of merging phenotype, genotype and pedigree data }
  \item{covar}{ a character vector for covariates in \code{test.dat} }
  \item{chr}{ chromosome number }
  \item{time}{ the character string of variable named for survival time }
}
\details{
   The \code{coxph.EC} function fits a Cox proportional hazards regression model with 
shared frailty (random effect) in each pedigree to test association between a survival 
phenotype and a single/pooled genetic variant with additive model. The trait-SNP 
association test is carried out by the \code{coxph} function from package \code{survival}. 
P-value from likelihood ratio test (LRT) is reported. This function is called in 
\code{coxph.ped} function to test all single/pooled variants.
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
  Therneau T (2014). A Package for Survival Analysis in S. R package version 2.37-7, 
http://CRAN.R-project.org/package=survival. 

  Terry M. Therneau and Patricia M. Grambsch (2000). Modeling Survival Data: Extending 
the Cox Model. Springer, New York. ISBN 0-387-98784-3.
}

\author{ Ming-Huei Chen <mhchen@bu.edu> and Qiong Yang <qyang@bu.edu> }

\examples{
\dontrun{
coxph.EC(snp=rsnps.dat[,"snp1"],snp1,phen="trait2",test.dat=rsnps.dat,
covar=c("age","sex"),chr=1,time="survival_time")}
}
