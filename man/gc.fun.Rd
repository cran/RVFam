\name{gc.fun}
\alias{gc.fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{function that does genomic control correction to single SNP analysis, sum of 
square test and RData for survival trait analysis}
\description{
   When high genomic control (GC) parameter (lambda) estimate is observed, \code{gc.fun} 
applies GC correction to SNPs with minor allele counts (\code{MAC}) less than a user specified
threshold that may have inflated type I error rate for survival traits in particular, adjusts 
RData output accordingly, and recomputes sum of square statistic. 
}
\usage{
gc.fun(path,phen,snpinfoRdata,snp.cor,mac,aggregateBy="SKATgene",
maf.file,mafRange,ssq.beta.wts=c(1,25))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{ path to directory that saves all 23 tab delimited single SNP analysis result files }
  \item{phen}{ a character string for the phenotype name of a trait of interest  }
  \item{snpinfoRdata}{ a character string naming the RData containing SNP info to be loaded, this 
should at least include 'Name' (for SNP name), 'Chr', and aggregateBy (default='SKATgene') columns }
  \item{snp.cor}{ a character string naming the RData containing lists of SNP correlation matrix 
within each 'SKATgene' }
  \item{mac}{ user specified MAC threshold for applying GC correction to SNPs with MAC under it}
  \item{aggregateBy}{ the column of SNP info on which single SNPs are to be aggregated for burden 
tests, default is 'SKATgene' }
  \item{maf.file}{ a character string naming the comma delimited file containing 'snp.names' for SNP name 
and 'maf' for MAF }
  \item{mafRange}{ range of MAF to include SNPs for gene-based burden tests, default is c(0,0.05) }
  \item{ssq.beta.wts}{ a vector of parameters of beta weights used in proposed sum of squares test, 
default=c(1,25) as in SKAT }
}
\details{
   When high lambda is observed from survival trait single SNP analysis, the \code{gc.fun} function 
applies GC correction to SNPs with user defined MAC, adjusts RData output based on GC 
corrected single SNP analysis results, recomputes sum of squares statistic and then outputs 
corrected single SNP analysis results, SSQ analysis results and RData. Initial single SNP analysis 
result files are required and the input arguments should be identical to the ones used in initial 
analysis (except for \code{path}). 
}
\value{
  No value is returned. Instead, tab delimited result files and an RData are generated. A single 
SNP result file, named with \code{phen} and singleSNP, contains columns: \code{gene}, \code{Name}, 
\code{maf}, \code{ntotal}, \code{nmiss}, \code{maf_ntotal}, \code{beta}, \code{se}, \code{Z}, 
\code{remark}, \code{p} (p-value from LRT), \code{MAC}, \code{n0}, \code{n1}, and \code{n2}. 
A SSQ test result file, named with \code{phen} and SSQ, contains columns: \code{gene}, \code{SSQ}, 
\code{cmafTotal}, \code{cmafUsed}, \code{nsnpsTotal}, \code{nsnpsUsed}, \code{nmiss}, \code{df}, 
and \code{p}. A generated RData that is a list that contains \code{scores}, \code{cov}, \code{n}, 
\code{maf} and \code{sey} for each gene with gene names being the names of the list. Note maf in 
RData is MAF based on ntotal. 
}

\author{ Ming-Huei Chen <mhchen@bu.edu> and Qiong Yang <qyang@bu.edu> }

\examples{
\dontrun{
gc.fun(path="/home/mhchen/",phen="trait1",mafRange=c(0,0.01),
snpinfoRdata="SNPinfo_EC.RData",aggregateBy="SKATgene",
maf.file="EC_MAF.csv",snp.cor="EC_SNPcor.RData",ssq.beta.wts=c(1,25))
}
}
