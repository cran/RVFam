\name{lme.ped}
\alias{lme.ped}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{function of single SNP analysis and gene-based tests for 
continuous traits with family data using Linear Mixed Effects model}
\description{
   Fit linear mixed effects (LME) model for single SNP analysis that 
tests associations between a continuous phenotype and each genotyped 
SNP on a chromosome in a genotype file and for gene-based tests in 
family data. The association test is carried out by \code{lme.EC} 
function. In each test, the \code{lmekin} function from package 
\code{coxme} is used. 
}
\usage{
lme.ped(phenfile,genfile,pedfile,phen,covars=NULL,mafRange=c(0,0.05),chr,
snpinfoRdata,sep.ped=",",sep.phe=",",sep.gen=" ",aggregateBy="SKATgene",
maf.file,snp.cor,ssq.beta.wts=c(1,25),singleSNP.outfile=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phenfile}{ a character string naming the phenotype file for reading }
  \item{genfile}{ a character string naming the genotype file for reading }
  \item{pedfile}{ a character string naming the pedigree file for reading }
  \item{phen}{ a character string for the phenotype name of a continuous 
trait of interest in \code{test.dat}  }
  \item{covars}{ a character vector for covariates in \code{phenfile} }
  \item{mafRange}{ range of MAF to include SNPs for gene-based burden tests, 
default is c(0,0.05) }
  \item{chr}{ chromosome number that can be 1,2,...,22, and 'X' }
  \item{snpinfoRdata}{ a character string naming the RData containing SNP 
info to be loaded, this should at least include 'Name' (for SNP name), 'Chr', 
and aggregateBy (default='SKATgene') columns }
  \item{sep.ped}{ the field separator character for pedigree file }
  \item{sep.phe}{ the field separator character for phenotype file }
  \item{sep.gen}{ the field separator character for genotype file }
  \item{aggregateBy}{ the column of SNP info on which single SNPs are to be 
aggregated for burden tests, default is 'SKATgene' }
  \item{maf.file}{ a character string naming the comma delimited file 
containing 'Name' for SNP name and 'maf' for MAF }
  \item{snp.cor}{ a character string naming the RData containing lists of 
SNP correlation matrix within each 'SKATgene' }
  \item{ssq.beta.wts}{ a vector of parameters of beta weights used in proposed 
sum of squares test, default=c(1,25) as in SKAT }
  \item{singleSNP.outfile}{ a logical value, TRUE indicating single SNP 
analysis has been done and result files are available for computing SSQ 
using a different \code{mafRange} }
}
\details{
   The \code{lme.ped} function reads in and merges phenotype, genotype, 
and pedigree files, and creates a relationship coefficient matrix using 
\code{pedfile} and \code{kinship2} package to perform single SNP analysis, 
two burden tests (weight=1 for Li & Leal 2008; weight=1/(MAF)/(1-MAF) for 
Madsen & Browning 2009), one sum of squares (SSQ) test (Wei 2009) using 
a LME model as implemented in \code{lmekin} function in \code{coxme} R 
package and to output an RData that is computed based on single SNP results 
and that is compatible with \code{seqMeta} for conducting meta-analysis. 
For burden tests and SSQ test, SNPs genotypes/results are aggregated by 
\code{aggregateBy} (default = "SKATgene") using SNPs selected according 
to user specified \code{mafRange} within each gene (by default). 
\code{genfile} contains unique individual numerical id and genotype data 
on a chromosome, with the column names being "id" and SNP names. For each 
SNP, the genotype data should be coded as 0, 1, 2 indicating the numbers 
of the coded alleles. The SNP name in genotype file should not have any 
dash, '-' and other special characters(dots and underscores are OK). 
\code{phenfile} contains unique individual id, phenotype and covariates 
data, with the column names being "id" and phenotype and covaraite names. 
\code{pedfile} contains pedigree informaion, with the column names being 
"famid","id","fa","mo","sex". Wald chi-square test is used in all genetic 
association tests.
}
\value{
    No value is returned. Instead, tab delimited result files and an 
RData are generated. A single SNP result file, named with \code{phen} and 
singleSNP, contains columns: \code{gene}, \code{Name}, \code{maf}, 
\code{ntotal}, \code{nmiss}, \code{maf_ntotal}, \code{beta}, \code{se}, 
\code{Z}, \code{remark}, \code{p} (p-value from LRT), \code{MAC}, \code{n0}, 
\code{n1}, and \code{n2}. A burden test result file, named with \code{phen} 
and T/MB for Li & Leal 2008/Madsen & Browning 2009 respectively, contains 
columns: \code{gene}, \code{beta}, \code{se}, \code{Z}, \code{cmafTotal}, 
\code{cmafUsed}, \code{nsnpsTotal}, \code{nsnpsUsed}, \code{nmiss}, 
\code{remark}, and \code{p}. A SSQ test result file, named with \code{phen} 
and SSQ, contains columns: \code{gene}, \code{SSQ}, \code{cmafTotal}, 
\code{cmafUsed}, \code{nsnpsTotal}, \code{nsnpsUsed}, \code{nmiss}, \code{df}, 
and \code{p}. A generated RData that is a list that contains \code{scores}, 
\code{cov}, \code{n}, \code{maf} and \code{sey} for each gene with gene names 
being the names of the list. Note \code{maf} in RData is MAF based on ntotal. 
 
  \item{gene }{gene name}
  \item{Name }{SNP name}
  \item{maf }{minor allele frequency based on genotyped sample}
  \item{ntotal }{number of individuals with genotype, phenotype and covariates}
  \item{nmiss }{number of individuals with missing genotype among \code{ntotal}}
  \item{maf_ntotal }{minor allele frequency based on \code{ntotal}}
  \item{beta }{regression coefficient of single SNP test or burden test}
  \item{se }{standard error of \code{beta}}
  \item{Z }{Wald Z statistic}
  \item{remark }{additional information of the analysis }
  \item{p }{p-value of single SNP test or burden test}
  \item{camfTotal }{sum of \code{maf_ntotal} of SNPs in a gene}
  \item{cmafUsed  }{sum of \code{maf_ntotal} of SNPs selected with \code{mafRange} 
in a gene for burden tests or SSQ test}
  \item{nsnpsTotal }{total number of SNPs in a gene}
  \item{nsnpsUsed }{number of SNPs selected and used in burden tests and SSQ test}
  \item{SSQ }{sum of squares statistics}
  \item{df }{degree of freedom of SSQ}
  \item{MAC }{minor allele count}
  \item{n0 }{the number of individuals with 0 copy of coded alleles}
  \item{n1 }{the number of individuals with 1 copy of coded alleles}
  \item{n2 }{the number of individuals with 2 copies of coded alleles}
  \item{scores }{\code{beta}/\code{se}^2 in output RData, where \code{beta} and 
\code{se} are vectors}
  \item{cov }{diag(1/se)*LD matrix*diag(1/se) in output RData}
  \item{n }{maximum \code{ntotal} in a gene in output RData}
  \item{sey }{residual standard error in output RData}

}
\references{ 
 coxme package: mixed-effects Cox models, sparse matrices, and modeling data from 
large pedigrees. Beth Atkinson (atkinson@mayo.edu) for pedigree functions.Terry 
Therneau (therneau@mayo.edu) for all other functions. 2007. Ref Type: Computer 
Program. http://cran.r-project.org/web/packages/coxme/. 

 Abecasis, G. R., Cardon, L. R., Cookson, W. O., Sham, P. C., & Cherny, S. S (2001). 
Association analysis in a variance components framework. \emph{Genet Epidemiol}, 
\bold{21} Suppl 1, S341-S346.

 Li, B. and Leal, S. M (2008). Methods for Detecting Associations with Rare Variants for Common Diseases: 
Application to Analysis of Sequence Data. \emph{Am J Hum Genet}, \bold{83(3)}, 311-321.

 Madsen, B. E. and Browning, S. R (2009). A Groupwise Association Test for Rare Mutations Using a Weighted 
Sum Statistic. \emph{PLoS Genet}, \bold{5(2)} e1000384.

 Wei P (2009). Asymptotic Tests of Association with Multiple SNPs in Linkage Disequilibrium. 
\emph{Genet Epidemiol}, \bold{33(6)}, 497-507. 
}

\author{ Ming-Huei Chen <mhchen@bu.edu> and Qiong Yang <qyang@bu.edu> }

\examples{
\dontrun{
lme.ped(genfile="EC_chr1.txt",phenfile="trait1.csv",pedfile="ped.csv",
phen="trait1",covars=NULL,sep.ped=",",sep.phe=",",sep.gen=" ",mafRange=c(0,0.01),
chr=1,snpinfoRdata="SNPinfo_EC.RData",aggregateBy="SKATgene",maf.file="EC_MAF.csv",
snp.cor="EC_SNPcor.RData",ssq.beta.wts=c(1,25))
}
}
