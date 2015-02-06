\name{RVFam-package}
\alias{RVFam-package}
\alias{RVFam}
\docType{package}
\title{
Rare Variants Association Analyses with Family Data
}
\description{
\code{RVFam} package provides functions to perform single SNP association analyses and gene-based tests for continuos, binary and 
survival traits against sequencing variant genotypes (e.g. exome chip and whole genome sequencing data) using family data. The 
gene-based tests include two burden tests, most powerful when effects are in the same direction across all included variants 
(Li and Leal 2008 and Madsen and Browning 2009), and one directional insensitive test (Wei 2009). For single SNP association 
analyses of continuous traits, \code{RVFam} fits linear mixed effects (LME) model with relationship coefficcient matrix as within 
pedigree correlation matrix to account for familial correlation; for binary traits, \code{RVFam} fits generalized linear mixed 
effects (GLMM) model that treats each pedigree as a cluster; while for survival traits, \code{RVFam} fits Cox proportional 
hazards regression model (COXPH) with frailty that adds a random effect for family clusters. 

}
\details{
\tabular{ll}{
Package: \tab RVFam\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-02-06\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Qiong Yang <qyang@bu.edu> and Ming-Huei Chen <mhchen@bu.edu>

Maintainer: Ming-Huei Chen <mhchen@bu.edu>
}