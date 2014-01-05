\name{patternCNV-package}
\alias{patternCNV-package}
\alias{patternCNV}
\docType{package}
\title{
patternCNV: a versatile tool for detecting copy number changes from exome-seq data
}
\description{
PatternCNV is a versatile tool to facilitate CNV detection and interpretation from either Germline or Tumor/Normal pairs of Exome-seq samples. The idea behind proposed patternCNV approach is to summarize common coverage trends of a same exome capture-kit, and then facilitate accurate CNV estimation. Across multiple control-samples (i.e. germline), we proposed to summarize average and variability of exon-coverage, defined as "average-pattern" and "variability-pattern", respectively. Then, exon-level CNV is reliably estimated by computing difference between observed read-depth and expected coverage
}
\details{
\tabular{ll}{
Package: \tab patternCNV\cr
Type: \tab Package\cr
Version: \tab 0.99.2\cr
Date: \tab 2013-12-18\cr
License: \tab Artistic-2.0\cr
Depends: \tab methods, GenomicRanges, IRanges\cr
}
 
}
\author{
Chen Wang

Maintainer: Chen Wang <Wang.Chen@mayo.edu> or <topsoil.ustc@gmail.com>
}
\references{
PatternCNV: a versatile tool for detecting copy number changes from exome sequencing data (submitted)
}
\keyword{ package }
\examples{

#=== load a simulation example
config.filename <- 'sim1.ini'
makeSimulation(config.filename)
sim.session <- createSession(config.filename)

#=== print summary information of created patternCNV session
summary(sim.session)
}
