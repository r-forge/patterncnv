\name{avgExonCoverage}
\alias{avgExonCoverage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
 Compute average coverage of exons across samples
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
 This function returns across-sample average coverage of all the exons or selected exons from given genes. The average calculation could be \code{mean} or \code{median}, while the later is default method considering its robustness towards outliers. 
}
\usage{
avgExonCoverage(obj, gene.name = NULL, type = c("RPKM"), method = "median")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{  an object of \code{\link{PatCNVData-class}}
%%     ~~Describe \code{obj} here~~
}
  \item{gene.name}{ a character vector of gene.name. A value of NULL instructs the function to use all the genes 
%%     ~~Describe \code{gene.name} here~~
}
  \item{type}{ "RPKM" or "raw" coverage. Default value "RPKM" refers to normalized coverage
%%     ~~Describe \code{type} here~~
}
  \item{method}{ "mean" or "median"; it specificies how average coverage is computed
%%     ~~Describe \code{method} here~~
}
}

\value{
a numerical vector of average coverages of exons
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
 
\author{
%%  ~~who you are~~
  Chen Wang
}
 
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
sample-level coverage summary \code{\link{sampleCoverage}}
}
\examples{

#=== load a simulation example
config.filename <- "sim1.ini"
makeSimulation(config.filename)
sim.session <- createSession(config.filename)

#=== scan coverages of multiple samples
germline.data <- scanMultiCovg(session.name=sim.session)

#=== print median/mean exon-coverage of selected genes
avgExonCoverage(germline.data,gene.name=c("BIK"))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}