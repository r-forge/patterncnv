\name{chrCoverageQC}
\alias{chrCoverageQC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Quality-check exome coverage across chromosomes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function displays coverage plot for chromosome 1-22, to facilitate quality-check. 
}
\usage{
chrCoverageQC(data.name, legend.layout = "topright", 
                     plot.cex = 2, legend.cex = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.name}{  an object of \code{\link{PatCNVData-class}} }
  \item{legend.layout}{  a character vector of legend location. It should be "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center" or "none". When "none" is used, no legend will be plotted for clarity}
  \item{plot.cex}{ plot point expansion factor}
  \item{legend.cex}{  legend character expansion factor relative to current par("cex") }
}
 
\author{
  Chen Wang
}
 
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
exon-level QC: \code{\link{geneCoverageQC}}
}
\examples{
 
#=== load a simulation example
config.filename <- 'sim1.ini'
makeSimulation(config.filename)
sim.session <- createSession(config.filename)

#=== scan coverages of multiple samples
germline.data <- scanMultiCovg(session.name=sim.session)

#=== only chr22 is supposed to have coverage
chrCoverageQC(data.name=germline.data,legend.layout='topleft')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{QC}
\keyword{visualization}

