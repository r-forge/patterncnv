\name{geneCoverageQC}
\alias{geneCoverageQC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Quality-check exon-coverage of given genes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function displays coverage plot for given genes, to facilitate quality-check if samples having reasonable coverage in particular gene regions.
}
\usage{

geneCoverageQC(data.name, sel.gene.name, 
                      legend.layout = "topleft", legend.cex = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.name}{  an object of \code{\link{PatCNVData-class}} }
  \item{sel.gene.name}{  a character vector of genes of interest}
  \item{legend.layout}{  a character vector of legend location. It should be "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center" or "none". When "none" is used, no legend will be plotted for clarity}
  \item{legend.cex}{  legend character expansion factor relative to current par("cex")}

}
\author{
%%  ~~who you are~~
  Chen Wang
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{chrCoverageQC}}
}
\examples{

#=== load a simulation example
config.filename <- "sim1.ini"
makeSimulation(config.filename)
sim.session <- createSession(config.filename)

#=== scan coverages of multiple samples
germline.data <- scanMultiCovg(session.name=sim.session)

#=== generate coverage plot for two genes
geneCoverageQC(germline.data, c("BIK","TTLL1"))

#=== generate the same coverage plot but hiding sample-legend for clarity
geneCoverageQC(germline.data, c("BIK","TTLL1"), legend.layout="none")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{QC}
\keyword{visualization}


