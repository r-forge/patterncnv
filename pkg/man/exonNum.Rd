\name{exonNum}
\alias{exonNum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Number of exons
}
\description{
This function gets number of exons in an object of \code{\link{PatCNVSession-class}} or \code{\link{PatCNVData-class}}
}
\usage{
exonNum(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{  an object of \code{\link{PatCNVSession-class}} or \code{\link{PatCNVData-class}}
%%     ~~Describe \code{obj} here~~
}
}
 
\value{
Number of exons in corresponding object
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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
retrieve exon information \code{\link{exonInfo}}; summarize average coverage of exons \code{\link{avgExonCoverage}}
}
\examples{
#=== load a simulation example
config.filename <- "sim1.ini"
makeSimulation(config.filename)
sim.session <- createSession(config.filename)

#=== print number of exons
cat("total number of exons:",exonNum(sim.session),"\n")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

