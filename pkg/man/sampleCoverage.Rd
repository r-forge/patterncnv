\name{sampleCoverage}
\alias{sampleCoverage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Summarize sample coverage
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function gets total coverage of samples from an object of \code{\link{PatCNVData-class}}
}
\usage{
sampleCoverage(obj, sample.type = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ 
  an object of \code{\link{PatCNVSession-class}} or \code{\link{PatCNVData-class}}
%%     ~~Describe \code{obj} here~~
  }
  \item{sample.type}{ 
   a character vector containing sample types. A value of NULL instructs the function to select all the types of samples
%%     ~~Describe \code{sample.type} here~~
    }
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
a numerical vector of total coverage per sample
}
\author{
%%  ~~who you are~~
Chen Wang
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
compute exon average coverage \code{\link{avgExonCoverage}}
}
\examples{

#=== load a simulation example
config.filename <- 'sim1.ini'
makeSimulation(config.filename)
sim.session <- createSession(config.filename)

#=== scan coverages of multiple samples
germline.data <- scanMultiCovg(session.name=sim.session)

#=== print BP-level coverage per sample
sampleCoverage(germline.data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

