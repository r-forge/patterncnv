\name{sampleInfo}
\alias{sampleInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Retrieve sample information
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function extracts specified sample information from an object of \code{\link{PatCNVData-class}}. The typical sample information could be "sample.name", "sample.type", and/or "subject.ID".
}
\usage{
sampleInfo(obj, sample.name = NULL, attrib = c("sample.name", "sample.type"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{  an object of \code{\link{PatCNVSession-class}} or \code{\link{PatCNVData-class}}
%%     ~~Describe \code{obj} here~~
}
  \item{sample.name}{  a character vector containing sample names. A value of NULL instructs the function to select all the samples
%%     ~~Describe \code{sample.name} here~~
}
  \item{attrib}{  a character vector containing selected sample attributes.
  The typical sample information could be "sample.name", "sample.type", and/or "subject.ID".
%%     ~~Describe \code{attrib} here~~
}
}
 
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
a data frame or a vector with given sample information attributes specified by \code{attrib}
}

\author{
%%  ~~who you are~~
Chen Wang
}
 

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
#=== load a simulation example
config.filename <- 'sim1.ini'
makeSimulation(config.filename)
sim.session <- createSession(config.filename)

#=== print sample information
print(sampleInfo(sim.session))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

