\name{sampleNum}
\alias{sampleNum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Number of samples
}
\description{
This function gets number of samples from an object of \code{\link{PatCNVSession-class}} or \code{\link{PatCNVData-class}}
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sampleNum(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  an object of \code{\link{PatCNVSession-class}} or \code{\link{PatCNVData-class}}
%%     ~~Describe \code{obj} here~~
}
}
 
\value{

Number of samples in corresponding object
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
retrieve sample information \code{\link{sampleInfo}}
}
\examples{
#=== load a simulation example
config.filename <- 'sim1.ini'
makeSimulation(config.filename)
sim.session <- createSession(config.filename)

#=== print number of samples
cat("total number of samples:",sampleNum(sim.session),"\n")
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

