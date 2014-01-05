\name{fitNullModel}
\alias{fitNullModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Fit Null distribution of copy number variation signals
}
\description{
This function fits null distribution of copy number variation signals, according to given reference/normal samples. The underlying assumption is reference/normal samples have very few number of CNV events, comparing to case/disease samples.
}
\usage{
fitNullModel(refdata.name, type = "autosome", is.plot = FALSE, capture.only = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{refdata.name}{  an object of \code{\link{PatCNVData-class}} of reference/normal samples
%%     ~~Describe \code{refdata.name} here~~
}
  \item{type}{  a character of either "autosome" or "all" to indicate if CNV signals of chr1-22 or all of the chromosomes are used for fitting Null distribution
%%     ~~Describe \code{type} here~~
}
  \item{is.plot}{  a logical value indicating if the null distribution is plotted 
%%     ~~Describe \code{is.plot} here~~
}
  \item{capture.only}{  a logical value indicating if only captured exons should be considered
%%     ~~Describe \code{capture.only} here~~
}
}
 
\value{
This function returns a list of null distribution model that is used by \code{\link{estimateFDR}}
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
\code{\link{estimateFDR}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
#=== load simulated data with pre-computed coverage and CNV
data(sim.cnvData)

#=== fit null distribution of CNV signals and plot distribution of CNV signals (log2ratio)
null.model <- fitNullModel(sim.cnvData,"all",TRUE)
}
 
\keyword{methods}

