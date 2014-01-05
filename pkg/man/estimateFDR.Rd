\name{estimateFDR}
\alias{estimateFDR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Estimate false discovery rate of copy number variation signals
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function estimates false discovery rate (FDR) of CNV log2-ratio values, with given Null distribution model returned by \code{\link{fitNullModel}}
}
\usage{
estimateFDR(data.name, null.model, 
                   FDR.type = "q.value", is.plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.name}{  an object of \code{\link{PatCNVData-class}} 
%%     ~~Describe \code{data.name} here~~
}
  \item{null.model}{  a list of Null distribution model returned by \code{\link{fitNullModel}}
%%     ~~Describe \code{null.model} here~~
}
  \item{FDR.type}{  "q.value" or "localFDR"
%%     ~~Describe \code{FDR.type} here~~
}
  \item{is.plot}{  logical, whether intermediate plots are generated
%%     ~~Describe \code{is.plot} here~~
}
}
 
\value{  an object of \code{\link{PatCNVData-class}}, with FDR results computed
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
compute CNV log2-ratio values of multiple samples \code{\link{computeMultiCNV}} 
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{

#=== load simulated data with pre-computed coverage and CNV
data(sim.cnvData)

#=== fit null distribution of CNV signals and plot distribution
null.model <- fitNullModel(sim.cnvData, "all")
germline.data <- estimateFDR(sim.cnvData, null.model)

#=== plot distribution of FDR of sample "n2"
hist(fdrMatrix(germline.data,"n2"), xlab="FDR", 
     main="FDR distribution of sample n2")

}
\keyword{methods}

