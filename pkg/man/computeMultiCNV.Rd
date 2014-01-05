\name{computeMultiCNV}
\alias{computeMultiCNV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Compute copy number variation (CNV) for multiple samples
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function returns an object of \code{\link{PatCNVData-class}}. 
}
\usage{
computeMultiCNV(session.name, data.name = NULL, sample.type = NULL,
                       ref.type = "average.pattern", episl = 1, 
                       small.delta = 1e-05, bin.size = 10, 
                       zero.median.adjust = TRUE, is.verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{session.name}{
  an object of \code{\link{PatCNVSession-class}}
%%     ~~Describe \code{session.name} here~~
}
  \item{data.name}{
  an object of \code{\link{PatCNVData-class}}. A NULL value instructs the function to create a new \code{\link{PatCNVData-class}} object
%%     ~~Describe \code{data.name} here~~
}
  \item{sample.type}{ 
  a character vector containing sample types. A value of NULL instructs the function to select all the types
%%     ~~Describe \code{sample.type} here~~
}
  \item{ref.type}{  "average.pattern" or "basic.paired". "average.pattern" instructs function make CNV calls using trained pattern files.  "basic.paired" instructs function make CNV calls using paired normal/Germline sample.
%%     ~~Describe \code{ref.type} here~~
}
  \item{episl}{  a small value adding to coverage to avoid doing log-transform on zero
%%     ~~Describe \code{episl} here~~
}
  \item{small.delta}{  a small value
%%     ~~Describe \code{small.delta} here~~
}
  \item{bin.size}{   numeric value of exon-level bin-size
%%     ~~Describe \code{bin.size} here~~
}
  \item{zero.median.adjust}{  logical value indicating if computed CNV should be adjusted to have zero median
%%     ~~Describe \code{zero.median.adjust} here~~
}
  \item{is.verbose}{   logical value indicating if the progress information is printed
%%     ~~Describe \code{is.verbose} here~~
}
}
 
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
 an object of \code{\link{PatCNVSession-class}}, with CNV log2-ratio computed
}
\author{
%%  ~~who you are~~
Chen Wang
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
Scan exome coverage of multiple samples \code{\link{scanMultiCovg}}
}
\examples{

#==========================================================
#=== load a simulation example with generated pattern files
#==========================================================
config.filename <- "sim2.ini"
makeSimulation(config.filename,with.pattern=TRUE)
sim.session <- createSession(config.filename)

#=== print session information
summary(sim.session)

#=== compute CNV and generate heatmap of selected genes
germline.data <- computeMultiCNV(sim.session)
plotGeneHeatmap(data.name=germline.data,
                       sel.gene.name=c("NAGA","ARFGAP3","CYB5R3","POLR3H","FAM109B"))

 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
 
