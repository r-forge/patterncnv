\name{exportCNVtable}
\alias{exportCNVtable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Exports exon-level CNV information as per-sample tables
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function outputs exon-level CNV information for each sample into a table, including basic exon information, raw and normalized coverage and CNV information.

}
\usage{
exportCNVtable(session.name, data.name, 
                      output_suffix = "_CNV_table.txt", 
                      ref.avg.method = "median", min.ref.avgRPKM = 2, 
                      capture.only = TRUE, FDR.threshold = 0.05, 
                      logR.cut = 0.1, is.verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{session.name}{  an object of \code{\link{PatCNVSession-class}}
%%     ~~Describe \code{session.name} here~~
}
  \item{data.name}{   an object of \code{\link{PatCNVData-class}}
%%     ~~Describe \code{data.name} here~~
}
  \item{output_suffix}{  filename suffix of output sample tables
%%     ~~Describe \code{output_suffix} here~~
}
  \item{ref.avg.method}{  a character indicating if exon-level "median" or "mean" coverage information should be outputted in the table
%%     ~~Describe \code{ref.avg.method} here~~
}
  \item{min.ref.avgRPKM}{  minimum normalized coverage of output exons  
%%     ~~Describe \code{min.ref.avgRPKM} here~~
}
  \item{capture.only}{  a logical value indicating if only captured exons should be outputted
%%     ~~Describe \code{capture.only} here~~
}
  \item{FDR.threshold}{  only exon with FDR below the threshold will be outputted
%%     ~~Describe \code{FDR.threshold} here~~
}
  \item{logR.cut}{  only exon with absolute CNV values larger or equal to the threshold will be outputted
%%     ~~Describe \code{logR.cut} here~~
}
  \item{is.verbose}{  logical value indicating if the progress information is printed
%%     ~~Describe \code{is.verbose} here~~
}
}
 
 
 
\author{
%%  ~~who you are~~
Chen Wang
}
 

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{computeMultiCNV}} 
}
\examples{
#=== load a simulation example
config.filename <- "sim2.ini"
makeSimulation(config.filename,with.pattern=TRUE)
sim.session <- createSession(config.filename)

#=== compute CNV and generate heatmap of selected genes
germline.data <- computeMultiCNV(sim.session)

#=== export per-sample CNV tables for exons with absolute CNV values larger or equal to 0.5
exportCNVtable(sim.session,germline.data,logR.cut=0.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}
 
