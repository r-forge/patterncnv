\name{segmentCNV}
\alias{segmentCNV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Output CNV segmentation results 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function outputs CNV segmentation results as .bed files and PDF files. It calls CBS segmentation method in \code{\link{DNAcopy}}
}
\usage{
segmentCNV(session.name, data.name, output_suffix = "CNV_seg",
                  ref.avg.method = "median", min.ref.avgRPKM = 2, 
                  capture.only = TRUE, is.plot = TRUE, 
                  ylim = c(-3, 3), cex = 0.6, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{session.name}{  an object of \code{\link{PatCNVSession-class}} 
%%     ~~Describe \code{session.name} here~~
}
  \item{data.name}{  an object of \code{\link{PatCNVData-class}} 
%%     ~~Describe \code{data.name} here~~
}
  \item{output_suffix}{  filename suffix of output CNV segmentation results
%%     ~~Describe \code{output_suffix} here~~
}
  \item{ref.avg.method}{  a character indicating if exon-level "median" or "mean" coverage information should be considered to filter out exons
%%     ~~Describe \code{ref.avg.method} here~~
}
  \item{min.ref.avgRPKM}{  minimum normalized exon-coverage threshold according to \code{ref.avg.method}
%%     ~~Describe \code{min.ref.avgRPKM} here~~
}
  \item{capture.only}{  a logical value indicating if only captured exons should be used
%%     ~~Describe \code{capture.only} here~~
}
  \item{is.plot}{  a logical value indicating if CNV segmentation plots need to be generated
%%     ~~Describe \code{is.plot} here~~
}
  \item{ylim}{  the y limits of the plot
%%     ~~Describe \code{ylim} here~~
}
  \item{cex}{  plot expansion factor relative to current \code{\link{par}}("cex")
%%     ~~Describe \code{cex} here~~
}
  \item{\dots}{  other graphical parameters (see \code{\link{par}})
%%     ~~Describe \code{\dots} here~~
}
}
 
 
\author{
%%  ~~who you are~~
Chen Wang
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
#=== load a simulation example
config.filename <- "sim2.ini"
makeSimulation(config.filename,with.pattern=TRUE)
sim.session <- createSession(config.filename)

#=== compute CNV 
germline.data <- computeMultiCNV(sim.session)

#=== generate CNV segmentation results
segmentCNV( sim.session, germline.data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}
\keyword{visualization}

