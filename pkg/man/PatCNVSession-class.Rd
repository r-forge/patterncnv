\name{PatCNVSession-class}
\Rdversion{1.1}
\docType{class}
\alias{PatCNVSession-class}
\alias{exonInfo,PatCNVSession-method}
\alias{exonNum,PatCNVSession-method}
\alias{sampleInfo,PatCNVSession-method}
\alias{sampleNum,PatCNVSession-method}


\title{Class \code{"PatCNVSession"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
This class of objects is returned by the \code{\link{createSession}} function to represent a session of patternCNV project, specifically, exon information, sample information and sample WIG files. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PatCNVSession", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{exon.info}:}{Object of class \code{"GRanges"}:  the corresponding genomic ranges for the observed exons }
    \item{\code{sample.info}:}{Object of class \code{"data.frame"}:  sample information, containing "sample.name", "sample.type", etc. }
    \item{\code{output.DIR}:}{Object of class \code{"list"}:  output directory information }
    \item{\code{pattern.file}:}{Object of class \code{"list"}:  average- and variability filenames }
    \item{\code{Misc}:}{Object of class \code{"list"}: other information }
    \item{\code{N.exon}:}{Object of class \code{"numeric"}: number of exons }
    \item{\code{N.sample}:}{Object of class \code{"numeric"}: number of samples }
  }
}
\section{Methods}{
  \describe{
    \item{exonInfo}{\code{signature(obj = "PatCNVSession")}:  retrieve exon information }
    \item{exonNum}{\code{signature(obj = "PatCNVSession")}:  number of exons}
    \item{sampleInfo}{\code{signature(obj = "PatCNVSession")}:  retrieve sample information}
    \item{sampleNum}{\code{signature(obj = "PatCNVSession")}:  number of samples}
    \item{show}{\code{signature(object = "PatCNVSession")}:  summarize patternCNV session information}
    \item{summary}{\code{signature(object = "PatCNVSession")}:  summarize patternCNV session information}
	 }
}
\author{
%%  ~~who you are~~
Chen Wang
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{PatCNVData-class}}
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
#=== load a simulation example
config.filename <- "sim1.ini"
makeSimulation(config.filename)
sim.session <- createSession(config.filename)

#=== test if session is created
summary(sim.session)
}
\keyword{classes}
