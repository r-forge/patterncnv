\name{PatCNVData-class}
\Rdversion{1.1}
\docType{class}
\alias{PatCNVData-class}
\alias{avgExonCoverage,PatCNVData-method}
\alias{cnvMatrix,PatCNVData-method}
\alias{coverageMatrix,PatCNVData-method}
\alias{exonInfo,PatCNVData-method}
\alias{exonNum,PatCNVData-method}
\alias{fdrMatrix,PatCNVData-method}
\alias{sampleCoverage,PatCNVData-method}
\alias{sampleInfo,PatCNVData-method}
\alias{sampleNum,PatCNVData-method}
 

\title{Class \code{"PatCNVData"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
This class of objects can be returned by the \code{\link{scanMultiCovg}} or \code{\link{computeMultiCNV}} functions to represent patternCNV results, specifically, coverage, CNV log2-ratio and FDR results. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PatCNVData", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{eSet}:}{Object of class \code{\link{SummarizedExperiment}}: a set container where rows represent exons and columns represent samples }
    \item{\code{method.info}:}{Object of class \code{"list"}: a list of method information }
	\item{\code{datatype.vec}:}{Object of class \code{"character"}: character vector of different data types}
	\item{\code{sample.covg.vec}:}{Object of class \code{"numeric"}: numerical vector of sample-level coverage}
    \item{\code{N.exon}:}{Object of class \code{"numeric"}: number of exons}
    \item{\code{N.sample}:}{Object of class \code{"numeric"}: number of samples }
  }
}
\section{Methods}{
  \describe{
    \item{avgExonCoverage}{\code{signature(obj = "PatCNVData")}: compute average coverage of exons across samples }
    \item{cnvMatrix}{\code{signature(obj = "PatCNVData")}: retrieve log2ratio CNV matrix }
    \item{coverageMatrix}{\code{signature(obj = "PatCNVData")}: retrieve coverage matrix }
    \item{exonInfo}{\code{signature(obj = "PatCNVData")}: retrieve exon information }
    \item{exonNum}{\code{signature(obj = "PatCNVData")}: number of exons }
    \item{fdrMatrix}{\code{signature(obj = "PatCNVData")}: retrieve FDR matrix }
    \item{sampleCoverage}{\code{signature(obj = "PatCNVData")}: summarize sample-level coverage }
    \item{sampleInfo}{\code{signature(obj = "PatCNVData")}: retrieve sample information }
    \item{sampleNum}{\code{signature(obj = "PatCNVData")}: number of samples }
    \item{show}{\code{signature(object = "PatCNVData")}: summarize patternCNV data information}
    \item{summary}{\code{signature(object = "PatCNVData")}: summarize patternCNV data information}
	 }
}
 
\author{
%%  ~~who you are~~
Chen Wang
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{PatCNVSession-class}}
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
#=== load simulated data with pre-computed coverage and CNV
data(sim.cnvData)
 
#=== display class of simulated data
typeof(sim.cnvData)

#=== print basic information of the created object
summary(sim.cnvData)

#=== list coverage information of a few top exons of sample "n3" and "n5"
head(coverageMatrix(sim.cnvData, c("n3","n5")))

#=== display CNV log2ratio of exons of gene "SLC25A17" of samples "n2" and "n5"
cnvMatrix(sim.cnvData, c("n2","n5"),"SLC25A17")

}
\keyword{classes}
