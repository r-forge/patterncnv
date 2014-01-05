\name{coverageMatrix}
\alias{coverageMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Retrieve coverage matrix
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Retrieve a raw or normalized coverage data matrix from a \code{\link{PatCNVData-class}} object
}
\usage{
coverageMatrix(obj, sample.name = NULL, gene.name = NULL,
               chr = NULL, pos.range = NULL, exon.score.vec = NULL,
               min.exon.score = NULL, capture.only = FALSE, 
               type = c("RPKM"), exon.ID.format = "gene+pos")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{  an object of \code{\link{PatCNVData-class}}
%%     ~~Describe \code{obj} here~~
}
  \item{sample.name}{  a character vector containing sample names. A value of NULL instructs the function to select all the samples
%%     ~~Describe \code{sample.name} here~~
}
  \item{gene.name}{  a character vector containing gene names. A value of NULL instructs the function to select all the genes
%%     ~~Describe \code{gene.name} here~~
}
  \item{chr}{  selected chromosome
%%     ~~Describe \code{chr} here~~
}
  \item{pos.range}{  a numeric vector of two values indicating start and end positions of selected chromosome
%%     ~~Describe \code{pos.range} here~~
}
  \item{exon.score.vec}{   a numeric vector of per-exon confidence score
%%     ~~Describe \code{exon.score.vec} here~~
}
  \item{min.exon.score}{  minimum exon confidence score
%%     ~~Describe \code{min.exon.score} here~~
}
  \item{capture.only}{  a logical value indicating if only captured exons should be used
%%     ~~Describe \code{capture.only} here~~
}
  \item{type}{ "RPKM" or "raw". default "RPKM" refers to normalized coverage, and "raw" refers to original BP-level coverage
%%     ~~Describe \code{type} here~~
}
  \item{exon.ID.format}{  "gene+pos" or "pos" indicating how to format exon IDs
%%     ~~Describe \code{exon.ID.format} here~~
}
}
 
\value{

a data matrix of exon by sample coverage 
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
retrieve CNV data matrix \code{\link{cnvMatrix}}
}

\examples{

#=== load simulated data with pre-computed coverage information
data(sim.covgData)

#=== list coverage information of a few top exons of sample "n3" and "n5"
head(coverageMatrix(sim.covgData, c("n3","n5")))

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

