\name{cnvMatrix}
\alias{cnvMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Retrieve CNV matrix
}
\description{
Retrieve CNV matrix from a \code{\link{PatCNVData-class}} object
}
\usage{

cnvMatrix(obj, sample.name = NULL, gene.name = NULL, 
          chr = NULL, pos.range = NULL, 
          exon.score.vec = NULL, min.exon.score = NULL,
          capture.only = FALSE, exon.ID.format = "gene+pos")
          

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
  \item{exon.score.vec}{  a numeric vector of per-exon confidence score
%%     ~~Describe \code{exon.score.vec} here~~
}
  \item{min.exon.score}{  minimum exon confidence score
%%     ~~Describe \code{min.exon.score} here~~
}
  \item{capture.only}{  a logical value indicating if only captured exons should be used 
%%     ~~Describe \code{capture.only} here~~
}
  \item{exon.ID.format}{  "gene+pos" or "pos" indicating how to format exon IDs
%%     ~~Describe \code{exon.ID.format} here~~
}
}
 
\value{
a data matrix of exon by sample log2-ratio CNV values
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
retrive coverage matrix \code{\link{coverageMatrix}}
}

\examples{

#=== load simulated data with pre-computed coverage and CNV
data(sim.cnvData)

#=== display CNV log2ratio of exons of gene "SLC25A17" of samples "n2" and "n5"
cnvMatrix(sim.cnvData, c("n2","n5"),"SLC25A17")

} 


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
 
