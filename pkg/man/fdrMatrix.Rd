\name{fdrMatrix}
\alias{fdrMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Retrieve FDR matrix
}
\description{

This function retrieves a data matrix of p-value or FDR from an object of \code{\link{PatCNVData-class}}. 

}
\usage{
fdrMatrix(obj, sample.name = NULL, gene.name = NULL, 
          chr = NULL, pos.range = NULL, exon.score.vec = NULL,
          min.exon.score = NULL, capture.only = FALSE, 
          type = "FDR", exon.ID.format = "gene+pos")
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
  \item{type}{  "FDR" or "pval". Default value "FDR" refers to false discovery rate and "pval" refers to nominal p-value
%%     ~~Describe \code{type} here~~
}
  \item{exon.ID.format}{  "gene+pos" or "pos" indicating how to format exon IDs
%%     ~~Describe \code{exon.ID.format} here~~
}
}
 
\value{

an FDR/p-value matrix of exon by sample  
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
\code{\link{estimateFDR}}, \code{\link{fitNullModel}}
}

\examples{
#=== load simulated data with pre-computed coverage and CNV
data(sim.cnvData)

#=== fit null distribution of CNV signals and plot distribution
null.model <- fitNullModel(sim.cnvData, "all")
germline.data <- estimateFDR(sim.cnvData, null.model)

#=== list FDR of a few top exons of samples "n2" and "n5"
head(fdrMatrix(germline.data,c("n2","n5")))

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

