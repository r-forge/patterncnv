\name{plotGeneHeatmap}
\alias{plotGeneHeatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generate CNV heatmap of given genes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function generates CNV heatmap of given genes. 
}
\usage{

plotGeneHeatmap(data.name, sel.gene.name, capture.only = TRUE,
                       min.ref.avgRPKM = 3, ref.avg.method = "median",
                       logR.cut = 1.5, font.cex = 0.6, heatmap.margin = c(5, 15), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.name}{ an object of \code{\link{PatCNVData-class}}
%%     ~~Describe \code{data.name} here~~
}
  \item{sel.gene.name}{ a character vector containing gene names. A value of NULL instructs the function to select all the genes
%%     ~~Describe \code{sel.gene.name} here~~
}
  \item{capture.only}{  a logical value indicating if only captured exons should be used
%%     ~~Describe \code{capture.only} here~~
}
  \item{min.ref.avgRPKM}{  minimum normalized exon-coverage threshold according to \code{ref.avg.method}
%%     ~~Describe \code{min.ref.avgRPKM} here~~
}
  \item{ref.avg.method}{  a character indicating if exon-level "median" or "mean" coverage information should be considered to filter out exons
%%     ~~Describe \code{ref.avg.method} here~~
}
  \item{logR.cut}{  only exon with absolute CNV values larger or equal to the threshold will be outputted
%%     ~~Describe \code{logR.cut} here~~
}
  \item{font.cex}{  expansion factor for font size
%%     ~~Describe \code{font.cex} here~~
}
  \item{heatmap.margin}{  numeric vector of length 2 containing the margins (see \code{\link{par}}(mar= *)) for column and row names, respectively
%%     ~~Describe \code{heatmap.margin} here~~
}
  \item{\dots}{  additional arguments passed on to \code{\link{image}}
%%     ~~Describe \code{\dots} here~~
}
}
 
\author{
%%  ~~who you are~~
Chen Wang
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
plot CNV results of given chromosome \code{\link{plotChrCNV}}; 
plot autosome CNV results \code{\link{plotAutosomeCNV}}
}
\examples{
#=== load simulated data with pre-computed coverage and CNV
data(sim.cnvData)

#===  plot gene CNV heatmap
plotGeneHeatmap(sim.cnvData,
                       c("NAGA","ARFGAP3","CYB5R3","POLR3H","FAM109B"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{visualization}

