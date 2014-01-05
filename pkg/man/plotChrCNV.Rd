\name{plotChrCNV}
\alias{plotChrCNV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generate chromosome-level CNV plot of given sample
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function generates chromosome-level CNV plot of given sample.
}
\usage{
plotChrCNV(data.name, sample.name, sel.chr = "chr1", pos.range = NULL, 
                  capture.only = TRUE, min.ref.avgRPKM = 2, ref.avg.method = "median", 
                  cex = 0.65, col = "steelblue3", ylim = c(-3, 3), 
                  xlab = "Mb", ylab = "log2-ratio", 
                  ideogram = FALSE, ideogram.cex = 0.75, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.name}{  an object of \code{\link{PatCNVData-class}}
%%     ~~Describe \code{data.name} here~~
}
  \item{sample.name}{  character of selected sample name
%%     ~~Describe \code{sample.name} here~~
}
  \item{sel.chr}{  selected chromosome
%%     ~~Describe \code{sel.chr} here~~
}
  \item{pos.range}{  numeric vector of length 2 containing start and stop positions of chromosome range
%%     ~~Describe \code{pos.range} here~~
}
  \item{capture.only}{  a logical value indicating if only captured exons should be outputted
%%     ~~Describe \code{capture.only} here~~
}
  \item{min.ref.avgRPKM}{  minimum normalized coverage of output exons 
%%     ~~Describe \code{min.ref.avgRPKM} here~~
}
  \item{ref.avg.method}{  a character indicating if exon-level "median" or "mean" coverage information should be outputted in the table
%%     ~~Describe \code{ref.avg.method} here~~
}
  \item{cex}{  plot expansion factor relative to current \code{\link{par}}("cex")
%%     ~~Describe \code{cex} here~~
}
  \item{col}{  color of the plot
%%     ~~Describe \code{col} here~~
}
  \item{ylim}{  the y limits of the plot
%%     ~~Describe \code{ylim} here~~
}
  \item{xlab}{  a title for the x axis
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{  a title for the y axis
%%     ~~Describe \code{ylab} here~~
}
  \item{ideogram}{  a logical value indicating if ideogram of corresponding chromsome should be plotted
%%     ~~Describe \code{ideogram} here~~
}
  \item{ideogram.cex}{  plot expansion factor relative to current \code{\link{par}}("cex") for plotting ideogram
%%     ~~Describe \code{ideogram.cex} here~~
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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
plot autosome CNV results \code{\link{plotAutosomeCNV}}
}
\examples{

#=== load simulated data with pre-computed coverage and CNV
data(sim.cnvData)

#=== plot CNV of sample "c2" in selected region in "chr22"
plotChrCNV(sim.cnvData,"c2","chr22",c(41e6,44e6),cex=1.2)

#=== plot CNV of sample "n3" in selected region in "chr22"
plotChrCNV(sim.cnvData,"n3","chr22",c(41e6,44e6),cex=1.2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{visualization}

