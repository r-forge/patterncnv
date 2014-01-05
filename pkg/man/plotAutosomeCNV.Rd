\name{plotAutosomeCNV}
\alias{plotAutosomeCNV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generate autosome (chr1-22) CNV plot of given sample
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function generates autosome (chr1-22) CNV plot of given sample, by connecting CNV values of all the chromosomes
}
\usage{
plotAutosomeCNV(data.name, sample.name, chr.Len.vec, 
                       capture.only = TRUE, min.ref.avgRPKM = 2, 
                       ref.avg.method = "median", ylim = c(-3, 3), cex = 0.6, 
                       xlab = "Mb", ylab = "log2-ratio", 
                       color.vec = c("red", "blue", "green", "orange", "brown",
                       "purple", "black"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.name}{  an object of \code{\link{PatCNVData-class}}
%%     ~~Describe \code{data.name} here~~
}
  \item{sample.name}{ character of selected sample name
%%     ~~Describe \code{sample.name} here~~
}
  \item{chr.Len.vec}{  a numerical vector recording length of 22 chromosomes
%%     ~~Describe \code{chr.Len.vec} here~~
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
  \item{ylim}{  the y limits of the plot
%%     ~~Describe \code{ylim} here~~
}
  \item{cex}{  plot expansion factor relative to current \code{\link{par}}("cex")
%%     ~~Describe \code{cex} here~~
}
  \item{xlab}{  a title for the x axis
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{  a title for the y axis
%%     ~~Describe \code{ylab} here~~
}
  \item{color.vec}{  a character vector of colors used to draw CNV of different chromosomes
%%     ~~Describe \code{color.vec} here~~
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
plot CNV results of given chromosome \code{\link{plotChrCNV}}
}

\examples{
#=== load simulated data with pre-computed coverage and CNV
data(sim.cnvData)

#=== load a numeric vector "hg19.Chr_length_vec" recording length of 22 autosome
data("hg19.Chr_length_vec",package="patternCNV")

#=== plot whole-exome CNV of sample "n3"
#=== only part of chr-22 is expected to be seen as this is according to the simulation
plotAutosomeCNV(sim.cnvData,"n3",hg19.Chr_length_vec)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{visualization}
