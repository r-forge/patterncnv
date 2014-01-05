\name{scanSingleCovg}
\alias{scanSingleCovg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Scan exome coverage of single sample
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is not supposed to be directly called by users. Rather, it should be called by \code{scanMultiCovg}
}
\usage{
scanSingleCovg(wig.filename, sample.ID, exon.bin.vec, is.capture.vec, 
                      is.plot = FALSE, log.for.plot = TRUE, plot.output.DIR, 
                      bin.size = 10, ylab = "# of exons", 
                      xlab = "log2(coverage per base)")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wig.filename}{the name of input WIG file }
  \item{sample.ID}{character name of input sample }
  \item{exon.bin.vec}{a numeric vector indicating number of genomic bins per exon }
  \item{is.capture.vec}{a logical vector indicating if corresponding is captured}
  \item{is.plot}{ a logical value indicating if generating coverage plot}
  \item{log.for.plot}{ a logical value indicating if the coverage plot is according to log-scale}
  \item{plot.output.DIR}{ name of directory for generating plot }
  \item{bin.size}{ genomic bin-size}
  \item{xlab}{a title for the x axis}
  \item{ylab}{a title for the y axis}
 
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Not for users
}
\value{
coverage information of a single sample
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
\code{\link{scanMultiCovg}}
}

\examples{

#=== locate an example wig file, which is average pattern of simulation data
wigfile <- system.file("extdata","avg_pattern.wig",package="patternCNV")

#=== load a simulation session with generated pattern files
config.filename <- 'sim2.ini'
makeSimulation(config.filename,with.pattern=TRUE)
sim.session <- createSession(config.filename)
is.capture <- exonInfo(sim.session,attri="is.captured")
exon.Nbin <- exonInfo(sim.session,attri="N.exon.bin")

#=== return coverage information from specificed wig file
wig.out <- scanSingleCovg(wigfile,"avg",exon.Nbin,is.capture)

#=== plot coverage versus position (Mb)
exon.startpos <- exonInfo(sim.session,attri="start")
plot(exon.startpos/1e6,wig.out,log="y",xlab="position (Mb)",ylab="BP-level coverage")

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}

