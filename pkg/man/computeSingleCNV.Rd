\name{computeSingleCNV}
\alias{computeSingleCNV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Compute copy number variation (CNV) for single sample
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This internal function computes copy number variation (CNV) for single sample. Not intended for users
}
\usage{
computeSingleCNV(session.name, sample.name, ref.type = "average.pattern",
                        episl = 1, small.delta = 1e-05, bin.size = 10, 
                        is.verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{session.name}{ an object of \code{\link{PatCNVSession-class}}
%%     ~~Describe \code{session.name} here~~
}
  \item{sample.name}{  a character vector of sample name defined in sample information file
%%     ~~Describe \code{sample.name} here~~
}
  \item{ref.type}{ "average.pattern" or "basic.paired". "average.pattern" instructs function make CNV calls using trained pattern files.  "basic.paired" instructs function make CNV calls using paired normal/Germline sample.
%%     ~~Describe \code{ref.type} here~~
}
  \item{episl}{  a small value adding to coverage to avoid doing log-transform on zero
%%     ~~Describe \code{episl} here~~
}
  \item{small.delta}{  a small value
%%     ~~Describe \code{small.delta} here~~
}
  \item{bin.size}{  numeric value of exon-level bin-size
%%     ~~Describe \code{bin.size} here~~
}
  \item{is.verbose}{  logical value indicating if the progress information is printed
%%     ~~Describe \code{is.verbose} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Not for users
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
numeric CNV log2-ratio value of single sample
}

\author{
%%  ~~who you are~~
Chen Wang
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
User should directly call \code{\link{computeMultiCNV}} 
}

\examples{

#==========================================================
#=== load a simulation example with generated pattern files
#==========================================================
config.filename <- 'sim2.ini'
makeSimulation(config.filename,with.pattern=TRUE)
sim.session <- createSession(config.filename)

#=== compute CNV of sample "c1"
c1.cnv <- computeSingleCNV(sim.session,"c1")

#=== plot CNV versus position (Mb)
exon.startpos <- exonInfo(sim.session,attri="start")
plot(exon.startpos/1e6, c1.cnv$CNV,
    xlab="position (Mb)",ylab="CNV log2ratio")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}

