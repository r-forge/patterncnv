\name{learnPatterns}
\alias{learnPatterns}

\title{
%%  ~~function to do ... ~~
  Learn average and variability patterns
}
\description{
  This function is used to learn patterns from normal/reference samples
}
\usage{
learnPatterns(session.name, refdata.name, sample.type = NULL, 
                     episl = 1, bin.size = 10, excluded.sample.name = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{session.name}{ an object of \code{\link{PatCNVSession-class} }
%%     ~~Describe \code{session.name} here~~
}
  \item{refdata.name}{ an object of \code{\link{PatCNVData-class} } of reference samples, which should typically be normal or Germline samples
%%     ~~Describe \code{refdata.name} here~~
}
  \item{sample.type}{
a character vector of sample type defined in sample information file: e.g. "Germline". A value of NULL instructs the function to select all the samples regardless sample types 
}
  \item{episl}{
%%     ~~Describe \code{episl} here~~
a small value adding to coverage to avoid doing log-transform on zero
}
  \item{bin.size}{
%%     ~~Describe \code{bin.size} here~~
numeric value of exon-level bin-size
}
  \item{excluded.sample.name}{
%%     ~~Describe \code{excluded.sample.name} here~~
a character vector containing sample names that are excluded for pattern training
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Average- and variability- pattern WIG files are expected to be generated under txt_output_DIF defined in configuration file, along with two other baseline coverage files.
}

 
\author{
%%  ~~who you are~~
  Chen Wang
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{scanMultiCovg}}; \code{\link{computeMultiCNV}}
}
\examples{ 
#=== load a simulation example
config.filename <- "sim1.ini"
makeSimulation(config.filename)
sim.session <- createSession(config.filename)

#=== print session information
summary(sim.session)

#=== scan coverages of multiple samples
germline.data <- scanMultiCovg(session.name=sim.session)

#=== learn average- and variability-patterns
learnPatterns(session.name=sim.session,refdata.name=germline.data)

#=== print session information; notice that pattern files have been generated
summary(sim.session)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

