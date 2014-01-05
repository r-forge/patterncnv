\name{makeSimulation}
\alias{makeSimulation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Make patternCNV simulation files for easy testing
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is designed to faciliate users to test patternCNV, with generated WIG and/or pattern files.
}
\usage{
makeSimulation(ini.file = "simu.ini",
                      sample.info.file = "sample_info.txt",
                      with.pattern = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ini.file}{  name of patternCNV configuration file, which is usually loaded by \code{\link{createSession}} }
  \item{sample.info.file}{  name of sample information file which records sample name/type and path to sample WIG file}
  \item{with.pattern}{  logical, whether trained pattern files need to be included in configuration file}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
With succesful run of this function, the user should be able to automatically generate a configuration file and a sample information file under current working directory. The \code{with.pattern} option instructs whether generated configuration file includes already trained pattern files for simulation study.
}
 
\author{
%%  ~~who you are~~
Chen Wang
}
 
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{createSession}}
}
\examples{

#=============================================================
#=== load a simulation example without generated pattern files
#=============================================================
config.filename <- "sim1.ini"
makeSimulation(config.filename)
sim.session <- createSession(config.filename)

#=== print session information
summary(sim.session)


#==========================================================
#=== load a simulation example with generated pattern files
#==========================================================
config.filename <- "sim2.ini"
makeSimulation(config.filename,with.pattern=TRUE)
sim.session <- createSession(config.filename)

#=== print session information
summary(sim.session)

#=== compute CNV and generate heatmap of selected genes
germline.data <- computeMultiCNV(sim.session)
plotGeneHeatmap(data.name=germline.data,
                       sel.gene.name=c("NAGA","ARFGAP3","CYB5R3","POLR3H","FAM109B"))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}

