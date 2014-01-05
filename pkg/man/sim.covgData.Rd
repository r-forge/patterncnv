\name{sim.covgData}
\alias{sim.covgData}
\docType{data}
\title{
simulated data of "PatCNVData" class with coverage information
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
This is a simulated data of "PatCNVData" class with coverage information for testing purpose.
}
\usage{data(sim.covgData)}
\format{
  sim.covgData of "PatCNVData" class 
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~

The purpose of this data is to faciliate users to test examples.
If interested, this simulation data can be easily generated as follows:

#=== load a simulation example

config.filename <- "sim1.ini"

makeSimulation(config.filename)

#=== create patternCNV session

sim.session <- createSession(config.filename)

#=== summarize coverage information

sim.covgData <- scanMultiCovg(sim.session)

}
\examples{
#=== load simulated data with coverage information
data(sim.covgData)

#=== display coverage of a few top exons in sample "n2" and "n3"
head(coverageMatrix(sim.covgData,c("n2","n3")))
}
\keyword{datasets}
