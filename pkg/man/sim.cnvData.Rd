\name{sim.cnvData}
\alias{sim.cnvData}
\docType{data}
\title{
simulated data of "PatCNVData" class with coverage and CNV information
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
This is a simulated data of "PatCNVData" class with coverage and CNV information for testing purpose.
}
\usage{data(sim.cnvData)}
\format{
  sim.cnvData of "PatCNVData" class 
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~

The purpose of this dataset is to faciliate users to test examples quickly.
If interested, this simulation dataset can be easily re-generated as follows:

#=== load a simulation example

config.filename <- "sim1.ini"

makeSimulation(config.filename)

#=== create patternCNV session

sim.session <- createSession(config.filename)

#=== compute CNV

sim.cnvData <- computeMultiCNV(sim.session)

}
\examples{
#=== load simulated data with pre-computed coverage and CNV
data(sim.cnvData)

#=== list CNV log2ratio of exons of gene "SLC25A17" of samples "n2" and "n5"
cnvMatrix(sim.cnvData,c("n2","n5"),"SLC25A17")
}
\keyword{datasets}
