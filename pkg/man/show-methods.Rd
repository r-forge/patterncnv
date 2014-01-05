\name{show-methods}
\docType{methods}
\alias{show-methods}
\alias{show,PatCNVData-method}
\alias{show,PatCNVSession-method}
\title{ Generic Methods for Function \code{show} }
\description{
 print basic information of an object of \code{\link{PatCNVData-class}} or \code{\link{PatCNVSession-class}} 
}
\section{Methods}{
\describe{

\item{\code{signature(object = "PatCNVData")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(object = "PatCNVSession")}}{
%%  ~~describe this method here~~
}
}}

\examples{

#=== load a simulation example
config.filename <- 'sim1.ini'
makeSimulation(config.filename)
sim.session <- createSession(config.filename)

#=== print basic information of created object of "PatCNVSession" class
show(sim.session)

#=== scan coverages of multiple samples
germline.data <- scanMultiCovg(session.name=sim.session)

#=== print basic information of created object of "PatCNVData" class
show(germline.data)

}

\keyword{methods}

