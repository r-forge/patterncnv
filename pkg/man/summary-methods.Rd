\name{summary-methods}
\docType{methods}
\alias{summary-methods}
\alias{summary,ANY-method}
\alias{summary,PatCNVData-method}
\alias{summary,PatCNVSession-method}
\title{ Generic Methods for Function \code{summary} }
\description{
print basic information of an object of \code{\link{PatCNVData-class}} or \code{\link{PatCNVSession-class}} 
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ANY")}}{
%%  ~~describe this method here~~
}

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
summary(sim.session)

#=== scan coverages of multiple samples
germline.data <- scanMultiCovg(session.name=sim.session)

#=== print basic information of created object of "PatCNVData" class
summary(germline.data)

}


\keyword{methods}

