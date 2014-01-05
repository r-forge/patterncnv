\name{createSession}
\alias{createSession}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create patternCNV session object
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function loads a given configuration file, and initialize a patternCNV session object with input sample information, exon information and output directories. 
}
\usage{
createSession(config.file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{config.file}{  the name of the configuration file which the patternCNV session information is to be read from. If it does not contain an absolute path, the file name is relative to the current working directory. The name could also be a URL.
%%     ~~Describe \code{config.file} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Configuratio file defines paths used by patternCNV functions. These include paths to where the results will be written and paths to the sample info config file along with the Exon Key. 

Example configuration file:

plot_output_DIR = '/path/to/plot_output/'

txt_output_DIR = '/path/to/txt_output/'

exon_key_file = '/path/to/PatternCNV.Exon.Key.txt'

sample_info_file ='/path/to/sample_info.txt'

}
\value{
An object of \code{\link{PatCNVSession-class} }
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

}
 
\author{
Chen Wang
}
 

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
Scan exome coverage of multiple samples \code{\link{scanMultiCovg}}; 
compute CNV log2-ratio values of multiple samples \code{\link{computeMultiCNV}} 
}
\examples{

##  load a simulation example
config.filename <- 'sim1.ini'
makeSimulation(config.filename)
sim.session <- createSession(config.filename)
  
##-- print basic session information
summary(sim.session)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}

