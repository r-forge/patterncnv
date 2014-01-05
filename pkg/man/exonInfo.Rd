\name{exonInfo}
\alias{exonInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Retrieve exon information
}
\description{
This function retrieves exon information from an object of \code{\link{PatCNVSession-class}} or \code{\link{PatCNVData-class}}
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
exonInfo(obj, gene.name = NULL, capture.only = FALSE, attri = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{  an object of \code{\link{PatCNVData-class}}
%%     ~~Describe \code{obj} here~~
}
  \item{gene.name}{  a character vector containing gene names. A value of NULL instructs the function to select all the genes
%%     ~~Describe \code{gene.name} here~~
}
  \item{capture.only}{  a logical value indicating if only captured exons should be outputted 
%%     ~~Describe \code{capture.only} here~~
}
  \item{attri}{  NULL or a character vector containing specified exon attributes
%%     ~~Describe \code{attri} here~~
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  If attri=NULL, a data.frame of exon information is outputted (each row is an exon):
  \item{gene.name}{  gene name or identifier of corresponding exon}
  \item{chr}{  exon chromosome}
  \item{start}{  exon start position}
  \item{start}{  exon end position}
  \item{is.captured}{  a binary value, whether this exon is captured}
  
  otherwise, 
  a data.frame of specified exon attributes is outputted
}

\author{
%%  ~~who you are~~
Chen Wang
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{exonNum}}, \code{\link{sampleInfo}}
}
\examples{
#=== load a simulation example
config.filename <- "sim1.ini"
makeSimulation(config.filename)
sim.session <- createSession(config.filename)

#=== print first parts of exon information
head(exonInfo(sim.session))

#=== print first parts of exon information with selected attributes
head(exonInfo(sim.session,attri=c("gene.name","chr","start"))) 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

