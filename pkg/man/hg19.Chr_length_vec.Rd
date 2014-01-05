\name{hg19.Chr_length_vec}
\alias{hg19.Chr_length_vec}
\docType{data}
\title{
 Length of chr1-22 according to hg19
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
load a numeric vector "hg19.Chr_length_vec" recording length of 22 autosome
}
\usage{data(hg19.Chr_length_vec)}
\format{
  hg19.Chr_length_vec
  \describe{
    An integer vector containing BP-level length of chromosome 1-22
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
This is mostly used by \code{\link{plotAutosomeCNV}} function to draw whole-exome CNV plot
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
Hg19 Genome size statistics
 %% \url{http://genomewiki.ucsc.edu/index.php/Hg19_Genome_size_statistics}
}
\examples{
data(hg19.Chr_length_vec)
#=== display the size(BP-level) of 22 chromosomes
print(hg19.Chr_length_vec)
}
\keyword{datasets}
