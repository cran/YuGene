\name{ascorbate}
\alias{ascorbate}
\docType{data}
\title{Ascorbate Experiment}
\description{log2 transformed samples using Illumina HumanWG-6 chips, 3 of which were controls, and three of which were sampled after the addition of ascorbate to the medium. Details and data available by searching `ascorbate' at www.stemformatics.org. This dataset is a random subset of 5000 genes for smaller package size and faster example times}
\usage{data(ascorbate)}
\format{

A list containing the following components:
  \describe{
    \item{\code{gene}}{data frame with 48803 rows and 6 columns. 
  	The expression levels of 48803 transcripts for the 6 subjects.}
    \item{\code{condition}}{a vector of 6 elements indicating the condition of each subject ('4ng.ml' or '100ng.ml')}
    }
}

\source{
The data were downloaded from \url{www.stemformatics.org} datasetID 5006.
}
\references{
Chung TL, Brena RM, Kolle G, Grimmond SM, Berman BP, Laird PW, Pera MF, Wolvetang EJ (2010). Vitamin C Promotes Widespread Yet Specific DNA Demethylation of the Epigenome in Human Embryonic Stem Cells; Stem Cells, 28 (10) 1848-1855, 
}
\keyword{datasets}
