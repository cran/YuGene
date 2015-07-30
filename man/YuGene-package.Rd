\name{YuGene-package}
\alias{YuGene-package}
\encoding{latin1}
\docType{package}
\title{
Transforms expression datasets to cumulative proportion for comparison
using YuGene transform
}
\description{
YuGene is a simple transform that answers the question: where is my gene of interest? across a range of datasets without the need to re-normalise datasets under consideration. 
}
\details{
\tabular{ll}{
Package: \tab YuGene\cr
Type: \tab Package\cr
Version: \tab 1.1.3\cr
Date: \tab 2015-07-30\cr
License: \tab GPL >= 2\cr
}
This package provides a single function (YuGene). It takes a log transformed dataset (ie multiple microarray samples in an experiment) and converts the values to a cumulative proportion. Values close to zero have the lowest expression, and values close to 1 have the highest expression. When many datasets have been YuGene transformed, relative expression levels (YuGene values) can be directly compared across experiments without re-normalization without significant loss of sensitivity when compared to quantile normalized data. 
}
\author{
KA Le Cao, Translational Research Institute, The University of Queensland Diamantina Institute, Australia \cr
Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, St Lucia, Australia \cr
Leo McHugh, Queensland Facility for Advanced Bioinformatics \cr
Othmar Korn, Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, St Lucia, Australia \cr
Christine A. Wells, Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, St Lucia, Australia \cr
}
\references{
Kim-Anh Lê Cao, Florian Rohart, Leo McHugh, Othmar Korn, Christine A. Wells. YuGene: A simple approach to scale gene expression data derived from different platforms for integrated analyses. Genomics. http://dx.doi.org/10.1016/j.ygeno.2014.03.001.
}

\keyword{ package }
\keyword{ manip }
