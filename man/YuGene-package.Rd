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
Version: \tab 1.0\cr
Date: \tab 2012-11-27\cr
License: \tab GPL >= 2\cr
}
This package provides a single function (YuGene). It takes a log transformed dataset (ie multiple microarray samples in an experiment) and converts the values to a cumulative proportion. Values close to zero have the lowest expression, and values close to 1 have the highest expression. When many datasets have been YuGene transformed, relative expression levels (YuGene values) can be directly compared across experiments without re-normalization without significant loss of sensitivity when compared to quantile normalized data. 
}
\author{
KA Le Cao, Queensland Facility for Advanced Bioinformatics, Institute for Molecular Biology, The University of Queensland, St Lucia, Australia \cr
Leo McHugh, Queensland Facility for Advanced Bioinformatics \cr
Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, St Lucia, Australia 
}
\references{
Le Cao K. A., McHugh L., Rohart F. and Wells C. A. YuGene: A simple approach to compare gene expression data derived from different platforms.  Submitted July 2013.
}

\keyword{ package }
\keyword{ manip }
