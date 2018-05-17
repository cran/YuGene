\name{YuGene-package}
\alias{YuGene-package}
\docType{package}
\title{
\packageTitle{YuGene}
}
\description{
\packageDescription{YuGene}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{YuGene}
\packageIndices{YuGene}
This package provides a single function (YuGene). It takes a log transformed dataset (ie multiple microarray samples in an experiment) and converts the values to a cumulative proportion. Values close to zero have the lowest expression, and values close to 1 have the highest expression. When many datasets have been YuGene transformed, relative expression levels (YuGene values) can be directly compared across experiments without re-normalization without significant loss of sensitivity when compared to quantile normalized data.
}
\author{
\packageAuthor{YuGene}

Maintainer: \packageMaintainer{YuGene}
}
\references{
Kim-Anh L Cao, Florian Rohart, Leo McHugh, Othmar Korn, Christine A. Wells. YuGene: A simple approach to scale gene expression data derived from different platforms for integrated analyses. Genomics. http://dx.doi.org/10.1016/j.ygeno.2014.03.001.
}

