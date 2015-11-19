\name{pca.default}
\encoding{latin1}
\alias{pca.default}

\title{Principal Components Analysis from the mixOmics package}

\description{ 
Performs a principal components analysis from the \code{\link[mixOmics]{pca}} function of the \code{mixOmics} package.
}

\usage{
    \method{pca}{default}(X, ncomp = 2, center = TRUE, scale = FALSE,
    max.iter = 500, tol = 1e-09,\dots)
}

\arguments{
  \item{X}{a numeric matrix (or data frame) which provides the 
    data for the principal components analysis. It can contain missing values.}
  \item{ncomp}{integer, if data is complete \code{ncomp} decides the number of components and associated 
    eigenvalues to display from the \code{pcasvd} algorithm and if the data has missing values, 
	\code{ncomp} gives the number of components to keep to perform the reconstitution of the data using 
	the NIPALS algorithm. If \code{NULL}, function sets \code{ncomp = min(nrow(X), ncol(X))}}	
  \item{center}{a logical value indicating whether the variables should be shifted to be zero centered. 
    Alternately, a vector of length equal the number of columns of \code{X} can be supplied. 
    The value is passed to \code{\link{scale}}.} 
  \item{scale}{a logical value indicating whether the variables should be scaled to have 
    unit variance before the analysis takes place. The default is \code{FALSE} for consistency with \code{prcomp} 
    function, but in general scaling is advisable. Alternatively, a vector of length equal the number of 
    columns of \code{X} can be supplied. The value is passed to \code{\link{scale}}.} 
  \item{max.iter}{integer, the maximum number of iterations in the NIPALS algorithm.}
  \item{tol}{a positive real, the tolerance used in the NIPALS algorithm.}
\item{\dots}{not used.}
}

\details{
see \code{\link[mixOmics]{pca}}
}
