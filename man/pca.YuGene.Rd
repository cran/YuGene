\name{pca.YuGene}
\alias{pca}
\alias{pca.YuGene}
\encoding{latin1}
\title{
Principal component analysis for the `YuGene' class.
}
\description{
Performs a principal components analysis thanks to the \code{\link[mixOmics]{pca}} function of the \code{mixOmics} package. The data are centered by study before performing the analysis, if the argument \code{study} is given.
}

\usage{
    \method{pca}{YuGene}(X, study, ncomp = 2, center = TRUE, scale = FALSE,
     max.iter = 500, tol = 1e-09,\dots)
}

\arguments{
  \item{X}{a numeric matrix (or data frame) which provides the
    data for the principal components analysis. It can contain missing values.}
  \item{study}{Factor of the study effect.}
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
If the argument \code{study} is given, the data are centered per study prior to performing the PCA with the \code{\link[mixOmics]{pca}} function of the \code{mixOmics} package.
Otherwise, the PCA is performed on the input data \code{X}.
}


\value{
Same outputs as the \code{pca} function from the \code{mixOmics} package.\cr
\code{pca} returns a list with class \code{"pca"} and \code{"prcomp"} containing the following components:
  \item{ncomp}{the number of principal components used.}
  \item{sdev}{the eigenvalues of the covariance/correlation matrix, though the calculation is actually 
    done with the singular values of the data matrix or by using NIPALS.} 
  \item{rotation}{the matrix of variable loadings (i.e., a matrix whose columns contain the eigenvectors).}
  \item{X}{if \code{retx} is true the value of the rotated data (the centred (and scaled if requested) data 
    multiplied by the rotation matrix) is returned.}  
  \item{center, scale}{the centering and scaling used, or \code{FALSE}.} 
}

\references{
Kim-Anh Lê Cao, Florian Rohart, Leo McHugh, Othmar Korn, Christine A. Wells. YuGene: A simple approach to scale gene expression data derived from different platforms for integrated analyses. Genomics. http://dx.doi.org/10.1016/j.ygeno.2014.03.001.
}
\author{
KA Le Cao, Translational Research Institute, The University of Queensland Diamantina Institute, Australia \cr
Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, St Lucia, Australia \cr
Leo McHugh, Queensland Facility for Advanced Bioinformatics \cr
Othmar Korn, Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, St Lucia, Australia \cr
Christine A. Wells, Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, St Lucia, Australia \cr
}
\examples{
#load data
data(array)

YuGene.data=t(YuGene(t(array$data.all))) # transpose the data to get the samples in columns

#PCA on YuGene data, centered by study
res.pca.yugene.center = pca(YuGene.data, ncomp = 3, scale = TRUE,
            center = TRUE, study = array$experiment.all)
expl.var = round(res.pca.yugene.center$sdev/sum(res.pca.yugene.center$sdev),4)*100

#plot of the results, one color per cell-type, one shape per study
plot(res.pca.yugene.center$x[,1],res.pca.yugene.center$x[,2],
            pch = as.numeric(array$experiment.all),
            col = as.numeric(array$type.all)+1, lwd = 2,
            cex = 1.5, cex.lab = 1.5,xlab=paste("PC1:",expl.var[1],"\%"),
            ylab=paste("PC2:",expl.var[2],"\%"))
title(paste('YuGene multi group data'), cex.main = 1.5)

#PCA on YuGene data, not centered by study
res.pca.yugene = pca(YuGene.data, ncomp = 3, scale = TRUE, center = TRUE)
expl.var = round(res.pca.yugene$sdev/sum(res.pca.yugene$sdev),4)*100

#plot of the results, one color per cell-type, one shape per study
plot(res.pca.yugene$x[,1],res.pca.yugene$x[,2],
            pch = as.numeric(array$experiment.all),
            col = as.numeric(array$type.all)+1, lwd = 2,
            cex = 1.5, cex.lab = 1.5,X.label=paste("PC1:",expl.var[1],"\%"),
            Y.label=paste("PC2:",expl.var[1],"\%"))
title(paste('YuGene data'), cex.main = 1.5)
}

