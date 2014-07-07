\name{YuGene}
\alias{YuGene}
\encoding{latin1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
YuGene: A simple method for comparing gene expression across platforms using a cumulative
proportion approach.
}
\description{
YuGene is a simple method for comparison of gene expression generated across different experiments, and on different platforms; that does not require global renormalization, and is
not restricted to comparison of identical probes. YuGene works on a range of microarray dataset distributions, such as between manufacturers. The resulting output allows direct comparisons of gene expression between experiments and experimental platforms.
}
\usage{
YuGene(data.prop, progressBar = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.prop}{ a matrix or data.frame of log intensity values, with samples in columns and expression levels in rows. Can be probe or transcript level. Can be raw or previously (i.e. quantile) normalized data. 
}
  \item{progressBar}{set to FALSE to supress progress bar}
}
\value{
returns an object of class `YuGene': a matrix of the same dimensions with each sample transformed to the cumulative proportion (YuGene) metric.
}

\seealso{\code{\link{pca}}}


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
\note{
Support for missing values not yet implemented. Will implement if requested. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

data(ascorbate) # gene expression data available in YuGene package
# apply the transform to the data
YuGene.transformed <- YuGene(ascorbate$gene)

# show distributions before and after YuGene
opar <- par()      # make a copy of current settings
par(mfrow=c(1,2))
plot(density(ascorbate$gene[,1]),main='Expression values', xlab='log2 expr.'); 
plot(density(YuGene.transformed[,1]),main='YuGene values',xlab='YuGene value');
par(opar)     # restore original settings

# unadjusted pvals from the quantile normalized data
quant.pvals <- apply(ascorbate$gene,1,function(row){return(t.test(row[1:3],row[4:6])$p.value)})
YuGene.pvals <- apply(YuGene.transformed,1,function(row){return(t.test(row[1:3],row[4:6])$p.value)})
plot(quant.pvals,YuGene.pvals,pch='.',main='comparison of pvals before and after YuGene Transform')
text(0.8,0.2,paste("Pearson cor: ",round(cor(quant.pvals,YuGene.pvals,method='pearson'),digits=3)))

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ package } %package     &   Package Summaries
\keyword{ manip } % manip       &   Data Manipulation

