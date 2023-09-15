\name{AutoTransQF}
\alias{AutoTransQF}
\title{
Transforms Datasets into Normality
}
\description{
This function helps to transform each vector of the matrix into normality based on the optimal test statistic of transformed vectors.
}
\usage{
## The function tries to transform each vector of mdata into normality
AutoTransQF(mdata, paramstruct = list(istat, iscreenwrite, FeatureNames))
}

\arguments{
  \item{mdata}{the matrix needs to be transformed.}
  \item{paramstruct}{A list with three entries istat, iscreenwrite and FeatureNames respectively. Missing entries will be set to default.}
  \item{istat}{a value representing the type of test statistic for evaluation of normality of the transformed vector with default to be \code{istat = 1}. If \code{istat = 1}, Anderson-Darling test statistic is chosen; if \code{istat = 2}, standard skewness statistic is chosen.}
  \item{iscreenwrite}{Whether there is screenwrite with default to be \code{iscreenwrite = 0}. If \code{iscreenwrite = 1}, to write progress to screen; if \code{iscreenwrite = 0}, no screenwrite.}
  \item{FeatureNames}{Contains feature names of each vector with default to be 'Feature1'}
}

\value{
Returns a list with three elements: 
  \item{data}{the transformed matrix}
  \item{beta}{a list of all shift parameters beta}
  \item{alpha}{a list of all shift parameters alpha}
}
\note{
When a vector of the original matrix is not transformed, its corresponding alpha and beta are both -1.
}
\references{
Feng, Q. , Hannig J. , Marron, J. S. (2016). A Note on Automatic Data Transformation. STAT, 5, 82-87. doi: 10.1002/sta4.104
}
\author{
Yue Hu, Hyeon Lee, J. S. Marron
}
\seealso{
\code{\link{ADStatQF}}, \code{\link{autotransfuncQF}}
}
\examples{
## Create a random matrix x.
x = matrix(rgamma(40, shape = 1, scale = 2), nrow = 4)

## Transform matrix x in default setting and 
## output transformed data
AutoTransQF(x)$data

## Transform matrix x in default setting and 
## output a list of shift parameter beta
AutoTransQF(x)$beta

## Transform matrix x with feature names and
## output a list of shift parameter alpha
Names = c('Feature1', 'Feature2', 'Feature3', 'Feature4')
AutoTransQF(x, paramstruct = list(FeatureNames = Names))$alpha

## Transform matrix x with feature names, progress to screen,
## and apply standard skewness statistic to transformed vectors
AutoTransQF(x, paramstruct = list(istat = 2, iscreenwrite = 1, FeatureNames = Names))

## Transform matrix x with progress to screen and 
## apply standard skewness statistic to transformed vectors
AutoTransQF(x, paramstruct = list(istat = 2, iscreenwrite = 1))
}

