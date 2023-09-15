\name{autotransfuncQF}
\alias{autotransfuncQF}
\title{
Transforms Individual Vectors into Normality
}
\description{
This function transforms individual vectors into normality. Paper from Feng et al. (2016) includes more details about the transformation mechanism.
}
\usage{
autotransfuncQF(vari, istat, paraindex)
}
\arguments{
  \item{vari}{a vector needs to be transformed}
  \item{istat}{a value representing the type of test statistic for evaluation of normality of the transformed vector. If \code{istat = 1}, Anderson-Darling test statistic is chosen; if \code{istat = 2}, standard skewness statistic is chosen.}
  \item{paraindex}{a value delivered to calculate the specific shift parameter beta}
}

\value{
Returns a transformed vector with the shift parameter calculated by the input paraindex
}

\references{
Feng, Q. , Hannig J. , Marron, J. S. (2016). A Note on Automatic Data Transformation. STAT, 5, 82-87. doi: 10.1002/sta4.104
}
\author{
Yue Hu, Hyeon Lee, J. S. Marron
}
\seealso{
\code{\link{ADStatQF}}, \code{\link{AutoTransQF}}
}
\examples{
## Generate a vector vec from Gamma distribution with parameters shape 1 and scale 2.
vec = rgamma(50, shape = 1, scale = 2)

## Choose Anderson-Darling test statistic for transformed vector.
## Assign paraindex to be 0.9.
autotransfuncQF(vec, istat = 1, paraindex = 0.9)

## Choose skewness for transformed vector
## Assign paraindex to be 0.9
autotransfuncQF(vec, istat = 2, paraindex = 0.9)
}
