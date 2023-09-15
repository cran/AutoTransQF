\name{ADStatQF}
\alias{ADStatQF}
\title{
Calculates Anderson-Darling Test Statistic
}
\description{
The function is used to calculate the Anderson-Darling test statistic of standard normal distribution. The input needs to be vectors with seven or more entries.
}
\usage{
## Calculate Anderson-Darling test statistic for vector x:
ADStatQF(x)
}

\arguments{
\item{x}{a vector with seven or more entries.}
}

\value{
Returns the Anderson-Darling test statistic for the given vector.
}
\references{
Feng, Q. , Hannig J. , Marron, J. S. (2016). A Note on Automatic Data Transformation. STAT, 5, 82-87. doi: 10.1002/sta4.104
}
\author{
Yue Hu, Hyeon Lee, J. S. Marron
}
\seealso{
\code{\link{autotransfuncQF}}, \code{\link{AutoTransQF}}
}
\examples{
## Generate a vector z from normal distribution with mean 0 and standard deviation 6.
z = rnorm(10, mean = 0, sd = 6)

## Calculate the Anderson-Darling test statistic for z
ADStatQF(z)
}
