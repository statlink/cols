\name{Positively constrained least squares with a multivariate response}
\alias{mvpls}
\title{
Positively constrained least squares with a multivariate response
}

\description{
Positively constrained least squares with a multivariate response.
}

\usage{
mvpls(y, x)
}

\arguments{
\item{y}{
The response variables, a numerical matrix with observations.
}
\item{x}{
A matrix with independent variables, the design matrix.
}
}

\details{
The constraint is that all beta coefficients (including the constant) are positive.
}

\value{
A list including:
\item{be}{
The positively constrained beta coefficients.
}
\item{mse}{
The mean squared error.
}
}

%\references{
%Hansen, B. E. (2022). Econometrics, Princeton University Press.
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{cls}
}
}

\examples{
y <- as.matrix( iris[, 1:2] )
x <- as.matrix( iris[, 3:4] )
mvpls(y, x)
}
