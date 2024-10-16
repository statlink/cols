\name{Positive and unit sum constrained least squares}
\alias{pcls}
\title{
Positive and unit sum constrained least squares
}

\description{
Positive and unit sum constrained least squares.
}

\usage{
pcls(y, x)
}

\arguments{
\item{y}{
The response variables, a numerical vector with observations.
}
\item{x}{
A matrix with independent variables, the design matrix.
}
}

\details{
The constraint is that all beta coefficients are positive and sum to 1.
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

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{pls}
}
}

\examples{
x <- as.matrix( iris[1:50, 1:4] )
y <- rnorm(50)
pcls(y, x)
}
