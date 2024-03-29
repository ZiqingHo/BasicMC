% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BasicMC.R
\name{propertrMat}
\alias{propertrMat}
\title{Check and initialise the transition matrix used}
\usage{
propertrMat(trMatrix)
}
\arguments{
\item{trMatrix}{Transition matrix of the Markov chain}
}
\value{
Transition matrix of the Markov chain
}
\description{
This method is used to check and initialise the transition matrix used
}
\examples{
egstate <- c("a", "b", "c")
MC <- matrix(data = c(1, 0, 0, 0, 1, 0, 0, 0, 1),
         nrow = 3, byrow = TRUE,
         dimnames = list(egstate, egstate))

propertrMat(MC)
}
\references{
Grimmett, G., & Stirzaker, D. (2001). Probability and Random Processes (3rd ed.). New York: Oxford University Press.
}
