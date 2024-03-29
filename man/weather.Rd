% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BasicMC.R
\name{weather}
\alias{weather}
\title{Convert the data into a transition matrix for weather forecasting}
\usage{
weather(y, weathercolumn, sunny, cloudy, rainy, other = NULL)
}
\arguments{
\item{y}{Data of the weather}

\item{weathercolumn}{The index of the column containing weather types}

\item{sunny}{The character vector for 'sunny' state}

\item{cloudy}{The character vector for 'cloudy' state}

\item{rainy}{The character vector for 'rainy' state}

\item{other}{The character vector for other weather types. The default value = NULL.}
}
\value{
The transition matrix of the Markov chain
}
\description{
This method is used to create a transition matrix from the data of weather
}
\examples{
data(wsoton)
rainnystate<-c("Light rain", "Moderate rain", "Patchy rain possible",
             "Light drizzle", "Light rain shower")
weather(wsoton,2,c("Sunny"),c("Cloudy","Partly cloudy"),rainnystate,"Overcast")
}
