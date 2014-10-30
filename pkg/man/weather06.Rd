\name{weather05}
\alias{weather06}
\alias{weather04}
\docType{data}
\title{
Weather data
}
\description{
    Weather data as produced by the \code{\link{weach}} function.
  These are for 2004, 2005 and 2006.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(weather06)}
\format{
  A data frame with 8760 observations on the following 8 variables.
  \describe{
    \item{\code{year}}{year}
    \item{\code{doy}}{day of the year}
    \item{\code{hour}}{hour}
    \item{\code{SolarR}}{solar radiation (micro mol/m2/s)}
    \item{\code{Temp}}{temperature (Celsius)}
    \item{\code{RH}}{relative humidity (0-1)}
    \item{\code{WS}}{wind speed (m s-1)}
    \item{\code{precip}}{precipitation (mm)}
  }
}
\source{
  These data are generated from the \code{\link{cmi0506}} object using the
  \code{weach} function.
}
\examples{
data(weather06)
summary(weather06)
## maybe str(weather06) ; plot(weather06) ...
}
\keyword{datasets}
