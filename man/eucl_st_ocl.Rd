\name{eucl_st_ocl}
\alias{eucl_st_ocl}
\encoding{UTF-8}
\title{Checking Correlation Model}
\description{
  The procedure controls if the correlation model
  inserted is correct.
}
\usage{
eucl_st_ocl(theta,fixed,coordx,coordy,ncoords,times,ntime,cc,datos,type_dist,maxdist,maxtime,
                      winc_s,winstp_s,winc_t,winstp_t,
                      weighted,type_sub,local,GPU,kernel)
}
\arguments{
  \item{theta}{Numeric; a vector which has fixed parameters info.}
  \item{fixed}{Numeric; a vector which has to be estimated parameters info.}
  \item{coordx}{Numeric; it indicates the correlation model}
  \item{coordy}{Numeric; it indicates the correlation model}
  \item{ncoords}{Numeric; it indicates the correlation model}
  \item{times}{Numeric; it indicates the correlation model}
  \item{ntime}{Numeric; it indicates the correlation model}
  \item{cc}{Numeric; it indicates the correlation model}
  \item{datos}{Numeric; it indicates the correlation model}
  \item{type_dist}{Numeric; it indicates the correlation model}
  \item{maxdist}{Numeric; it indicates the correlation model}
  \item{maxtime}{Numeric; it indicates the correlation model}
  \item{winc_s}{Numeric; it indicates the correlation model}
  \item{winstp_s}{Numeric; it indicates the correlation model}
  \item{winc_t}{Numeric; it indicates the correlation model}
  \item{winstp_t}{Numeric; it indicates the correlation model}
  \item{weighted}{Numeric; it indicates the correlation model}
  \item{type_sub}{Numeric; it indicates the correlation model}
  \item{local}{Numeric; it indicates the correlation model}
  \item{GPU}{Numeric; it indicates the correlation model}
  \item{kernel}{Character; It has the function name for binary convertion}
  
  
}

\value{Return the \code{theta} parameters with correct names.} 




\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
}

\keyword{Composite}
