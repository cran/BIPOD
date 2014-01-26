\name{BBSim}
\alias{BBSim}
\title{
Function for simulation of p dimensional Brownian bridge
}
\description{
Simulation of p-dimensional driftless SDE with constant
diffusion, conditional on end points:
dV_t=Sigma dW_t, conditional on V_0 and V_T. This function makes a
call to C++ and it is therefore relatively fast.
}
\usage{
BBSim(start, end, n, Sigma=diag(2), T, t0 = 0,seed = 1)
}
\arguments{
  \item{start}{
Numerical vector of length p: Starting point for the process
}
  \item{end}{
Numerical vector of length p: Ending point for the process
}
  \item{n}{
Positive integer: Number of time points where the process is simulated
}
  \item{Sigma}{
p*p matrix: The diffusion matrix for the process
}
  \item{T}{
Positive number: End of time interval.
}
  \item{t0}{
Non negative number, defaults to 0. Start of time interval.
}
  \item{seed}{
Integer, defaults to 1. Specifies seed for random generator. If <=0 it
is set randomly.
}
}
\details{
An n*p matrix with columns representing simulations for each coordinate.}
\value{
  An n*p matrix
}
\examples{
(tmp <- BBSim(start = c(1,2),
              end   = c(3,5),
              n     = 10,
              Sigma = diag(2),
              T     = 2,
              t0    = 0,
              seed  = 1))
matplot(tmp,type="l")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
