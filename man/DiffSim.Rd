\name{DiffSim}
\alias{DiffSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of a 2-dimensional diffusion process. See the Model
argument for options.
}
\description{
Function for simulation of 2-dimensional diffusion processes, using
the Euler-maruyama scheme.
}
\usage{
DiffSim(n, start, Delta, driftpar, Sigma, seed, thin=1, Model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
positive integer: Length of simulation.
}
  \item{start}{
Numerical vector: Starting point for the simulation.
}
  \item{Delta}{
Numerical: Time interval between observations.
}
\item{driftpar}{
  Numerical vector. Parameters of the FitzHugh-Nagumo model.
  }
  \item{Sigma}{
    2*2 diffusion matrix.
}
\item{seed}{
  Integer: Gives the seed for the random number generator. If seed <= 0
  the seed is 'randomly' chosen.
}
\item{thin}{
  Integer: Output only every 'thin' simulation.
}
\item{Model}{
  Character specifying the model. Currently one of 'OU', 'FHN', 'FHN5'
  and 'CIR'.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
An (n/thin) by 2 matrix.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
FH <- DiffSim(n        = 10000,
              start    = c(1,1),
              Delta    = .001,
              driftpar = c(10,0.6,1.5,0.0),
              Sigma    = diag(c(.5,.3)),
              seed     = 1,
              thin     = 100,
              Model    = "FHN")
matplot(FH,type="l")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
