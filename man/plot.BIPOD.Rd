\name{plot.BIPOD}
\alias{plot.BIPOD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot function for class BIPOD
}
\description{
Graphical summaries of output from Gibbs sampler
}
\usage{
\method{plot}{BIPOD}(x, theta = NULL, subset = NULL, type, lag = 20,
       interval = NULL, treshold = 0.1, speed = 0.1, truepath = NULL,
       BY = 1, prop = c(.05,.95), diffPriorMean = NULL,
       diffPriorCovar = NULL, log = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
x: An object of class BIPOD.
}
  \item{theta}{
Optional numerical vector with true parameter values.
}
  \item{subset}{
    Numeric vector. Which parameters should be used for plotting.
    Defaults to all.
}
  \item{type}{
    Character choosing plotting type: Either "trace", "hist", "acp",
    "pairs","SDtrace", "accept", "movie" or "cover". See details.
}
  \item{lag}{
Positive integer: Number of lags used in autocorrelation plot
}
  \item{interval}{
Positive integer or numericla vector: If integer, used as burn in for
the Gibbs sampler. If vector, used to subsample Gibbs output. Defaults
to no subsampling and no burn in.
}
  \item{treshold}{
Positive numeric: Cut off for acceptance rate.
}
  \item{speed}{
Positive number used for type="movie": How much time to pause between
each frame of the movie?
}
\item{truepath}{
  Vector or matrix with latent data. Optional.
}
\item{BY}{
integer: Only relevant for type="movie". How many frames to skip for
each iteration?
}
\item{prop}{
Numeric vector with values between 0 and 1: Only relevant for
type="cover". Specifies quantiles for the sampled paths
}
\item{diffPriorMean}{
  Numeric of length 2, only for type=="hist". The prior mean for the
  diffusion coefficients.
}
\item{diffPriorCovar}{
  2*2 matrix, only for type=="hist". The prior covariance for the
  diffusion coefficients.
}
\item{log}{
  Boolean, only for type=="hist". If TRUE, the prior density and the
  estimate of the diffusion coefficients are log transformed before plotting.
}
\item{...}{
Additional arguments to be passed to matplot, density or acf.
}
}
\details{
  Different 'type'-argument gives different plots. More details to come...
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
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
% \examples{
% theta     <- c(5,1,-2,1,.5,.3)
% delta.tmp <- .001
% thin      <- 100
% delta     <- thin*delta.tmp
% diffPriorMean <- c(0,0)
% diffPriorCovar <- diag(2)*2
% Model     <- "OU"

% Data      <- DiffSim(n        = 10000,
%                      start    = c(1,1),
%                      Delta    = delta.tmp,
%                      driftpar = theta[1:4],
%                      Sigma    = diag(theta[5:6]),
%                      seed     = 3,
%                      thin     = thin,
%                      Model    = Model)

% Res <- BIPOD(data    = Data[,1],
%            Delta   = .001*100,
%            ImputeN = 5,
%            seed    = 1,
%            GibbsN  = 300,
%            parKnown= list("gam"=theta[3],"beta"=theta[4]),
%            Start   = theta,
%            driftPriorMean = driftPriorMean,
%            driftPriorCovar     = driftPriorCovar,
%            SigmaRWV= diag(c(.01,.02)),
%            LatentPathStart = 0,
%            Model   = Model)

% plot(A,type="trace",interval=1,theta=theta,subset=c(1:2,5:6))
% plot(A,type="hist",theta=theta,interval=100,subset=c(1:2,5:6),diffPriorMean=diffPriorMean,diffPriorCovar=diffPriorCovar,log=T)
% plot(A,type="acp",lag=1000,subset=c(1:2,5:6),interval=100)
% plot(A,type="SDtrace",interval=50,subset=c(1:2,5:6))
% plot(A,type="accept",interval=50,treshold=.1)
% plot(A,type="movie",truepath=Data[,2],BY=1,speed=.01)
% plot(A,type="cover",truepath=Data[,2],interval=100)
% }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
