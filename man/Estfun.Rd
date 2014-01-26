\name{Estfun}
\alias{Estfun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameter estimation for some two dimensional diffusions.
}
\description{
Applies a Gibbs sampler to parameters and augmented data for
two-dimensional stochastic differential equations. Currently the Ornstein-Uhlenbeck, the
stochastic FitzHugh-Nagumo model and the extended FitzHugh-Nagumo
model are implemented.
}
\usage{
Estfun(data, Delta, ImputeN = 5, seed, GibbsN = 1000,
       parKnown = list(), Start = c(0,0,0,0,1,1), diffPriorMean,
       diffPriorCovar, diffRW=diag(2), LatentPathStart,Model = NULL,
       driftPriorMean, driftPriorCovar, driftRW, LatentMeanY0 = 0,
       LatentVarY0 = 1, RWrhoPaths = 1, RWrho2PathPoints = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    Data to estimate parameters from. Matrix or numeric. Dimensions must
    be n*1 or n*2 depending on whether second coordinate is observed.
  }
  \item{Delta}{
    Positive numeric: Time between observations.
  }
  \item{ImputeN}{
    Positive integer>=3: M-2 is the number of imputed data points between consecutive observed data.
  }
  \item{seed}{
    Positive integer giving the seed for the random number generator. Defaults to random.
  }
  \item{GibbsN}{
    Positive integer: Number of iterations of the Gibbs sampler.
  }
  \item{parKnown}{
    List of named values for the known drift and diffusion parameters.
  }
  \item{Start}{
    Numerical vector with starting values for the drift and diffusion
    parameters in the Gibbs sampler.
  }
  \item{diffPriorMean}{
    numerical vector of length 2. Prior mean for diffusion coefficients
  }
  \item{diffPriorCovar}{
    2*2 matrix. Prior variance for diffusion coefficients.
  }
  \item{diffRW}{
    Random walk variance for the MH step for the diffusion ocefficients.
  }
  \item{LatentPathStart}{
    Numeric of length one or same length as Data. Starting value for the latent
    path. If LatentPathStart is a single number then all starting values
    take this value.
  }
  \item{Model}{
    Charater, specifying the model. Currently the only options are 'OU','FHN'
    and 'FHN5'.
  }
  \item{driftPriorMean}{
    prior mean for the drift parameters
  }

  \item{driftPriorCovar}{
    Prior covariance for the drift parameters
  }
  \item{driftRW}{
    Covariance matrix for the RW update of the drift parameters
  }
  \item{LatentMeanY0}{
    Prior mean for the first data point of the unobserved coordinate.
  }
  \item{LatentVarY0}{
    Prior variance for the first data point of the unobserved
    coordinate.  If 0, the point is fixed at first value of LatentPathStart.
  }
  \item{RWrhoPaths}{
    Numeric in [0,1]. Parameter for random walk update of the latent path between
    observation times. The value 0 samples a BB, the value 1 keeps the
    current value of the (skeleton) path
  }
  \item{RWrho2PathPoints}{
    Parameter for random walk update of the latent coordinate at
    observation times. The value 0 samples a middle point of a BB, the
    value 1 keeps the current value of the points
  }
}
\details{
 More details for the help page will be added soon.
 }
\value{
  An object of class BIPOD.
  \item{Drift}{
    Output of the Gibbs sampler for the drift parameters.
  }
  \item{Diff}{
    Output of the Gibbs sampler for the diffusion parameters.
  }
  \item{AccRate1}{
    Accept/reject (1/0) for each path interval and each iteration of the
    sampler.
  }
  \item{AccRate2}{
    Accept/reject (1/0) for each path endpoint of the latent coordinate
    and each iteration of the sampler. Only valid if second coordinate
    is latent.
  }
  \item{LatentPath}{
    Output of the Gibbs sampler for the endpoints of the latent
    path. Only valid when one coordinate is observed.
  }
  \item{diffAcc}{
    Accept/reject (1/0) for the MH step of the diffusion coefficient.
    }
  \item{Info}{
   List with information about the estimated model.
    }
    \item{driftPriormu}{
      Prior mean of the drift parameters.
    }
  \item{driftPriorOmega}{
      Prior variance in the drift parameters.
    }
  \item{driftRW}{
Random Walk variance for updating drift parameters.
    }
}
%references{
%% ~put references to the literature/web site here ~
%
\author{
Anders Chr. Jensen
}
%note{
%%  ~~further notes~~
%

%% ~Make other sections like Warning with \section{Warning }{....} ~

%seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%
\examples{
Data <- DiffSim(n=5000,
                start=c(0,0),
                Delta=.001,
                driftpar=c(10,5,1.5,.6),
                Sigma=diag(c(.5,.3)),
                seed=1,
                thin=100,
                Model="FHN")

A <- Estfun(data    = Data[,1],
            Delta   = .001*100,
            ImputeN = 5,
            seed    = 2,
            GibbsN  = 500,
            parKnown = list("drift3"=1.5,"drift4"=.6,"diff2"=.3),
            Start=c(10,10,10,10,1,1),
            diffPriorMean= c(0,0),
            diffPriorCovar= diag(2),
            diffRW = diag(c(.01,.02)),
            LatentPathStart = .5,
            Model="FHN",
            driftPriorMean = NULL,
            driftPriorCovar = NULL,
            driftRW = diag(4),
            LatentMeanY0 = 0,
            LatentVarY0 = 1,
            RWrhoPaths = 0,
            RWrho2PathPoints = 0)

class(A);names(A)
plot(A,type="trace",interval=1,theta=c(10,5,1.5,.6,.5,.3),subset=c(1,2,5))
### plot(A,type="movie",truepath=Data[,2],speed=.01,BY=10,interval=1)
}
%
%
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
