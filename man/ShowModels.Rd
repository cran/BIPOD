\name{ShowModels}
\alias{ShowModels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prints form of supported stochastic differential equations}
\description{
Print function displaying the model structures currently
supported. Used to fix the parametrization of the parameters.
}
\usage{
ShowModels(Model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Model}{
Character specifying the model. Current options are "OU" for the
Ornstein Uhlenbeck process, "FHN" for the stochastic FitzHugh-Nagumo
process "FHN5" for the extended FitzHugh-Nagumo model and "CIR" for the
Cox-Ingersoll-Ross model.
}
}
\details{
This function is used to identity the parameter names in the supported
models. This is necessary when specifying the "parKnown" argument in the
"Estfun" function.
}
\value{
  List with three entries:
  \item{Model}{
    A 2*1 Matrix with character entries.
  }
  \item{Ndrift}{
    Numeric giving the number of drift parameters
  }
  \item{Ndiff}{
    Numeric giving the number of diffusion parameters
    }
}
\examples{
ShowModels(Model="FHN")
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
