\name{angina}
\alias{angina}
\non_function{}
\title{ Dose Response Data Set }
\usage{data(angina)}
\description{
 Dose response study of an angina drug; data set taken from 
 Westfall et al. (1999, p. 164). Different trend tests in a balanced 
 one-way layout.
}
\format{
  This data frame contains the following variables
  \describe{
    \item{dose}{Treatments at 5 levels: \code{0}, \code{1}, \code{2},
                \code{3}, \code{4}}
    \item{response}{Response variable: change from pretreatment 
                    as measured in minutes of pain-free walking.}
  }
}
\details{
  See Westfall et al. (1999, p. 164)
}
\source{
  P. H. Westfall, R. D. Tobias, D. Rom, R. D. Wolfinger, Y. Hochberg (1999).
  \emph{Multiple Comparisons and Multiple Tests Using the SAS System}.
  Cary, NC: SAS Institute Inc.  
}
\examples{
data(angina)

# perform a dose-response analysis using simultaneous confidence 
# intervals for Williams' contrasts
summary(simint(response~dose, data=angina, alternative="greater",
               type="Williams"))

# consider the same data set as a changepoint problem
summary(simint(response~dose, data=angina, alternative="greater",
               type="Changepoint"))

# compute now adjusted p-values for McDermott's test on trend
summary(simtest(response~dose, data=angina, type="McDermott",
                alternative="greater",ttype="logical"))
}
\keyword{datasets}
