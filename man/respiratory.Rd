\name{respiratory}
\alias{respiratory}
\non_function{}
\title{ Respiratory Health Data Set }
\usage{data(respiratory)}
\description{
  Sub-group and whole-group analysis of a respiratory therapy drug; data set taken from Westfall et al. (1999, p. 195). Three-way ANOVA with weighted contrasts.
}
\format{
  This data frame contains the following variables
  \describe{
    \item{Treatment}{Treatment at 2 levels: \code{Active} and
\code{Placebo}.}
    \item{AgeGroup}{Age group at 2 levels: \code{Younger} and \code{Older}.}
    \item{InitHealth}{Inititial health condition at 2 levels: \code{Good}
                      and \code{Poor}.}
    \item{Score}{Response variable: time-weighted score comprised 
                 of four post-treatment evaluations.}
  }
}
\details{
  See Westfall et al. (1999, p. 195)
}
\source{
  P. H. Westfall, R. D. Tobias, D. Rom, R. D. Wolfinger, Y. Hochberg (1999).
  \emph{Multiple Comparisons and Multiple Tests Using the SAS System}.
  Cary, NC: SAS Institute Inc.  
}
\examples{
data(respiratory)

# compute the contrast matrix in several steps 
# overall active vs. placebo
CA  <- c(0, 13, 0, 11, 0, 13, 0, 17, 0)
CP  <- c(0, 0, 14, 0, 12,  0, 19, 0, 12)
CA  <- CA/sum(CA)
CP  <- CP/sum(CP)
C1  <- CP-CA

# for older subgroup only
CAO <- c( 0, 13,  0,  0,  0, 13,  0,  0,  0 ) 
CPO <- c( 0,  0, 14,  0,  0,  0, 19,  0,  0 ) 
CAO <- CAO/sum(CAO)
CPO <- CPO/sum(CPO)
C2  <- CPO - CAO

# for younger subgroup only 
CAY <- c( 0,  0,  0, 11,  0,  0,  0, 17,  0 ) 
CPY <- c( 0,  0,  0,  0, 12,  0,  0,  0, 12 ) 
CAY <- CAY/sum(CAY)
CPY <- CPY/sum(CPY)
C3  <- CPY - CAY

# subgroup with inital good health
CAG <- c( 0, 13,  0, 11,  0,  0,  0,  0,  0 ) 
CPG <- c( 0,  0, 14,  0, 12,  0,  0,  0,  0 ) 
CAG <- CAG/sum(CAG)
CPG <- CPG/sum(CPG)
C4  <- CPG - CAG

# subgroup with inital poor health
CAP <- c( 0,  0,  0,  0,  0, 13,  0, 17,  0 ) 
CPP <- c( 0,  0,  0,  0,  0,  0, 19,  0, 12 ) 
CAP <- CAP/sum(CAP)
CPP <- CPP/sum(CPP)
C5  <- CPP - CAP

# all 4 subgroup combinations of age and initial health condition 
C6  <- c( 0, -1,  1,  0,  0,  0,  0,  0,  0 )
C7  <- c( 0,  0,  0,  0,  0, -1,  1,  0,  0 )
C8  <- c( 0,  0,  0, -1,  1,  0,  0,  0,  0 )
C9  <- c( 0,  0,  0,  0,  0,  0,  0, -1,  1 )

C   <- rbind(C1, C2, C3, C4, C5, C6, C7, C8, C9)   
# numerate the contrasts
colnames(C) <- NULL
rownames(C) <- c("Overall", "Older", "Younger", "Good Init", "Poor Init",
                 "Old x Good", "Old x Poor", "Young x Good", "Young x Poor") 

summary(simtest(Score ~ Treatment:AgeGroup:InitHealth,
                data=respiratory, ttype="logical",
                alternative="greater", cmatrix=C))
}
\keyword{datasets}
