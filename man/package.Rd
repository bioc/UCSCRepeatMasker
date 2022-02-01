\name{UCSCRepeatMasker}
\docType{package}

\alias{UCSCRepeatMasker-package}
\alias{UCSCRepeatMasker}

\title{UCSC RepeatMasker annotations}

\description{
  This package stores AnnotationHub resource metadata for
  RepeatMasker annotations from UCSC.
}

\note{
  This data package was made from resources at UCSC
  based on version 20Jan2022.
}

\author{Robert Castelo <robert.castelo@upf.edu> [cre]}

\seealso{
  \code{\link[AnnotationHub]{AnnotationHub}}
}

\examples{
library(AnnotationHub)

ah <- AnnotationHub()
query(ah, "RepeatMasker")
}

\keyword{package}
\keyword{data}
