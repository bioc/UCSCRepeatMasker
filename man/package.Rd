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
\dontrun{
library(AnnotationHub)

ah <- AnnotationHub()
rmsk <- query(ah, "RepeatMasker")
rmsk
}
}

\keyword{package}
\keyword{data}
