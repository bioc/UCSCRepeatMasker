## adapted from AnnotationForge:::.getSubDirs()
#' @importFrom RCurl getURL
#' @importFrom XML xmlGetAttr htmlTreeParse
#' @import GenomeInfoDb
#' @import AnnotationHub
.getSubDirs <- function(dname)
{
    getLinks <- function() {
        links <- character(0)
        list(a = function(node, ...) {
                   links <<- c(links, xmlGetAttr(node, "href"))
                   node
                 },
             links = function() links)
    }
    h1 <- getLinks()
    htmlContent <- getURL(dname) ## required for https protocol
    htmlTreeParse(htmlContent, handlers = h1)
    res <- h1$links()
    res <- res[!(res %in% c("?C=N;O=D", "?C=M;O=A", "?C=S;O=A", "?C=D;O=A",
                            "/download/current/", "/annotationhub/"))]
    res
}
