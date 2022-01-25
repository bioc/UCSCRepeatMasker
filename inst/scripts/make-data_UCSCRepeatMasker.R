suppressPackageStartupMessages(library(RCurl))
suppressPackageStartupMessages(library(httr))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomeInfoDb))
suppressPackageStartupMessages(library(GenomicRanges))

options(timeout=120)

## adapted from https://stackoverflow.com/questions/39388499/r-curl-download-remote-file-only-when-changed
#' Download a file only if it has a more recent modification date than
#' \code{last_modified} and, in such a case, set the same last modification date
#' and time as file had in the server from where it has been downloaded
#'
#' @param URL url of file
#' @param fil path to write file
#' @param last_modified \code{POSIXct}.
#' @param overwrite overwrite the file if it exists?
#' @param .verbose output a message if the file was unchanged?
get_file <- function(URL, fil, last_modified=NULL, overwrite=TRUE, .verbose=TRUE) {

  if ((!file.exists(fil)) || is.null(last_modified)) {
    res <- GET(URL, write_disk(fil, overwrite))
    lastmod <- httr::parse_http_date(res$headers$`last-modified`)
    if (overwrite)
      Sys.setFileTime(fil, lastmod)
    return(lastmod)
  } else if (inherits(last_modified, "POSIXct")) {
    res <- HEAD(URL)
    cur_last_mod <- httr::parse_http_date(res$headers$`last-modified`)
    if (cur_last_mod > last_modified) {
      res <- GET(URL, write_disk(fil, overwrite))
      lastmod <- httr::parse_http_date(res$headers$`last-modified`)
      if (overwrite)
        Sys.setFileTime(fil, lastmod)
      return(lastmod)
    }
    if (.verbose) message(sprintf("'%s' no more recent than %s", URL, last_modified))
    return(last_modified)
  } 
}

citdata <- bibentry(bibtype="Manual",
                    author=c(person("A. Smit"), person("R. Hubley"),
                             person("P. Green")),
                    title="RepeatMasker Open-3.0",
                    year="1996-2010",
                    url="https://www.repeatmasker.org")

#' Process a rmsk.txt.gz file from UCSC and store its
#' data into a GRanges object
#'
#' @param fname path to the rmsk.txt.gz file
processRMSKfile <- function(fname) {
  message(sprintf("processing %s\n",fname))
  rmsktbl <- read.table(gzfile(file.path(g, "rmsk.txt.gz")), header=FALSE,
                        sep="\t", colClasses=c("NULL",      ## bin
                                               "integer",   ## swScore
                                               "NULL",      ## milliDiv
                                               "NULL",      ## milliDel
                                               "NULL",      ## milliIns
                                               "character", ## genoName
                                               "integer",   ## genoStart
                                               "integer",   ## genoEnd
                                               "NULL",      ## genoLeft
                                               "character", ## strand
                                               "character", ## repName
                                               "character", ## repClass
                                               "character", ## repFamily
                                               "NULL",      ## repStart
                                               "NULL",      ## repEnd
                                               "NULL"       ## id
                                               ))
  colnames(rmsktbl) <- c("swScore", "genoName", "genoStart", "genoEnd",
                         "strand", "repName", "repClass", "repFamily")
  rmskGR <- GRanges(seqnames=rmsktbl$genoName,
                    ranges=IRanges(rmsktbl$genoStart+1, rmsktbl$genoEnd),
                    strand=rmsktbl$strand,
                    swScore=rmsktbl$swScore,
                    repName=rmsktbl$repName,
                    repClass=rmsktbl$repClass,
                    repFamily=rmsktbl$repFamily)
  si <- getChromInfoFromUCSC(g, as.Seqinfo=TRUE)
  seqlevels(rmskGR) <- seqlevels(si)
  seqinfo(rmskGR) <- si
  rmskGR <- sort(rmskGR)
  rmskGR
}

## fetch all current genome versions at UCSC
allg <- ucscGenomes(TRUE)

## for some strange reason the release date for hg19,
## hg38, mm10 and mm39 is missing
mt <- match(c("hg19", "hg38", "mm10", "mm39"), allg$db)
stopifnot(all(!is.na(mt))) ## QC
allg$date[mt] <- c("Feb. 2009", "Jun. 2013", "Dec. 2011", "Jun. 2020")
## fix release date for wuhCor1
mt <- match("wuhCor1", allg$db)
allg$date[mt] <- "Jan. 2020"
## fix release date for eboVir3
mt <- match("eboVir3", allg$db)
allg$date[mt] <- "Jun. 2014"

## split genome versions by organism
gbyo <- split(allg$db, allg$organism)

## order genome version within organism
gbyo <- lapply(gbyo, function(g) {
                 mt <- gregexpr("[0-9]+", g)
                 vst <- unlist(mt, use.names=FALSE)
                 vw <- sapply(mt, attr, "match.length")
                 vn <- as.integer(substr(g, vst, vst+vw-1))
                 g[order(vn, decreasing=TRUE)]
               })

twoversionsgenomes <- c("Homo sapiens", "Mus musculus")
for (i in 1:length(gbyo)) {
  maxg <- 1
  if (names(gbyo)[i] %in% twoversionsgenomes)
    maxg <- 2

  theseg <- gbyo[[i]][1:min(c(length(gbyo[[i]]), maxg))]
  for (g in theseg) {
    rmskurl <- sprintf("https://hgdownload.soe.ucsc.edu/goldenPath/%s/database/rmsk.txt.gz", g)
    if (url.exists(rmskurl)) {
      message(sprintf("create %s\n", g))
      dir.create(g, showWarnings=FALSE)
      dldexit <- 0
      if (!file.exists(file.path(g, basename(rmskurl)))) {
        moddate <- get_file(rmskurl, file.path(g, basename(rmskurl)))
        if (is.null(moddate)) {
          unlink(file.path(g, basename(rmskurl)))
          dldexit <- 1
        }
      } else {
        moddate <- file.info(file.path(g, "rmsk.txt.gz"))$mtime
        moddate <- get_file(rmskurl, file.path(g, basename(rmskurl)), moddate)
        if (is.null(moddate)) {
          unlink(file.path(g, basename(rmskurl)))
          dldexit <- 1
        }
      }
      if (dldexit != 0) {
        message(sprintf("couldn't download rmsk.txt.gz for %s\n", g))
      } else {
        rmskGRfname <- file.path(g, sprintf("rmsk.%s.%s.rds", g,
                                            format(moddate, "%b%Y")))
        if (!file.exists(rmskGRfname)) {
          rmskGR <- processRMSKfile(file.path(g, basename(rmskurl)))
          mt <- match(g, allg$db)
          refgenomeGD <- GenomeDescription(organism=allg[mt, "organism"],
                                           common_name=allg[mt, "species"],
                                           provider="UCSC",
                                           provider_version=g,
                                           release_date=allg[mt, "date"],
                                           release_name=g,
                                           seqinfo=seqinfo(rmskGR))
          metadata(rmskGR) <- list(srcurl=sprintf("https://hgdownload.soe.ucsc.edu/goldenPath/%s/database/rmsk.txt.gz", g),
                                   srcVersion=format(moddate, "%b%Y"),
                                   citation=citdata,
                                   gdesc=refgenomeGD)
          saveRDS(rmskGR, file.path(g, sprintf("rmsk.%s.%s.rds", g,
                                               format(moddate, "%b%Y"))))
          cat(sprintf("%s\n", rmskGRfname))
        } else
          message(sprintf("rmsk.txt.gz for %s already processed, skipping it\n", g))
      }
    } else
      message(sprintf("couldn't download rmsk.txt.gz for %s\n", g))
  }
}
