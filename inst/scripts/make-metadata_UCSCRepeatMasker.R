.UCSCRepeatMaskerMetadataFromUrl <- function(baseUrl) {
  ## genus+species to NCBI taxon identifier from GenomeInfoDb
  load(system.file("extdata", "assembly_accessions.rda", package="GenomeInfoDb"))
  assembly_accessions <- get("assembly_accessions") ## just to avoid a NOTE during R CMD check

  genomes <- .getSubDirs(baseUrl)
  genomes <- gsub("/", "", genomes)
  metadf <- data.frame(title=character(0),
                       speciies=character(0),
                       taxonomyId=integer(0),
                       genome=character(0),
                       sourceUrl=character(0),
                       description=character(0),
                       rDataPath=character(0))
  for (g in genomes) {
    rmskobjs <- .getSubDirs(paste0(baseUrl, g, "/"))
    rmskobjs <- rmskobjs[grep(g, rmskobjs)]
    for (objfname in rmskobjs) {
      obj <- readRDS(gzcon(url(sprintf("%s%s/%s", baseUrl, g, objfname), open="rb")))
      gd <- metadata(obj)$gdesc
      taxonId <- as.integer(subset(assembly_accessions,
                                   organism_name == organism(gd))$taxid)[1]
      rDataPath <- sprintf("%s/%s", g, objfname)
      tmpdf <- data.frame(title=sprintf("UCSC RepeatMasker annotations (%s) for %s (%s)",
                                        metadata(obj)$srcVersion, commonName(gd), providerVersion(gd)),
                          species=organism(gd),
                          taxonomyId=taxonId,
                          genome=providerVersion(gd),
                          sourceUrl=sprintf("%s%s/%s", baseUrl, g, objfname),
                          sourceVersion=metadata(obj)$srcVersion,
                          description=sprintf("UCSC RepeatMasker annotations (%s) for %s -%s- (%s)",
                                              metadata(obj)$srcVersion, commonName(gd), organism(gd),
                                              providerVersion(gd)),
                          rDataPath=rDataPath)
      metadf <- rbind(metadf, tmpdf)
    }
  }
  rownames(metadf) <- NULL
  metadf
}

makeMetadata_UCSCRepeatMasker <- function()
{
  biocver <- "3.15"
  baseUrl <- "https://functionalgenomics.upf.edu/annotationhub/repeatmasker/"
  meta <- .UCSCRepeatMaskerMetadataFromUrl(baseUrl)
  n <- nrow(meta)
  data.frame(
    BiocVersion=rep(biocver, n),
    Description=meta$description,
    Genome=meta$genome,
    SourceUrl=meta$sourceUrl,
    SourceType=rep("UCSC track", n),
    SourceVersion=meta$sourceVersion,
    Species=meta$species,
    TaxonomyId=meta$taxonomyId,
    Title=meta$title,
    ResourceName=meta$title,
    RDataPath=meta$rDataPath,
    Coordinate_1_based=rep(TRUE, n),
    DataProvider=rep("UCSC", n),
    Maintainer=rep("Robert Castelo <robert.castelo@upf.edu>", n),
    RDataClass=rep("GRanges", n),
    DispatchClass=rep("RDS", n),
    Location_Prefix=baseUrl,
    Tags=rep(paste("RepeatMasker", "UCSC", sep=":"), n))
}

library(XML)
library(RCurl)
library(GenomeInfoDb)

source("../../R/utils.R") ## for .getSubDirs()

metadata <- makeMetadata_UCSCRepeatMasker()
write.csv(metadata, file="../extdata/metadata_UCSCRepeatMasker.csv", row.names=FALSE)
