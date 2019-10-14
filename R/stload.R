#' STLoad
#'
#' Load ST data file, extract coordinates, and annotate genes
#'
#' @param f path to ST data file
#' @param species either "hs" or "mm", used for gene annotation
#' @import EnsDb.Hsapiens.v86
#' @import EnsDb.Mmusculus.v79
#' @importFrom foreach %do%
#' @importFrom AnnotationDbi select
#' @importFrom utils read.delim
#' @export
#'
STLoad <- function(f, species = c("hs", "mm")) {
  species <- match.arg(species)
  dat <- read.delim(f, sep="\t",
                    header = TRUE,
                    as.is = TRUE)
  x <- as.numeric(gsub("x.*", "", dat[,1]))
  y <- as.numeric(gsub(".*x", "", dat[,1]))
  if(species == "hs") {
    edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  } else {
    edb <- EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79
  }
  dat <- dat[,-1]
  gg <- AnnotationDbi::select(edb, keys = colnames(dat), keytype = "GENEID",
                 columns = c("ENTREZID", "SYMBOL"))
  ix <- match(colnames(dat), gg$GENEID)
  if(sum(is.na(gg$ENTREZID)) > .5 * length(ix))
    warning("Less than 50% of genes annotated...is species correct?")
  return(list(x = x, y = y, exp = dat, entrez = gg$ENTREZ[ix], symbol = gg$SYMBOL[ix]))
}


