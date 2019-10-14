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

#' STPlotGene
#'
#' Lookup gene data and create 2D heatmap of expression
#'
#' @param x An ST data object such as returned from STLoad
#' @param gene gene identifier. Identifier type is identified
#'             automatically.
#'
#' @import ggplot2
#' @importFrom foreach %do% foreach
#' @export
#'
STPlotGene <- function(x, genes, size = 1) {
  ix <- foreach::foreach(g = genes, .combine=c) %do% {
    if(g %in% colnames(x$exp)) {
      ix <- which(colnames(x$exp) == g)[1]
    } else if(g %in% x$symbol) {
      ix <- which(x$symbol == g)[1]
    } else if(g %in% x$entrez) {
      ix <- which(x$entrez == g)[1]
    } else {
      ix <- NA
    }
    ix
  }
  if(length(which(is.na(ix))) == length(ix)) {
    stop("Did not identify any matching gene identifiers.")
  }
  genes <- genes[which(!is.na(ix))]
  ix <- ix[which(!is.na(ix))]
  gd <- as.data.frame(x$exp[,ix])
  gd$x <- x$x
  gd$y <- x$y
  gd <- reshape2::melt(gd, id.vars=c("x", "y"))
  gd$variable <- factor(gd$variable, labels=genes)
  ggplot2::ggplot(gd, aes(x=x, y=y, color=value)) +
    geom_point(size = size) +
    scale_color_gradientn(colors=c("#3488C1", "#F4EA93", "#DC4052"), na.value="transparent") +
    facet_wrap(~variable, scales = "free") +
    theme_bw()
}

