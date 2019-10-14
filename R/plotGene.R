#' STPlotGene
#'
#' Lookup gene data and create 2D heatmap of expression
#'
#' @param x An ST data object such as returned from STLoad
#' @param genes gene identifier(s). Identifier type is identified
#'             automatically.
#'
#' @param size Size of data in plot
#' @param normalize Normalize expression to 0:1?
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom foreach %do% foreach
#' @export
#'
STPlotGene <- function(x, genes, size = 1, normalize = FALSE) {
  g <- NULL # avoid no visible binding warning on check
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
  if(normalize) {
    gd <- as.data.frame(apply(gd, 2, function(x) { x <- x-min(x); x / max(x)}))
    legname = "Normalized \nCounts"
  } else {
    legname = "Counts"
  }
  gd$x <- x$x
  gd$y <- x$y
  gd <- melt(gd, id.vars=c("x", "y"))
  gd$variable <- factor(gd$variable, labels=genes)
  ggplot(gd, aes(x=gd$x, y=gd$y, color=gd$value)) +
    geom_point(size = size ,) +
    scale_color_gradientn(colors=c("#3488C1", "#F4EA93", "#DC4052"),
                          na.value="transparent",
                          name = legname) +
    facet_wrap(~variable, scales = "free") +
    xlab("X") +
    ylab("Y") +
    theme_bw()
}
