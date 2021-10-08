library(ggplot2)


#' Plot cells on UMAP projection
#'
#' @param obs (data.frame) cell metadata.
#' @param x (str) variable to plot along x-axis.
#' @param y (str) variable to plot along y-axis.
#' @param color (str) variable to color individual cells by.
#'
#' @return (ggplot) scatter plot where each cell is a dot.
#'
#' @examples
#' p <- plot_cells(df, 'umap1', 'umap2', 'cell.type')
plot_cells <- function(obs, x, y, color) {
    obs[ , color] <- as.factor(obs[ ,color])
    p <- ggplot(obs, aes_string(x=x, y=y, color=color)) +
        geom_point()
    return(p)
}

if (exists('snakemake')) {
    # read in count matrix
    X <- read.csv(, header=FALSE, row.names=NULL)
    # read in cell metadata
    obs <- read.csv(, header=TRUE, row.names=1)
    # read in gene metadata
    var <- read.csv(, header=TRUE, row.names=1)
    # plot cells on UMAP projection
    p <- plot_cells(, x='umap1', y='umap2',)
    ggsave(, p)
}