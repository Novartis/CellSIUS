#' 2D visualization of CellSIUS sub clusters
#'
#'
#'
#' @description
#'
#'
#' \code{CellSIUS_plot} function returns a ggplot visualization of the 2D distribution
#' of cells highlighting the sub clusters identified by CellSIUS.
#'
#' @param CellSIUS.out data.table: Output of \code{\link[CellSIUS]{CellSIUS}}
#' @param coord matrix: 2D cell distribution coordinates, \emph{e.g.} tSNE or PCA,
#' where its \code{rownames} correspond to cell indices.
#'
#' @seealso
#' Main CellSIUS function: \code{\link[CellSIUS]{CellSIUS}}
#' @import RColorBrewer
#' @import data.table
#' @import ggplot2
#' @export


CellSIUS_plot = function(coord,CellSIUS.out){

  coord_dt = data.table(coord1 = coord[,1], coord2 = coord[,2], cell_idx = rownames(coord))
  coord_dt = merge(coord_dt, CellSIUS.out[,c('cell_idx','main_cluster','sub_cluster')],
                  by = c('cell_idx'), all = T)
  coord_dt[is.na(main_cluster),main_cluster:='Other']
  coord_dt[main_cluster == 'Other',sub_cluster:='none']
  coord_dt[grepl('_0$',sub_cluster),sub_cluster:= 'none']

  setkey(coord_dt, 'cell_idx')
  coord_dt = unique(coord_dt)

  rc_cols = brewer.pal(10,"Spectral")[rep(c(1,9,7,2,6,10,3,8),3)]

  p = ggplot(coord_dt, aes(x = coord1, y= coord2)) +
    geom_point(color = "darkgray", alpha = 0.5, size = 1.5)+
    theme_bw() + theme(text = element_text(size = 15))
  p = p + geom_point(data = coord_dt[sub_cluster!='none'], aes(x=coord1, y=coord2, color = sub_cluster))+
    scale_color_manual(values = rc_cols) + guides(color = guide_legend(title = 'Subcluster'))
  return(p)
}
