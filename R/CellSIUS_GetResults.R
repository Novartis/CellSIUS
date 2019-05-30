#' Summary of CellSIUS results
#'
#'
#'
#' @description
#'
#'
#' \code{CellSIUS_GetResults} function prints a summary of the CellSIUS results.
#' It returns a nested list where each entry of this list corresponds to one subcluster.
#' For each subcluster, the cell index of its members and the ID and log2FC of its marker
#' genes are returned.
#'
#' @param CellSIUS.out data.table: Output of \code{\link[CellSIUS]{CellSIUS}}
#'
#' @seealso
#' Main CellSIUS function: \code{\link[CellSIUS]{CellSIUS}}
#' @import RColorBrewer
#' @import data.table
#' @export


CellSIUS_GetResults = function(CellSIUS.out){
   type = "summary"
  out = switch(type,
               "summary" = cellsius_print_summary(CellSIUS.out),
               "cells_per_clust" = cellsius_get_cells_per_subclust(CellSIUS.out))
}

#' Print Summary
#' @description
#' Accessory function of \code{\link[CellSIUS]{CellSIUS_GetResults}}
#' @keywords internal



cellsius_print_summary = function(CellSIUS.out){

  summary_list = list()

  cat('--------------------------------------------------------\n',
      'Summary of CellSIUS output\n',
      '--------------------------------------------------------\n\n')
  for(clust in unique(CellSIUS.out$main_cluster)){

    if(!any(CellSIUS.out[main_cluster==clust]$N_cells!=0)){
      next
    }

    cat('Main cluster: ', clust,  '\n', '---------------\n')
    subclusts = unique(CellSIUS.out[main_cluster==clust & gene_cluster!=0][order(gene_cluster)]$gene_cluster)

    for(subclust in subclusts){

      cells = unique(CellSIUS.out[main_cluster==clust & sub_cluster == paste(clust,subclust,1,sep="_")]$cell_idx)
      genes = unique(CellSIUS.out[main_cluster==clust & gene_cluster == subclust][,c("gene_id","log2FC")])
      n_cells = length(cells)
      subclust_list = list(cells = cells, genes = genes)

      cat('Subcluster: ', subclust, '\n',
          'Number of cells: ', n_cells,
          '\n Marker genes: \n')

      print(genes)
      cat('\n\n')

      subclust_name = paste(clust, subclust, sep='_')
      summary_list[[subclust_name]] = subclust_list

    }
  }
  return(summary_list)
}

#' Get cells per subcluster
#' @description
#' Accessory function of \code{\link[CellSIUS]{CellSIUS_GetResults}}
#' @keywords internal

cellsius_get_cells_per_subclust = function(CellSIUS.out){
  print('Number of cells per subcluster:')
  print(CellSIUS.out[,length(unique(cell_idx)),by="sub_cluster"])
  out = unique(CellSIUS.out[,length(unique(cell_idx)),by="sub_cluster"])
}
