############################################
# OPTIONAL: Final assignment to unique clusters
# Note: This is different from the previous subcluster asignemnt, where a cell can potentially be
# a member of multiple subgroups.
############################################

#' Final assignment to unique clusters
#'
#'
#'
#' @description
#'
#'
#' \code{CellSIUS_final_cluster_assignment} merges the main clusters with the sub clusters
#' identified by \code{\link[CellSIUS]{CellSIUS}} and returns the final cluster assignment.
#'
#'
#' @param CellSIUS.out data.table: Output of \code{\link[CellSIUS]{CellSIUS}}.
#' @param group_id Character: vector with cluster cell assignment.
#' Make sure that the order of cell cluster
#' assignment reflects the order of the columns of norm.mat.
#' @param min_n_genes Integer: Minimum number of genes in a signature.
#'
#' @details
#'
#' \code{CellSIUS_final_cluster_assignment} returns a \code{character}.
#' Note that a cell can potentially be a member of multiple subgroups.
#'
#' @seealso
#'
#' Main CellSIUS function: \code{\link[CellSIUS]{CellSIUS}}.
#' @import data.table
#' @importFrom data.table ":="
#' @export

CellSIUS_final_cluster_assignment = function(CellSIUS.out, group_id, min_n_genes = 3){

  CellSIUS.out[,n_genes:=length(unique(gene_id)),by='sub_cluster']

  assignments = data.table::data.table(cell_idx = names(group_id), group=as.character(group_id))
  names(assignments) = c('cell_idx', 'group')
  assignments$group = as.character(assignments$group)
  assignments = merge(assignments, CellSIUS.out[n_genes>=min_n_genes,c('cell_idx','main_cluster','sub_cluster')],by='cell_idx',all=T)
  assignments = unique(assignments)

  final_assignment = function(main_cluster,sub_cluster){

    if(length(sub_cluster)==1){
      if(is.na(sub_cluster) || grepl("0$",sub_cluster)){
        out = main_cluster
      } else {
        out = gsub('_\\d$','',sub_cluster)
      }
    } else {
      subclusts = gsub('_\\d$', '',sub_cluster[grepl("1$",sub_cluster)])
      out = paste(subclusts,collapse='-')
      if(out == ''){out = main_cluster}
    }
    return(out)
  }

  assignments[,final:=final_assignment(group,sub_cluster),by="cell_idx"]
  assignments = unique(assignments[,c('cell_idx','final')])
  out = as.character(assignments$final)
  names(out)=assignments$cell_idx
  out = out[names(group_id)]
  return(out)
}
