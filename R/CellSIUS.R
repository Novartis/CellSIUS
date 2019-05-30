#####################################################
#
# Cell subtype identification using CellSIUS
# ___________________________________________________
#
# This implements the CellSIUS method
#
# Authors:
#   Rebekka Wegmann (wegmann@imsb.biol.ethz.ch)
#   Marilisa Neri (marilisa.neri@novartis.com)
####################################################

#########################################
# Identifying rare cell types
#########################################

#' CellSIUS: Cell Subtype Identification from Upregulated gene Sets
#'
#'
#'
#' @description
#'
#' \strong{Identification of cluster-specific gene sets and cell classification according to their expression.}
#'
#' CellSIUS enables the identification and characterization
#' of (rare) cell sub-populations from complex scRNA-seq datasets: it takes as input expression values of N cells grouped into M(>1) clusters.
#' Within each cluster, genes with a bimodal distribution are selected and only genes with cluster-specific
#' expression are retained. Among these candidate marker genes, sets with correlated expression patterns
#' are identified by graph-based clustering. Finally, cells are assigned to subgroups based on their average
#' expression of each gene set. The CellSIUS algorithm output provides the rare/ sub cell types by cell indices and
#' their transcriptomic signatures.
#'
#' @section Warning:
#' The execution of CellSIUS requires
#' the external GNU Markov Cluster Algorithm (MCL) software for UNIX. It needs to be installed by the user [\href{https://micans.org/mcl/}{link}].
#'
#'
#'
#' @param mat.norm Numeric Matrix: normalized gene expression matrix where \code{rownames} and \code{colnames} are gene and cell ID, respectively.
#' @param group_id Character: vector with cluster cell assignment. Make sure that the order of cell cluster assignment reflects the order of the columns of norm.mat.
#' @param min_n_cells Integer: when identifying bimodal gene distributions, this specifies the minimum number of cells per mode. Clusters with a total number of cells below \code{min_n_cells} will be entirely ignored. Defaults to 10.
#' @param min_fc Numeric: minimum difference in mean [log2] between the two modes of the gene expression distribution. Defaults to 2.
#' @param corr_cutoff Numeric: correlation cutoff for MCL clustering of candidate marker genes. If \code{NULL}, it will be set automatically for each cluster, which in general works much better than forcing a fixed value. Therefore, leave this at the default unless having a good reason to change it. Defaults to \code{NULL}.
#' @param iter 0 or 1: relevant for the final step (assigning cells to subgroups). If set to 1, the first mode of gene expression will be discarded. Cells are then assigned based on the 2nd and 3rd mode, which results in a more stringent assignment. The default is 0 and is usually stringent enough. Defaults to 0.
#' @param max_perc_cells Numeric: maximum percentage of cells that can be part of a subcluster. Defaults to 50, implying that a “subgroup” cannot contain more than half of the total observations.
#' @param fc_between_cutoff Numeric: minimum difference [log2] in gene expression between cells in the subcluster and all other cells. The higher, the more cluster-specific is the gene signature. Note that this should not be set higher than min_fc.
#' @param mcl_path Character: path to the MCL executable. The external GNU MCL software for UNIX needs to be installed by the user.
#'
#' @details
#'
#' \code{CellSIUS} returns a \code{data.table} of a collection of key/ value pairs:
#' \itemize{
#'   \item \strong{cell_idx} Cell indices correspond to the \code{colnames} of \code{norm.mat}
#'   \item \strong{gene_id} Gene IDs correspond to the \code{rownames} of \code{norm.mat}
#'   \item \strong{main_cluster} Name of the main cluster defined in \code{group_id} input
#'   \item \strong{expr} Average gene expression of ‘gene_id’ across the cells in the correspondent ‘main_cluster’
#'   \item \strong{sub_cluster} identifies if the cell is \emph{not} a member of a CellSIUS subcluster. Subclusters are named as follows: [Name-of-main-cluster]_[\emph{number of subcluster}]_[0 or 1], where in the last position 0 means the cell is not a member, 1 means it is a member.
#'   \item \strong{N_cells} Number of cells in the correspondent CellSIUS sub cluster
#'   \item\strong{log2FC}   log2 fold change of between the two modes of the bimodal gene expression distribution across cells in the correspondent ‘main_cluster’
#'}
#'
#' @seealso
#' Accessory functions are provided to help the user to summarize and visualize the CellSIUS results:
#' \itemize{
#'   \item \code{\link[CellSIUS]{CellSIUS_GetResults}}
#'   \item \code{\link[CellSIUS]{CellSIUS_plot}}
#'   \item \code{\link[CellSIUS]{CellSIUS_final_cluster_assignment}}
#'}
#' @references
#'
#' Wegmann, R.\emph{et Al.}, 2018.
#' CellSIUS for sensitive and specific identification and characterization
#' of rare cell populations from complex single-cell RNA sequencing data.
#' \emph{Nature Comunications Submission}
#'
#' @examples
#' require(CellSIUS)
#' require(ggplot2)
#' require(data.table)
#'
#' # 2D cell projection
#' qplot(my_tsne[,1],my_tsne[,2],xlab="tSNE1",ylab="tSNE1",main='2D scRNAseq cell projection')
#' qplot(my_tsne[,1],my_tsne[,2],xlab="tSNE1",ylab="tSNE1",color=as.factor(clusters),main="Cells colored by clusters")
#'
#' # Run CellSIUS
#' # !WARNING: before to run the CellSIUS function, take care to correctly set the mcl_path to point
#' # to the executable of your local installation of MCL algorithm
#'
#' CellSIUS.out<-CellSIUS(mat.norm = norm.counts,group_id = clusters,min_n_cells=10, min_fc = 2,
#'              corr_cutoff = NULL, iter=0, max_perc_cells = 50,
#'              fc_between_cutoff = 1,mcl_path = "~/local/bin/mcl")
#'
#' #____________________________
#' # EXPLORE CellSIUS RESULTS:
#' #____________________________
#'
#' # Summary of sub-populations identified by CellSIUS
#'
#' Result_List=CellSIUS_GetResults(CellSIUS.out=CellSIUS.out)
#'
#' # 2D visualization of of sub-populations identified by CellSIUS
#' require(RColorBrewer)
#' CellSIUS_plot(coord = my_tsne,CellSIUS.out = CellSIUS.out)
#'
#' # Final clustering assignment
#'
#' Final_Clusters = CellSIUS_final_cluster_assignment(CellSIUS.out=CellSIUS.out, group_id=clusters, min_n_genes = 3)
#' table(Final_Clusters)
#' qplot(my_tsne[,1],my_tsne[,2],xlab="tSNE1",ylab="tSNE1",color=as.factor(Final_Clusters))
#'
#' @importFrom data.table ":="
#' @import Ckmeans.1d.dp
#' @import igraph
#' @export

CellSIUS = function(mat.norm,group_id,min_n_cells=10, min_fc = 2,
                         corr_cutoff = NULL, iter=0, max_perc_cells = 50,
                         fc_between_cutoff = 1,mcl_path = "~/local/bin/mcl"){
  options(warn=-1)
  if(class(mat.norm)!="matrix") stop("mat.norm is not of 'matrix' class")
  if(length(unique(group_id))<2) stop("The number of clusters in 'group_id' must be > 1")
  if(!file.exists(mcl_path)) stop("I don't find mcl executable. External UNIX installation of MCL is required: https://micans.org/mcl/")
  if(is.null(rownames(mat.norm))| is.null(colnames(mat.norm))) stop("Column and row names of mat.norm matrix must be non-emptys")
  group_id = data.table::data.table(cell_idx = colnames(mat.norm),main_cluster=as.character(group_id))
  expr_dt = data.table::data.table(gene_id = rownames(mat.norm),mat.norm)
  expr_dt_melt = data.table::melt(expr_dt,id.vars="gene_id",val="expr",var="cell_idx")
  expr_dt_melt = merge(expr_dt_melt,group_id, by="cell_idx")

  #Identify genes with significant bimodal distribution

  expr_dt_melt[,c("N_cells","within_p","pos0","pos1","Dpos"):=cellsius_find_bimodal_genes(expr,min_n_cells = min_n_cells, max_perc_cells = max_perc_cells),by=c('gene_id','main_cluster')]
  expr_dt_melt[,sig := within_p<100 & Dpos > min_fc]
  expr_dt_melt[sig==T, within_adj_p:=p.adjust(within_p),by=c('cell_idx')] #correct for multiple testing, only consider genes where test has actually been run
  expr_dt_melt[,sig:=within_adj_p<0.1]
  expr_dt_melt = expr_dt_melt[gene_id %in% expr_dt_melt[!is.na(sig) & sig==T]$gene_id]

  # If no bimodal gene were found, exit and return NA
  if(dim(expr_dt_melt)[1] == 0){
    print("No genes with bimodal distribution found, returning NA.")
    return(NA)
  }
  # Check whether these genes are specific to the subcluster

  for(clust in unique(expr_dt_melt$main_cluster)){
    expr_dt_melt = expr_dt_melt[,paste0(clust,"_",c("p_between","fc")):=cellsius_test_cluster_specificity(expr,main_cluster,clust, fc_between_cutoff = fc_between_cutoff),by="gene_id"]
    expr_dt_melt[main_cluster==clust,keep:=(expr_dt_melt[main_cluster==clust][[paste0(clust,"_p_between")]] < 0.1)]
  }

  expr_dt_melt = expr_dt_melt[keep==TRUE & !is.na(sig)]

  # If there are still non-specific genes, discard them (this can happen for
  # very high expressed genes like mitochondrial genes)
  expr_dt_melt[,n_clust_per_gene:=length(unique(main_cluster)),by='gene_id']
  expr_dt_melt = expr_dt_melt[n_clust_per_gene==1]
  expr_dt_melt[,n_clust_per_gene:=NULL]

  # Identify correlated gene sets with MCL
  expr_dt_melt = expr_dt_melt[,gene_cluster:=0]
  expr_dt_melt = cellsius_find_gene_sets(expr_dt_melt, corr_cutoff = corr_cutoff,mcl_path=mcl_path)

  # discard gene sets that only contain one gene (those are assigned to cluster 0)
  expr_dt_melt = expr_dt_melt[gene_cluster !=0 ]

  if(dim(expr_dt_melt)[1] == 0){
    print("No subclusters found, returning NA.")
    return(NA)
  }

  # Extract cell subclusters
  expr_dt_melt[,sub_cluster:=main_cluster]
  expr_dt_melt[,mean_expr := mean(expr), by = c('main_cluster','gene_cluster','cell_idx')]
  expr_dt_melt[,sub_cluster:=cellsius_sub_cluster(mean_expr,sub_cluster,gene_cluster, iter=iter),by=c('main_cluster','gene_cluster')]

  # Check how many cells belong to the subgroup relative to the total cluster size.
  # If a sub cluster contains more than max_perc_cells cells, discard it.
  clust_list = expr_dt_melt[,list(sub = length(unique(cell_idx))) ,by=c('sub_cluster','main_cluster')]
  clust_list[,tot := sum(sub)/(length(sub_cluster)/2), by= 'main_cluster']
  clust_list = clust_list[grep('_1$',sub_cluster)]
  clust_list[,perc:=sub/tot*100]
  discard_sub_clust = clust_list[perc > max_perc_cells]$sub_cluster
  discard_sub_clust = append(discard_sub_clust,gsub('_1$','_0',discard_sub_clust))

  expr_dt_melt = expr_dt_melt[!sub_cluster%in%discard_sub_clust]

  keep.columns = c("cell_idx", "gene_id", "expr", "main_cluster", "N_cells", "Dpos", "sub_cluster","gene_cluster")
  expr_dt_melt = expr_dt_melt[,..keep.columns]
  setnames(expr_dt_melt, 'Dpos', 'log2FC')
  return(expr_dt_melt)
}

#' Identify genes with bimodal distribution
#' @description
#' Accessory function of \code{\link[CellSIUS]{CellSIUS}}
#' @keywords internal

##################################################
# STEP 1: Identify genes with bimodal distribution
##################################################
cellsius_find_bimodal_genes = function(expr, min_n_cells, max_perc_cells){

  #skip genes with 0 expression
  if(sum(expr)==0){
    return(list(-1,100,-1,-1,-1))
  }
  # run k-means
  k1d = Ckmeans.1d.dp::Ckmeans.1d.dp(expr,k=2)
  # check if a cluster with more than n cells exists
  indx = which(k1d$size>min_n_cells)
  if(length(indx)>1 ){

    # do statistic only if in pos2 cells are less than max_perc_cells% of the total cells in the cluster
    if(k1d$size[2] < round(length(expr)*max_perc_cells/100)){

      t1=tryCatch({t.test(expr[which(k1d$cluster==2)],y=expr[which(k1d$cluster==1)])},
                  error = function(cond){return(0)},
                  finally={}
      )

      if(!is.numeric(t1)){

        p1=t1$p.value
        N0=k1d$size[1] # number of cells where the gene is downregulated
        N1=k1d$size[2] # number of cells  where the gene is upregulated
        pos0=k1d$centers[1]
        pos1=k1d$centers[2]
        Dpos=pos1-pos0
        return(list(N1,p1,pos0,pos1,Dpos))
      } #else {print(paste("ttest failed, dpos = ",pos1-pos0))} # for testing
    }
  }
  # if no cluster was found, return a list of dummy values
  return(list(-1,100,-1,-1,-1))
}

#' Check whether the subset of genes are specific to one cell subgroup
#' @description
#' Accessory function of \code{\link[CellSIUS]{CellSIUS}}
#' @keywords internal


#############################################################################
# STEP 2: Check whether these genes are specific to one cell subgroup
############################################################################
cellsius_test_cluster_specificity = function(exprs, cluster, current_cluster, fc_between_cutoff){

  in_clust = which(cluster == current_cluster)
  k1d = Ckmeans.1d.dp::Ckmeans.1d.dp(exprs[in_clust],k=2)
  in_subclust = in_clust[which(k1d$cluster==2)]

  mean_in = mean(exprs[in_subclust])
  mean_out = mean(exprs[-in_subclust])
  mean_out_nozero = mean(exprs[-in_subclust][exprs[-in_subclust]>0])

  # If there are subclusters, but all cells outside the subcluster express 0,
  # set mean_out_nozero to 0
  if(length(in_subclust>0) && !any(exprs[-in_subclust]>0)){mean_out_nozero=0}

  fc = mean_in - mean_out

  ts = tryCatch({t.test(exprs[in_subclust],exprs[-in_clust])},
                error = function(cond){ return(0)})

  if(!is.numeric(ts)){pv = ts$p.value} else {
    #print(paste("ttest failed, fc = ",mean_in-mean_out_nozero)) #for testing only
    pv=999}

  if(!is.nan(mean_out_nozero) && (mean_in-mean_out_nozero < fc_between_cutoff)) pv = 999
  return(list(pv,fc))
}

#' MCL clustering to find correlated gene sets
#' @description
#' Accessory function of \code{\link[CellSIUS]{CellSIUS}}
#' @keywords internal

#####################################################
# STEP 3: MCL clustering to find correlated gene sets
#####################################################
cellsius_find_gene_sets = function(expr_dt_melt, corr_cutoff = NULL, min_corr = 0.35, max_corr = 0.5,mcl_path = mcl_path){
  for(clust in unique(expr_dt_melt$main_cluster)){

    if(length(unique(expr_dt_melt[main_cluster == clust]$gene_id))==1) { next }

    mat = data.table::dcast.data.table(expr_dt_melt[main_cluster==clust], gene_id ~ cell_idx,
                                       value.var = 'expr')
    mat = mat[rowSums(mat[,-1,with=F])!=0,]
    corr.mat = cor(t(mat[,-1,with=F]))
    dimnames(corr.mat) = list(mat$gene_id,mat$gene_id)

    if(is.null(corr_cutoff)){
      corr_cutoff = max(quantile(corr.mat[corr.mat!=1],0.95),min_corr)
      corr_cutoff = min(corr_cutoff, max_corr)}
    adj.corr = corr.mat
    adj.corr[adj.corr<corr_cutoff] = 0
    adj.corr[adj.corr>=corr_cutoff] = 1
    diag(adj.corr) = 0 # no self-loop for MCL

    graphs = igraph::get.data.frame( igraph::graph_from_adjacency_matrix(adj.corr), what = "edges") # gene connection for graphs

    # if a graph has no edges (i.e. all genes are uncorrelated),
    # assign all genes to cluster "0" and go to next iteration
    if(dim(graphs)[1]==0){
      expr_dt_melt = expr_dt_melt[main_cluster == clust, gene_cluster := 0]
      next
    }

    graphs = data.frame(graphs,CORR=sapply(seq(dim(graphs)[1]), function(i) corr.mat[graphs$from[i],graphs$to[i]] -corr_cutoff))
    write.table(graphs, file = "tmp.mcl.inp",row.names=F,col.names=F,sep = " ")
    system(paste0(mcl_path, " tmp.mcl.inp --abc -o tmp.mcl.out"))
    x = scan("tmp.mcl.out", what="", sep="\n")
    y = strsplit(x, "[[:space:]]+")
    y = lapply(seq(length(y)), function(i){
      tmp = sapply(seq(length(y[[i]])),function(j){
        gsub('\"','',y[[i]][j])
      })
    })

    for(i in seq(length(y))){
      if(length(y[[i]] > 1)){
        expr_dt_melt = expr_dt_melt[main_cluster==clust & gene_id %in% y[[i]],gene_cluster:=i]
      }
    }
  }

  return(expr_dt_melt)
}

#' Assign cells to subgroups
#' @description
#' Accessory function of \code{\link[CellSIUS]{CellSIUS}}
#' @keywords internal

############################################
# Step 4: Assign cells to subgroups
############################################
cellsius_sub_cluster = function(mean_expr,sub_cluster,gene_cluster, iter = iter){

  k1d = Ckmeans.1d.dp::Ckmeans.1d.dp(mean_expr,k=2)$cluster
  cells_sub = (k1d==2)

  if(iter == 0){return(paste0(sub_cluster,"_",gene_cluster,"_",as.numeric(cells_sub)))}

  # if iter is set higher than 0, a second step of kmeans clustering.
  # This will remove the lowest peak and can sometimes help to get a more
  # accurate classification.

  k1d = Ckmeans.1d.dp::Ckmeans.1d.dp(mean_expr[cells_sub],k=2)$cluster

  if (max(k1d)>1) {
    cells_sub[cells_sub] = (k1d==2)
    return(paste0(sub_cluster,"_",gene_cluster,"_",as.numeric(cells_sub)))
  }
  return(paste0(sub_cluster,"_",gene_cluster,"_",0))
}

