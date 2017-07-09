### R code from vignette source 'monocle-vignette.Rnw'

###################################################
### code chunk number 1: package_loads
###################################################
library(Biobase)
library(knitr)
library(reshape2)
library(ggplot2)

knitr::opts_chunk$set(autodep=TRUE, cache=FALSE, warning=FALSE)
set.seed(0)


###################################################
### code chunk number 2: init_monocle
###################################################
library(HSMMSingleCell)
library(monocle)
data(HSMM_expr_matrix)
data(HSMM_gene_annotation)
data(HSMM_sample_sheet)


###################################################
### code chunk number 3: cluster_cells_unsup_gene_pick
###################################################
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)


###################################################
### code chunk number 4: cluster_cells_unsup_no_covariate
###################################################
# HSMM@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(HSMM, return_all = F) # norm_method = 'log',
HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 6, 
                        reduction_method = 'tSNE', verbose = T) 
HSMM <- clusterCells(HSMM,
                     num_clusters=2)
plot_cell_clusters(HSMM, 1, 2, color="CellType", markers=c("MYF5", "ANPEP"))


###################################################
### code chunk number 5: cluster_cells_unsup_plot_by_media
###################################################
plot_cell_clusters(HSMM, 1, 2, color="Media")


###################################################
### code chunk number 6: cluster_cells_unsup_control_for_media
###################################################
HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 2, reduction_method = 'tSNE', 
                       residualModelFormulaStr="~Media + num_genes_expressed", verbose = T) #
HSMM <- clusterCells(HSMM, num_clusters=2)
plot_cell_clusters(HSMM, 1, 2, color="CellType")


###################################################
### code chunk number 7: cluster_cells_unsup_plot_by_cell_type
###################################################
HSMM <- clusterCells(HSMM, num_clusters=2)
plot_cell_clusters(HSMM, 1, 2, color="Cluster") + facet_wrap(~CellType)


###################################################
### code chunk number 8: cluster_cells_diff_table
###################################################
marker_diff <- markerDiffTable(HSMM[expressed_genes,], 
                                 cth, 
                                 residualModelFormulaStr="~Media + num_genes_expressed",
                                 cores=1)


###################################################
### code chunk number 9: cluster_cells_semisup_show_marker_spec
###################################################
candidate_clustering_genes <- row.names(subset(marker_diff, qval < 0.01))
marker_spec <- calculateMarkerSpecificity(HSMM[candidate_clustering_genes,], cth)
head(selectTopMarkers(marker_spec, 3))


###################################################
### code chunk number 10: cluster_cells_semisup_pick_genes
###################################################
semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)
HSMM <- setOrderingFilter(HSMM, semisup_clustering_genes)
plot_ordering_genes(HSMM)


###################################################
### code chunk number 11: cluster_cells_semisup_clustering_no_impute
###################################################
plot_pc_variance_explained(HSMM, return_all = F) # norm_method = 'log',
HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 2, reduction_method = 'tSNE', 
                       residualModelFormulaStr="~Media + num_genes_expressed", verbose = T) 
HSMM <- clusterCells(HSMM, num_clusters=2) 

plot_cell_clusters(HSMM, 1, 2, color="CellType")
 @ %def
 
 \subsection{Imputing cell type}
 
 Note that we've reduce the number of ``contaminating'' fibroblasts in the myoblast cluster, and vice versa. But what about the ``Unknown'' cells? If you provide \Rfunction{clusterCells} with a the CellTypeHierarcy, Monocle will use it classify \emph{whole clusters}, rather than just individual cells. Essentially, cluserCells works exactly as before, except after the clusters are built, it counts the frequency of each cell type in each cluster. When a cluster is composed of more than a certain percentage (in this case, 10\%) of a certain type, all the cells in the cluster are set to that type. If a cluster is composed of more than one cell type, the whole thing is marked ``Ambiguous''. If there's no cell type thats above the threshold, the cluster is marked ``Unknown''. Thus, Monocle helps you impute the type of each cell even in the presence of missing marker data.
 
 <<cluster_cells_semisup_clustering_with_impute, fig.width = 7, fig.height = 4, fig.align="center", eval=TRUE>>=
HSMM <- clusterCells(HSMM,
                      num_clusters=2, 
                      frequency_thresh=0.1,
                      cell_type_hierarchy=cth)
 plot_cell_clusters(HSMM, 1, 2, color="CellType", markers = c("MYF5", "ANPEP"))
 @ %def
 
As you can see, the clusters are fairly pure in terms of \emph{MYF5} expression. There are some cells expressing \emph{ANPEP} in both clusters, but those in the myoblast cluster also express \emph{MYF5}. This is not surprising, as ANPEP isn't a very specific marker of fibroblasts. Overall, we've successfully classified all the cells:
 
 <<count_cells_semisup_pie, fig.width = 5, fig.height = 5, fig.align="center", eval=TRUE>>=
 pie <- ggplot(pData(HSMM), aes(x = factor(1), fill = factor(CellType))) +
   geom_bar(width = 1)
 pie + coord_polar(theta = "y") + 
   theme(axis.title.x=element_blank(), axis.title.y=element_blank())
 @ %def
 
 Finally, we subset the CellDataSet object to create \Robject{HSMM\_myo}, which includes only myoblasts. We'll use this in the rest of the analysis.
 <<select_myoblasts, eval=TRUE>>=
 HSMM_myo <- HSMM[,pData(HSMM)$CellType == "Myoblast"]	
 HSMM_myo <- estimateDispersions(HSMM_myo)
 @ %def
 
 \section{Constructing single cell trajectories}
 
 During development, in response to stimuli, and througout life, cells transition from one functional ``state'' to another. Cells in different states express different sets of genes, producing a dynamic repetoire of proteins and metabolites that carry out their work. As cells move between states, undergo a process of transcriptional re-configuration, with some genes being silenced and others newly activated. These transient states are often hard to characterize because purifying cells in between more stable endpoint states can be difficult or impossible. Single-cell RNA-Seq can enable you to see these states without the need for purification. However, to do so, we must determine where each cell is the range of possible states. 
 
 Monocle introduced the strategy of using RNA-Seq for \emph{single cell trajectory analysis}. Rather than purifying cells into discrete states experimentally, Monocle uses an algorithm to learn the sequence of gene expression changes each cell must go through as part of a dynamic biological process. Once it has learned the overall ``trajectory'' of gene expression changes, Monocle can place each cell at its proper position in the trajectory. You can then use Monocle's differential analysis toolkit to find genes regulated over the course of the trajectory, as described in section \ref{pseudotime_diff}. If there are multiple outcome for the process, Monocle will reconstruct a ``branched'' trajectory. These branches correspond to cellular ``decisions'', and Monocle provides powerful tools for identifying the genes affected by them and involved in making them. You can see how to analyze branches in section \ref{BEAM}. Monocle relies on a machine learning technique called \emph{reversed graph embedding} to construct single-cell trajectories. You can read more about the theoretical foundations of Monocle's approach in section \ref{theory}, or consult the references shown below in section \ref{references}.
 
 \subsection{``Pseudotime'': a measure of progress through a biological process }
 
 In many biological processes, cells do not progress in perfect synchrony.  In single-cell expression studies of processes such as cell differentiation, captured cells might be widely distributed in terms of progress.  That is, in a population of cells captured at exactly the same time, some cells might be far along, while others might not yet even have begun the process.  This asynchrony creates major problems when you want to understand the sequence of regulatory changes that occur as cells transition from one state to the next. Tracking the expression across cells captured at the same time produces a very compressed sense of a gene's kinetics, and the apparent variability of that gene's expression will be very high. 
 
 By ordering each cell according to its progress along a learned trajectory, Monocle alleviates the problems that arise due to asynchrony. Instead of tracking changes in expression as a function of time, Monocle tracks changes as a function of progress along the trajectory, which we term ``pseudotime''. Pseudotime is an abstract unit of progress: it's simply the distance between a cell and the start of the trajectory, measured along the shortest path. The trajectory's total length is defined in terms of the total amount of transcriptional change that a cell undergoes as it moves from the starting state to the end state. For further details, see section \ref{theory}.
 
 \subsection{The ordering algorithm}
 
 \subsubsection{Choosing genes for ordering}
 Inferring a single-cell trajectory is a hard machine learning problem. The first step is to select the genes Monocle will use as input for its machine learning approach. This is called \emph{feature selection}, and it has a major impact in the shape of the trajectory. In single-cell RNA-Seq, genes expressed at low levels are often very noisy, but some may contain important information regarding the state of the cell.  If we simply provide all the input data, Monocle might get confused by the noise, or fix on a feature of the data that isn't biologically meaningful, such as batch effects arising from collecting data on different days. Monocle provides you with a variety of tools to select genes that will yield a robust, accurate, and biologically meaningful trajectory. You can use these tools to either perform a completely ``unsupervised'' analysis, in which Monocle has no forehand knowledge of which gene you consider important. Alternatively, you can make use of expert knowledge in the form of genes that are already known to define biolgical progress to shape Monocle's trajectory. We consider this mode ``semi-supervised'', because Monocle will augment the markers you provide with other, related genes. Monocle then uses these genes to produce trajectories consistent with known biology but that often reveal new regulatory structure. We return to the muscle data to illustrate both of these modes.
 
 \subsubsection{Reducing the dimensionality of the data}
 Once we have selected the genes we will use to order the cells, Monocle applies a \emph{dimensionality reduction} to the data, which will drastically improves the quality of the trajectory. Monocle reduces the dimensionality of the data with the \Rfunction{reduceDimension} function. This function has a number of options, so you should familiarize yourself with them by consulting the its manual page. You can choose from two algorithms in \Rfunction{reduceDimension}. The first, termed Independent Component Analysis, is a classic linear technique for decomposing data that powered the original version of Monocle. The second, called DDRTree, is a much more powerful nonlinear technique that is the default for Monocle 2.  For more on how these both work, see section \ref{theory}.   
 
 \subsubsection{Ordering the cells in pseudotime}
 With the expression data projected into a lower dimensional space, Monocle is ready to learn the trajectory that describes how cells transition from one state into another. Monocle assumes that the trajectory has a tree structure, with one end of it the ``root'', and the others the ``leaves''. A cell at the beginning of the biological process starts at the root and progresses along the trunk until it reaches the first branch, if there is one. It then chooses a branch, and moves further and further along the tree until it reaches a leaf. These mathematical assumptions translate into some important biological ones. First, that the data includes all the major stages of the biological process. If your experiment failed to capture any cells at a key developmental transition, Monocle won't know its there. Second, that gene expression changes are smooth as a cell moves from one stage to the next. This assumption is realistic: major discontinuities in the trajectory would amount to a cell almost instantaneously turning over its transcriptome, which probably doesn't happen in most biolgical processes.
 
 To order your cells, Monocle uses the \Rfunction{orderCells} function. This routine has one important argument, which allows you to set the root of the tree, and thus the beginning of the process. See the manual page for \Rfunction{orderCells} for more details.
 
 \subsection{Unsupervised ordering}
 
 In this section, we discuss ordering cells in a completely unsupervised fashion. First, we must decide which genes we will use to define a cell's progress through myogenesis. Monocle orders cells by examining the pattern of expression of these genes across the cell population. Monocle looks for genes that vary in ``interesting'' (i.e. not just noisy) ways, and uses these to structure the data. We ultimately want a set of genes that increase (or decrease) in expression as a function of progress through the process we're studying. 
 
 Ideally, we'd like to use as little prior knowledge of the biology of the system under study as possible. We'd like to discover the important ordering genes from the data, rather than relying on literature and textbooks, because that might introduce bias in the ordering. One effective way to isolate a set of ordering genes is to simply compare the cells collected at the beginning of the process to those at the end and find the differentially expressed genes, as described above. The command below will find all genes that are differentially expressed in response to the switch from growth medium to differentiation medium:
 
 <<ordering_not_run, eval=FALSE>>=
 diff_test_res <- differentialGeneTest(HSMM_myo[expressed_genes,], 
                                       fullModelFormulaStr="~Media")
 ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
 @ %def
 
 Choosing genes based on differential analysis of time points is often highly effective, but what if we don't have time series data? If the cells are asynchronously moving through our biological process (as is usually the case), Monocle can often reconstruct their trajectory from a single population captured all at the same time. Below are two methods to select genes that require no knowledge of the design of the experiment at all.
 
\subsubsection{Selecting genes with high dispersion across cells} 
Genes that vary a lot are often highly informative for identifying cell subpopulations or ordering cells along a trajectory. In RNA-Seq, a gene's variance typically depends on its mean, so we have to be a bit careful about how we select genes based on their variance.
 
 <<select_ordering_genes, eval=TRUE>>=
 disp_table <- dispersionTable(HSMM_myo)
 ordering_genes <- subset(disp_table, 
                          mean_expression >= 0.5 & 
                          dispersion_empirical >= 1 * dispersion_fit)$gene_id
 @ %def
 
Once we have a list of gene ids to be used for ordering, we need to set them in the \Robject{HSMM} object, because the next several functions will depend on them.
 


###################################################
### code chunk number 12: set_ordering_filter
###################################################
HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
plot_ordering_genes(HSMM_myo)


###################################################
### code chunk number 13: order_cells
###################################################
HSMM_myo <- orderCells(HSMM_myo)


###################################################
### code chunk number 14: plot_ordering_with_GM_state
###################################################
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$Hours)[,"0"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  }else { return (1) } 
}
HSMM_myo <- orderCells(HSMM_myo, root_state=GM_state(HSMM_myo))
plot_cell_trajectory(HSMM_myo, color_by="Pseudotime")


###################################################
### code chunk number 15: init_hsmm_facet_state
###################################################
plot_cell_trajectory(HSMM_myo, color_by="State") + facet_wrap(~State, nrow=1)


###################################################
### code chunk number 16: init_hsmm_jitter_state
###################################################
blast_genes <- row.names(subset(fData(HSMM_myo), gene_short_name %in% c("CCNB2", "MYOD1", "MYOG")))
plot_genes_jitter(HSMM_myo[blast_genes,], grouping="State", min_expr=0.1)


###################################################
### code chunk number 17: plot_markers_linear
###################################################
HSMM_expressed_genes <-  row.names(subset(fData(HSMM_myo), num_cells_expressed >= 10))
HSMM_filtered <- HSMM_myo[HSMM_expressed_genes,]
 
my_genes <- row.names(subset(fData(HSMM_filtered), 
                              gene_short_name %in% c("CDK1", "MEF2C", "MYH3"))) 
 
cds_subset <- HSMM_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by="Hours")


###################################################
### code chunk number 18: set_ordering_gene_dpfeature
###################################################
HSMM_myo <- detectGenes(HSMM_myo, min_expr=0.1)
fData(HSMM_myo)$use_for_ordering <- fData(HSMM_myo)$num_cells_expressed > 0.05 * ncol(HSMM_myo)


###################################################
### code chunk number 19: plot_pc_variance
###################################################
plot_pc_variance_explained(HSMM_myo, return_all = F) #look at the plot and decide how many dimensions you need. It is determined by a huge drop of variance at that dimension. pass that number to num_dim in the next function.


###################################################
### code chunk number 20: perform tSNE dimension reduction
###################################################
  HSMM_myo <- reduceDimension(HSMM_myo, max_components=2, norm_method = 'log', num_dim = 3, 
                              reduction_method = 'tSNE', verbose = T)


###################################################
### code chunk number 21: cluster cells
###################################################
HSMM_myo <- clusterCells(HSMM_myo, verbose = F)


###################################################
### code chunk number 22: check the clustering results
###################################################
plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Cluster)')
plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Hours)')


###################################################
### code chunk number 23: plot_rho_delta
###################################################
plot_rho_delta(HSMM_myo, rho_threshold = 2, delta_threshold = 4 )


###################################################
### code chunk number 24: recluster_cells_again
###################################################
HSMM_myo <- clusterCells(HSMM_myo,  
                         rho_threshold = 2, 
                         delta_threshold = 4, 
                         skip_rho_sigma = T, 
                         verbose = F)


###################################################
### code chunk number 25: check_clustering_again
###################################################
plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Cluster)')
plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Hours)')


###################################################
### code chunk number 26: perform_DEG_analysis
###################################################
clustering_DEG_genes <- differentialGeneTest(HSMM_myo[HSMM_expressed_genes,], 
                                             fullModelFormulaStr = '~Cluster', 
                                             cores = 1)


###################################################
### code chunk number 27: order_cells_dp_feature
###################################################
HSMM_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000] 

HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes = HSMM_ordering_genes)
HSMM_myo <- reduceDimension(HSMM_myo)
HSMM_myo <- orderCells(HSMM_myo)
HSMM_myo <- orderCells(HSMM_myo, root_state=GM_state(HSMM_myo))
plot_cell_trajectory(HSMM_myo, color_by="Hours")



###################################################
### code chunk number 28: semi_sup_ordering_myoblast_cth
###################################################
CCNB2_id <- row.names(subset(fData(HSMM_myo), gene_short_name == "CCNB2"))
MYH3_id <- row.names(subset(fData(HSMM_myo), gene_short_name == "MYH3"))

cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "Cycling myoblast", classify_func=function(x) {x[CCNB2_id,] >= 1})
cth <- addCellType(cth, "Myotube", classify_func=function(x) {x[MYH3_id,] >=1})
cth <- addCellType(cth, "Reserve cell", classify_func=function(x) {x[MYH3_id,] == 0 & x[CCNB2_id,] == 0})
HSMM_myo <- classifyCells(HSMM_myo, cth)


###################################################
### code chunk number 29: semi_sup_ordering_markers
###################################################
marker_diff <- markerDiffTable(HSMM_myo[HSMM_expressed_genes,], 
                               cth, 
                               cores=1)
#semisup_clustering_genes <- row.names(subset(marker_diff, qval < 0.05))
semisup_clustering_genes <- row.names(marker_diff)[order(marker_diff$qval)][1:500] 


###################################################
### code chunk number 30: semi_sup_ordering_trajectory
###################################################
HSMM_myo <- setOrderingFilter(HSMM_myo, semisup_clustering_genes)
#plot_ordering_genes(HSMM_myo)
HSMM_myo <- reduceDimension(HSMM_myo, max_components=2)
HSMM_myo <- orderCells(HSMM_myo)
HSMM_myo <- orderCells(HSMM_myo, root_state=GM_state(HSMM_myo))
plot_cell_trajectory(HSMM_myo, color_by="CellType") + theme(legend.position="right")


###################################################
### code chunk number 31: semi_sup_ordering_genes_in_pseudotime
###################################################
HSMM_filtered <- HSMM_myo[HSMM_expressed_genes,]

my_genes <- row.names(subset(fData(HSMM_filtered), 
                             gene_short_name %in% c("CDK1", "MEF2C", "MYH3"))) 

cds_subset <- HSMM_filtered[my_genes,]
plot_genes_branched_pseudotime(cds_subset, 
                               branch_point=1, 
                               color_by="Hours", 
                               ncol=1)


###################################################
### code chunk number 32: select_genes
###################################################
marker_genes <- row.names(subset(fData(HSMM_myo), 
                                 gene_short_name %in% c("MEF2C", "MEF2D", "MYF5", 
                                                        "ANPEP", "PDGFRA","MYOG", 
                                                        "TPM1",  "TPM2",  "MYH2", 
                                                        "MYH3",  "NCAM1", "TNNT1", 
                                                        "TNNT2", "TNNC1", "CDK1", 
                                                        "CDK2",  "CCNB1", "CCNB2", 
                                                        "CCND1", "CCNA1", "ID1")))


###################################################
### code chunk number 33: basic_diff
###################################################
diff_test_res <- differentialGeneTest(HSMM_myo[marker_genes,], 
                                      fullModelFormulaStr="~Media")

# Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res, qval < 0.1)

sig_genes[,c("gene_short_name", "pval", "qval")]


###################################################
### code chunk number 34: plot_myog_jitter
###################################################
MYOG_ID1 <- HSMM_myo[row.names(subset(fData(HSMM_myo), 
                                      gene_short_name %in% c("MYOG", "CCNB2"))),]
plot_genes_jitter(MYOG_ID1, grouping="Media", ncol=2)


###################################################
### code chunk number 35: setup_test_genes
###################################################
to_be_tested <- row.names(subset(fData(HSMM), 
                                 gene_short_name %in% c("UBC", "NCAM1", "ANPEP"))) 
cds_subset <- HSMM[to_be_tested,]


###################################################
### code chunk number 36: all_in_one_test
###################################################
diff_test_res <- differentialGeneTest(cds_subset, fullModelFormulaStr="~CellType")
diff_test_res[,c("gene_short_name", "pval", "qval")]


###################################################
### code chunk number 37: jitter_plot_diff_res
###################################################
plot_genes_jitter(cds_subset, grouping="CellType", color_by="CellType", 
                  nrow=1, ncol=NULL, plot_trend=TRUE)


###################################################
### code chunk number 38: piecewise_test (eval = FALSE)
###################################################
## full_model_fits <- fitModel(cds_subset,  modelFormulaStr="~CellType")
## reduced_model_fits <- fitModel(cds_subset, modelFormulaStr="~1")
## diff_test_res <- compareModels(full_model_fits, reduced_model_fits)
## diff_test_res


###################################################
### code chunk number 39: setup_test_genes_pt
###################################################
to_be_tested <- row.names(subset(fData(HSMM), 
                                 gene_short_name %in% c("MYH3", "MEF2C", "CCNB2", "TNNT1"))) 
cds_subset <- HSMM_myo[to_be_tested,]


###################################################
### code chunk number 40: piecewise_test_pt
###################################################
diff_test_res <- differentialGeneTest(cds_subset,  fullModelFormulaStr="~sm.ns(Pseudotime)")


###################################################
### code chunk number 41: all_in_one_test_pt
###################################################
diff_test_res[,c("gene_short_name", "pval", "qval")]


###################################################
### code chunk number 42: plot_diff_res_pt
###################################################
plot_genes_in_pseudotime(cds_subset, color_by="Hours")


###################################################
### code chunk number 43: plot_diff_res_pt_heatmap
###################################################
diff_test_res <- differentialGeneTest(HSMM_myo[marker_genes,], 
                                      fullModelFormulaStr="~sm.ns(Pseudotime)")

sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))


plot_pseudotime_heatmap(HSMM_myo[sig_gene_names,], 
                        num_clusters = 3, 
                        cores = 1, 
                        show_rownames = T)


###################################################
### code chunk number 44: plot_diff_res_multi
###################################################
to_be_tested <- row.names(subset(fData(HSMM), 
                                 gene_short_name %in% c("TPM1", "MYH3", "CCNB2", "GAPDH"))) 
cds_subset <- HSMM[to_be_tested,]

diff_test_res <- differentialGeneTest(cds_subset,  
                                      fullModelFormulaStr="~CellType + Hours", 
                                      reducedModelFormulaStr="~Hours")
diff_test_res[,c("gene_short_name", "pval", "qval")]
plot_genes_jitter(cds_subset, 
                  grouping="Hours", color_by="CellType", plot_trend=TRUE) + 
  facet_wrap( ~ feature_label, scales="free_y")


###################################################
### code chunk number 45: init_treutlein
###################################################
lung <- load_lung()

plot_cell_trajectory(lung, color_by="Time")



###################################################
### code chunk number 46: lung_beam_test
###################################################
BEAM_res <- BEAM(lung, branch_point=1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]


###################################################
### code chunk number 47: lung_beam_branched_heatmap
###################################################
plot_genes_branched_heatmap(lung[row.names(subset(BEAM_res, qval < 1e-4)),], 
                            branch_point = 1, 
                            num_clusters = 4, 
                            cores = 1, 
                            use_gene_short_name = T, 
                            show_rownames = T)


###################################################
### code chunk number 48: lung_beam_branched_pseudotime
###################################################
lung_genes <- row.names(subset(fData(lung), gene_short_name %in% c("Ccnd2", "Sftpb", "Pdpn")))
plot_genes_branched_pseudotime(lung[lung_genes,], 
                               branch_point=1, 
                               color_by="Time", 
                               ncol=1)


###################################################
### code chunk number 49: sessi
###################################################
sessionInfo()


