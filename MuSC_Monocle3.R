#' ---
#' title: "Monocle3 Pseudo trajectory analysis on MuSC Old mice C26 and Ctrl "
#' output: html_document
#' date: "2025-12-29"
#' ---
#' 
## ----setup, include=FALSE-------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
## ----Libraries------------------------------------------------------------------------------------------------------------
library(SingleCellExperiment);library(Seurat);library(tidyverse);library(writexl);library(paletteer);
library(Matrix);library(scales);library(cowplot);library(RCurl);library(openxlsx) ;library(knitr);library(monocle3);
library(SeuratWrappers);library(hdf5r);library(slingshot);library(BUSpaRse);library(tidyverse);library(tidymodels);
library(Seurat);library(scales);library(viridis);library(Matrix);library(grDevices); library(RColorBrewer);library(ggplot2);
library(patchwork);library(ClusterGVis);library(readxl);library(tidydr);library(dplyr);library(tibble);library(writexl)

#' 
## ----Loading and Pre-processing-------------------------------------------------------------------------------------------
load("MuSC/MuSC.RData")
DefaultAssay(Musc)<- "RNA"
Idents(Musc)<- "celltype1"

## Subset by Sample/Treatment
Idents(Musc)<-"Sample"
C26 <- subset(Musc, idents = "C26")
Ctrl<- subset(Musc, idents = c("Ctrl"))
DefaultAssay(C26)<- "RNA"
DefaultAssay(Ctrl)<- "RNA"

#' 
#' #--------------------------------------
#' CytoTRACE2 for developmental potential
#' #--------------------------------------
## -------------------------------------------------------------------------------------------------------------------------
library(CytoTRACE2) 
#--------------
# Cachexia sample
DefaultAssay(C26)<- "RNA"
C26 <- cytotrace2(C26, is_seurat = T, slot_type = "counts", species = "mouse")
annotation <- C26$celltype1
annotation <- data.frame(annotation)
rownames(annotation) <- colnames(C26)
colnames(annotation)[1] <- "celltype" 
plots <- plotData(cytotrace2_result = C26, annotation = annotation, expression_data = C26, is_seurat = T )
pdf("MuSC/CytoTRACE2/C26_cytotrace.pdf", height = 5, width = 7)
plots$CytoTRACE2_UMAP
plots$CytoTRACE2_Potency_UMAP
plots$CytoTRACE2_Relative_UMAP
plots$Phenotype_UMAP
plots$CytoTRACE2_Boxplot_byPheno
dev.off()
save(C26, file = "MuSC/CytoTRACE2/C26_cytotrace.RData")

#-------------
# Ctrl sample
DefaultAssay(Ctrl)<- "RNA"
Ctrl <- cytotrace2(Ctrl, is_seurat = T, slot_type = "counts", species = "mouse")
annotation <- Ctrl$celltype1
annotation <- data.frame(annotation)
rownames(annotation) <- colnames(Ctrl)
colnames(annotation)[1] <- "celltype" 
plots <- plotData(cytotrace2_result = Ctrl, annotation = annotation, expression_data = Ctrl, is_seurat = T )
pdf("MuSC/CytoTRACE2/Ctrl_cytotrace.pdf", height = 5, width = 7)
plots$CytoTRACE2_UMAP
plots$CytoTRACE2_Potency_UMAP
plots$CytoTRACE2_Relative_UMAP
plots$Phenotype_UMAP
plots$CytoTRACE2_Boxplot_byPheno
dev.off()
save(Ctrl, file = "MuSC/CytoTRACE2/Ctrl_cytotrace.RData")

#' #-----------------------------------------
#' # Pseudotrajectory analysis using Monocle3
#' #-------------------------------------------
library(monocle3)
library(SingleCellExperiment)
library(Seurat)

#---------
# Monocle3 in C26 
DefaultAssay(C26)<- "RNA"
Idents(C26)<- "celltype1"

counts_mat <- GetAssayData(C26, assay="RNA", layer="counts")
C26.cds <- new_cell_data_set(expression_data = counts_mat,cell_metadata = data.frame( my_cluster = C26$celltype1, row.names = colnames(C26)),gene_metadata = data.frame( gene_short_name = rownames(counts_mat), row.names = rownames(counts_mat)))
# preprocess + reduce 
C26.cds <- preprocess_cds(C26.cds, num_dim = 50, scaling = TRUE, method = "PCA")
C26.cds <- reduce_dimension(C26.cds, preprocess_method = "PCA", reduction_method = "UMAP")
# overwrite monocle UMAP with Seurat UMAP 
reducedDims(C26.cds)$UMAP <- C26@reductions$umap@cell.embeddings
# monocle recluster 
C26.cds <- cluster_cells(C26.cds)  
# learn graph 
C26.cds <- learn_graph(C26.cds, use_partition = FALSE, learn_graph_control = list(ncenter = 100))
# Ordering cells
C26.cds <- order_cells(C26.cds)

#-------
# Viz
p <- plot_cells(cds = C26.cds,color_cells_by = "pseudotime",show_trajectory_graph = TRUE,
  label_branch_points = TRUE,label_roots = T,label_leaves = TRUE, cell_size = 0.9,                 
  graph_label_size = 3.5,           trajectory_graph_color = "grey20", trajectory_graph_segment_size = 0.6) +
  guides(color = guide_colorbar( title = "Pseudotime", barheight = unit(55, "mm"), barwidth  = unit(5, "mm"), ticks = TRUE)) +
  coord_equal() + theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 14, color = "black", face = "bold"), axis.text  = element_text(size = 12, color = "black"),
    axis.line  = element_line(linewidth = 0.8, color = "black"), legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 10), legend.key.height = unit(6, "mm"), plot.title = element_text(face = "bold", hjust = 0.5))
# save cds object
save(C26.cds, file = "MuSC/Monocle3/C26.cds.RData")

# Plot genes 
plot_cells(C26.cds, genes= c("Cdkn1a", "Myod1", "Mt2", "Tnfrsf12a"),label_cell_groups=FALSE, show_trajectory_graph=TRUE, min_expr=1)

#---
# Gene Correlation with pseudotime
C26_pt_res <- graph_test(C26.cds, neighbor_graph = "principal_graph", cores = 8)
C26_pt_res<- na.omit(C26_pt_res)

write_xlsx(C26_pt_res, "MuSC/Monocle3/C26_gene_predicted_pseudo.xlsx")
C26_pt_res<- C26_pt_res[C26_pt_res$p_value < 0.05 & C26_pt_res$status == "OK",]

#' 
## -------------------------------------------------------------------------------------------------------------------------
#---------
# Monocle3 in Control 
Control<- Ctrl # tmp variable
DefaultAssay(Control) <- "RNA"

counts_mat <- GetAssayData(Control, assay="RNA", layer="counts")
Control.cds <- new_cell_data_set(
  expression_data = counts_mat,
  cell_metadata = data.frame(my_cluster = Control$celltype1,row.names = colnames(Control)), gene_metadata = data.frame(gene_short_name = rownames(counts_mat),row.names = rownames(counts_mat)))

# preprocess
Control.cds <- preprocess_cds(Control.cds, num_dim = 50, scaling = TRUE, method = "PCA")
Control.cds <- reduce_dimension(Control.cds, preprocess_method = "PCA", reduction_method = "UMAP")
# overwrite monocle UMAP with Seurat UMAP
reducedDims(Control.cds)$UMAP <- Control@reductions$umap@cell.embeddings
# monocle recluster 
Control.cds <- cluster_cells(Control.cds)  
# learn graph 
Control.cds <- learn_graph(Control.cds, use_partition = FALSE, learn_graph_control = list(ncenter = 80)) 

# Order cells
Control.cds <- order_cells(Control.cds)
save(Control.cds, file = "MuSC/Monocle3/Control.cds.RData")
#----
# Viz
q <- plot_cells(cds = Control.cds,color_cells_by = "pseudotime",show_trajectory_graph = TRUE,
  label_branch_points = TRUE,label_roots = TRUE,label_leaves = TRUE,
  cell_size = 0.9,                 graph_label_size = 3.5,         trajectory_graph_color = "grey20",trajectory_graph_segment_size = 0.6) +guides(color = guide_colorbar( title = "Pseudotime", barheight = unit(55, "mm"), barwidth  = unit(5, "mm"), ticks = TRUE)) + coord_equal() +theme_classic(base_size = 12) +
  theme( axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text  = element_text(size = 12, color = "black"),axis.line  = element_line(linewidth = 0.8, color = "black"),
    legend.title = element_text(size = 12, face = "bold"),legend.text  = element_text(size = 10),
    legend.key.height = unit(6, "mm"), plot.title = element_text(face = "bold", hjust = 0.5))

#---
# Gene Correlation with pseudotime
Control_pt_res <- graph_test(Control.cds, neighbor_graph = "principal_graph", cores = 8)
Control_pt_res<- na.omit(Control_pt_res)
write_xlsx(Control_pt_res, "MuSC/Monocle3/Control_gene_predicted_pseudo.xlsx")
Control_pt_res<- Control_pt_res[Control_pt_res$p_value < 0.05 & Control_pt_res$status == "OK",]

#' 
#' # Final Vizualization
#' #-------
## -------------------------------------------------------------------------------------------------------------------------
# cols
ct_cols <- c(
  "Activated MuSC"   = as.character(paletteer_d("rcartocolor::Teal")[5]),
  "Quiescent MuSC"   = as.character(paletteer_d("rcartocolor::Vivid")[5]),
  "Cachexia MuSC"   = as.character(paletteer_d("rcartocolor::Teal")[7]))

# Dimplot/UMAP of each condition
a <- DimPlot(C26, group.by = "celltype1", label = F, cols = ct_cols) + ggtitle("C26") +
  theme_classic(base_size = 12) + coord_equal() +
  theme(axis.title = element_text(size = 14, color = "black", face = "bold"),axis.text  = element_text(size = 12, color = "black"), axis.line  = element_line(linewidth = 0.8, color = "black"),legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 10),legend.key.height = unit(6, "mm"), plot.title = element_text(face = "bold", hjust = 0.5, size = 15))

c <- DimPlot(Control, group.by = "celltype1", label = F, cols = ct_cols) +ggtitle("Ctrl") +coord_equal() +
  theme_classic(base_size = 12) +theme(axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text  = element_text(size = 12, color = "black"),axis.line  = element_line(linewidth = 0.8, color = "black"),
    legend.title = element_text(size = 12, face = "bold"),legend.text  = element_text(size = 10), legend.key.height = unit(6, "mm"), plot.title = element_text(face = "bold", hjust = 0.5, size = 15) )

pdf("MuSC/Monocle3/C26_Pseudotime.pdf", width = 10, height = 5)
a | p
dev.off()

pdf("MuSC/Monocle3/Ctrl_Pseudotime.pdf", width = 10, height = 5)
c | q
dev.off()

## ----Heatmaps of pseudotime trajectory analysis---------------------------------------------------------------------------
#-----
# C26 output 
load("MuSC/Monocle3/C26.cds.RData")
C26_pt_res<- read_xlsx("MuSC/Monocle3/C26_gene_predicted_pseudo.xlsx")
C26_pt_res<- column_to_rownames(C26_pt_res, var = "gene_short_name")

# processing
genes <- row.names(subset(C26_pt_res, q_value < 0.01 & morans_I > 0.2))
pre_pseudotime_matrix <- getFromNamespace("pre_pseudotime_matrix","ClusterGVis")
mat <- pre_pseudotime_matrix(cds_obj = C26.cds, gene_list = genes)
ck <- clusterData(obj = mat,clusterMethod = "kmeans", clusterNum = 4)
pdf('MuSC/Monocle3/C26_heatmap_monocle3.pdf',height = 10,width = 8,onefile = F)
visCluster(object = ck, plotType = "both",addSampleAnno = F, markGenes = sample(rownames(mat),50,replace = F))
dev.off()

#-----
# Control output 
load("MuSC/Monocle3/Control.cds.RData")
Control_pt_res<- read_xlsx("MuSC/Monocle3/Control_gene_predicted_pseudo.xlsx")
Control_pt_res<- column_to_rownames(Control_pt_res, var = "gene_short_name")

# processing
genes <- row.names(subset(Control_pt_res, q_value < 0.01 & morans_I > 0.2))
pre_pseudotime_matrix <- getFromNamespace("pre_pseudotime_matrix","ClusterGVis")
mat <- pre_pseudotime_matrix(cds_obj = Control.cds,gene_list = genes)
ck <- clusterData(obj = mat, clusterMethod = "kmeans", clusterNum = 4)
pdf('MuSC/Monocle3/Control_heatmap_monocle3.pdf',height = 10,width = 8,onefile = F)
visCluster(object = ck, plotType = "both", addSampleAnno = F, markGenes = sample(rownames(mat),50,replace = F))
dev.off()

## ----session-info---------------------------------------------------------------------------------------------------------
sessioninfo::session_info()

# Session Information

# R version 4.5.1 (2025-06-13 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 26200)
# 
# Matrix products: default
#   LAPACK version 3.12.1
# 
# locale:
# [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
# [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] CytoTRACE2_1.1.0            RSpectra_0.16-2             Rfast_2.1.5.2               RcppParallel_5.1.11-2      
#  [5] zigg_0.0.2                  Rcpp_1.1.1-1.1              plyr_1.8.9                  magrittr_2.0.5             
#  [9] HiClimR_2.2.1               doParallel_1.0.17           iterators_1.0.14            foreach_1.5.2              
# [13] data.table_1.18.4           tidydr_0.0.6                readxl_1.5.0                ClusterGVis_0.99.9         
# [17] patchwork_1.3.2             RColorBrewer_1.1-3          viridis_0.6.5               viridisLite_0.4.3          
# [21] yardstick_1.4.0             workflowsets_1.1.1          workflows_1.3.0             tune_2.1.0                 
# [25] tailor_0.1.0                rsample_1.3.2               recipes_1.3.3               parsnip_1.6.0              
# [29] modeldata_1.5.1             infer_1.1.0                 dials_1.4.3                 broom_1.0.13               
# [33] tidymodels_1.5.0            BUSpaRse_1.22.1             slingshot_2.16.0            TrajectoryUtils_1.16.1     
# [37] princurve_2.1.6             hdf5r_1.3.12                SeuratWrappers_0.4.0        monocle3_1.4.26            
# [41] knitr_1.51                  openxlsx_4.2.8.1            RCurl_1.98-1.19             cowplot_1.2.0              
# [45] scales_1.4.0                Matrix_1.7-3                paletteer_1.7.0             writexl_1.5.4              
# [49] lubridate_1.9.5             forcats_1.0.1               stringr_1.6.0               dplyr_1.2.1                
# [53] purrr_1.2.2                 readr_2.2.0                 tidyr_1.3.2                 tibble_3.3.1               
# [57] ggplot2_4.0.3               tidyverse_2.0.0             Seurat_5.5.0                SeuratObject_5.4.0         
# [61] sp_2.2-1                    SingleCellExperiment_1.30.1 SummarizedExperiment_1.38.1 Biobase_2.68.0             
# [65] GenomicRanges_1.60.0        GenomeInfoDb_1.44.3         IRanges_2.44.0              S4Vectors_0.48.0           
# [69] BiocGenerics_0.54.1         generics_0.1.4              MatrixGenerics_1.20.0       matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] R.methodsS3_1.8.2        dichromat_2.0-0.1        progress_1.2.3           nnet_7.3-20             
#   [5] goftest_1.2-3            Biostrings_2.78.0        vctrs_0.7.3              spatstat.random_3.5-0   
#   [9] digest_0.6.39            png_0.1-9                plyranges_1.28.0         ggrepel_0.9.8           
#  [13] deldir_2.0-4             parallelly_1.47.0        MASS_7.3-65              reshape2_1.4.5          
#  [17] httpuv_1.6.17            withr_3.0.2              ggfun_0.2.0              xfun_0.55               
#  [21] survival_3.8-3           memoise_2.0.1            Seqinfo_1.0.0            zoo_1.8-15              
#  [25] pbapply_1.7-4            R.oo_1.27.1              prettyunits_1.2.0        rematch2_2.1.2          
#  [29] KEGGREST_1.50.0          promises_1.5.0           otel_0.2.0               httr_1.4.8              
#  [33] restfulr_0.0.16          globals_0.19.1           fitdistrplus_1.2-6       rstudioapi_0.18.0       
#  [37] UCSC.utils_1.4.0         miniUI_0.1.2             ncdf4_1.24               curl_7.1.0              
#  [41] polyclip_1.10-7          GenomeInfoDbData_1.2.14  SparseArray_1.10.8       xtable_1.8-8            
#  [45] evaluate_1.0.5           S4Arrays_1.10.1          BiocFileCache_3.0.0      hms_1.1.4               
#  [49] irlba_2.3.7              filelock_1.0.3           ROCR_1.0-12              reticulate_1.46.0       
#  [53] spatstat.data_3.1-9      lmtest_0.9-40            later_1.4.8              lattice_0.22-7          
#  [57] spatstat.geom_3.8-1      future.apply_1.20.2      scattermore_1.2          XML_3.99-0.23           
#  [61] scuttle_1.18.0           RcppAnnoy_0.0.23         class_7.3-23             pillar_1.11.1           
#  [65] nlme_3.1-168             compiler_4.5.1           beachmat_2.26.0          stringi_1.8.7           
#  [69] gower_1.0.2              tensor_1.5.1             minqa_1.2.8              GenomicAlignments_1.44.0
#  [73] crayon_1.5.3             abind_1.4-8              BiocIO_1.18.0            bit_4.6.0               
#  [77] codetools_0.2-20         plotly_4.12.0            mime_0.13                splines_4.5.1           
#  [81] fastDummies_1.7.6        dbplyr_2.5.2             DiceDesign_1.10          cellranger_1.1.0        
#  [85] blob_1.3.0               AnnotationFilter_1.32.0  lme4_2.0-1               fs_2.1.0                
#  [89] listenv_0.10.1           Rdpack_2.6.6             tzdb_0.5.0               pkgconfig_2.0.3         
#  [93] tools_4.5.1              cachem_1.1.0             rbibutils_2.4.1          RSQLite_3.53.1          
#  [97] DBI_1.3.0                fastmap_1.2.0            grid_4.5.1               ica_1.0-3               
# [101] Rsamtools_2.24.1         BiocManager_1.30.27      dotCall64_1.2            RANN_2.6.2              
# [105] rpart_4.1.24             farver_2.1.2             reformulas_0.4.4         yaml_2.3.12             
# [109] VGAM_1.1-14              rtracklayer_1.68.0       cli_3.6.6                lifecycle_1.0.5         
# [113] uwot_0.2.4               sessioninfo_1.2.4        lava_1.9.1               backports_1.5.1         
# [117] BiocParallel_1.44.0      timechange_0.4.0         gtable_0.3.6             rjson_0.2.23            
# [121] ggridges_0.5.7           progressr_0.19.0         jsonlite_2.0.0           RcppHNSW_0.7.0          
# [125] bitops_1.0-9             bit64_4.8.2              Rtsne_0.17               yulab.utils_0.2.4       
# [129] spatstat.utils_3.2-3     zip_2.3.3                zeallot_0.2.0            spatstat.univar_3.2-0   
# [133] R.utils_2.13.0           timeDate_4052.112        lazyeval_0.2.3           shiny_1.13.0            
# [137] htmltools_0.5.9          sctransform_0.4.3        rappdirs_0.3.4           ensembldb_2.32.0        
# [141] glue_1.8.1               spam_2.11-4              httr2_1.2.2              XVector_0.50.0          
# [145] BSgenome_1.76.0          gridExtra_2.3            boot_1.3-31              igraph_2.3.2            
# [149] R6_2.6.1                 GenomicFeatures_1.60.0   cluster_2.1.8.2          ipred_0.9-15            
# [153] nloptr_2.2.1             DelayedArray_0.36.0      tidyselect_1.2.1         ProtGenerics_1.40.0     
# [157] xml2_1.5.2               AnnotationDbi_1.72.0     future_1.70.0            rsvd_1.0.5              
# [161] KernSmooth_2.23-26       S7_0.2.2                 furrr_0.4.0              htmlwidgets_1.6.4       
# [165] biomaRt_2.64.0           rlang_1.2.0              spatstat.sparse_3.2-0    spatstat.explore_3.8-1  
# [169] remotes_2.5.0            hardhat_1.4.3            prodlim_2026.03.11      
