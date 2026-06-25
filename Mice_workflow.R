#' ---
#' title: "Mice dataset integrated workflow for C26 vs Control mice and Muscle Stem cells characterization"
#' author: "Carlos Alfaro"
#' date: "2025-05-01"
#' output: html_document
#' ---
#' 
## ----setup, include=FALSE-------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
## ----libraries------------------------------------------------------------------------------------------------------------
library(SingleCellExperiment);library(Seurat);library(tidyverse);library(Matrix);library(scales);
library(cowplot);library(RCurl);library(openxlsx);library(knitr);library(SeuratWrappers);library(hdf5r);library(paletteer);
library(ggrepel);library(ggrastr);library(ggprism);library(dplyr);library(ggplot2);library(viridis);
library(clustree);library(presto);library(pheatmap);library(RColorBrewer)

#' 
#' # Loading Cellbender files from Cachexia and Normal mice 
## ----Pre processing-------------------------------------------------------------------------------------------------------
file_paths <- list(
  Ctrl = "data/Partek_Miller_scRNAseq_tdtomato_8457-MR-0106_S1_filtered_feature_bc_matrix.h5",
  C26 = "data/Partek_Miller_scRNAseq_tdtomato_8457-MR-0107_S1_filtered_feature_bc_matrix.h5")

# Cutoffs are relaxed for initial exploration
min.cells <- 0  
min.features <- 0
seurat_list<-list()
cat("Creating Mice Seurat Objects ..")
for (name in names(file_paths)) {
  df <- Read10X_h5(file_paths[[name]])
  seurat_list[[name]] <- CreateSeuratObject(counts = df,project = name,min.cells = min.cells,min.features = min.features)}
cat("done...!!!")

# Merging the files into a single object.
cat("Merging all samples into a single one Seurat...")
GT<-merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list),project = "Muscle_SC")
cat("done!...")
View(GT@meta.data)

#' 
## ----Addidng metadata-----------------------------------------------------------------------------------------------------
metadata_mapping <- data.frame(SampleID = c("Ctrl", "C26"),Treatment = c("Control","Cachexia"),Age = c("Unknown", "Unknown"))
# Sample name modification 
cell_sample_ids <- sapply(strsplit(colnames(GT), "_"), function(x) paste(x[-length(x)], collapse = "_"))
cell_metadata <- data.frame(Treatment = metadata_mapping$Treatment[match(cell_sample_ids, metadata_mapping$SampleID)],Age = metadata_mapping$Age[match(cell_sample_ids, metadata_mapping$SampleID)])
rownames(cell_metadata) <- colnames(GT)
# Merge metadata to original one 
GT <- AddMetaData(GT, metadata = cell_metadata)

#' 
## ----pMito, pRibo, Complexity scores--------------------------------------------------------------------------------------
# Calculating the pMito per cell.
GT[['pMito']] <- PercentageFeatureSet(GT, pattern = '^mt-')
if(sum(GT[['pMito']], na.rm = TRUE) == 0) {
  GT[['pMito']] <- PercentageFeatureSet(GT, pattern = '^MT-')}
if(sum(GT[['pMito']], na.rm = TRUE) == 0) {
  GT[['pMito']] <- PercentageFeatureSet(GT, pattern = '^Mt-')}

# Calculating the pRibo per cell.
GT[["pRibo"]] <- PercentageFeatureSet(GT, pattern = "^Rp[sl]")
# Calculating cell complexity per cell per cell.
GT@meta.data<-GT@meta.data |>dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA) |>dplyr::mutate(log10GenesPerUMI = log10(nGene) / log10(nUMI))
# Adding one column for sample name
GT@meta.data$Sample<- sapply(strsplit(colnames(GT), "_"), function(x) paste(x[-length(x)], collapse = "_"))

# Save unfiltered Seurat object for later purposes.
cat("Saving unfiltered object...")
saveRDS(GT, file = "Mouse_combined_unfiltered_GT.rds")
cat('done!\n')

## ----QC viz---------------------------------------------------------------------------------------------------------------
pdf("QC/QC_control_features_before.pdf", width = 15, height = 8)
VlnPlot(GT, features = c("nUMI", "nGene", "pMito", "pRibo", "log10GenesPerUMI"), pt.size = 0,group.by = 'Sample', ncol = 4) + theme(legend.position = "none")
invisible(dev.off())

#' 
## ----Subset quality cells based on QC stats-------------------------------------------------------------------------------
# STATS
quantile(GT$nUMI,  probs = c(0.001, 0.01, 0.05, 0.5, 0.95, 0.99, 0.999), na.rm = TRUE)
quantile(GT$nGene, probs = c(0.001, 0.01, 0.05, 0.5, 0.95, 0.99, 0.999), na.rm = TRUE)
quantile(GT$pMito, probs = c(0.50, 0.80, 0.90, 0.95, 0.99), na.rm = TRUE)
quantile(GT$log10GenesPerUMI, probs = c(0.01, 0.05, 0.1, 0.5), na.rm = TRUE)

# Cutoffs
GT_filtered<-subset(GT, nUMI > 500 & nUMI < 50000 &  nGene > 200 & pMito <= 20 &log10GenesPerUMI > 0.8) 

## Visualization QC after filtering 
pdf("QC/QC_after_control_features.pdf", width = 15, height = 8)
VlnPlot(GT_filtered, features = c("nUMI", "nGene", "pMito", "pRibo", "log10GenesPerUMI"), pt.size = 0, group.by = 'Sample', ncol = 4) + theme(legend.position = "none")
invisible(dev.off())
nrow(GT_filtered) # 22544 cells

#' 
## ----Integration workflow-------------------------------------------------------------------------------------------------
# Joining Layers
cat("Joining the layers...")
GT_filtered<-JoinLayers(GT_filtered)
cat("done...")

# SCTransfrom by Sample
cat("Splitting the object...")
seuObject_split <- SplitObject(GT_filtered, split.by = "Sample")
cat("done happily...")

# Running SC Transform and regressing out variables
options(future.globals.maxSize = 8000 * 1024^2) 
cat("SC_transfomr is processing now...")
for (i in 1:length(seuObject_split)) {
  message("Running SC Transform on : ", names(seuObject_split)[i])
  seuObject_split[[i]] <- SCTransform(seuObject_split[[i]],vars.to.regress = c("nUMI", "pMito", "pRibo"),verbose = FALSE) 
  gc()}
cat("done SC transform! ready for integration now!...")

# Calculating CellCycle Scores
load("MouseCellCycleGenes.rda")
seuObject_split <- lapply(seuObject_split, function(x) {
  x<-CellCycleScoring(x,s.features = s_genes,g2m.features = g2m_genes,set.ident = TRUE,nbin =12)
  return(x)})
# Create a Difference score bettwen S - G2M per cell
seuObject_split <- lapply(seuObject_split, function(x) {x$CC.Difference <- x$S.Score - x$G2M.Score
return(x)})

# Select integration features = 3000, across the groups
integ_features <- SelectIntegrationFeatures(object.list = seuObject_split, nfeatures = 3000)
head(integ_features,10)

# Preparing the SCT object for integration
cat("Preparing the SCT object...")
seuObject_split <- PrepSCTIntegration(object.list = seuObject_split,anchor.features = integ_features)
cat("Done")

# Finding integration anchors
cat("Finding integration anchors....")
start_time <- Sys.time()
integ_anchors <- FindIntegrationAnchors(object.list = seuObject_split,normalization.method = "SCT",anchor.features = integ_features)
end_time <- Sys.time()
cat("Finished...!")

# Checking integration features 
length(integ_features)
table(integ_anchors@anchors[, "dataset1"])
table(integ_anchors@anchors[, "dataset2"])

# Integrate the data set into a single Seurat object
library(future)
cat("Integrating dataset...")
seuObject_integrated <- IntegrateData(anchorset = integ_anchors,new.assay.name = "integrated",
  normalization.method = "SCT",dims = 1:50, k.weight = 100,sd.weight = 1, eps = 0.5,verbose = TRUE)
cat('done...')

DefaultAssay(seuObject_integrated) <- "integrated"
# Dimensional reduction (PCA)
cat("Starting dimensionality reduction...")
seuObject_integrated <- RunPCA(seuObject_integrated,features = NULL, weight.by.var = TRUE, ndims.print = 1:5,
                               nfeatures.print = 30,npcs = 50,reduction.name = "pca")
ElbowPlot(seuObject_integrated, ndims = 50)

# Using the first 30 PCs to find neighbors and clusters
seuObject_integrated <- FindNeighbors(object=seuObject_integrated,reduction = "pca",dims = 1:30,nn.eps = 0.5)
seuObject_integrated <- FindClusters(seuObject_integrated,resolution = seq(0.1, 1.2, by = 0.1),algorithm = 1, n.iter = 1000) 

## ----UMAP viz-------------------------------------------------------------------------------------------------------------
pdf("Cluster_tree_reductions.pdf", width = 12, height = 10)
clustree(seuObject_integrated@meta.data, prefix = "integrated_snn_res.", node_colour = "sc3_stability")
invisible(dev.off())

set.seed(7081998)
seuObject_integrated <- RunUMAP(seuObject_integrated, dims = 1:30,reduction = "pca")
# Show resolutions created with UMAP
resolutions <- seq(0.1, 1.2, by = 0.1)
for (res in resolutions) {
  res_col <- paste0("integrated_snn_res.", res)
  Idents(seuObject_integrated) <- seuObject_integrated[[res_col]][,1]
  p <- DimPlot(seuObject_integrated, reduction = "umap", label = TRUE) +ggtitle(paste("Resolution", res)) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(filename = paste0("UMAP_by_resolution/UMAP_res_", res, ".png"),plot = p,width = 8, height = 6, dpi = 300)
  print(p)}

# Select the RNA counts slot to be the default assay and normalize
DefaultAssay(seuObject_integrated) <- "RNA"
seuObject_integrated <- NormalizeData(object = seuObject_integrated,normalization.method = "LogNormalize",scale.factor = 10000)

#' 
#' # scDblFinder for doublets identification
## -------------------------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
library(tidyverse);library(ggrepel);library(emmeans);library(SingleCellExperiment)
library(scater);library(BiocParallel);library(ggpubr);library(speckle);library(magrittr)
;library(broom);library(muscat);library(Seurat);library(clustree);library(leiden);library(data.table);library(cowplot)
;library(scDblFinder);library(BiocSingular);library(scds)})

#' 
## -------------------------------------------------------------------------------------------------------------------------
# Join layers 
seuObject_integrated_slim_Joined<-JoinLayers(seuObject_integrated)
# Transforming into a SingleCellExperiment 
sce <- as.SingleCellExperiment(seuObject_integrated_slim_Joined)
# Finding Doublets per sample
cat("Start running the doubleting detection.....")
set.seed(7081998)
sce <- scDblFinder(sce,samples="Sample", BPPARAM=BiocParallel::SnowParam(workers = 10),nfeatures = 3000,
                   dims = 30,dbr.sd = 1)
cat("Finished..!")

# Adding doublet identification on the metadata
seuObject_integrated_slim_Joined@meta.data$Doublets <- sce$scDblFinder.class
# Subset doublets from singlets
seuObject_slim_nodoub <- subset(seuObject_integrated_slim_Joined, subset = Doublets == "singlet")
# Determing percentage of singlets and doublets
table(sce$scDblFinder.class)  
table(sce$scDblFinder.class, sce$Sample) 
round(100 * prop.table(table(sce$scDblFinder.class, sce$Sample), margin = 2), 2)

#' 
## ----Saving Objects-------------------------------------------------------------------------------------------------------
## Slimming Seurat
seuObject_integrated_slim <- DietSeurat(seuObject_integrated,counts = TRUE,data = TRUE,scale.data = FALSE,assays="RNA",dimreducs = c("pca","umap"))

save(seuObject_integrated_slim, file = "output/Slim_Patient_SeuratObj_Final.RData")
save(seuObject_integrated, file = "Seurat_object.RData")


#' 
## ----Presto Markers-------------------------------------------------------------------------------------------------------
# Normalizing the singlets Seurat Object 
DefaultAssay(seuObject_slim_nodoub) <- "RNA"
seuObject_slim_nodoub <- NormalizeData(object = seuObject_slim_nodoub,normalization.method = "LogNormalize",scale.factor = 10000)

# Resolution
resolution <- "integrated_snn_res.0.3"
Idents(seuObject_slim_nodoub) <- resolution
DefaultAssay(seuObject_slim_nodoub)<- "RNA"

# Use function below in case presto is updated
wilcoxauc.Seurat <- function(X,group_by = NULL,assay = "data",groups_use = NULL, seurat_assay = "RNA",
    ...
) {
    requireNamespace("Seurat")
    X_matrix <- Seurat::GetAssayData(X, assay = seurat_assay, layer = assay)
    if (is.null(group_by)) {
        y <- Seurat::Idents(X)
    } else {
        y <- Seurat::FetchData(X, group_by) %>% unlist %>% as.character()
    }
    wilcoxauc(X_matrix, y, groups_use)
}
all_markers_clustID <- wilcoxauc.Seurat(seuObject_slim_nodoub, group_by ='integrated_snn_res.0.3')
all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")
all_markers.Sign <- all_markers_clustID %>%dplyr::filter(padj < 0.05, logFC > 0)
top20 <- presto::top_markers(all_markers.Sign,n = 20,auc_min = 0.5, pval_max = 0.05)

# Saving 
openxlsx::write.xlsx(all_markers.Sign,
                     file = "0.3_PrestoByCluster_Filteredmarkers_padjLT05_logfcGT0.xlsx",
                     colNames = TRUE,rowNames = FALSE,borders = "columns",sheetName="Markers")
openxlsx::write.xlsx(top20,
                     file = "0.3_PrestoByCluster_Top20.xlsx",colNames = TRUE,
                     rowNames = FALSE,borders = "columns",sheetName="Markers")

#' 
## ----Testing--------------------------------------------------------------------------------------------------------------
FeaturePlot(seuObject_slim_nodoub, features = "Pdgfra", min.cutoff = 'q10', max.cutoff = "q95", 
            cols = c("gray90", "purple4"), label = T)
DotPlot(seuObject_slim_nodoub, features = c("Cd19", "Cd22",
                                            "Ctsg", "Mcpt4",
                                            "Cdh5", "Aqp1",
                                            "Pdgfra", "Fbn1", "Has1",
                                            "Adgre1", "Mrc1",
                                            "Pax7", "Myf5", "Myod1",
                                            "Plp1", "Kcna1", "Prx",
                                            "Cxcr2", "Mmp9",
                                            "Acta2", "Myl9",
                                            "Eno3", "Ampd1",
                                            "Skap1", "Nkg7", "Itk",
                                            "Scx", "Cilp2"
                                            ), cols = "RdYlBu")+ RotatedAxis() + CoordFlip

FeaturePlot(seuObject_slim_nodoub, features = "Tdtom", min.cutoff = 'q10', max.cutoff = "q95", 
            cols = c("gray90", "purple4"), label = T)
DimPlot(seuObject_slim_nodoub, label = T)

## ----Cluster annotation---------------------------------------------------------------------------------------------------
cluster_annotations <- c(
  "0" = "FAPs",
  "1" = "FAPs", 
  "2" = "FAPs", 
  "3" = "SKeletal Muscle",
  "4" = "Endothelial cells",
  "5" = "SKeletal Muscle", 
  "6" = "Tenocytes",
  "7" = "Macrophages",
  "8" = "Neutrophils",
  "9" = "T and NK cells",
  "10" = "Schwann cells",
  "11" = "FAPs",
  "12" = "MuSC Progenitors",
  "13" = "Peryctes",
  "14" = "Basophils",
  "15" = "B-cells",
  "16" = "Lymphatic endothelial cells",
  "17" = "Adipocytes")
Idents(seuObject_slim_nodoub)<-seuObject_slim_nodoub$integrated_snn_res.0.3
seuObject_slim_nodoub@meta.data$celltype<- cluster_annotations[as.character(Idents(seuObject_slim_nodoub))]

DimPlot(seuObject_slim_nodoub, group.by = "celltype", label = T)
View(seuObject_slim_nodoub@meta.data)

## ----Final UMAP viz-------------------------------------------------------------------------------------------------------
Idents(seuObject_slim_nodoub)<-"celltype"
#PLottheme to use
plotTheme <- theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 20, face = "plain"),   
    axis.text  = element_text(size = 15, face = "plain"),   
    axis.line  = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks = element_line(linewidth = 0.6, colour = "black"),
    legend.position = "right",
    plot.title = element_text(size = 16, face = "plain", hjust = 0.5),
    plot.margin = margin(8, 8, 8, 8, "pt"))

# Define colors (one per celltype)
ct_cols <- c(
  "FAPs"                          = as.character(paletteer_d("rcartocolor::Teal")[5]),
  "SKeletal Muscle"               = as.character(paletteer_d("rcartocolor::Vivid")[5]),
  "Endothelial cells"             = as.character(paletteer_d("rcartocolor::Teal")[2]),
  "Tenocytes"                     = as.character(paletteer_d("rcartocolor::Safe")[1]),
  "Macrophages"                   = as.character(paletteer_d("rcartocolor::Antique")[1]),
  "Neutrophils"                   = as.character(paletteer_d("rcartocolor::Antique")[3]),
  "T and NK cells"                = as.character(paletteer_d("rcartocolor::Vivid")[1]),
  "Schwann cells"                 = as.character(paletteer_d("rcartocolor::Safe")[2]),
  "MuSC Progenitors"              = as.character(paletteer_d("rcartocolor::Vivid")[2]),
  "Peryctes"                      = as.character(paletteer_d("rcartocolor::Safe")[3]),
  "Basophils"                     = as.character(paletteer_d("rcartocolor::Antique")[6]),
  "B-cells"                       = as.character(paletteer_d("rcartocolor::Pastel")[2]),
  "Lymphatic endothelial cells"   = as.character(paletteer_d("rcartocolor::Teal")[3]),
  "Adipocytes"                    = as.character(paletteer_d("rcartocolor::Pastel")[5])
)
pdf("New_cellytype_General_umap.pdf", width = 8, height = 6)
DimPlot(seuObject_slim_nodoub, cols = ct_cols) + plotTheme
dev.off()

#' 
## ----Dotplot Viz----------------------------------------------------------------------------------------------------------
features = c("Cd19", "Cd22", "Ctsg", "Mcpt4", "Cdh5", "Aqp1", "Pdgfra", "Fbn1", "Has1",
             "Adgre1", "Mrc1","Pax7", "Myf5", 
             "Myod1", "Plp1", "Kcna1", "Prx",
             "Cxcr2", "Mmp9","Acta2", "Myl9",
             "Eno3", "Ampd1", "Skap1", "Nkg7", "Itk",
             "Scx", "Cilp2","Adipoq", "Plin1" )

celltype_levels <- c(
  "FAPs",
  "SKeletal Muscle",
  "Endothelial cells",
  "Tenocytes",
  "Macrophages",
  "Neutrophils",
  "T and NK cells",
  "Schwann cells",
  "MuSC Progenitors",
  "Peryctes",
  "Basophils",
  "B-cells",
  "Lymphatic endothelial cells",
  "Adipocytes"
)

# First Viz
celltype_levels <- names(features)
seuObject_slim_nodoub$celltype <- factor(seuObject_slim_nodoub$celltype, levels = celltype_levels)
features_vec <- unname(unlist(features))
features_vec <- features_vec[!duplicated(features_vec)]
pdf("Modified_Macrophages_dotplo.pdf", width = 6, height = 10)
DotPlot(seuObject_slim_nodoub, features = features_vec, scale = T, group.by = "celltype") +  
        theme_bw() + theme(axis.text.x  = element_text(size = 14),  axis.text.y  = element_text(size = 14),axis.title.x = element_blank(),axis.title.y = element_blank()) +  RotatedAxis() +
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +scale_colour_viridis(option="magma") +
       guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + coord_flip()
dev.off()

# Second Viz
p<-DotPlot(seuObject_slim_nodoub, features = features_vec, scale = TRUE, group.by = "celltype") +  
  theme_bw() +
  theme(axis.text.x  = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),axis.text.y  = element_text(size = 14),
    axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "top",legend.direction = "horizontal",legend.box = "horizontal") +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
  scale_colour_gradient2(low = "#5CACDB",mid = "white",high = "#EA7FA3",midpoint = 0) +
  guides(
    size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"),nrow = 1, byrow = TRUE),
    colour = guide_colourbar(direction = "horizontal")) + coord_flip()

leg <- get_legend(
  p + theme(legend.position = "top",legend.direction = "horizontal",legend.box = "horizontal",legend.text = element_text(size = 9),legend.title = element_text(size = 9)))

p_noleg <- p + theme(legend.position = "none")
pdf("Modified_Macrophages_dotplo.pdf", width = 4, height = 10)
p_noleg
dev.off()

pdf("Legend_Modified_Macrophages_dotplo.pdf", width = 5, height = 2)
grid::grid.newpage()
grid::grid.draw(leg)
dev.off()

#' 
## ----FeaturePlots Viz-----------------------------------------------------------------------------------------------------
fp <- list()
featlist <- c("C1qa", "C1qc", "Ccr2", "Arg1")
for (feat in featlist) {
    fp[[feat]] <- FeaturePlot(seuObject_slim_nodoub, features = feat,
cols = c("#CFD1D3", "#8077A3"), min.cutoff = 0, max.cutoff = 6, order = T) + theme_prism() + NoAxes()
    if (feat != "Arg1") {
        fp[[feat]] <- fp[[feat]] + NoLegend()}}

p <- cowplot::plot_grid(plotlist = fp, nrow = 1, rel_widths = c(1, 1, 1, 1.4))
pdf("Modified_FeaturePlot.pdf", height = 4, width = 10)
p
dev.off()

# save object
save(seuObject_slim_nodoub, file = "output/Seurat_noboud.RData")

#' 
## ----Cell Proportions-----------------------------------------------------------------------------------------------------
Idents(seuObject_slim_nodoub)<-  seuObject_slim_nodoub$celltype
nClust <- length(unique(Idents(seuObject_slim_nodoub)))
colCls <- colorRampPalette(brewer.pal(n = 10, name = "Paired"))(nClust)

plotTheme <- theme(
  text = element_text(size = 18),
  axis.text.x = element_text(angle = 45, hjust = 1),
  axis.text.y = element_text(size = 16),
  legend.title = element_blank(),
  panel.grid = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(color = "black")
)
ggData = data.frame(prop.table(table(seuObject_slim_nodoub$celltype, seuObject_slim_nodoub$Sample), margin = 2))
colnames(ggData) = c("cluster", "sample", "value")
p1 <- ggplot(ggData, aes(sample, value, fill = cluster)) +
  geom_col() + xlab("Sample") + ylab("Proportion of Cells (%)") +
  scale_fill_manual(values = colCls)+ plotTheme + coord_flip()
p1
ggsave(p1, 
       width = 10, height = 5, filename = "Treatment_proportions.jpg")

save(seuObject_slim_nodoub, file = "Seurat_noboud.RData")




#' 
#' #-------------------------
#' ## Subset MuSC Progenitors
#' #-------------------------
## -------------------------------------------------------------------------------------------------------------------------
Idents(seuObject_slim_nodoub)<- seuObject_slim_nodoub$celltype
DefaultAssay(seuObject_slim_nodoub)<- "RNA"

# Subset
Musc <- subset(seuObject_slim_nodoub, idents = "MuSC Progenitors")
# Pre-process
source("Utils_v2.R")
DefaultAssay(Musc)<-"RNA"
Musc<-  processing_seurat_sctransform(Musc, vars_to_regress = c("nUMI", "pMito", "pRibo"), npcs = 30, res = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))

dir.create("MuSC/SCT_UMAP_resolutions/")
resolutions <- seq(0.1, 1, by = 0.1)
for (res in resolutions) {
  res_col <- paste0("SCT_snn_res.", res)
  Idents(Musc) <- Musc[[res_col]][,1]
  p <- DimPlot(Musc, reduction = "umap", label = TRUE) +ggtitle(paste("Resolution", res)) +theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(filename = paste0("MuSC/SCT_UMAP_resolutions/UMAP_res_", res, ".png"),plot = p,width = 8, height = 6, dpi = 300)
  print(p)}

DefaultAssay(Musc) <- "RNA"
Musc <- NormalizeData(object = Musc,normalization.method = "LogNormalize",scale.factor = 10000)

# Presto Markers
#-------------------
wilcoxauc.Seurat <- function(X,group_by = NULL,assay = "data",groups_use = NULL, seurat_assay = "RNA",
    ...
) {
    requireNamespace("Seurat")
    X_matrix <- Seurat::GetAssayData(X, assay = seurat_assay, layer = assay)
    if (is.null(group_by)) {
        y <- Seurat::Idents(X)
    } else {
        y <- Seurat::FetchData(X, group_by) %>% unlist %>% as.character()
    }
    wilcoxauc(X_matrix, y, groups_use)
}
#all_markers_clustID <- wilcoxauc.Seurat(Musc, group_by ='SCT_snn_res.0.5')
all_markers_clustID <- wilcoxauc.Seurat(Musc, group_by ='celltype1')
all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")
all_markers.Sign <- all_markers_clustID %>%dplyr::filter(padj < 0.05, logFC > 0)
top20 <- presto::top_markers(all_markers.Sign,n = 20,auc_min = 0.5, pval_max = 0.05)

# Saving 
openxlsx::write.xlsx(all_markers.Sign,
                     file = "MuSC/celtype_PrestoByCluster_Filteredmarkers_padjLT05_logfcGT0.xlsx",
                     colNames = TRUE,rowNames = FALSE,borders = "columns",sheetName="Markers")
openxlsx::write.xlsx(top20,
                     file = "MuSC/0.3_PrestoByCluster_Top20.xlsx",colNames = TRUE,
                     rowNames = FALSE,borders = "columns",sheetName="Markers")


Idents(Musc)<- Musc$SCT_snn_res.0.5
DimPlot(Musc, split.by = "Sample", group.by = "SCT_snn_res.0.5", label = T)
FeaturePlot(Musc, features = "Cxcl1", min.cutoff = "q10", max.cutoff = "q95", cols = c("gray90", "purple4"),
            split.by = "Sample", label = T)

## Cluster Main Annotation 
cluster_annotations <- c(
  "0" = "Activated MuSC",
  "1" = "Quiescent MuSC", 
  "2" = "Cachexia MuSC")
Idents(Musc)<-Musc$SCT_snn_res.0.5
Musc@meta.data$celltype1<- cluster_annotations[as.character(Idents(Musc))]
table(Musc$celltype1, Musc$Sample)
# factor
Musc$celltype1 <- factor(Musc$celltype1, levels = c("Quiescent MuSC","Activated MuSC", "Cachexia MuSC"))

## Cluster Second main  for publication purposes 
cluster_annotations <- c(
  "0" = "Cluster 2",
  "1" = "Cluster 1", 
  "2" = "Cluster 3")
Idents(Musc)<-Musc$SCT_snn_res.0.5
Musc@meta.data$Pubcelltype1<- cluster_annotations[as.character(Idents(Musc))]
table(Musc$Pubcelltype1, Musc$Sample)

#-----------------------
## Making better UMAP
Idents(Musc)<- Musc$Pubcelltype1
plotTheme <- theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 20, face = "plain"),   
    axis.text  = element_text(size = 15, face = "plain"),   
    axis.line  = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks = element_line(linewidth = 0.6, colour = "black"),
    legend.position = "right",
    plot.title = element_text(size = 16, face = "plain", hjust = 0.5),
    plot.margin = margin(8, 8, 8, 8, "pt"))

# Define colors (one per celltype)
ct_cols <- c(
  "Cluster 2"   = as.character(paletteer_d("rcartocolor::Teal")[5]),
  "Cluster 1"   = as.character(paletteer_d("rcartocolor::Vivid")[5]),
  "Cluster 3"   = as.character(paletteer_d("rcartocolor::Teal")[7]))

pdf("MuSC/V2New_cellytype_General_umap.pdf", width = 8, height = 6)
DimPlot(Musc, cols = ct_cols, pt.size = 2) + plotTheme
dev.off()

pdf("MuSC/V2Sample_cellytype_General_umap.pdf", width = 10, height = 6)
DimPlot(Musc, cols = ct_cols, pt.size = 1.2, split.by = "Sample") + plotTheme +
  theme(strip.text = element_text(size = 15, face = "bold"))
dev.off()

#---------
# Dotplot 
# First Viz
features <- list(
  `Cluster 1` = c("Spry1", "Hes1", "Ckm", "Jun", "Fos", "Zfp36l2", "Txnip"),
  `Cluster 2` = c("Ifrd1", "Hspa5", "Cebpb", "Errfi1", "Tnfrsf12a", "Cdkn1a"),
  `Cluster 3` = c("Slc39a14", "Mt1", "Mt2", "Cebpd", "Stat3")
)
celltype_levels <- names(features)
Musc$Pubcelltype1 <- factor(Musc$Pubcelltype1, levels = celltype_levels)
features_vec <- unname(unlist(features))
features_vec <- features_vec[!duplicated(features_vec)]
pdf("MuSC/V2Modified_Macrophages_dotplo.pdf", width = 4.5, height = 8)
DotPlot(Musc, features = features_vec, scale = T, group.by = "Pubcelltype1") +  
        theme_bw() + theme(axis.text.x  = element_text(size = 14),  axis.text.y  = element_text(size = 14),axis.title.x = element_blank(),axis.title.y = element_blank()) +  RotatedAxis() +
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +scale_colour_viridis(option="magma") +
       guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + coord_flip()
dev.off()

pdf("MuSC/V3Modified_Macrophages_dotplo.pdf", width = 4.5, height = 8)
DotPlot(Musc, features = features_vec, scale = T, group.by = "Pubcelltype1", cols = "RdYlBu") +  
        theme_bw() + theme(axis.text.x  = element_text(size = 14),  axis.text.y  = element_text(size = 14),axis.title.x = element_blank(),axis.title.y = element_blank()) +  RotatedAxis() +
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
       guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + coord_flip()
dev.off()

#--------------------------------
## Modified dotplot for publication paper:
featlist <- c(
  "Hes1","Pax7", "Myf5", "Spry1", "Junb", "Lgals1", 
  "Tdtom", "Ifrd1", "Cdkn1a", "Tnfrsf12a","Kdm6b", "Tubb6","Myog", 
  "C4b", "Tnfrsf1a", "Myod1", "Ccl2", "Slc39a14","Stat3", "Cxcl1", "Cfh", "Serping1","Gpx3","Igfbp5")

pdf("MuSC/V4Modified_Macrophages_dotplo.pdf", width = 10, height = 3.5)
DotPlot(Musc, features = featlist, scale = T, group.by = "Pubcelltype1", cols = "RdYlBu") +  
        theme_bw() + theme(axis.text.x  = element_text(size = 14),  axis.text.y  = element_text(size = 14),axis.title.x = element_blank(),axis.title.y = element_blank(), legend.position = "bottom") +  RotatedAxis() +
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
       guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

#' 
## ----MuSC GO--------------------------------------------------------------------------------------------------------------
Idents(Musc)<- "celltype1"
DefaultAssay(Musc)<- "RNA"
## DEGS
DEGs_MUSC <- Reduce("rbind",lapply(unique(Musc$celltype1), function(x) {
    Markers <- FindMarkers(Musc, ident.1 = x, ident.2 = NULL, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.3, pseudocount.use = 0.1)
    Markers <- Markers[which(Markers$p_val_adj < 0.05),]
    Markers$gene <- rownames(Markers)
    Markers$Cluster <- rep(paste("Cluster", x),nrow(Markers))
    return(Markers)}))

library(scToppR)
toppData <-  toppFun(DEGs_MUSC,topp_categories = NULL, cluster_col = "Cluster", gene_col = "gene", pval_cutoff = 0.05,p_val_col = "p_val_adj",logFC_col = "avg_log2FC", min_genes = 10, max_genes = 500, max_results = 50 )

toppPlot(toppData, category = "GeneOntologyMolecularFunction", num_terms = 10, p_val_adj = "BH", p_val_display = "log", save = TRUE, save_dir = "MuSC/pseudobulk_moderate_GO", width = 5, height = 6)

toppPlot(toppData, category = "GeneOntologyBiologicalProcess", num_terms = 10, p_val_adj = "BH", p_val_display = "log", save = TRUE, save_dir = "MuSC/pseudobulk_moderate_GO", width = 5, height = 6)
# Dataframe of output
GO <- as.data.frame(toppData)
write.xlsx(GO, file = "MuSC/pseudobulk_moderate_GO/GO.xlsx")
save(toppData, file = "MuSC/pseudobulk_moderate_GO/GO.RData")

#' 
## ----MuSC Feature Plots---------------------------------------------------------------------------------------------------
DefaultAssay(Musc)<- "RNA"
FeaturePlot(Musc, features = c("Pax7", "Myf5", "Tdtom","Myod1",
                               "Myog","Cfh", "Ier3", "Slc39a14",
                               "Gpx3", "Lgals1"), min.cutoff = "q10", max.cutoff = "q95", 
            cols = c("gray90", "purple4"), label = T, ncol = 5)

## -------------------------------------------------------------------------------------------------------------------------
fp <- list()
featlist <- c("Pax7", "Myf5", "Tdtom","Myod1", "Myog","Cfh", "Ier3", "Slc39a14", "Gpx3", "Lgals1")

# Cutoff for gene expression
#-----------------------------
for (gene in featlist) {
  cat("\n", gene, "\n")
  print(
    quantile(
      FetchData(Musc, vars = gene)[,1],
      probs = c(0, 0.9, 0.95, 0.99, 1),
      na.rm = TRUE
    )
  )
}

for (feat in featlist) {
    fp[[feat]] <- FeaturePlot(Musc, features = feat,
cols = c("gray90", "purple4"), min.cutoff = "q5", max.cutoff = "q99", order = T) + theme_prism() + NoAxes()
    #if (feat != "Lgals1") {
    #    fp[[feat]] <- fp[[feat]] + NoLegend()}}
    if (feat != tail(featlist, 1)) {
        fp[[feat]] <- fp[[feat]] + NoLegend()}}

#p <- cowplot::plot_grid(plotlist = fp, nrow = 1, rel_widths = c(1, 1, 1, 1.4))
#p <- cowplot::plot_grid(plotlist = fp,nrow = 1,rel_widths = rep(1, length(fp)))
p <- cowplot::plot_grid(plotlist = fp, ncol = 5)
pdf("Modified_FeaturePlot.pdf", height = 4, width = 10)
p
dev.off()

#' 
## ----modified featureplots for publication--------------------------------------------------------------------------------
library(cowplot)

fp <- list()

featlist <- c(
  "Pax7", "Myf5", "Tdtom", "Myod1", "Myog","Spry1", "Hes1", "Ifrd1", "Cdkn1a", "Tnfrsf12a",
  "Kdm6b", "Tubb6", "Lgals1", "Cfh", "Serping1","C4b", "Tnfrsf1a", "Ccl2", "Slc39a14", "Gpx3","Igfbp5",
  "Stat3", "Cxcl1", "Junb")

fp <- list()
for (feat in featlist) {
  fp[[feat]] <- FeaturePlot(Musc,features = feat,cols = c("gray90", "red4"),
    min.cutoff = "q05",max.cutoff = "q99",order = TRUE, pt.size = 0.6) +
    theme_prism() +NoAxes() +NoLegend()}

# Create legend
legend_plot <- FeaturePlot(Musc, features = "Lgals1",cols = c("gray90", "red4"),
  min.cutoff = "q05",max.cutoff = "q99",order = TRUE) +
  theme_prism() +NoAxes()
legend <- cowplot::get_legend(legend_plot)

# Main grid
p_main <- cowplot::plot_grid(plotlist = fp,ncol = 6)
# Combine grid + legend
p <- cowplot::plot_grid(p_main,legend,ncol = 2,rel_widths = c(1, 0.08))

pdf("MuSC/v3FeaturePlots.pdf", height = 9, width = 15)
p
dev.off()

#' 
## ----MuSC Violin Plots----------------------------------------------------------------------------------------------------
featlist <- c(
  "Pax7", "Myf5", "Tdtom", "Myod1", "Myog","Spry1", "Hes1", "Ifrd1", "Cdkn1a", "Tnfrsf12a",
  "Kdm6b", "Tubb6", "Lgals1", "Cfh", "Serping1","C4b", "Tnfrsf1a", "Ccl2", "Slc39a14", "Gpx3","Igfbp5",
  "Stat3", "Cxcl1", "Junb")

cols_2 <- c(
  as.character(paletteer::paletteer_d("rcartocolor::Vivid")[1:12]),
  as.character(paletteer::paletteer_d("rcartocolor::Antique")[1:8]),
  as.character(paletteer::paletteer_d("rcartocolor::Safe")[1:4]))

#Viz
gene_cols <- setNames(cols_2, featlist)
p <- VlnPlot(Musc, features = featlist, stack = TRUE, flip = TRUE, pt.size = 0.05) +theme_classic(base_size = 15) + geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.15,linewidth = 0.25, colour = "black") +
  scale_fill_manual(values = gene_cols) + NoLegend() +theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),  panel.spacing = unit(0.15, "lines")) +  theme(strip.text.y = element_text(angle = 0)) + RotatedAxis()

pdf("MuSC/V3_vlnplot.pdf", width = 4.5, height = 10)
p
dev.off()

#' 
## -------------------------------------------------------------------------------------------------------------------------
save(Musc, file = "MuSC/MuSC.RData")

#' 
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
# [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] scds_1.24.0                 BiocSingular_1.26.1         scDblFinder_1.22.0          leiden_0.4.3.1             
#  [5] muscat_1.22.0               broom_1.0.13                magrittr_2.0.5              speckle_1.8.0              
#  [9] ggpubr_0.6.3                BiocParallel_1.44.0         scater_1.36.0               scuttle_1.18.0             
# [13] emmeans_2.0.3               RColorBrewer_1.1-3          pheatmap_1.0.13             presto_1.0.0               
# [17] data.table_1.18.4           Rcpp_1.1.1-1.1              clustree_0.5.1              ggraph_2.2.2               
# [21] viridis_0.6.5               viridisLite_0.4.3           ggprism_1.0.7               ggrastr_1.0.2              
# [25] ggrepel_0.9.8               paletteer_1.7.0             hdf5r_1.3.12                SeuratWrappers_0.4.0       
# [29] knitr_1.51                  openxlsx_4.2.8.1            RCurl_1.98-1.19             cowplot_1.2.0              
# [33] scales_1.4.0                Matrix_1.7-3                lubridate_1.9.5             forcats_1.0.1              
# [37] stringr_1.6.0               dplyr_1.2.1                 purrr_1.2.2                 readr_2.2.0                
# [41] tidyr_1.3.2                 tibble_3.3.1                ggplot2_4.0.3               tidyverse_2.0.0            
# [45] Seurat_5.5.0                SeuratObject_5.4.0          sp_2.2-1                    SingleCellExperiment_1.30.1
# [49] SummarizedExperiment_1.38.1 Biobase_2.68.0              GenomicRanges_1.60.0        GenomeInfoDb_1.44.3        
# [53] IRanges_2.44.0              S4Vectors_0.48.0            BiocGenerics_0.54.1         generics_0.1.4             
# [57] MatrixGenerics_1.20.0       matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] R.methodsS3_1.8.2        dichromat_2.0-0.1        progress_1.2.3           goftest_1.2-3           
#   [5] Biostrings_2.78.0        TH.data_1.1-5            vctrs_0.7.3              spatstat.random_3.5-0   
#   [9] digest_0.6.39            png_0.1-9                corpcor_1.6.10           shape_1.4.6.1           
#  [13] deldir_2.0-4             parallelly_1.47.0        MASS_7.3-65              reshape2_1.4.5          
#  [17] httpuv_1.6.17            foreach_1.5.2            withr_3.0.2              xfun_0.55               
#  [21] survival_3.8-3           memoise_2.0.1            ggbeeswarm_0.7.3         Seqinfo_1.0.0           
#  [25] gtools_3.9.5             zoo_1.8-15               GlobalOptions_0.1.4      pbapply_1.7-4           
#  [29] R.oo_1.27.1              Formula_1.2-5            prettyunits_1.2.0        rematch2_2.1.2          
#  [33] promises_1.5.0           otel_0.2.0               httr_1.4.8               restfulr_0.0.16         
#  [37] rstatix_0.7.3            globals_0.19.1           fitdistrplus_1.2-6       rstudioapi_0.18.0       
#  [41] UCSC.utils_1.4.0         miniUI_0.1.2             curl_7.1.0               ScaledMatrix_1.18.0     
#  [45] polyclip_1.10-7          GenomeInfoDbData_1.2.14  SparseArray_1.10.8       xtable_1.8-8            
#  [49] doParallel_1.0.17        evaluate_1.0.5           S4Arrays_1.10.1          hms_1.1.4               
#  [53] irlba_2.3.7              colorspace_2.1-2         ROCR_1.0-12              reticulate_1.46.0       
#  [57] spatstat.data_3.1-9      lmtest_0.9-40            later_1.4.8              lattice_0.22-7          
#  [61] spatstat.geom_3.8-1      future.apply_1.20.2      XML_3.99-0.23            scattermore_1.2         
#  [65] RcppAnnoy_0.0.23         pillar_1.11.1            nlme_3.1-168             iterators_1.0.14        
#  [69] caTools_1.18.3           compiler_4.5.1           beachmat_2.26.0          RSpectra_0.16-2         
#  [73] stringi_1.8.7            tensor_1.5.1             minqa_1.2.8              GenomicAlignments_1.44.0
#  [77] plyr_1.8.9               BiocIO_1.18.0            crayon_1.5.3             abind_1.4-8             
#  [81] blme_1.0-7               locfit_1.5-9.12          graphlayouts_1.2.3       bit_4.6.0               
#  [85] sandwich_3.1-1           codetools_0.2-20         multcomp_1.4-30          GetoptLong_1.1.1        
#  [89] plotly_4.12.0            remaCor_0.0.20           mime_0.13                splines_4.5.1           
#  [93] circlize_0.4.18          fastDummies_1.7.6        here_1.0.2               clue_0.3-68             
#  [97] lme4_2.0-1               listenv_0.10.1           Rdpack_2.6.6             ggsignif_0.6.4          
# [101] estimability_1.5.1       statmod_1.5.2            tzdb_0.5.0               fANCOVA_0.6-1           
# [105] tweenr_2.0.3             pkgconfig_2.0.3          tools_4.5.1              cachem_1.1.0            
# [109] RhpcBLASctl_0.23-42      rbibutils_2.4.1          numDeriv_2016.8-1.1      fastmap_1.2.0           
# [113] rmarkdown_2.31           grid_4.5.1               ica_1.0-3                Rsamtools_2.24.1        
# [117] patchwork_1.3.2          coda_0.19-4.1            BiocManager_1.30.27      dotCall64_1.2           
# [121] carData_3.0-6            RANN_2.6.2               farver_2.1.2             reformulas_0.4.4        
# [125] aod_1.3.3                tidygraph_1.3.1          mgcv_1.9-3               yaml_2.3.12             
# [129] rtracklayer_1.68.0       cli_3.6.6                lifecycle_1.0.5          uwot_0.2.4              
# [133] glmmTMB_1.1.14           mvtnorm_1.4-1            sessioninfo_1.2.4        bluster_1.18.0          
# [137] backports_1.5.1          timechange_0.4.0         gtable_0.3.6             rjson_0.2.23            
# [141] ggridges_0.5.7           progressr_0.19.0         pROC_1.19.0.1            parallel_4.5.1          
# [145] limma_3.66.0             jsonlite_2.0.0           edgeR_4.8.2              RcppHNSW_0.7.0          
# [149] bitops_1.0-9             xgboost_3.2.1.1          bit64_4.8.2              Rtsne_0.17              
# [153] spatstat.utils_3.2-3     BiocNeighbors_2.2.0      zip_2.3.3                metapod_1.16.0          
# [157] dqrng_0.4.1              spatstat.univar_3.2-0    R.utils_2.13.0           pbkrtest_0.5.5          
# [161] lazyeval_0.2.3           shiny_1.13.0             htmltools_0.5.9          sctransform_0.4.3       
# [165] rappdirs_0.3.4           glue_1.8.1               spam_2.11-4              XVector_0.50.0          
# [169] rprojroot_2.1.1          scran_1.36.0             gridExtra_2.3            EnvStats_3.1.0          
# [173] boot_1.3-31              igraph_2.3.2             variancePartition_1.38.1 TMB_1.9.21              
# [177] R6_2.6.1                 DESeq2_1.48.2            gplots_3.3.0             cluster_2.1.8.2         
# [181] nloptr_2.2.1             DelayedArray_0.36.0      tidyselect_1.2.1         vipor_0.4.7             
# [185] ggforce_0.5.0            car_3.1-5                future_1.70.0            rsvd_1.0.5              
# [189] KernSmooth_2.23-26       S7_0.2.2                 htmlwidgets_1.6.4        ComplexHeatmap_2.26.0   
# [193] rlang_1.2.0              spatstat.sparse_3.2-0    spatstat.explore_3.8-1   lmerTest_3.2-1          
# [197] remotes_2.5.0            beeswarm_0.4.0          
