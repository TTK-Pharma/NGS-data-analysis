library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)


# step-1 load the 10xgnomics file
glio <- Read10X_h5("data\\Parent_SC3v3_Human_Glioblastoma_raw_feature_bc_matrix.h5")

#step-2create a seurat object to create a centralized storage
glio_object <- CreateSeuratObject(counts = glio, min.cells = 3, min.features = 200)

View(glio_object)
#accessing the information about the features and samples
glio_object

#step-3 first check the mitochondrial gene content in the cells
#in this step we will use thr regex pattern to find the gene with stating MT in it
glio_object[["percent.mt"]] <- PercentageFeatureSet(glio_object, pattern = "^MT-")
View(glio_object@meta.data)
glio_object


#step-4 Plot the qulaity of the data seeing the mess visually
raw_plot <- VlnPlot(glio_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
VlnPlot(glio_object, features = c("nFeature_RNA", "nCount_RNA"))
ggsave(filename = "raw_plot.png", plot = raw_plot, path = "plots")

#step-5 Apply some parameters to clean the data with false positives 
glio_object <- subset(glio_object, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 &
                        nCount_RNA > 1000 & nCount_RNA < 40000 &
                        percent.mt < 5)
glio_object

#normalise the data to make the gene that is comparable across cell (eliminate technical bias)
glio_object <- NormalizeData(glio_object)
norm_data <- glio_object[["RNA"]]$data
norm_data[5:15, 5:15]

#Find the most varaible gene across the cells with the default value 2000

glio_object <- FindVariableFeatures(glio_object,selection.method = "vst", nfeatures = 2000)
# to visualize it use 
variablePlot <- VariableFeaturePlot(glio_object)
ggsave(filename = "variableplot.png", plot = variablePlot, path = "plots")

#before performing the PCA analysis there is a main stp called scaling
#where the normalization is done for difference between cells (sequencing depth)
#scaling is to normalize the difference between genes (highly expressed genes will
#have high no of variations the PCA will onlly care about higher valuse nor for biology)
#only access the genes so create a variable with only gene name

all_genes <- rownames(glio_object)
glio_object <- ScaleData(glio_object, features = all_genes)
dim(glio_object[["RNA"]]$scale.data)

features_present <- intersect(VariableFeatures(glio_object), rownames(glio_object[["RNA"]]$scale.data))
length(features_present)
#here comes the linear dimentionality where the PCA will organize the gene with the 
# name of PC1 which means principle component one the genes that have more variations 
#across cells. This gene makes the cells distinct

glio_object <- RunPCA(glio_object, features = features_present)
PCA_plot<-ElbowPlot(glio_object)
ggsave(filename = "elbowplot.png", plot = PCA_plot, path = "NGS-data-analysis\\plots")

#now find the neighors from the PCA data like if cell A and cell B will get many of the
#genes similar then they will have strong connection
#first do for 1st 10 PCs
glio_object <- FindNeighbors(glio_object, dims = 1:10)

#here comes the final part the identification of cell types by grouping the cells into
#clusters we can later find which cluster belongs to which cell type
#important concept called resolution parameter
#low resolution ranges between 0.1-0.4 results in few clusters but the size is huge
#use this only we want broad categories like malignant cells and non-malignant
#High resolution rage 0.8-1.2 more cluster with small size use this if we want
#specific sub-type of T cells
glio_object <- FindClusters(glio_object, resolution = 0.8)
#next step will show how many cells are assigned in a specific cluster
table(Idents(glio_object))
#here comes the part of non-linear dimentionality UMAP
#used to view the cluster in 2D structure keeping the similar cells nearby 
#and keeping differnet cells far away fro top 10 PCs
glio_object <- RunUMAP(glio_object, dims = 1:10)
#use dimplot to visualize it 
dim_plot <- DimPlot(glio_object, reduction = "umap")
ggsave(filename = "dimplot.png", plot = dim_plot, path = "NGS-data-analysis\\plots")

#here we will find the markers acroos clusters based on the expression of the specific gene


Idents(glio_object) <- "seurat_clusters"
markers <- FindAllMarkers(glio_object, only.pos = TRUE, min.pct = 0.25, 
                          logfc.threshold = 0.25, test.use = "wilcox")
write.csv(markers, file = "markersByclusters.csv")

#use the following commands to view the markers that are upregulated across clusters
head(markers, 20)
markers[15:25, ]
markers[ , 1:4]
slice_min(markers, order_by = avg_log2FC)
slice_max(markers, order_by = avg_log2FC)

markers_List <- read.csv("NGS-data-analysis\\markersByclusters.csv")
head(markers_List)
markers_final <- as_tibble(markers_List, .rows = NULL)
View(markers_final)
#now lests find the genes that are expressed highNULL#now lests find the genes that are expressed highly from each cluster
#use top_n() function
#top_5 <- markers_List %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top_5 <- markers_List %>% filter(avg_log2FC > 1 & p_val_adj < 0.05) %>%
                  group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 5)
head(top_5)
DotPlot(glio_object, features = unique(top_5$gene)) + theme(axis.text.x = element_text(size = 5)) + RotatedAxis()
expressed_plot <- DotPlot(glio_object, features = unique(top_5$gene)) + theme_bw(ink = "white") + theme(axis.text.x = element_text(size = 5)) + RotatedAxis()
ggsave(filename = "expressed_gene5.png", plot = expressed_plot, path = "NGS-data-analysis\\plots")

#now come to the part of cluster annotations
top_5 %>% filter(cluster == 12)

VlnPlot(glio_object, features = c("nFeature_RNA", "nCount_RNA"), group.by="seurat_clusters")
#now visualize it in feature plot for confirmation
FeaturePlot(glio_object, features = c(
  "CX3CR1",
  "TMEM119",
  "P2RY12",
  "CD14"
))

VlnPlot(glio_object, features = c("MBP","PLP1","OLIG2"), group.by = "seurat_clusters")

#Assign the cluster types to specific clusters
levels(Idents(glio_object))
#Rename it 
glio_object <- RenameIdents(glio_object, 
                            "0" = "OPC-like malignant cells",
                            "1" = "mesenchymal-like malignant cells",
                            "2" = "Peripheral macrophages", 
                            "3" = "Hypoxia-associated malignant Cells", 
                            "4" = "Tumor-Associated Macrophages",
                            "5" = "Mature oligodendrocytes",
                            "6" = "OPC-like malignant cells",
                            "7" = "Tumor-associated macrophages",
                            "8" = "Neuronal progenitor-like cells",
                            "9" = "Microglia",
                            "10" = "Proliferating tumor cells", 
                            "11" = "T cells",
                            "12" = "Cancer-associated fibroblasts")
DimPlot(glio_object, label = TRUE, repel = TRUE)
#save this in a separate column 
glio_object$celltypes <- Idents(glio_object)
head(glio_object@meta.data)
Annotated<-DimPlot(glio_object, group.by = "celltypes", label = TRUE)
ggsave(filename = "Annotated-UMAp.png", plot = Annotated, path = "NGS-data-analysis\\plots",
       height = 8, width = 19)
