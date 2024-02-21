
#Load libraries
library(Seurat)
library(tidyverse)
library(dplyr)

#Colon sample
#Read in results/matrix from cell ranger alignment
data_dir <- '/Users/andson/Desktop/gdTcells/COL/filtered_feature_bc_matrix/'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
COL.object = CreateSeuratObject(counts = expression_matrix)

DefaultAssay(COL.object)<- "RNA"

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
COL.object[["percent.mt"]] <- PercentageFeatureSet(COL.object, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(COL.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.0001)

COL.object <- subset(COL.object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 4)

COL.object <- NormalizeData(COL.object, normalization.method = "LogNormalize", scale.factor = 10000)


COL.object <- FindVariableFeatures(COL.object, selection.method = "vst", nfeatures = 2000)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(COL.object), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(COL.object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(COL.object)
COL.object <- ScaleData(COL.object, features = all.genes)

COL.object <- RunPCA(COL.object, features = VariableFeatures(object = COL.object))

ElbowPlot(COL.object)

COL.object <- FindNeighbors(COL.object, dims = 1:12)
COL.object <- FindClusters(COL.object, resolution = c(0.1, 0.5))

COL.object <- RunUMAP(COL.object, dims = 1:12)


Idents(COL.object)<- "RNA_snn_res.0.1"

DimPlot(COL.object, reduction = "umap")


VlnPlot(object = COL.object, features = c("TRDV2"))



# run sctransform
COL.object <- SCTransform(COL.object, vars.to.regress = "percent.mt", verbose = FALSE)

#COL.object <- NormalizeData(COL.object, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

# Visualization and clustering
COL.object <- RunPCA(COL.object, verbose = FALSE)
ElbowPlot(COL.object)
COL.object <- RunUMAP(COL.object, dims = 1:15, verbose = FALSE)

COL.object <- FindNeighbors(COL.object, dims = 1:15, verbose = FALSE)
COL.object <- FindClusters(COL.object, resolution =  c(0.3, 0.4, 0.5, 1))

library(clustree)

clustree(COL.object, prefix = "SCT_snn_res.")

#Idents(COL.object)<- "SCT_snn_res.0.3"
Idents(COL.object)<- "SCT_snn_res.0.3"

DimPlot(COL.object, label = TRUE)



#Which are T cells

#Which are Nk cells
FeaturePlot(COL.object, 
            reduction = "umap", 
            features = c("GNLY", "NK"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

#Which are Nk cells
FeaturePlot(COL.object, 
            reduction = "umap", 
            features = c("TRDV2", "TRGV9"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(COL.object, 
            reduction = "umap", 
            features = c("TRAC", "TRBC1", "TRBC2"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(COL.object, 
            reduction = "umap", 
            features = c("KLRB1", "TRAV1-2"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
FeaturePlot(COL.object, 
            reduction = "umap", 
            features = c("FCGR3A", "EOMES","TBX21", "CD8A", "CD8B"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(COL.object, 
            reduction = "umap", 
            features = c("FCGR3A", "EOMES","TBX21", "CD8A", "CD8B"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

