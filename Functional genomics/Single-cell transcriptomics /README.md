# Contributors
    Primary contributor: Maria Nikoghosyan
    Contributing authors: Melina Tamazyan, Ekaterina Kostiuk
  # Contents  
  ### Lecture 
  [Lecture slides](https://docs.google.com/presentation/d/1kXkgaH0Ym0tvMoHXpu6AIuAgSQb9s06g/edit?slide=id.g37061a333be_1_0#slide=id.g37061a333be_1_0) 
  
### Practice 

The tutorials are based on Satija lab tutorials on [scRNA-seq pipeline](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) and on [integration](https://satijalab.org/seurat/articles/integration_introduction.html)

Alternative single cell processing pipeline with scanpy for python enthusiasts (for future reference):
[tutorial](https://scanpy.scverse.org/en/stable/tutorials/basics/clustering.html)

#### Set up the Seurat object

For this tutorial, we will be analyzing a dataset of Peripheral Blood Mononuclear Cells (PBMC) from a healthy donor freely available from 10X Genomics. The raw data can be found [here]([https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz](https://www.10xgenomics.com/datasets/Flexv2_16k_Human_PBMCs_TotalseqC_4plex)), samples (A01, B01).

Install via browser direct download, or

In linux: 
```
curl -O https://cf.10xgenomics.com/samples/cell-exp/10.0.0/16k_Human_PBMCs_TotalseqC_4plex_PBMC_TotalSeqC_C-A01/16k_Human_PBMCs_TotalseqC_4plex_PBMC_TotalSeqC_C-A01_sample_raw_feature_bc_matrix.h5
```
```
wget https://cf.10xgenomics.com/samples/cell-exp/10.0.0/16k_Human_PBMCs_TotalseqC_4plex_PBMC_TotalSeqC_C-A01/16k_Human_PBMCs_TotalseqC_4plex_PBMC_TotalSeqC_C-A01_sample_raw_feature_bc_matrix.h5 . 
```
In windows: 
```
Invoke-WebRequest -Uri "https://cf.10xgenomics.com/samples/cell-exp/10.0.0/16k_Human_PBMCs_TotalseqC_4plex_PBMC_TotalSeqC_C-A01/16k_Human_PBMCs_TotalseqC_4plex_PBMC_TotalSeqC_C-A01_sample_raw_feature_bc_matrix.h5" -OutFile "16k_Human_PBMCs_TotalseqC_4plex_PBMC_TotalSeqC_C-A01_sample_raw_feature_bc_matrix.h5"
```
This will download the file into your current directory.

We start by reading in the data. This more recent versions of cellranger outputs the [h5 file format](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices), which can be read in using the `[Read10X_h5()](https://satijalab.org/seurat/reference/read10x_h5)` function in Seurat.

We next use the count matrix to create a `Seurat` object. The object serves as a container that contains both data (like the count matrix) and analysis (like PCA, or clustering results) for a single-cell dataset. For more information, check out the [docs for SeuratObject](https://satijalab.github.io/seurat-object/index.html) or the [section on object interaction](https://satijalab.org/seurat/articles/essential_commands#seurat-object-data-access) in our list of essential commands. For example, in Seurat v5, the count matrix is stored in `seurat[["RNA"]]$counts`.

```
library(dplyr)
library(Seurat)
library(patchwork)
library(hdf5r)
library(ggplot2)

# Load the downloaded PBMC datasets
counts <- Read10X_h5("~/Documents/Rstudio/ASD_Glial/first_sc_tutorial/16k_Human_PBMCs_TotalseqC_4plex_PBMC_TotalSeqC_C-A01_sample_raw_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
seurat_obj <- CreateSeuratObject(counts = counts, project = "PBMC_10X")
# Inspect the file
head(rownames(seurat_obj))
head(colnames(seurat_obj))
meta <- seurat_obj@meta.data
head(meta)

```

```
print(seurat_obj)
An object of class Seurat 
18714 features across 422298 samples within 1 assay 
Active assay: RNA (18714 features, 0 variable features)
 2 layers present: counts.Gene Expression, counts.Antibody Capture
```
```
head(rownames(seurat_obj))

   [1] "TSPAN6"    "TNMD"      "DPM1"      "SCYL3"     "C1orf112"

head(colnames(seurat_obj))
   [1] "AAACCAATCAACCGGAAGCCCGAACT-1" "AAACCAATCAACCTTGAGCCCGAACT-1"
   [3] "AAACCAATCAACTGAAAGCCCGAACT-1" "AAACCAATCAAGCCAGAGCCCGAACT-1"
```
The genes in the matrix are rows and cells columns.
```
# Add sample_id as a metadata column to differentiate the samples later
seurat_obj$sample_id <- "A01"
```
Then since our object contains two layers in the same RNA assay, it is best to split them into individual assays.
```
seurat_obj[["ADT"]] <- CreateAssayObject(counts = seurat_obj[["RNA"]]$`counts.Antibody Capture`)
seurat_obj[["RNA"]] <- CreateAssayObject(counts = seurat_obj[["RNA"]]$`counts.Gene Expression`)
DefaultAssay(seurat_obj) <- "RNA"
print(seurat_obj[["RNA"]])
```
```
> print(seurat_obj[["RNA"]])Assay data with 18530 features for 422298 cells 
First 10 features: TSPAN6, TNMD, DPM1, SCYL3, C1orf112, FGR, CFH, FUCA2, GCLC, NFYA
```
Now repeat the process for the second sample. Don't forget to give it a different name and assign `sample_id`. You can also visually inspect the counts matrix.
```
> seurat_obj2[["RNA"]]$counts[1:5, 1:5]
```
```
5 x 5 sparse Matrix of class "dgCMatrix"
         `AACCAATCAACCGTCAGCCCTGCAT-1 AAACCAATCAACTATCAGCCCTGCAT-1 AAACCAATCAACTGAAAGCCCTGCAT-1
TSPAN6                              .                            .                            .
TNMD                                .                            .                            .
DPM1                                .                            .                            .
SCYL3                               .                            .                            .
C1orf112                            .                            .                            .
         AAACCAATCAAGCCAGAGCCCTGCAT-1 AAACCAATCAAGCGTGAGCCCTGCAT-1
TSPAN6                              .                            .
TNMD                                .                            .
DPM1                                .                            .
SCYL3                               .                            .
C1orf112                            .                            .

```
Now we continue with the workflow onto the QC stage. In general the entire pre-processing pipeline looks like this: 
```
Individual Datasets 
          │
          ▼
QC & Filtering (Mito %, nFeature_RNA, nCount_RNA on each sample)
          │
          ▼
Doublet Removal (Per sample)
          │
          ▼
Merging Objects 
          │
          
Normalize & Scale (NormalizeData / SCTransform)
          │
          ▼
Integration (Harmony, CCA, etc.)
```

#### QC

A few QC metrics [commonly used](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4758103/) by the community include:
- The number of unique genes detected in each cell.
    - Low-quality cells or empty droplets will often have very few genes
    - Cell doublets or multiplets may exhibit an aberrantly high gene count
- Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
- The percentage of reads that map to the mitochondrial genome
    - Low-quality / dying cells often exhibit extensive mitochondrial contamination
    - We calculate mitochondrial QC metrics with the `[PercentageFeatureSet()](https://satijalab.org/seurat/reference/percentagefeatureset)` function, which calculates the percentage of counts originating from a set of features
    - We use the set of all genes starting with `MT-` as a set of mitochondrial genes
```
#Calculate MT%
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
#Visualize QC metrics as a violin plot
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Plot a feature scatter plot
plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p <- plot1 + plot2
p
```
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)`

VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)`

Do the same for the second Seurat object.

There are specific tools that are used to classify and remove doublets from single cell datasets. One of them is [`DoubletFinder`](https://github.com/chris-mcginnis-ucsf/DoubletFinder). It's an R package that integrates with Seurat, imulates arftificial doublets based on average gene expression and lastly classifies real cells into doublets and singlets based on proportions of artificial doublets in their vicinity in PCA-reduced space. Due to time constraints, we are not going to use DoubletFinder in this tutorial, but its a nice tool nonetheless. 

#### Merging

Before proceeding to normalisation, scaling and variable features extraction, we need to merge the two Seurat objects into one.
```
seurat_merged <- merge(
  x = seurat_obj,
  y = seurat_obj2,
  add.cell.ids = c("A01", "B01"),
  project = "PBMC_4plex"
)
```
```
> head(seurat_merged@meta.data)
> seurat_merged
An object of class Seurat
18714 features across 752 samples within 2 assays
Active assay: RNA (18530 features, 0 variable features)
2 layers present: counts, data
1 other assay present: ADT
 ```
 Merging does not equal integration! Integration methods such as Harmony,(which we are using in this tutorial), iteratively cluster cells from different datasets with respect to user-defined penalties (depend on the datasets) essentially regressing batch effects and ensuring that cells group together by biology rather than technical artefacts.

For more information on Harmony algorithm refer to this [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6884693/).
#### Normalization
By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. In Seurat v5, normalized values are stored in `seurat[["RNA"]]$data`.

```
seurat_merged <- NormalizeData(seurat_merged, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_merged[["RNA"]]$data
```


# Ensure metadata contains a sample/batch identifier column (e.g., seurat_merged$sample)
table(seurat_merged$sample)

# Since the data is pre-merged, check if layers are split by sample.
# For integration to work accurately, variable features need to be calculated 
# per batch/sample.

# If layers are currently joined, split the RNA assay by sample:
seurat_merged[["RNA"]] <- split(seurat_merged[["RNA"]], f = seurat_merged$sample)

# Identify variable features across split layers (uses existing normalized data)
seurat_merged <- FindVariableFeatures(seurat_merged, selection.method = "vst", nfeatures = 2000)

# Optional: Visualize top features
top10 <- head(VariableFeatures(seurat_merged), 10)
plot1 <- VariableFeaturePlot(seurat_merged)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale data across the identified variable features prior to PCA
seurat_merged <- ScaleData(seurat_merged, features = VariableFeatures(seurat_merged))

# Run unintegrated PCA
seurat_merged <- RunPCA(seurat_merged, features = VariableFeatures(object = seurat_merged), reduction.name = "pca")

# Check unintegrated structure (to confirm if batch effects exist)
seurat_merged <- RunUMAP(seurat_merged, dims = 1:20, reduction.name = "umap.unintegrated")
DimPlot(seurat_merged, reduction = "umap.unintegrated", group.by = "sample") + 
  ggtitle("Unintegrated UMAP (Check Batch Effects)")


### Option A: Seurat v5 CCA Anchor-based Integration (Default)
seurat_merged <- IntegrateLayers(
  object = seurat_merged,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)

### Option B: Harmony Integration (Alternative fast approach)
# library(harmony)
# seurat_merged <- IntegrateLayers(
#   object = seurat_merged,
#   method = HarmonyIntegration,
#   orig.reduction = "pca",
#   new.reduction = "integrated.harmony",
#   verbose = FALSE
# )

# Re-join layers post-integration for unified downstream differential expression
seurat_merged[["RNA"]] <- JoinLayers(seurat_merged[["RNA"]])

# Target reduction choice (switch to "integrated.harmony" if using Harmony)
target_reduction <- "integrated.cca" 

# Graph building and clustering on integrated space
seurat_merged <- FindNeighbors(seurat_merged, reduction = target_reduction, dims = 1:20)
seurat_merged <- FindClusters(seurat_merged, resolution = 0.5)

# Run UMAP on integrated space
seurat_merged <- RunUMAP(seurat_merged, reduction = target_reduction, dims = 1:20, reduction.name = "umap")

# Visualizations
DimPlot(seurat_merged, reduction = "umap", group.by = "sample") + ggtitle("Integrated by Sample")
DimPlot(seurat_merged, reduction = "umap", label = TRUE) + ggtitle("Integrated Clusters")
DimPlot(seurat_merged, reduction = "umap", split.by = "sample", label = TRUE)

# CRITICAL RULE: Perform marker checks and DE on original unintegrated RNA expression values!
DefaultAssay(seurat_merged) <- "RNA"

# 5.1 Evaluate canonical markers
markers.to.check <- c("CD3D", "MS4A1", "NKG7", "CD14", "LYZ", "FCGR3A", "PPBP")
FeaturePlot(seurat_merged, features = markers.to.check, reduction = "umap")
DotPlot(seurat_merged, features = markers.to.check) + RotatedAxis()

# 5.2 Find cluster markers across all clusters
seurat_merged.markers <- FindAllMarkers(
  seurat_merged, 
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25
)

# Save processed output
saveRDS(seurat_merged, file = "data/seurat_merged_normalized_integrated.rds")

#### Paper reading 

Homework 1 - to read the paper about the general guidelines for scRNA-seq analysis, prepare a short summary report \
 ``` Practice/hw1.pdf  ``` \
Homework 2 - to read the paper about best practices in scRNA-seq analysis, prepare a short summary report \
``` Practice/hw2.pdf ``` 

### Useful links
[scRNA-seq analysis book](https://www.sc-best-practices.org/introduction/prior_art.html)

