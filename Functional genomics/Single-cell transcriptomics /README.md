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

Install Seurat - https://satijalab.org/seurat/articles/install_v5
#### Set up the Seurat object

For this tutorial, we will be analyzing a dataset of Peripheral Blood Mononuclear Cells (PBMC) before and after interferon (IFN-b) stimulation. The raw data can be found [here](https://www.ncbi.nlm.nih.gov/gds/302560248) and [here](https://www.ncbi.nlm.nih.gov/gds/302560249), control and ifn, respectively.

Install via browser direct download, or

In linux: 
```
## Control
# 1. Download the Matrix file (.mtx.gz)
wget https://ftp.ncbi.nlm.nih.gov/geo/conditions/GSM2560nnn/GSM2560248/suppl/GSM2560248_2.1.mtx.gz

# 2. Download the Barcodes file (.tsv.gz)
wget https://ftp.ncbi.nlm.nih.gov/geo/conditions/GSM2560nnn/GSM2560248/suppl/GSM2560248_barcodes.tsv.gz

# 3. Download the Genes/Features file (from Series level GSE96583)
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96583/suppl/GSE96583_batch2.genes.tsv.gz

```
```
# 1. Download the Matrix file (.mtx.gz)
wget https://ftp.ncbi.nlm.nih.gov/geo/conditions/GSM2560nnn/GSM2560249/suppl/GSM2560248_2.2.mtx.gz

# 2. Download the Barcodes file (.tsv.gz)
wget https://ftp.ncbi.nlm.nih.gov/geo/conditions/GSM2560nnn/GSM2560249/suppl/GSM2560249_barcodes.tsv.gz

# 3. Genes/Features is common
```
Rename files for easier conversion

```
mv GSM2560248_2.1.mtx.gz matrix.mtx.gz
mv GSM2560248_barcodes.tsv.gz barcodes.tsv.gz
mv GSE96583_batch2.genes.tsv.gz features.tsv.gz
```
In windows: 
```
Invoke-WebRequest -Uri "same/link/to/file" -OutFile "matrix.mtx.gz"
Invoke-WebRequest -Uri "same/link/to/file" -OutFile "barcodes.tsv.gz"
Invoke-WebRequest -Uri "same/link/to/file" -OutFile "features.tsv.gz"
```

The alternative links for files
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2560248
[https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2560249)

```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2560nnn/GSM2560248/suppl/GSM2560248_2.1.mtx.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2560nnn/GSM2560248/suppl/GSM2560248_barcodes.tsv.gz

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96583/suppl/GSE96583_batch2.genes.tsv.gz

wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2560nnn/GSM2560249/suppl/GSM2560249_2.2.mtx.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2560nnn/GSM2560249/suppl/GSM2560249_barcodes.tsv.gz
```
We start by reading in the data. 
We next use the count matrix to create a `Seurat` object. The object serves as a container that contains both data (like the count matrix) and analysis (like PCA, or clustering results) for a single-cell dataset. For more information, check out the [docs for SeuratObject](https://satijalab.github.io/seurat-object/index.html) or the [section on object interaction](https://satijalab.org/seurat/articles/essential_commands#seurat-object-data-access) in our list of essential commands. For example, in Seurat v5, the count matrix is stored in `seurat[["RNA"]]$counts`.

```
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Matrix)

# Load the downloaded PBMC datasets
input_dir <- "~/Documents/Rstudio/ASD_Glial/first_sc_tutorial"
counts <- readMM(file.path(input_dir, "matrix.mtx.gz"))

# Read genes and barcodes
genes <- read.table(file.path(input_dir, "features.tsv.gz"), header = FALSE, sep = "\t")
barcodes <- read.table(file.path(input_dir, "barcodes.tsv.gz"), header = FALSE, sep = "\t")
```
Now, genes file contains duplicated human gene symbols (HGNC), so we first have to make them unique
```
> length((genes$V2))
[1] 35635
> length(unique(genes$V2))
[1] 32938
> genes |>
+   janitor::get_dupes(V2)
              V2 dupe_count              V1
1    AC012314.20         10 ENSG00000223660
2    AC012314.20         10 ENSG00000273386
```
```
genes$V2 <- make.unique(genes$V2)
rownames(counts) <- genes$V2  # or genes$V1 for Ensembl IDs
colnames(counts) <- barcodes$V1
```
```
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = counts, project = "GSM2560248")

# Inspect the file
head(rownames(pbmc))
head(colnames(pbmc))
meta <- pbmc@meta.data
head(meta)

```
Repeat the steps for the interferon-stimulated dataset (`pbmc_ifn`).
```
#### Load the downloaded PBMC datasets (stim/interferon condition)
counts_ifn <- readMM(file.path(input_dir, "matrix.mtx.gz"))

#### Read genes and barcodes
genes_ifn <- read.table(file.path(input_dir, "features.tsv.gz"), header = FALSE, sep = "\t")
barcodes_ifn <- read.table(file.path(input_dir, "barcodes.tsv.gz"), header = FALSE, sep = "\t")

#### Make gene symbols unique (same duplicate HGNC issue as control)
genes_ifn$V2 <- make.unique(genes_ifn$V2)
rownames(counts_ifn) <- genes_ifn$V2  # or genes_ifn$V1 for Ensembl IDs
colnames(counts_ifn) <- barcodes_ifn$V1

#### Initialize the Seurat object with the raw (non-normalized) data
pbmc_ifn <- CreateSeuratObject(counts = counts_ifn, project = "GSM2560249")
```
```
print(pbmc)
An object of class Seurat 
35635 features across 14619 conditions within 1 assay 
Active assay: RNA (35635 features, 0 variable features)
 1 layer present: counts
print(pbmc_inf)
An object of class Seurat 
35635 features across 14446 conditions within 1 assay 
Active assay: RNA (35635 features, 0 variable features)
 1 layer present: counts

```
The genes in the matrix are rows and cells columns.
```
# Add `condition` as a metadata column to differentiate the conditions later
pbmc$condition <- "Control"
pbmc_ifn$condition <- "Stim"
```
You can also visually inspect the counts matrix.
```
> pbmc[["RNA"]]$counts[1:5, 1:5]
```
```
5 x 5 sparse Matrix of class "dgCMatrix"
             AAACATACAATGCC-1 AAACATACATTTCC-1 AAACATACCAGAAA-1 AAACATACCAGCTA-1 AAACATACCATGCA-1
MIR1302-10                  .                .                .                .                .
FAM138A                     .                .                .                .                .
OR4F5                       .                .                .                .                .
RP11-34P13.7                .                .                .                .                .
RP11-34P13.8                .                .                .                .                .
```
Now we continue with the workflow onto the QC stage. In general the entire pre-processing pipeline looks like this: 
```
Individual Datasets 
          │
          ▼
QC & Filtering (Mito %, nFeature_RNA, nCount_RNA on each condition)
          │
          ▼
Doublet Removal (Per condition)
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
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Plot a feature scatter plot
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p <- plot1 + plot2
p

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)`

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)`
```
Do the same for the second Seurat object.

There are specific tools that are used to classify and remove doublets from single cell datasets. One of them is [`DoubletFinder`](https://github.com/chris-mcginnis-ucsf/DoubletFinder). It's an R package that integrates with Seurat, imulates arftificial doublets based on average gene expression and lastly classifies real cells into doublets and singlets based on proportions of artificial doublets in their vicinity in PCA-reduced space. Due to time constraints, we are not going to use DoubletFinder in this tutorial, but its a nice tool nonetheless. 

#### Merging

Before proceeding to normalisation, scaling and variable features extraction, we need to merge the two Seurat objects into one.
```
seurat_merged <- merge(
  x = pbmc,
  y = pbmc_ifn,
  add.cell.ids = c("Control", "STIM"),
  project = "PBMC"
)
```
```
> head(seurat_merged@meta.data)
> seurat_merged

 ```
 Merging does not equal integration! Integration methods such as Harmony,(which we are using in this tutorial), iteratively cluster cells from different datasets with respect to user-defined penalties (depend on the datasets) essentially regressing batch effects and ensuring that cells group together by biology rather than technical artefacts.

For more information on Harmony algorithm refer to this [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6884693/).
#### Normalization
By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. In Seurat v5, normalized values are stored in `seurat[["RNA"]]$data`.

```
seurat_merged <- NormalizeData(seurat_merged, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_merged[["RNA"]]$data
```

#### Ensure metadata contains a condition/batch identifier column (e.g., seurat_merged$condition)
```
table(seurat_merged$condition)
```

#### Since the data is pre-merged, check if layers are split by condition.
```
seurat_merged[["RNA"]] <- split(seurat_merged[["RNA"]], f = seurat_merged$condition)
```
#### Identify variable features across split layers (uses existing normalized data)
```
seurat_merged <- FindVariableFeatures(seurat_merged, selection.method = "vst", nfeatures = 2000)
```
####  Optional: Visualize top features
```
top10 <- head(VariableFeatures(seurat_merged), 10)
top10

plot1 <- VariableFeaturePlot(seurat_merged)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
```
#### Scale data across the identified variable features prior to PCA
```
seurat_merged <- ScaleData(seurat_merged, features = VariableFeatures(seurat_merged))
```
#### Run unintegrated PCA
```
seurat_merged <- RunPCA(seurat_merged, features = VariableFeatures(object = seurat_merged), reduction.name = "pca")
```
#### Check unintegrated structure (to confirm if batch effects exist)
```
seurat_merged <- RunUMAP(seurat_merged, dims = 1:20, reduction.name = "umap.unintegrated")
DimPlot(seurat_merged, reduction = "umap.unintegrated", group.by = "condition") + 
  ggtitle("Unintegrated UMAP (Check Batch Effects)")
```
#### Seurat v5 CCA Anchor-based Integration (Default)
```
seurat_merged <- IntegrateLayers(
  object = seurat_merged,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)
```
#### Re-join layers post-integration for unified downstream differential expression
```
seurat_merged[["RNA"]] <- JoinLayers(seurat_merged[["RNA"]])
```
#### Target reduction choice (switch to "integrated.harmony" if using Harmony)
```
target_reduction <- "integrated.cca" 
```
#### Graph building and clustering on integrated space
```
seurat_merged <- FindNeighbors(seurat_merged, reduction = target_reduction, dims = 1:20)
seurat_merged <- FindClusters(seurat_merged, resolution = 0.5)
```
#### Run UMAP on integrated space
```
seurat_merged <- RunUMAP(seurat_merged, reduction = target_reduction, dims = 1:20, reduction.name = "umap")
```
#### Visualizations
```
DimPlot(seurat_merged, reduction = "umap", group.by = "condition") + ggtitle("Integrated by condition")
DimPlot(seurat_merged, reduction = "umap", label = TRUE) + ggtitle("Integrated Clusters")
DimPlot(seurat_merged, reduction = "umap", split.by = "condition", label = TRUE)
```
#### Save processed output
```
saveRDS(seurat_merged, file = "data/seurat_merged_normalized_integrated.rds")
```
#### Paper reading 

Homework 1 - to read the paper about the general guidelines for scRNA-seq analysis, prepare a short summary report \
 ``` Practice/hw1.pdf  ``` \
Homework 2 - to read the paper about best practices in scRNA-seq analysis, prepare a short summary report \
``` Practice/hw2.pdf ``` 

### Useful links
[scRNA-seq analysis book](https://www.sc-best-practices.org/introduction/prior_art.html)

