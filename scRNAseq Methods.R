# Single-cell RNA sequencing methods
# David Coffey (davidcoffey@miami.edu)
# March 1, 2023

# Import R libraries
library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(Libra)

# Step 1: Load data, perform QC, integrate datasets, and classify cells

# Import gene matrix files and convert to list of Seurat objects
paths = list.files("/Path/To/GeneMatrix/", full.names = TRUE, all.files = TRUE, recursive = FALSE, include.dirs = FALSE)
names = basename(paths)

data.list = list()
i = 1
for(i in 1:length(names)){ 
  mtx = Read10X(data.dir = paths[i])
  mtx = CreateSeuratObject(counts = mtx, project = names[i])
  data.list = c(data.list, list(mtx))
}

# Calculate nCount_RNA and nFeature_RNA
data.list = lapply(X = data.list, FUN = function(x) {
  calcn = as.data.frame(x = Seurat:::CalcN(object = x))
  colnames(x = calcn) = paste(colnames(x = calcn), "RNA", sep = '_')
  x = AddMetaData(object = x, metadata = calcn)
})

# Calculate percent mitochondrial genes 
data.list = lapply(X = data.list, FUN = function(x) {
  x = PercentageFeatureSet(object = x, pattern = '^MT-', col.name = "Percent_mt", assay = "RNA")
})

# Remove cells with < 500 UMI, < 250 genes, or > 20% mitochondrial genes
data.list.filtered = lapply(X = data.list, FUN = function(x) {
  cells.use = x[["nCount_RNA", drop = TRUE]] > 500 | x[["nFeature_RNA", drop = TRUE]] > 250 | x[["Percent_mt", drop = TRUE]] < 20
  x = x[, cells.use]
})

# Normalize
data.list.normalized = lapply(X = data.list.filtered, FUN = function(x) {
  x = SCTransform(x, verbose = FALSE, conserve.memory = TRUE)
})

# Perform feature selection
features = SelectIntegrationFeatures(object.list = data.list.normalized)

# Run PCA
data.list.normalized = lapply(X = data.list.normalized, FUN = function(x) {
  x = RunPCA(x, features = features, verbose = FALSE)
})

# Integrate
anchors.integrate = FindIntegrationAnchors(object.list = data.list.normalized, reduction = "rpca", dims = 1:50)
integrated = IntegrateData(anchorset = anchors.integrate, dims = 1:50)

# Automate cell classification using a labeled reference PBMC scRNAseq dataset
reference = LoadH5Seurat("pbmc_multimodal.h5seurat") # Available for download here: https://azimuth.hubmapconsortium.org

# Load raw, aggregated data output from Cell Ranger pipeline
data = Read10X_h5(filename = "raw_feature_bc_matrix.h5")
data = CreateSeuratObject(counts = data)

# Normalize the data 
data = SCTransform(data, verbose = TRUE)

# Find anchors between reference and query using a precomputed supervised PCA (spca) transformation
anchors = FindTransferAnchors(
  reference = reference,
  query = data,
  normalization.method = "SCT",
  reference.reduction = "spca",
  recompute.residuals = FALSE,
  dims = 1:50)

# MapQuery is a wrapper around three functions: TransferData, IntegrateEmbeddings, and ProjectUMAP. TransferData is used to transfer cell type labels and impute the ADT values. IntegrateEmbeddings and ProjectUMAP are used to project the query data onto the UMAP structure of the reference.
data = MapQuery(
  anchorset = anchors,
  query = data,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)
  
# Step 2: Perform differential abundance using EdgeR

# Extract cell counts
cell.counts = seurat.sce@meta.data[,c("SampleID", "celltype.l2", "Barcode")]
cell.counts = aggregate(data = cell.counts, Barcode~SampleID+celltype.l2, length)

# Create matrix with samples as rows and cell types as columns
count.matrix = pivot_wider(data = cell.counts, names_from = celltype.l2, values_from = Barcode, values_fill = 0)
count.matrix = column_to_rownames(count.matrix, var = "SampleID")

# Create design matrix and contrasts
features = read.csv("scRNAseq features.csv") # Data frame of sample IDs, patient IDs, transplant history, time point, and MRD status
design = model.matrix(~ 0 + Transplant:Baseline + NoTransplant:Baseline, data = features) # Compares baseline samples from patients with and without transplant
colnames(design) = gsub(colnames(design), pattern = ":", replacement = "")
contrasts = makeContrasts(TransplantBaseline-BaselineNoTransplant, levels = design)

# EdgeR - fit a negative binomial generalized log-linear model to cell abundance
y = DGEList(count.matrix)
y = estimateDisp(y, design, trend.method = "none")
fit = edgeR::glmFit(y, design)
lrt = glmLRT(fit, contrast = contrasts)
differential.abundance = as.data.frame(edgeR::topTags(lrt, n = Inf, adjust.method = "fdr", sort.by = "PValue"))

# Step 3: Perform differential gene expression using EdgeR likelihood ratio testing
differential.expression = Libra::run_de(subset(x = seurat.sce, subset = TimePoint == "Baseline"), replicate_col = "SampleID", cell_type_col = "celltype.l2", label_col = "MRDNegative") # Compares baseline samples from patients with and without sustained MRD negativity

# Step 4: Define exhausted T cells

# Subset terminal effector memory T cells
Idents(seurat.sce) = "celltype.l2"
tem.cells = subset(x = seurat.sce, idents = c("CD4 TEM", "CD8 TEM"))

# Define markers
markers = c("TIGIT", "LAG3", "PDCD1", "CXCL13", "LAYN", "TOX", "CTLA4", "HAVCR2", "ITGAE")

# Convert gene expression to binary (1 = cell expression marker, 0 = cell does not express marker)
counts = GetAssayData(object = tem.cells)
counts@x[counts@x > 0] = 1
counts = counts[rownames(counts) %in% markers, ]
counts = data.frame(t(as.data.frame.matrix(counts)))
counts$Barcode = rownames(counts)