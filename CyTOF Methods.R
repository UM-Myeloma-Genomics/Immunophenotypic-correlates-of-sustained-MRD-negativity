# CyTOF methods
# David Coffey (davidcoffey@miami.edu)
# March 1, 2023

# Import R libraries
library(orloj)
library(tidyverse)

# Step 1: Load data and perform QC

# Load astrolabe data
experiment = orloj::loadExperiment("Astrolabe/")

# Extract cell counts
cell.counts = orloj::experimentCellSubsetCounts(experiment = experiment, level = "Assignment")

# Create matrix with samples as rows and cell types as columns
count.matrix = pivot_wider(data = cell.counts[,c("CellSubset", "Name", "N")], names_from = CellSubset, values_from = N, values_fill = 0)
count.matrix = column_to_rownames(count.matrix, var = "Name")

# Remove cell types where there are less than 3 cells in at least half of samples
min_samples = ncol(count.matrix)/2
tf = count.matrix >= 3
ix_keep = apply(tf, 1, function(r) sum(r) >= min_samples)
count.matrix = count.matrix[ix_keep, , drop = FALSE]

# Step 2: Perform differential abundance using EdgeR

# Create design matrix and contrasts
features = read.csv("CyTOF features.csv") # Data frame of sample IDs, patient IDs, transplant history, time point, and MRD status
design = model.matrix(~ 0 + Transplant:Baseline + NoTransplant:Baseline, data = features) # Compares baseline samples from patients with and without transplant
colnames(design) = gsub(colnames(design), pattern = ":", replacement = "")
contrasts = makeContrasts(TransplantBaseline-BaselineNoTransplant, levels = design)

# EdgeR - fit a negative binomial generalized log-linear model to cell abundance
y = DGEList(count.matrix)
y = estimateDisp(y, design, trend.method = "none")
fit = edgeR::glmFit(y, design)
lrt = glmLRT(fit, contrast = contrasts)
results = as.data.frame(edgeR::topTags(lrt, n = Inf, adjust.method = "fdr", sort.by = "PValue"))