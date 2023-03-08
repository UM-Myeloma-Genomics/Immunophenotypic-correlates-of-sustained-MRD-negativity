# Single-cell VDJ sequencing methods
# David Coffey (davidcoffey@miami.edu)
# March 1, 2023

# Import R libraries
library(data.table)
library(stringr)
library(tidyverse)
library(plyr)

# Merge all contig annotations
paths = list.files("Path/To/VDJ/", full.names = TRUE, pattern = "all_contig_annotations.csv", recursive = TRUE)
all.contigs = plyr::llply(paths, data.table::fread, data.table = FALSE, .progress = "text")
all.contigs = plyr::ldply(all.contigs, data.frame, .id = "SampleID")

# Create contig file of clonotypes
clonotype.contigs = all.contigs[all.contigs$high_confidence == TRUE & 
                                  all.contigs$is_cell == TRUE & 
                                  all.contigs$raw_clonotype_id != "None" &
                                  all.contigs$productive != "None" & 
                                  all.contigs$raw_consensus_id != "None" &
                                  all.contigs$chain != "None" &
                                  all.contigs$cdr3 != "", ]

# Create a matrix with cells as rows and VDJ chains as columns
aggregate.chains = aggregate(data = clonotype.contigs, cdr3~SampleID+chain+barcode, paste, collapse = ",")
vdj.matrix = pivot_wider(data = aggregate.chains, names_from = chain, values_from = cdr3, values_fill = "")

# Mark multi-chain receptors (defined as a TCR or BCR with more than one chain of the same type)
vdj.matrix$multi_chain = ifelse(grepl(vdj.matrix$TRB, pattern = ",") | grepl(vdj.matrix$TRA, pattern = ","), TRUE, FALSE)

# Mark dual receptors
vdj.matrix$dualAlpha = ifelse(grepl(vdj.matrix$TRA, pattern = ","), TRUE, FALSE)
vdj.matrix$dualBeta = ifelse(grepl(vdj.matrix$TRB, pattern = ","), TRUE, FALSE)
vdj.matrix$dualIGH = ifelse(grepl(vdj.matrix$IGH, pattern = ","), TRUE, FALSE)
vdj.matrix$dualIGKL = ifelse(grepl(vdj.matrix$IGK, pattern = ",") | grepl(vdj.matrix$IGL, pattern = ","), TRUE, FALSE)

# Mark paired receptors
vdj.matrix$pairedTCR = ifelse(vdj.matrix$TRA != "" & vdj.matrix$TRB != "", TRUE, FALSE)
vdj.matrix$pairedBCR = ifelse((vdj.matrix$IGK != "" | vdj.matrix$IGL != "") & vdj.matrix$IGH != "", TRUE, FALSE)

# Mark doublets (defined as a single cell with a BCR and TCR chain)
vdj.matrix$doublet = ifelse((vdj.matrix$SampleID %in% doublets$Var1 & vdj.matrix$barcode %in% doublets$Var2) |
                                     (vdj.matrix$IGH != "" | vdj.matrix$IGK != "" | vdj.matrix$IGL != "") &
                                     (vdj.matrix$TRB != "" | vdj.matrix$TRA != ""), TRUE, FALSE)

# Filter clones (exclude doublets, unpaired, and multi-chain receptors)
vdj.matrix.filter = vdj.matrix[(vdj.matrix$pairedTCR == TRUE | vdj.matrix$pairedBCR == TRUE) &
                                       vdj.matrix$multi_chain == FALSE & vdj.matrix$doublet == FALSE,]

# Create unique identifier for each clonotype
vdj.matrix.filter$CDR3 = ifelse(vdj.matrix.filter$pairedBCR == TRUE,
                            paste("BCR:", vdj.matrix.filter$IGH, "-", vdj.matrix.filter$IGK, vdj.matrix.filter$IGL, sep = ""),
                            paste("TCR:", vdj.matrix.filter$TRB, "-", vdj.matrix.filter$TRA, sep = ""))