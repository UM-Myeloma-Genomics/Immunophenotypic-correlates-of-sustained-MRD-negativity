# T cell receptor beta sequencing methods
# David Coffey (davidcoffey@miami.edu)
# March 1, 2023

# Import R libraries
library(LymphoSeq)

# Step 1: Load data
trb = readImmunoSeq(path = "Path/To/ImmunoSeqV2/")
features = read.csv("scRNAseq features.csv") # Data frame of sample IDs, patient IDs, transplant history, time point, and MRD status

# Compute reportoire statistics
metrics = clonality(trb)

# Step 2: Perform differential abundance between sustained and not sustained MRD- patients

# Merge samples by sustained MRD- status
Sustained.mrd = mergeFiles(file.list = trb, samples = as.character(features[features$SustainedMRD == "Yes", "SampleID"]))
NotSustained.mrd = mergeFiles(file.list = trb, samples = as.character(features[features$SustainedMRD == "No", "SampleID"]))
mrd = c(Sustained = list(Sustained.mrd), NotSustainedMRD = list(NotSustained.mrd))
mrd.aa = productiveSeq(mrd, prevalence = TRUE)

common.mrd = commonSeqs(productive.aa = mrd.aa, samples = names(mrd.aa))
Sustained.mrd.only = mrd.aa$Sustained[!(mrd.aa$sustained$aminoAcid %in% common.mrd$aminoAcid),]
NotSustained.mrd.only = mrd.aa$NotSustainedMRD[!(mrd.aa$NotSustainedMRD$aminoAcid %in% common.mrd$aminoAcid),]

# Merge samples from the same patient
patients = unique(features$PatientID)
trb.patients = list()
i = 1
for(i in 1:length(patients)){
  samples = features[features$PatientID == patients[i], "SampleID"]
  if(length(files) > 0){
    patient = list(mergeFiles(file.list = trb, samples))
    names(patient) = patients[i]
    trb.patients = c(trb.patients, patient)
  }
}

# Create top frequency matrix from merged samples
trb.patients.aa = productiveSeq(trb.patients, prevalence = TRUE)
unique.seqs = uniqueSeqs(trb.patients.aa)
sequence.matrix.patients = seqMatrix(productive.aa = trb.patients.aa, sequences = unique.seqs$aminoAcid)
top.freq.patients = topFreq(productive.aa = trb.patients.aa, percent = 0)
top.freq.matrix.patients = merge(top.freq.patients, sequence.matrix.patients)

# Identify recurrent TCR β sequences (> 2 out of 4 sustained MRD- patients, 75%) present in patients with sustained MRD- only (fisher.test(matrix(c(3,1,0,9), nrow = 2)) p = 0.014)
top.freq.matrix.sustained = top.freq.matrix.patients[top.freq.matrix.patients$aminoAcid %in% Sustained.mrd.only$aminoAcid & top.freq.matrix.patients$numberSamples > 2,]

# Identify recurrent TCR β sequences (> 6 out of 9 not sustained MRD- patients, 78%) present in patients without sustained MRD- only (fisher.test(matrix(c(0,4,7,2), nrow = 2)) p = 0.021)
top.freq.matrix.NotSustained = top.freq.matrix.patients[top.freq.matrix.patients$aminoAcid %in% NotSustained.mrd.only$aminoAcid & top.freq.matrix.patients$numberSamples > 6,]