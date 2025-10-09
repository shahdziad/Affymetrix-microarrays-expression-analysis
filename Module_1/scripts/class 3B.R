library(affy)                 # Pre-processing of Affymetrix microarray data (RMA normalization)
library(arrayQualityMetrics)  # QC reports for microarray data
library(dplyr)                # Data manipulation
library(affyio)
library(oligo)

getwd()
setwd("C:/Bioinformatics/AI_Omics_Internship_2025/Module_1/raw_data/E-MTAB-9637-RawData")

# Read sdrf file for phenotype data


sdrf_file <- read.delim(file.choose()) 
sum(is.na(sdrf_file$Factor.Value.compound.)) # check if there are any missing names as it's important for labelling groups next
rownames(sdrf_file) <- sdrf_file$Array.Data.File
SDRF <- AnnotatedDataFrame(sdrf_file)


# Read CEL files into R and create an expression set

raw_data_dir <- "C:\\Bioinformatics\\AI_Omics_Internship_2025\\Module_1\\raw_data\\E-MTAB-9637-RawData"

raw_data <- oligo::read.celfiles(filenames = file.path(raw_data_dir, 
                                                       SDRF$Array.Data.File),
                                 verbose = FALSE, phenoData = SDRF)
stopifnot(validObject(raw_data))



head(Biobase::pData(raw_data))
Biobase::pData(raw_data) <- Biobase::pData(raw_data)[, c("Source.Name",
                                                         "Characteristics.individual.",
                                                         "Array.Data.File",
                                                         "Factor.Value.compound.")]

# Check QC before normalization

arrayQualityMetrics(
  expressionset = raw_data,
  outdir = "C:/Bioinformatics/AI_Omics_Internship_2025/Module_1/Results/QC_Raw_Data",
  force = TRUE,
  do.logtransform = TRUE
)


# Normalize raw data

eset <- oligo::rma(raw_data, target = "core")

# Check QC before normalization

arrayQualityMetrics(
  expressionset = eset,
  outdir = "C:/Bioinformatics/AI_Omics_Internship_2025/Module_1/Results/QC_Normalized_Data",
  force = TRUE)

# Extract normalized expression values into a data frame

processed_data <- as.data.frame(exprs(eset))
dim(processed_data)


# Filter Low-Variance Transcripts
# Calculate median intensity per probe across samples

row_median <- rowMedians(as.matrix(processed_data))

# Visualize distribution of probe median intensities
hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")

# Set a threshold to remove low variance probes (dataset-specific, adjust accordingly)
threshold <- 5 
abline(v = threshold, col = "black", lwd = 2)

# Select probes above threshold
indx <- row_median > threshold 
filtered_data <- processed_data[indx, ] 

dim(filtered_data)

# Overwrite processed data with filtered dataset
processed_data <- filtered_data 

# Check class of sdrf column before 
class(sdrf_file$Factor.Value.compound.)

# Define experimental groups (normal vs cancer)

table(sdrf_file$Factor.Value.compound.)
groups <- factor(sdrf_file$Factor.Value.compound.,
                 levels = c("atorvastatin", "none"),
                 label = c("treated", "diseased"))

class(groups)
levels(groups)

# Rename column names in processed data
colnames(processed_data) <- (sdrf_file$Source.Name)

# Save results

save.image(file = "ShahdZiadMohamedRashad_Class_3B_Assignment.RData")
