args <- commandArgs(trailingOnly = TRUE)

#
# Set script arguments
#
output_dir <- "" 
mclust_model <- NULL
minjoin <- NULL
ntrial <- NULL
description <- "No description provided by user."
if (length(args) >= 1){
	output_dir <- args[1]
}
if (length(args) >= 2){
	mclust_model <- args[2]	
}
if (length(args) >= 3){
	minjoin <- as.numeric(args[3])
}
if(length(args) >= 4){
	ntrial <- as.numeric(args[4])
}
if (length(args) >= 5){
	description <- args[5]
}

print(paste("PARAMETERS:", "output_dir=", output_dir, "mclust_model=", mclust_model, "minjoin=", minjoin, "ntrial=", ntrial, "description=", description))
#
# Load source libraries
# TODO: Organize dependencies
#
setwd("~/code/cnprep_clustering/scripts")
source("helperFunctions.R")
source("segmentClusteringLibrary2.R")

# TODO: Do not include segments with lower than 5K bp (see paper)

reference <- "hN31"
sample_dir <- "./resources/FACETS_Reference_hN31_7_28_18_2/"

#
# Load input
#
cd_cnprep()
normal_samples <- load_samples(classes = c("N"), sampleList = "./resources/sampleList.csv")
normal_samples <- normal_samples[normal_samples != reference] # Remove reference organoid

cytobands <- retrieveCytobands(dir = "./resources/cytoBand.txt")
chromosomeSizes <- readRDS("./resources/chromosomeSizes.rds")
normalSegments <- selectSegmentsWithEvents(events = c("A", "D", "N"), samples = normal_samples, chromosomeSizes = chromosomeSizes, 
                                           dir = sample_dir, sample_subdir="/", reference=reference, extension = "cnv.facets.v0.5.2.txt", inSampleFolder = TRUE, 
                                           rescaleInput = TRUE, ampCall = 0.2, delCall = -0.235)
target_samples <- load_samples(classes = c("N", "T", "F", "M"), sampleList = "./resources/sampleList.csv")
target_samples <- target_samples[target_samples != reference] # Remove reference organoid

# Generate norminput argument
norminput <- retrieveNormInput(normalSegments)
norminput <- filterNormInput(norminput, length_threshold=10000000)

# Create folder with output
dir.create(file.path("./output/", output_dir))

for(target_samples.i in 1:length(target_samples)) {
  sample <- target_samples[target_samples.i]
  
  print(paste("Analyzing sample", sample))
  
  facets_segment_data <- retrieveFacetsSegments(sample, sample_subdir = "/", reference = reference, dir = sample_dir)
  facets_snp_data <- retrieveFacetsSnps(sample, sample_subdir = "/", reference = reference, dir = sample_dir)
  
  # Generate seginput argument
  seginput <- retrieveSegInput(facets_segment_data, sample, chromosomeSizes, cytobands)
  print(paste("Retrieved segment input for sample", sample))
  
  # Generate ratinput argument
  ratinput <- retrieveRatInput(facets_snp_data, sample)
  print(paste("Retrieved ratio input for sample", sample))
  
  # Run CNprep:CNpreprocessing
  try({
  	segtable <- runCNpreprocessing(seginput = seginput, ratinput = ratinput, norminput = norminput, modelNames = mclust_model, minjoin = minjoin, ntrial = ntrial) #TODO: Is there a distrib="Grid"?
	  print(paste("Produced segtable for sample", sample))
  
	  write.table(segtable, paste("./output/", output_dir,"/", sample, "_segtable.tsv", sep = ""), row.names = F, sep = "\t", quote = FALSE)
	  print(paste("Wrote output for sample", sample))
  }, silent = TRUE)
}

# TODO: Do this for all HPC scripts. May need to make function in helperFunctions.R to reduce duplicate code
info_filename <- paste("./output/", output_dir, "/JobInformation.txt", sep = "")
file.create(info_filename)
fileConn <- file(info_filename)
writeLines(	c("UGE JOB SUBMISSION NOTES FOR segmentClusteringHPC.R",
	 	paste("User description of job: ", description),
		"The script centers the input segments and clusters the segments using GMM",
		paste("The output files wrote to: ", output_dir),
		paste("This input parameters used mclust_model =", mclust_model, "and minjoin =", minjoin, "and ntrial=", ntrial)),
	 fileConn)
close(fileConn)
