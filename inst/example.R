library('neuronorm')
library('oro.nifti')


# Get general folder
folder <- system.file("extdata", package = "neuronorm")

# Get covariates
covariates <- system.file("covariates.txt", package = "neuronorm")

# Read covariates information
clinical_info <- read.csv(file = covariates, sep = ';')

# Preprocess MRI scans: 'N4' inhomogeneity correction, 'SyN' non-linear transformation to MNI152 atlas template
# Brain extraction, Spatial informed MRI scan , a.k.a., brain segmentation and RAVEL intensity normalization.
paths_preprocess_patients <- preprocess_patients(folder, clinical_info)

# Outputs paths of the preprocessed MRI scans.
paths_preprocess_patients

# visualize a preprocessed MRI scan for a patient.
img <- readNIfTI(file.path(paths_preprocess_patients$patient04$ravel))
orthographic(img)
