
## ----------------------- Installation  --------------------------

# install.packages("devtools")
devtools::install_github("DavidPayares/neuronorm")

## ----------------------- Data Loading  --------------------------

# Get general folder
folder <- system.file("extdata", package = "neuronorm")

# Get covariates
covariates <- system.file("covariates.txt", package = "neuronorm")

# Read covariates information
clinical_info <- read.csv(file = covariates, sep = ';')

## ----------------------- Preprocessing  --------------------------

require('neuronorm')
# Preprocess MRI scans: 'N4' inhomogeneity correction, 'SyN' non-linear transformation to MNI152 atlas template
# Brain extraction, Spatial informed MRI scan , a.k.a., brain segmentation and RAVEL intensity normalization.
paths_preprocess_patients <- preprocess_patients(folder, clinical_info)

# Outputs paths of the preprocessed MRI scans.
paths_preprocess_patients

## ----------------------- Image Visualization  --------------------------

require('oro.nifti')
# visualize a preprocessed MRI scan for a patient.
img <- readNIfTI(file.path(paths_preprocess_patients$patient04$ravel))
orthographic(img)


#---------------------------------------------------------------------
folder <- '/mnt/d/M.Sc. Gesopatial Tecnologies/Thesis/Resources/Code/neuronorm/inst/extdata'

patients_mri <- load_mri_group(folder)

folder.patient <- '/mnt/d/M.Sc. Gesopatial Tecnologies/Thesis/Resources/Code/neuronorm/inst/extdata/patient01/'
modalities <- c('T1','T2','FLAIR')
mri.patient <- patients_mri$patient01
transformation = "SyN"
inhomogeneity = "N4"
atlas <- MNITemplate::getMNIPath()
mask <- MNITemplate::readMNI("Brain_Mask")

preprocess_modalities(mri.patient, folder.patient, atlas, mask, inhomogeneity = "N4", transformation = "SyN")

preprocess_modalities <- function(mri.patient, folder.patient, atlas, mask, inhomogeneity = "N4", transformation = "SyN"){

  # empty list to save paths
  mri_paths <- list()

  # Inhomogeneity Correction: N4
  cat(paste0('*********************************************\n****** Inhomogeneity Correction: ', inhomogeneity ,' *********\n*********************************************\n--Running...\n'))
  bias_files <- lapply(modalities, function(x) file.path(paste0(folder.patient, x, '_bias.nii.gz')))
  bias_mri <- mapply(extrantsr::bias_correct, file = mri.patient, correction = inhomogeneity, outfile = bias_files, verbose = FALSE)
  mri_paths[['bias']]<- bias_files
  cat('--Complete.\n')

  # Coregistration to T1-weighted image : Rigid Transformation
  cat(paste0('*********************************************\n****** Coregistration to T1 sequence *********\n*********************************************\n--Running...\n'))
  coregisteredImg <- extrantsr::within_visit_registration(fixed = mri.patient$T1, moving = bias_files[2:length(mri.patient)], typeofTransform = "Rigid", interpolator = "Linear", verbose = FALSE)
  bias_mri_comp <- coregistration_images(coregisteredImg)
  cat('--Complete.\n')

  # Registration to Template (SyN: Non-linear)
  cat(paste0('*********************************************\n******* Spatial Registration : ',transformation , ' **********\n*********************************************\n--Running...\n'))
  syn_files <- lapply(modalities, function(x) file.path(paste0(folder.patient, x, '_SyN_MNI152.nii.gz')))
  bias_mris <- create_bias_list(modalities, bias_mri$T1, bias_mri_comp)
  syn_mri <- mapply(extrantsr::ants_regwrite, filename = bias_mris, outfile = syn_files, template.file = atlas, typeofTransform = transformation,  verbose = FALSE)
  mri_paths[['syn']] <- syn_files
  cat('--Complete.\n')

  # Brain Mask
  cat('*********************************************\n*************** Brain Mask ******************\n*********************************************\n--Running...\n')
  mask_files <- lapply(modalities, function(x) file.path(paste0(folder.patient, x, '_masked')))
  mask_mri <- lapply(syn_files, neurobase::mask_img, mask)
  mapply( oro.nifti::writeNIfTI, nim = mask_mri, filename = mask_files)
  mri_paths[['syn_masked']] <- paste0(mask_files,'.nii.gz')
  cat('--Complete.\n')

  # Spatially Informed layer
  cat('*********************************************\n******** Spatially Informed Layer ***********\n*********************************************\n--Running...\n')
  spatial_file <- file.path(folder.patient, 'T1_spatially_informed.nii.gz')
  spatial_mri = fslr::fast(file = mask_mri[[1]], outfile = spatial_file, opts = "--nobias", verbose = FALSE, type = 'T1')
  mri_paths[['spatial']] <- file.path(folder.patient, 'T1_spatially_informed_pveseg.nii.gz')
  cat('--Complete.\n')

  # CSF tissue mask for RAVEL algorithm
  cat('*********************************************\n********* CSF Tissue Mask (RAVEL) ***********\n*********************************************\n--Running...\n')
  tissue_mask <- spatial_mri
  tissue_mask[tissue_mask != 1] <- 0
  tissue_mask_file = file.path(folder.patient, 'T1_CSF_tissue')
  oro.nifti::writeNIfTI(tissue_mask, tissue_mask_file)
  mri_paths[['csf_mask']] <- paste0(tissue_mask_file,'.nii.gz')
  cat('--Complete.\n\n')

  # Path for RAVEL
  ravel_files <- lapply(modalities, function(x) file.path(paste0(folder.patient, x, '_RAVEL')))
  mri_paths[['ravel']] <- paste0(ravel_files ,'.nii.gz')

  return(mri_paths)
}


bias_files <- lapply(modalities, function(x) file.path(paste0(folder.patient, x, '_bias.nii.gz')))
bias_mri <- mapply(extrantsr::bias_correct, file = mri.patient, correction = inhomogeneity, outfile = bias_files, verbose = FALSE)
mri_paths[['bias']]<- bias_files

coregisteredImg <- extrantsr::within_visit_registration(fixed = mri.patient$T1, moving = bias_files[2:length(mri.patient)], typeofTransform = "Rigid", interpolator = "Linear", verbose = FALSE)
bias_mri_comp <- coregistration_images(coregisteredImg)

syn_files <- lapply(modalities, function(x) file.path(paste0(folder.patient, x, '_SyN_MNI152.nii.gz')))
bias_mris <- create_bias_list(modalities, bias_mri$T1, bias_mri_comp)


modalities[[1]]
syn_mri <- mapply(extrantsr::ants_regwrite, filename = bias_mris, outfile = syn_files, template.file = atlas, typeofTransform = transformation,  verbose = FALSE)



lala <- 'x'
lolo <- list('m','n')

nini <- list(lala, lolo)
new.pp <- unlist(nini,recursive=FALSE)

