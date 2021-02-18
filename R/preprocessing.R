#' Preprocess T1-w MRI scan for one patient
#'
#' This function preprocesses a raw T1-w MRI scan and generates a spatially informed MRI data using the fast algorithm.'
#' The preprocesising steps comprises imhomogeneity correction 'N4', registration to the MNI152 template with isotropic voxel size of 2mm
#' using the 'SyN' transformation, and skull stripping.
#'
#' @param mri_patient path of the T1-w scan.
#' @param folder_patient folder containing the T1-w scan. This folder usually refers as the patient.
#' @param atlas atlas template to spatially register the T1-w scans. By default the MNI152 atlas template is used.
#' @param mask brain mask of the atlas template to performed the skull stripping.
#' @param inhomogeneity inhomogeneity correction algorithm to be applied. The correction by default is the 'N4' bias correction.
#' @param trasnformation non-linear transformation for registering the T1-w MRI scan to the reference template. 'SyN' transformation is used by default.
#' @return paths of preprocessed MRI scans.
#' @export
preprocess_modality_t1 <- function(mri_patient, folder_patient, atlas, mask, inhomogeneity = "N4", transformation = "SyN"){

  # empty list to save paths
  mri_paths <- list()

  # Inhomogeneity Correction: N4
  cat(paste0('*********************************************\n****** Inhomogeneity Correction: ', inhomogeneity ,' *********\n*********************************************\n--Running...\n'))
  bias_file <- file.path(folder_patient, 'T1_bias.nii.gz')
  bias_mri <- extrantsr::bias_correct(mri_patient, correction = inhomogeneity, outfile = bias_file, verbose = FALSE)
  mri_paths[['bias']] <- bias_file
  cat('--Complete.\n')

  # Registration to Template (SyN: Non-linear)
  cat(paste0('*********************************************\n******* Spatial Registration : ',transformation , ' **********\n*********************************************\n--Running...\n'))
  syn_file <-file.path(folder_patient, 'T1_SyN.nii.gz')
  syn_mri <- extrantsr::ants_regwrite(filename = bias_mri, outfile = syn_file, template.file = atlas, typeofTransform = transformation, verbose = FALSE)
  mri_paths[['syn']] <- syn_file
  cat('--Complete.\n')

  # Brain Mask
  cat('*********************************************\n*************** Brain Mask ******************\n*********************************************\n--Running...\n')
  mask_file <- file.path(folder_patient, 'T1_masked_SyN')
  mask_mri = neurobase::mask_img(syn_file, mask)
  oro.nifti::writeNIfTI(mask_mri, mask_file)
  mri_paths[['syn_masked']] <- paste0(mask_file,'.nii.gz')
  cat('--Complete.\n')

  # Spatially Informed layer
  cat('*********************************************\n******** Spatially Informed Layer ***********\n*********************************************\n--Running...\n')
  spatial_file <- file.path(folder_patient, 'T1_spatially_informed.nii.gz')
  spatial_mri = fslr::fast(file = mask_mri, outfile = spatial_file, opts = "--nobias", verbose = FALSE)
  mri_paths[['spatial']] <- spatial_file
  cat('--Complete.\n')

  # CSF tissue mask for RAVEL algorithm
  cat('*********************************************\n********* CSF Tissue Mask (RAVEL) ***********\n*********************************************\n--Running...\n')
  tissue_mask <- spatial_mri
  tissue_mask[tissue_mask != 1] <- 0
  tissue_mask_file = file.path(folder_patient, 'T1_CSF_tissue')
  oro.nifti::writeNIfTI(tissue_mask, tissue_mask_file)
  mri_paths[['csf_mask']] <- paste0(tissue_mask_file,'.nii.gz')
  cat('--Complete.\n\n')

  # Path for RAVEL
  ravel_file = file.path(folder_patient, 'T1_Ravel_norm')
  mri_paths[['ravel']] <- paste0(ravel_file ,'.nii.gz')

  return(mri_paths)
}

#' Wrapper function for RAVEL normalization in T1-w images
#'
#' Ravel intensity normalization including control voxels and clinical covariates.'
#' @param masked_paths list or vector of paths of the input NIfTI images to be normalized.
#' @param csf_paths NIfTI image paths for the binary control region masks.
#' @param ravel_paths list or vector of paths of the output NIfTI images.
#' @param demographics table of covariates associated to the MRI scans. Number of rows should be equal to the number of images.
#' @param brain_mask NIfTI image path for the binary brain mask. Must have value 1 for the brain and 0 otherwise.
#' @param patients_folder folder to save the output control mask.
#' @return RAVEL-corrected images are saved in disk.
#' @export
### Ravel Normalization
image_normalization_ravel <- function(masked_paths, csf_paths, ravel_paths, demographics, brain_mask, patients_folder){

  ### Control region mask for all patients (masks intersect)
  mask_intersect_path <- file.path(patients_folder, 'CSF_control_mask.nii.gz')
  intersect_mask <- RAVEL::maskIntersect(csf_paths, output.file = mask_intersect_path , prob = 0.8)

  # Control for biological covariates
  mod_cov <- model.matrix(~., demographics[,-1])

  #RAVEL
  RAVEL::normalizeRAVEL(input.files = masked_paths, output.files = ravel_paths, brain.mask = brain_mask,
                 control.mask = intersect_mask, mod = mod_cov , WhiteStripe_Type	= 'T1',
                 k = 1, returnMatrix = FALSE , writeToDisk = TRUE)
}


#' Preprocess T1-w MRI scan for multiple patients
#'
#' This function preprocesses raw T1-w MRI scans and generates a spatially informed MRI scans using the fast algorithm.'
#' The preprocesising steps comprises imhomogeneity correction 'N4', registration to the MNI152 template with isotropic voxel size of 2mm
#' using the 'SyN' transformation, skull stripping, and RAVEL intensity normalization.
#'
#' @param patients_folder general folder containing folders per patient with raw T1-w images.
#' @param clinical_covariates table of covariates associated to the MRI scans. Number of rows should be equal to the number of images.
#' @return paths of preprocessed MRI scans.
#' @export
preprocess_patients <- function(patients_folder, clinical_covariates){

  # getting patients scans
  patients_mri <- load_mri_group(folder_patients)

  # empty list to save paths
  paths_mri <- list()

  # Atlas template (2mm MNI152)
  mniTemplate <- MNITemplate::getMNIPath("Brain", res="2mm")

  # Atlas brain mask (2mm MNI152)
  brainMask = MNITemplate::readMNI("Brain_Mask", res='2mm')

  # preprocess each patients T1 sequence scan
  for (patient in patients_mri){

    # selecting T1-w scan
    pathT1 <- patient$T1

    # selecting patient's folder
    patient_folder <- dirname(patient$T1)

    # getting folder name
    folder_name = unlist(strsplit(patient_folder, '/'))
    folder_name = folder_name[length(folder_name)]

    # Preprocess T1-w image
    cat(paste0('--------------------------------------------------\n Preprocesing T1-w image in ',folder_name ,'\n--------------------------------------------------\n\n'))
    patient_scans <- preprocess_modality_t1(pathT1, patient_folder, mniTemplate, brainMask)

    # Save preprocessed images paths
    paths_mri[[folder_name]] <- patient_scans
  }

  ### RAVEL
  cat('*********************************************\n*********** RAVEL normalization *************\n*********************************************\n--Running...\n')

  ## Group paths for ravel
  ravel_paths <- lapply(paths_mri, function(x){x$ravel})
  csf_paths <- lapply(paths_mri, function(x){x$csf_mask})
  masked_paths <- lapply(paths_mri, function(x){x$syn_masked})

  ## ravel algorithm
  image_normalization_ravel(masked_paths, csf_paths, ravel_paths, clinical_covariates, brainMask, patients_folder)

  cat(paste0('--------------------------------------------------\n Preprocess Complete \n--------------------------------------------------\n\n'))

  return(paths_mri)
}
