#' @title Preprocesses T1-weighted MRI scan for one patient
#'
#' @description  This function preprocesses a raw T1-w MRI scan and generates a segmentation MRI scan using the FAST algorithm.
#' The preprocesising steps comprises imhomogeneity correction 'N4', registration to the MNI152 template with isotropic voxel size of 1mm^3
#' using the 'SyN' transformation, and skull stripping.
#'
#' @param mri.patient path of the T1-weighted scan.
#' @param folder.patient folder containing the T1-weighted scan. This folder usually refers to the patient.
#' @param atlas atlas template in NifTI format to spatially register the T1-weighted scans. By default the MNI152 atlas template is used.
#' @param mask brain mask in NifTI format of the atlas template to performed the skull stripping.
#' @param inhomogeneity inhomogeneity correction algorithm to be applied. The correction by default is the 'N4' bias correction.
#' @param transformation non-linear transformation for registering the T1-w MRI scan to the reference template. 'SyN' transformation is used by default.
#' @return paths of preprocessed MRI scans.
#' @importFrom extrantsr bias_correct ants_regwrite
#' @importFrom neurobase mask_img
#' @importFrom fslr fast
#' @importFrom oro.nifti writeNIfTI
#' @export
preprocess_modality_t1 <- function(mri.patient, folder.patient, atlas, mask, inhomogeneity = "N4", transformation = "SyN"){

  # Empty list to save paths
  mri_paths <- list()

  # Inhomogeneity Correction: N4
  cat(paste0('*********************************************\n****** Inhomogeneity Correction: ', inhomogeneity ,' *********\n*********************************************\n--Running...\n'))
  bias_files <- file.path(folder.patient, 'T1_bias.nii.gz')
  bias_mri <- extrantsr::bias_correct(mri.patient, correction = inhomogeneity, outfile = bias_files[[1]], verbose = FALSE)
  mri_paths[['bias']] <- bias_file
  cat('--Complete.\n')

  # Registration to Template (SyN: Non-linear)
  cat(paste0('*********************************************\n******* Spatial Registration : ',transformation , ' **********\n*********************************************\n--Running...\n'))
  syn_file <-file.path(folder.patient, 'T1_SyN.nii.gz')
  syn_mri <- extrantsr::ants_regwrite(filename = bias_mri, outfile = syn_file, template.file = atlas, typeofTransform = transformation, verbose = FALSE)
  mri_paths[['registered']] <- syn_file
  cat('--Complete.\n')

  # Brain Mask
  cat('*********************************************\n*************** Brain Mask ******************\n*********************************************\n--Running...\n')
  mask_file <- file.path(folder.patient, 'T1_masked_SyN')
  mask_mri = neurobase::mask_img(syn_file, mask)
  oro.nifti::writeNIfTI(mask_mri, mask_file)
  mri_paths[['striped']] <- paste0(mask_file,'.nii.gz')
  cat('--Complete.\n')

  # Spatially Informed layer (Segmentation layer)
  cat('*********************************************\n******** Segmentation (HMRF) ***********\n*********************************************\n--Running...\n')
  spatial_file <- file.path(folder.patient, 'T1_spatially_informed.nii.gz')
  spatial_mri = fslr::fast(file = mask_mri, outfile = spatial_file, opts = "--nobias", verbose = FALSE)
  mri_paths[['segmented']] <- file.path(folder.patient, 'T1_spatially_informed_pveseg.nii.gz')
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
  ravel_file = file.path(folder.patient, paste0(folder.patient, '_T1_Ravel_norm'))
  mri_paths[['ravel']] <- paste0(ravel_file ,'.nii.gz')

  return(mri_paths)
}


#' @title Preprocesses group of MRI scan for one patient
#'
#' @description This function preprocesses raw T1-weighted, T2-weighted and FLAIR MRI scans and generates a segmentation MRI scan using the FAST algorithm.
#' The preprocesising steps comprises imhomogeneity correction 'N4', coregistration of other sequences to the T1-weighted scan,
#' non-linear registration to the MNI152 template with an isotropic voxel size of 1mm,
#' using the 'SyN' transformation, skull stripping, brain segmentation and intensity normalization using the RAVEL or White Stripe algorithms.
#'
#' @param mri.patient path of the MRI scans.
#' @param folder.patient folder containing the MRI scans. This folder usually refers to the patient.
#' @param atlas atlas template in NifTI format to spatially register the MRI scans. By default the MNI152 atlas template is used.
#' @param modalities vector of strings containing the modalities to be preprocessed. It must always contains the T1-weighted sequence scan.
#' @param mask brain mask in NifTI format of the atlas template to performed the skull stripping.
#' @param inhomogeneity inhomogeneity correction algorithm to be applied. The correction by default is the 'N4' bias correction.
#' @param transformation non-linear transformation for registering the T1-w MRI scan to the reference template. 'SyN' transformation is used by default.
#' @return paths of preprocessed MRI scans.
#' @importFrom extrantsr bias_correct ants_regwrite
#' @importFrom neurobase mask_img
#' @importFrom fslr fast
#' @importFrom oro.nifti writeNIfTI
#' @export
preprocess_modalities <- function(mri.patient, folder.patient, modalities, atlas, mask, inhomogeneity = "N4", transformation = "SyN"){

  # empty list to save paths
  mri_paths <- list()

  # Inhomogeneity Correction: N4
  cat(paste0('*********************************************\n****** Inhomogeneity Correction: ', inhomogeneity ,' *********\n*********************************************\n--Running...\n'))
  bias_files <- lapply(modalities, function(x) file.path(folder.patient, paste0( x, '_bias.nii.gz')))
  bias_mri <- mapply(extrantsr::bias_correct, file = mri.patient, correction = inhomogeneity, outfile = bias_files, verbose = FALSE)
  mri_paths[['bias']]<- bias_files
  cat('--Complete.\n')

  # Coregistration to T1-weighted image : Rigid Transformation
  if (length(modalities) > 1){
    cat(paste0('*********************************************\n****** Coregistration to T1 sequence ********\n*********************************************\n--Running...\n'))
    coregisteredImg <- extrantsr::within_visit_registration(fixed = mri.patient$T1, moving = bias_files[2:length(mri.patient)], typeofTransform = "Rigid", interpolator = "Linear", verbose = FALSE)
    bias_mri_comp <- coregistration_images(coregisteredImg)
    cat('--Complete.\n')
  }

  # Registration to Template (SyN: Non-linear)
  cat(paste0('*********************************************\n******* Spatial Registration : ',transformation , ' **********\n*********************************************\n--Running...\n'))
  syn_files <- lapply(modalities, function(x) file.path(folder.patient, paste0( x, '_SyN_MNI152.nii.gz')))
  if (length(modalities) > 1){
    bias_mris <- create_bias_list(modalities, bias_mri$T1, bias_mri_comp)
    syn_mri <- mapply(extrantsr::ants_regwrite, filename = bias_mris, outfile = syn_files, template.file = atlas, typeofTransform = transformation,  verbose = FALSE)
  }else{
    syn_mri <- extrantsr::ants_regwrite(filename = bias_files[[1]], outfile = syn_files[[1]], template.file = atlas, typeofTransform = transformation, verbose = FALSE)
  }
  mri_paths[['registered']] <- syn_files
  cat('--Complete.\n')

  # Brain Mask
  cat('*********************************************\n*************** Brain Mask ******************\n*********************************************\n--Running...\n')
  mask_files <- lapply(modalities, function(x) file.path(folder.patient, paste0( x, '_masked')))
  mask_mri <- lapply(syn_files, neurobase::mask_img, mask)
  mapply( oro.nifti::writeNIfTI, nim = mask_mri, filename = mask_files)
  mri_paths[['stripped']] <- paste0(mask_files,'.nii.gz')
  cat('--Complete.\n')

  # Spatially Informed layer
  cat('*********************************************\n******** Brain Segmentation ***********\n*********************************************\n--Running...\n')
  spatial_file <- file.path(folder.patient, 'T1_spatially_informed.nii.gz')
  spatial_mri = fslr::fast(file = mask_mri[[1]], outfile = spatial_file, opts = "--nobias", verbose = FALSE, type = 'T1')
  mri_paths[['segmented']] <- file.path(folder.patient, 'T1_spatially_informed_pveseg.nii.gz')
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
  ravel_file = file.path(folder.patient, 'T1_Ravel_norm')
  mri_paths[['ravel']] <- paste0(ravel_file ,'.nii.gz')

  return(mri_paths)
}


#' @title Wrapper function for RAVEL normalization of T1-weighted images
#'
#' @description Ravel intensity normalization using control voxels and clinical covariates.
#' @param masked.paths list or vector of paths of the preprocessed input NIfTI images to be normalized.
#' @param csf.paths NIfTI image paths for the binary control region masks.
#' @param ravel.paths list or vector of paths of the output NIfTI images.
#' @param demographics table of covariates associated to the MRI scans. Number of rows should be equal to the number of images.
#' @param brain_mask NIfTI image path for the binary brain mask. Must have value 1 for the brain tissue and 0 otherwise.
#' @param patients.folder folder to save the output control mask.
#' @return RAVEL-corrected images are saved in disk.
#' @importFrom RAVEL maskIntersect normalizeRAVEL
#' @export
### Ravel Normalization
image_normalization_ravel <- function(masked.paths, csf.paths, ravel.paths, demographics, brain.mask, patients.folder, modality = 'T1'){

  ### Control region mask for all patients (masks intersect)
  mask_intersect_path <- file.path(patients.folder, 'CSF_control_mask.nii.gz')
  intersect_mask <- RAVEL::maskIntersect(csf.paths, output.file = mask_intersect_path , prob = 0.8)

  # Control for biological covariates
  mod_cov <- model.matrix(~., demographics[,-1])

  #RAVEL
  RAVEL::normalizeRAVEL(input.files = masked.paths, output.files = ravel.paths, brain.mask = brain.mask,
                 control.mask = intersect_mask, mod = mod_cov , WhiteStripe_Type	= modality,
                 k = 1, returnMatrix = FALSE , writeToDisk = TRUE)
}


#' @title Wrapper function to preprocess T1-weighted, T2-weighted and/or FLAIR MRI scan for multiple patients
#'
#' @description This function preprocesses raw T1-weighted, T2-weighted and/or FLAIR MRI scans and generates a brain segmentation MRI scans using the FAST algorithm.
#' The preprocessing steps comprise imhomogeneity correction 'N4', linear coregistration of T2-weighted and/or FLAIR to the T1-weighted, registration of all available modalities to the MNI152 template with an isotropic voxel size of 1mm^3
#' using the 'SyN' transformation, skull stripping, and RAVEL intensity normalization.
#'
#' @param patients.general general folder containing sub-folders per patient with raw MRI images.
#' @param clinical.covariates data.frame of covariates associated to the MRI scans. Number of rows should be equal to the number of images.
#' @return paths of preprocessed MRI scans. MRI preprocessed images are stored in the patient's folder.
#' @importFrom MNITemplate getMNIPath readMNI
#' @export
preprocess_patients <- function(patients.folder, clinical.covariates){

  if(missing(clinical.covariates)) {
    message("No covariates provided. Intensity normalization step will use the White Stripe algorithm instead of RAVEL.\n")
  } else {
    if(nrow(clinical.covariates) < 1){
      stop('no covariates provided. File is empty.')
    }
  }

  # getting patients scans
  patients_mri <- load_mri_group(patients.folder)

  # empty list to save paths
  paths_mri <- list()

  # Atlas template (MNI152)
  mniTemplate <- MNITemplate::getMNIPath()

  # Atlas brain mask (MNI152)
  brainMask <- MNITemplate::readMNI("Brain_Mask")

  # preprocess each patients T1 sequence scan
  for (patient in patients_mri){

    # selecting patient's folder
    patient_folder <- dirname(patient$T1)

    # getting modalities
    modalities <- get_modalities(patient)

    # getting folder name
    folder_name = unlist(strsplit(patient_folder, '/'))
    folder_name = folder_name[length(folder_name)]

    # Preprocess images
    message(paste0('--------------------------------------------------\n Preprocesing images in ',folder_name ,'\n--------------------------------------------------\n\n'))
    patient_scans <- preprocess_modalities(patient, patient_folder, modalities, mniTemplate, brainMask)

    # Save preprocessed images paths
    paths_mri[[folder_name]] <- patient_scans
  }

  ### RAVEL

  ## Group paths for ravel
  ravel_paths <- lapply(paths_mri, function(x){x$ravel})
  csf_paths <- lapply(paths_mri, function(x){x$csf_mask})
  masked_paths <- lapply(paths_mri, function(x){x$syn_masked})
  masked_paths_T1 <- lapply(masked_paths, function(x) x[grepl("T1", x)])


  if(missing(clinical.covariates)) {
    cat('*********************************************\n******** White Stripe normalization **********\n*********************************************\n--Running...\n')
    RAVEL::normalizeWS(input.files = masked_paths_T1, output.files = ravel_paths, brain.mask = brainMask, WhiteStripe_Type = "T1", returnMatrix = FALSE , writeToDisk = TRUE)
  } else {
    if (length(ravel_paths) < nrow(clinical.covariates)){
      stop("more covariates information than images. Information per MRI scan should be provided.")
    }

    if (length(ravel_paths) > nrow(clinical.covariates)){
      stop("more images than covariates information. Information per MRI scan should be provided.")
    }

    if (length(ravel_paths) < ncol(clinical.covariates)){
      stop("more covariates than images. Number of covariates can not be greater than number of images.")
    }

    ## ravel algorithm
    image_normalization_ravel(masked_paths_T1, csf_paths, ravel_paths, clinical.covariates, brainMask, patients.folder)
  }

  cat(paste0('--------------------------------------------------\n Preprocessing Complete \n--------------------------------------------------\n\n'))
  return(paths_mri)
}


