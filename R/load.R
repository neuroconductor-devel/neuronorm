#' Load MRI scans per patient
#'
#' This function loads the MRI scans from one patient. It assumes that the MRI
#' scans are contained in the same folder and refer to MRI modalities T1, T2, Flair.
#' Any scan with duplicated modality keyword will be ommited with the first one being
#' kept.
#'
#' @param folder folder containing the MRI scans. The MRI scans should be in format NiFTI.
#' @param modalities the MRI modalities to be considered. Should be at least one of T1, T2 or FLAIR. By default, all modalities are searched within the folder.
#' @return paths of MRI scans for a patient if they exist.
#' @export
load_mri_patient <- function(folder, modalities = c('T1','T2','FLAIR')){
  folder = file.path(folder)
  mri_images <- list()
  for (modality in modalities){
    data <- list.files(folder, pattern = modality, full.names = TRUE)
    if (length(data) > 0){
      mri_images[[modality]] <- data[1]
    }else{
      cat(paste0('-- No MRI images found for ', modality, ' modality\n'))
    }
  }
  return(mri_images)
}


#' Load MRI per group or disease
#'
#' This function loads the MRI scans from multiple patients. It assumes that the patients
#' folder with the the MRI scans are contained in a general folder.
#'
#' @param folder folder containing the MRI scans.
#' @return paths of MRI scans per patient if they exist.
#' @export
load_mri_group <- function(folder){
  folders <- list.dirs(folder, full.names = TRUE, recursive = TRUE)
  mri_images <- list()
  cat('--------------------------------------------------\n')
  for (fold in folders){
    fold_name = unlist(strsplit(fold, '/'))
    fold_name = fold_name[length(fold_name)]
    message (paste0('Reading folder ', fold_name , '\n'))
    mri_patient <- load_mri_patient(fold)
    if (length(mri_patient) > 0){
      mri_images[[fold_name]] <- mri_patient
    }
  }
  if(length(mri_images) == 0){
    stop("T1-w scans must be provided.")
  }
  cat('--------------------------------------------------\n')
  return(mri_images)
}
