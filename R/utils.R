#' Create a vector of the coregister images based on the images modalities
#'
#' This function creates a vector of corregistered images for a patient.
#' The vector contains whether a vector of length one for only one modality (T2 or FLAIR)
#' or a vector of length two including both modalities (T2 and FLAIR).
#'
#' @param vector output object from the coregistration function.
#' @export
coregistration_images <- function(vector){
  if (length(vector)== 2){
    imgs <- lapply(vector, function(x) x$outfile)
  }else{
    imgs <- vector$outfile
  }
  return(imgs)
}


#' Create a vector of bias images based on the images modalities
#'
#' This function creates a vector of bias images for a patient.
#' The vector contains whether a vector of length one for only one modality (T1)
#' or a vector of length two/three including all modalities.
#'
#' @param vector output object from the coregistration function.
#' @export
create_bias_list <- function(modalities, bias_T1, list_corregister){
  bias_mris <- list()
  bias_mris$T1 <- bias_T1
  for (n in 1:length(list_corregister)){
    modality <- modalities[[n+1]]
    bias_mris[[modality]] <- list_corregister[[1]]
  }
  return(bias_mris)
}


#' Create a vector with the image modalities in each patient folder
#'
#' This function creates a vector with the name of the image modalitities fpr a patient.
#' @param patient paths of MRI scan per patient.
#' @export
get_modalities <- function(patient){
  names <- names(patient)
  if ('T1' %in% names){
    modalities <- names
  }else{
    stop('Preprocessing can not be performed without a T1-weighted scan. Please make sure your folder contains a T1-w image.')
  }
  return(modalities)
}
