#' Example Data for r2tfDiffusion Package
#'
#' A list containing example network gene expression data and transcription factor activity andclinical outcomes
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{network}{A dataframe representing the TLR4 signaling network (edges).}
#'   \item{feedbackNodes}{A character vector of feedback or negative regulator genes.} 
#'   \item{tfAct}{Transcription factor activity changes ex vivo.}
#'   \item{Mexvivo}{Ex vivo gene expression before stimulation.}
#'   \item{Mclinical}{Gene expression in clinical samples before treatment.}
#'   \item{patientData}{Data frame of clinical patient outcomes}
#' }
#'
#' @usage data(r2tfData)
#' @examples
#' data(r2tfData)
#' names(r2tfData)
#' 
#' @keywords datasets
"r2tfData"
