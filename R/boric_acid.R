#' Developmental toxicity of boric acid in mice
#'
#' A dataset containing the outcomes of a developmental toxicity experiment of boric acid by the National Toxicology
#' Program. Pregnant mice were provided feed containing boric acid. At the end of the study, each mouse's uterus was 
#' examined to evaluate the fetal outcomes.
#'
#' @format A data frame with 67 rows and 4 variables:
#' \describe{
#'   \item{Trt}{Exposure level to boric acid, %}
#'   \item{ClusterSize}{Number of fetuses in the litter}
#'   \item{NResp}{Number of malformed, dead, or resorped fetuses}
#'   \item{Freq}{Number of litters with the same exposure/response pattern}
#' }
#' @source Heindel, J J and Price, C J and Schwetz, B A. The developmental toxicity of boric acid in mice, rats, and rabbits.
#'	Environmental Health Perspectives, 102 Suppl 7,	1994, 107--12,	.
"boric_acid"