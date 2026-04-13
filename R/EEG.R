#' @title  ElectroEncephaloGraphy Data
#' 
#' @description This data arises from a large study to examine EEG
#' correlates of genetic predisposition to alcoholism.
#' See details.
#' 
#' @format `list(alcholic, control)`.
#' 
#' @details
#' The original data are fully open-access and available in UCI Machine
#' Learning Repository (<http://kdd.ics.uci.edu/databases/eeg/>).
#' It includes two groups of subjects: 77 alcoholic and 45 control.
#' In the original study, each subject was exposed to either a single
#' stimulus (S1) or two stimuli (S1 and S2), which were pictures chosen from
#' the 1980 Snodgrass and Vanderwart picture set.
#' Each subject underwent 120 trials under each condition.
#' 
#' Here we provide this preprocessed data of the averages of 120 trials under
#' S1 condition, which has been studied in several literature.
#' The dataset is structured as a list containing two arrays,
#' \itemize{
#'  \item `EEG$alcholic`: Array of dimensions `256 x 64 x 77`.
#'  \item `EEG$control`:  Array of dimensions `256 x 64 x 45`.
#' }
#' 
#' @references
#' If you use this dataset, we would be very grateful if you could cite both
#' the original data source and our work:
#' 
#' Zengchao Xu, Shan Luo, and Binyan Jiang. "Doubly Regularized Matrix-Variate Regression". *Submitted*.
#' 
#' @examples
#' data(EEG, package="drrglm")
#' 
"EEG"
#> [1] "EEG"
