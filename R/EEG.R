#' @title  ElectroEncephaloGraphy Data
#' 
#' @description This data arises from a large study to examine EEG
#' correlates of genetic predisposition to alcoholism.
#' See <http://kdd.ics.uci.edu/databases/eeg/> for details.
#' 
#' @format `list(alcholic, control)`.
#' 
#' @details There were two groups of subjects: 77 alcoholic and 45
#' control. In the original data, each subject was exposed to either 
#' a single stimulus (S1) or to two stimuli (S1 and S2) which were 
#' pictures of objects chosen from the 1980 Snodgrass and Vanderwart 
#' picture set. 
#' Under each condition, each subject was observed repeatly for 120 trials.
#' 
#' Here, this dataset include the averages of 120 trials under S1. 
#' It is a list including two arrays,
#' \itemize{
#'      \item `alcholic`: Array of `256 x 64 x 77`.
#'      \item `control`:  Array of `256 x 64 x 45`.
#' }
#' 
"EEG"
#> [1] "EEG"
