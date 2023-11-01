correctForSensSpec <- function(positiveResults,negativeResults,sensitivity,specificity) {
#' This computes the real number of positives and negatives given a number of positive and negative test results as well as sensitivity and specificity values
#' 
#' Derivation see picture of hand written notes
#' 
#' NOTE: this corresponds to the Rogan-Gladen estimator to adjust prevalence data for test performance!!!
#' 
#' TP True Positives
#' TN True Negatives
#' FP False Positives
#' FN False Negatives
  
  result = list()
  
  result$totalResults <- positiveResults + negativeResults
  
  result$TP <- (positiveResults + negativeResults - negativeResults/specificity) / (1/specificity + 1/sensitivity - 1/(sensitivity*specificity))
  result$TN <- (positiveResults + negativeResults - positiveResults/sensitivity) / (1/specificity + 1/sensitivity - 1/(sensitivity*specificity))
  result$FN <- negativeResults - result$TN
  result$FP <- positiveResults - result$TP
  
  result$realPositives <- result$TP + result$FN
  result$realNegatives <- result$TN + result$FP
  
  
  return(result)
  
}



#' rogan gladen estimator
#' inputs: apparent prevalence (positive tests / sample size), sensitivity, specificity (as proportions)
#' value: corrected prevalence
rogangladen <- function(aparent_prevalence, sensitivity, specificity) {
  return((aparent_prevalence + specificity - 1) / (sensitivity + specificity - 1))
}

dnorm_rescale <- function(x,scale,mean,sd,log = FALSE) {
  return(scale * dnorm(x,mean,sd,log))
}

dnorm_mixture <- function(x,scale1,mean1,sd1,mean2,sd2,log = FALSE) {
  return(scale1 * dnorm(x,mean1,sd1,log) + (1-scale1) * dnorm(x,mean2,sd2,log))
}
