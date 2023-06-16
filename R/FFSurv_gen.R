#' Generate Survival Data for Flexible Feed-Forward Estimation
#'
#' @description This function generates survival data according to specified parameters.
#' The generated data can be used for flexible feed-forward estimation in discrete-time panel
#' survival models, including models with right-censoring, unobserved heterogeneity,
#' and grouped responses.
#'
#' @param n_person Integer. The number of individuals to generate data for.
#' @param n_incident Integer or vector. The number of incidents per person. If a single
#' integer is given, the same number of incidents is generated for each person.
#' If a vector is given, it should have length equal to n_person, specifying the number of
#' incidents for each person individually.
#' @param n_features Integer. The number of covariate features to include. Default is 3.
#' @param beta Numeric vector. The coefficients for the covariate features. Default is c(-1, .4, -1, 1).
#' @param het Vector. Specifies the distribution of unobserved heterogeneity.
#' The first element is a character string indicating the distribution ('gamma', 'lognorm', 'dual', 'uniform'),
#' and the remaining elements are parameters for the distribution.
#'   * 'gamma': Gamma distribution, with shape and scale parameters. Example: c('gamma', 1, 1)
#'   * 'lognorm': Log-normal distribution, with mean and standard deviation parameters. Example: c('lognorm', 1, 1)
#'   * 'dual': Binary distribution, with the probability of success and two numeric values
#'             representing the two possible outcomes. Example: c('dual', 0.5, 1, 5)
#'   * 'uniform': Uniform distribution, with minimum and maximum parameters. Example: c('uniform', 1, 5)
#' Default is c('gamma',1,1).
#' @param X_features List. Specifies the distributions of the covariate features.
#' The list has three elements:
#'   * 'dist': A character vector indicating the distribution of each feature ('normal' or 'dummy').
#'   * 'mean': A numeric vector indicating the mean of each feature.
#'   * 'sd': A numeric vector indicating the standard deviation of each feature. This can be NA for 'dummy' features.
#' Default is list( 'dist' = c('normal','normal','dummy'), 'mean' = c(5,5,0.5), 'sd' = c(.2,.2,NA) ).
#' @param baseline_hazard List. The function defining the baseline hazard and its maximum value.
#' The list has two elements:
#'   * 'formula': A formula according to which the baseline hazard function is constructed. It has to be in the form of "y = f(x)"
#'                with f(.) being the desired function.
#'   * 'max': The maximum value for the baseline hazard.
#' Default is list('formula' = 'y = log(1+x/125)', 'max' = 100).
#' @param psi Numeric. The grouping interval. Default is 1.
#' @param censoring List. Specifies the censoring mechanism. The first element 'mode' can be 'fixed', 'random', or 'None'.
#'   * 'fixed': Fixed censoring time. Specified by 'val' (fixed censoring time at a certain value) or 'ratio' (quantile for censoring time).
#'   * 'random': Random censoring times. Specified by 'val' (mean censoring time) and 'sd' (standard deviation).
#'   * 'None': No censoring.
#' Default is list( 'mode' = 'fixed', 'val' = 80, 'ratio'=NA, 'sd' = NA).
#'
#' @return A data frame with columns for 'id', 'response', 'censor', 'intercept',
#' and the covariate features.
#'
#' @examples
#' \dontrun{
#'   data <- FFSurv_gen(n_person = 100, n_incident = 1, n_features = 3,
#'                      beta = c(-1, .4, -1, 1), het = c('gamma',1,1),
#'                      X_features = list( 'dist' = c('normal','normal','dummy'),
#'                                         'mean' = c(5,5,0.5), 'sd' = c(.2,.2,NA) ),
#'                      baseline_hazard = list('formula' = 'y = log(1+x/125)', 'max' = 100),
#'                      psi = 1,
#'                      censoring = list( 'mode' = 'fixed', 'val' = 80, 'ratio'=NA, 'sd' = NA))
#' }
#'
#' @references Bateni[2023]: Flexible Feed-Forward Estimation in Discrete-time Panel Survival Models.
#'
#' @export
FFSurv_gen = function(n_person, n_incident, n_features = 3, beta = c(-1, .4, -1, 1), het = c('gamma',1,1),
                      X_features = list( 'dist' = c('normal','normal','dummy'), 'mean' = c(5,5,0.5), 'sd' = c(.2,.2,NA) ),
                      baseline_hazard = list('formula' = 'y = log(1+x/125)', 'max' = 100),
                      psi = 1,
                      censoring = list( 'mode' = 'fixed', 'val' = 80, 'ratio'=NA, 'sd' = NA)){

  if(length(n_incident)==1){
    X = matrix(NA, nrow = n_person * n_incident, ncol = n_features + 4)
    X[,1] = rep(1:n_person, each = n_incident)
  } else{
    if(length(n_incident) != n_person){
      stop('the size of the number of incidents vector is not the same as the number of persons.')
    }
    X = matrix(NA, nrow = sum(n_incident), ncol = n_features + 4)
    X[,1] = unlist(mapply(rep, 1:n_person, each = n_incident))
  }

  colnames(X) = c('id','response','censor','intercept',paste('x',1:n_features,sep=''))
  X = as.data.frame(X)
  X$intercept = 1

  if(het[1] == 'gamma'){
    heterogeneity = rgamma(n_person,as.numeric(het[2]),as.numeric(het[3]))
  }
  if(het[1] == 'lognorm'){
    heterogeneity = rlnorm(n_person,as.numeric(het[2]),as.numeric(het[3]))
  }
  if(het[1] == 'dual'){
    heterogeneity = rbinom(n_person,size = 1, prob = as.numeric(het[2]))
    heterogeneity = ifelse(heterogeneity == 1, as.numeric(het[3]), as.numeric(het[4]))
  }
  if(het[1] == 'uniform'){
    heterogeneity = runif(n_person, as.numeric(het[2]) , as.numeric(het[3]) )
  }
  if(!(het[1] %in% c('gamma','lognorm','dual','uniform'))){
    stop('the heterogeneity is not from any of the specified forms.')
  }

  heterogeneity = as.vector(unlist(mapply(rep, heterogeneity, each = n_incident)))

  for(i in 1:n_features){
    if(X_features$dist[i] == 'normal'){
      features = rnorm(nrow(X), mean = X_features$mean[i], sd = X_features$sd[i])
    }
    if(X_features$dist[i] == 'dummy'){
      features = rbinom(nrow(X), 1, prob = X_features$mean[i])
    }
    if(!(X_features$dist[i] %in% c('normal','dummy'))){
      stop('some feature distribution is not from the specified forms')
    }
    X[,4+i] = features
  }

  if(length(beta) != (n_features +1) ){
    stop('the size of the coefficients does not match the number of features')
  }

  exp_beta = exp( as.matrix(X[,4:ncol(X)]) %*% beta )
  het_exp_beta = heterogeneity * exp_beta

  baseline_eval = function(formula, X) {
    y <- numeric(length(X))  # Pre-allocate a numeric vector for y
    for(i in seq_along(X)) {
      x <- X[i]  # Current x value
      y[i] <- eval(parse(text = gsub("y =", "", formula)))
    }
    return(y)
  }
  time_seq = seq(0, baseline_hazard$max, psi/10)
  hazard = baseline_eval(baseline_hazard$formula, time_seq)
  hazard = c(hazard, 1e12)
  time_seq = c(time_seq, baseline_hazard$max + psi/10)
  if(hazard[1]!=0){hazard = c(0,hazard)}
  cum_hazard = cumsum(hazard)
  for(i in 1:nrow(exp_beta)){
    CDF_vec = 1 - exp(- het_exp_beta[i,] * cum_hazard)
    PDF_vec = c(0,diff(CDF_vec))
    X$response[i] = sample(time_seq, size=1, prob = PDF_vec)+1
  }
  X$response = psi * floor(X$response/psi)
  if(censoring$mode == 'fixed'){
    if(!(is.na(censoring$val))){
      censor = rep(censoring$val, nrow(X))
    }else{
      cens_val = quantile(X$response, censoring$ratio)
      censor = rep(cens_val, nrow(X))
    }
    censor = psi * floor(censor / psi)
  }
  if(censoring$mode == 'random'){
    censor = rnorm(nrow(X), mean = censoring$val, sd = censoring$sd)
    censor = psi * floor(censor / psi)
  }
  X$censor = 0
  if(censoring$mode != 'None'){
    X$censor   = ifelse(X$response > censor,1,0)
    X$response = ifelse(X$censor == 0, X$response, censor)
  }
  return(X)
}
