#' Flexible Feed-Forward Estimation in Discrete-time Panel Survival Models
#'
#' @description This function implements the maximum likelihood method for discrete time panel
#' survival models that is proposed in the paper: "Bateni[2023]: Flexible Feed-Forward Estimation
#' in Discrete-time Panel Survival Models". The model allows for right-censoring, unobserved
#' heterogeneity and grouped responses. It estimates the effect of covariates in a panel setup
#' where individuals may undergo multiple duration spells.
#'
#' @param X A data frame that must contain the following columns:
#'  * 'id': a numeric or character column for unique identification of individuals.
#'  * 'response': a numeric column representing the response variable.
#'  * 'censor': a numeric binary column indicating whether the observation was right-censored or not (1 = censored, 0 = not censored).
#'  Other columns can be included and will be treated as additional covariates.
#'
#' @param intercept Logical. If TRUE, an intercept term is included in the model. Default is TRUE.
#'
#' @param psi Numeric. Parameter for the grouping interval. Default is 1.
#'
#' @param init_params Numeric vector. Initial values for the parameter estimates.
#' Not generally recommended to set manually as some initial values violate the
#' assumptions of the model and can result in bogus results. Default is NA.
#'
#' @param n_iter Numeric. The maximum number of iterations for the likelihood optimization. Default is 1e4.
#'
#' @param bootstrap Logical. If TRUE, bootstrap resampling is performed for statistical inference. Default is FALSE.
#'
#' @param bootstrap_n_iter Numeric. The maximum number of iterations for the bootstrap process. Default is 1e3.
#'
#' @param n_bootstrap Numeric. The number of bootstrap resamples to take. Default is 1e3.
#'
#' @return A list containing:
#' * Coefficients: a data frame containing the estimated coefficients and their associated standard errors, confidence intervals, and p-values.
#' * Frailty: a data frame showing the estimated shape and rate parameters for the frailty distribution.
#' * Baseline Hazard Plot: a plot object showing the estimated baseline hazard function.
#'
#' @examples
#' \dontrun{
#'   results <- FFSurv_est(X = data, intercept = TRUE, psi = 1, n_iter = 1e4)
#' }
#'
#' @references Bateni[2023]: Flexible Feed-Forward Estimation in Discrete-time Panel Survival Models.
#'
#' @export
FFSurv_est = function(X, intercept = T , psi = 1 , init_params = NA , n_iter = 1e4 , bootstrap = F, bootstrap_n_iter = 1e3, n_bootstrap=1e3){

  if(intercept){
    int_0 = matrix(rep(1, nrow(X)),ncol=1)
    colnames(int_0) = 'intercept'
    X = cbind(int_0,X)}
  X = as.data.frame(X)
  Array_func = function(X,intercept=F,psi=1){

    data_check = function(X,psi=1){

      Y = X$response
      Y_psi = trunc(Y/psi*10)/10
      if(!(all(is.na(Y_psi) | Y_psi == round(Y_psi)))){stop("Error: Some of the duration outcomes are not multiples of the grouping interval length.")}
      if(!(all(is.na(Y) | Y >= 0))){stop("Error: Some of the duration outcomes are negative.")}
      if(!all(sapply(Y, is.numeric))) {stop("Error: All duration outcomes must be numeric.")}

      C = X$censor
      if(!all(sapply(C, function(y) is.numeric(y) & y >= 0 | is.logical(y)))) {stop("Error: All censoring indicators must be non-negative numerics or boolean")}

      Z = X[, -match(c("id", "response","censor"), names(X))]
      if(!all(sapply(Z, is.numeric))) {stop("Error: All predictors must be numeric.")}
    }

    data_check(X,psi)

    preds = X[, -match(c("id", "response","censor"), names(X))]
    max_obs = max(table(X$id))
    num_predictors = ncol(preds)
    num_individuals = length(unique(X$id))

    Z = array(NA, dim = c(max_obs, num_predictors, num_individuals))
    Y = matrix(NA, nrow = num_individuals, ncol = max_obs)
    C = matrix(NA, nrow = num_individuals, ncol = max_obs)

    bool_v <- function(vec) {ifelse(vec > 0 | vec == TRUE, 1, ifelse(vec == 0 | vec == FALSE, 0, NA))}

    for (i in unique(X$id)) {
      preds_i = preds[X$id == i, ]
      X_i = X[X$id == i, ]
      j = which( unique(X$id) == i)

      Z[1:nrow(preds_i), , j ] = as.matrix(preds_i)
      Y[j, 1:nrow(preds_i)] = X_i$response
      C[j, 1:nrow(preds_i)] = bool_v(X_i$censor)
    }

    Y = Y/psi
    Y = trunc(Y*10)/10

    #identifiability assumption: assuming the highest observation is censored
    C[ which ( Y == max ( Y , na.rm = T ) ) ] = 1

    return(list(pred = Z, resp = Y, cens = C))
  }

  param_func = function(data){

    beta_0 = rep( 0 , ncol( data$pred ) )

    observed_durations = c( na.omit ( unique ( c (data$resp) ) ) )
    delta_0 = rep( -1e4 , max( observed_durations ) )
    #identifiability assumption: assuming every increment in which no outcome is observed is zero.
    delta_0[observed_durations] = 0

    delta_0 = c(-1e4,delta_0)
    #adding a zero first increment because the model specifies that each delta_k corresponds to the integral increment for [(k-1)*psi , k*psi). Hence, if we observe an outcome at y = k, it would correspond to delta_(k+1)

    nu_0 = kappa_0 = 0

    return( c( beta_0 , delta_0, nu_0 , kappa_0 ) )

  }

  Flex_lkl_1 = function(data,params){

    pred = data$pred
    resp = data$resp
    cens = data$cens

    length_beta = ncol(pred)

    beta = params[1:length_beta]
    exp_delta = exp( params[(length_beta + 1) : (length(params) - 2)] )
    nu_0 = exp( params[ length(params) - 1 ] )
    kappa_0 = exp( params[ length(params) ] )

    n.person = nrow(resp)
    n.incid  = ncol(resp)

    cum_baseline_hazard = cumsum(exp_delta)

    exp_pred_beta = exp(aperm(sapply(1:dim(pred)[3], function(i) pred[,,i] %*% beta), c(2,1)))

    cum_hazard = array( cum_baseline_hazard[ resp ] , dim= c(n.person , n.incid) )

    cum_hazard_ahead = array( cum_baseline_hazard[ resp+1 ] , dim= c(n.person , n.incid) )

    nu    = matrix( rep( nu_0    , n.person*n.incid ) , nrow = n.person )
    kappa = matrix( rep( kappa_0 , n.person*n.incid ) , nrow = n.person )

    for(i in 2:n.incid){
      nu[,i] = nu[,i-1] + (1-cens[,i-1])

      kappa[,i] = kappa[,(i-1)] +
        exp_pred_beta[ , (i-1)] *
        ( cum_hazard [ , (i-1)] + .5 * exp_delta[ ( resp[ , (i-1)] + 1 ) ] * ( 1 - cens[ , (i-1)] ) )
    }

    kappa_exp_pred_beta = kappa^(-1) * exp_pred_beta
    xi      = 1 + kappa_exp_pred_beta * cum_hazard
    xi_plus = 1 + kappa_exp_pred_beta * cum_hazard_ahead

    ind_likelihood =  (( xi )^(-nu) - ( xi_plus )^(-nu)*(1-cens))
    approx = (ind_likelihood<=0)
    approx = ifelse(is.na(approx),FALSE,approx)

    if(sum(approx)){

      xi_approx = exp_pred_beta * cum_hazard
      xi_plus_approx = exp_pred_beta * cum_hazard_ahead * (1-cens)

      #Laurent series
      ind_likelihood_approx = (nu * ( xi_plus_approx - xi_approx ) / kappa) + ( nu * (nu+1) * (xi_approx^2 - xi_plus_approx^2) ) / ( 2 * kappa^2 ) - ( nu * (nu+1) * (nu+2) * (xi_approx^3 - xi_plus_approx^3) ) / ( 6 * kappa^3 )

      ind_likelihood[approx] = ind_likelihood_approx[approx]

    }

    ind_log_likelihood = log(ind_likelihood)

    joint_log_likelihood  =  sum(  ind_log_likelihood ,na.rm = T )
    joint_log_likelihood = ifelse ( joint_log_likelihood==0, sum(ind_log_likelihood), joint_log_likelihood)

    return(-joint_log_likelihood)
  }

  Flex_derivative_1 = function(data,params,likelihood=F,raw = F){

    pred = data$pred
    resp = data$resp
    cens = data$cens

    length_beta = ncol(pred)

    beta = params[1:length_beta]
    exp_delta = exp( params[(length_beta + 1) : (length(params) - 2)] )
    nu_0 = exp( params[ length(params) - 1 ] )
    kappa_0 = exp( params[ length(params) ] )

    n.person = nrow(resp)
    n.incid  = ncol(resp)

    cum_baseline_hazard = cumsum(exp_delta)

    exp_pred_beta = exp(aperm(sapply(1:dim(pred)[3], function(i) pred[,,i] %*% beta), c(2,1)))

    cum_hazard = array( cum_baseline_hazard[ resp ] , dim= c(n.person , n.incid) )

    cum_hazard_ahead = array( cum_baseline_hazard[ resp+1 ] , dim= c(n.person , n.incid) )

    nu    = matrix( rep( nu_0    , n.person*n.incid ) , nrow = n.person )
    kappa = matrix( rep( kappa_0 , n.person*n.incid ) , nrow = n.person )

    for(i in 2:n.incid){
      nu[,i] = nu[,i-1] + (1-cens[,i-1])

      kappa[,i] = kappa[,(i-1)] +
        exp_pred_beta[ , (i-1)] *
        ( cum_hazard [ , (i-1)] + .5 * exp_delta[ ( resp[ , (i-1)] + 1 ) ] * ( 1 - cens[ , (i-1)] ) )
    }

    kappa_exp_pred_beta = kappa^(-1) * exp_pred_beta
    xi      = 1 + kappa_exp_pred_beta * cum_hazard
    xi_plus = 1 + kappa_exp_pred_beta * cum_hazard_ahead

    ind_likelihood =  (( xi )^(-nu) - ( xi_plus )^(-nu)*(1-cens))
    approx = (ind_likelihood<=0)
    approx = ifelse(is.na(approx),FALSE,approx)

    if(sum(approx)){

      xi_approx = exp_pred_beta * cum_hazard
      xi_plus_approx = exp_pred_beta * cum_hazard_ahead * (1-cens)

      #Laurent series
      ind_likelihood_approx = (nu * ( xi_plus_approx - xi_approx ) / kappa) + ( nu * (nu+1) * (xi_approx^2 - xi_plus_approx^2) ) / ( 2 * kappa^2 ) - ( nu * (nu+1) * (nu+2) * (xi_approx^3 - xi_plus_approx^3) ) / ( 6 * kappa^3 )

      ind_likelihood[approx] = ind_likelihood_approx[approx]

    }

    ind_log_likelihood = log(ind_likelihood)


    if(likelihood){
      joint_log_likelihood  =  sum(  ind_log_likelihood ,na.rm = T )
      joint_log_likelihood = ifelse ( joint_log_likelihood==0, sum(ind_log_likelihood), joint_log_likelihood)
    }

    pred_new = aperm( pred , c(3,1,2) )
    dimensions = dim( pred_new )

    ##d/d.beta

    d_xi_d_beta      =  pred_new * array( exp_pred_beta * cum_hazard       , dim = dimensions )
    d_xi_plus_d_beta =  pred_new * array( exp_pred_beta * cum_hazard_ahead , dim = dimensions )

    d_kappa_d_beta = ( exp_pred_beta * ( cum_hazard * (1 + cens) + ( cum_hazard_ahead * (1 - cens) ) / 2 ) )
    d_kappa_d_beta = pred_new * array( d_kappa_d_beta , dim = dimensions )
    d_kappa_d_beta = aperm(apply(d_kappa_d_beta, c(1, 3), cumsum),c(2,1,3))

    d_kappa_d_beta [ , 2:ncol(d_kappa_d_beta) , ] = d_kappa_d_beta [ , 1:(ncol(d_kappa_d_beta)-1) , ]
    d_kappa_d_beta[ , 1 , ] = matrix( 0 , nrow=n.person , ncol = dim(pred)[2] )

    ##d/d.delta

    dimensions_delta = c ( n.person , n.incid , length(exp_delta) )

    d_xi_d_delta = array( 0, dimensions_delta )
    d_xi_plus_d_delta = array( 0, dimensions_delta )

    for(i in 1:length(exp_delta)){
      d_xi_d_delta[,,i] = exp_delta[i]*(resp>=i)
      d_xi_plus_d_delta[,,i] = exp_delta[i]*((resp+1)>=i)
    }

    cens_expanded =   array(rep( cens , dimensions_delta[3] ), dimensions_delta )
    d_kappa_d_delta = array((d_xi_d_delta * (1 + cens_expanded) + d_xi_plus_d_delta * (1 - cens_expanded)) / 2,dim = dimensions_delta)

    exp_pred_beta_expanded = array( rep ( exp_pred_beta , dimensions_delta[3] ), dimensions_delta )
    d_xi_d_delta =           array( d_xi_d_delta      * exp_pred_beta_expanded , dim = dimensions_delta)
    d_xi_plus_d_delta =      array( d_xi_plus_d_delta * exp_pred_beta_expanded , dim = dimensions_delta)
    d_kappa_d_delta =        array( d_kappa_d_delta   * exp_pred_beta_expanded , dim = dimensions_delta)

    d_kappa_d_delta = aperm(apply(d_kappa_d_delta, c(1, 3), cumsum),c(2,1,3))

    d_kappa_d_delta[,2:ncol(d_kappa_d_delta),] = d_kappa_d_delta[,1:(ncol(d_kappa_d_delta)-1),]
    d_kappa_d_delta[,1,] = matrix(0,nrow=dimensions_delta[1], ncol = dimensions_delta[3])

    ##derivative of the log-likelihood
    fg.der = function(g,f,g.p,f.p,nu,configs){       ##Gives the derivative of the matrix (1+(f(x)^-1)g(x))^(-nu) over the vector x. The output is a 3D derivative array.

      f_inv   = f^(-1)

      g_f_der =     array( rep( -g/(f^2), dim(f.p)[3] ), dim = dim(f.p))
      f_inv_ext   = array( rep( f_inv   , dim(g.p)[3] ), dim = dim(g.p))

      term1   = array( f.p * g_f_der + g.p * f_inv_ext , dim = configs)

      term2   = (-nu) * ( 1 + (f_inv) * g ) ^ (-nu-1)

      term2   = array( rep( term2   , dim(term1)[3]), dim=dim(term1))

      return( list(  'org' = array( term1 * term2 , configs ), 'ext' = term1  )  )
    }


    ##beta:

    term1_d_beta = fg.der(exp_pred_beta * cum_hazard       , kappa , d_xi_d_beta      , d_kappa_d_beta , nu , dimensions)
    term2_d_beta = fg.der(exp_pred_beta * cum_hazard_ahead , kappa , d_xi_plus_d_beta , d_kappa_d_beta , nu , dimensions)

    d_term1_d_beta = term1_d_beta$org
    d_term2_d_beta = array( term2_d_beta$org * (1 - array(rep(cens,dimensions[3]),dimensions)) , dimensions )

    ##delta:
    term1_d_delta = fg.der(exp_pred_beta * cum_hazard       , kappa , d_xi_d_delta      , d_kappa_d_delta , nu , dimensions_delta)
    term2_d_delta = fg.der(exp_pred_beta * cum_hazard_ahead , kappa , d_xi_plus_d_delta , d_kappa_d_delta , nu , dimensions_delta)

    d_term1_d_delta = term1_d_delta$org
    d_term2_d_delta = array( term2_d_delta$org * (1 - array(rep(cens,dimensions_delta[3]),dimensions_delta)) , dimensions_delta )

    ## kappa.0:

    U.kappa.1 = xi^(-nu-1)      * ( ( exp_pred_beta * cum_hazard       * nu ) / (kappa^2) )

    U.kappa.2 = xi_plus^(-nu-1) * ( ( exp_pred_beta * cum_hazard_ahead * nu ) / (kappa^2) ) * (1-cens)

    ## nu.0:

    U.nu.1 = - ( xi      ^ (-nu) ) * log(xi)
    U.nu.2 = - ( xi_plus ^ (-nu) ) * log(xi_plus) * (1-cens)

    ##

    U.beta  = ( d_term1_d_beta  - d_term2_d_beta  ) / array(rep( ind_likelihood , dimensions[3]       ) , dimensions       )
    U.delta = ( d_term1_d_delta - d_term2_d_delta ) / array(rep( ind_likelihood , dimensions_delta[3] ) , dimensions_delta )
    U.kappa = ( (U.kappa.1  -  U.kappa.2) / ind_likelihood ) *  kappa_0
    U.nu    = ( (U.nu.1     -  U.nu.2)    / ind_likelihood ) *  nu_0

    U.beta =  apply(U.beta  , 3 , function(x) sum(x,na.rm = T) )
    U.delta = apply(U.delta , 3 , function(x) sum(x,na.rm = T) )
    U.kappa = sum(U.kappa , na.rm = T)
    U.nu =    sum(U.nu    , na.rm = T)

    if(likelihood){
      return(list("ell"     = - joint_log_likelihood ,
                  "U_beta"  = - U.beta               ,
                  "U_delta" = - U.delta              ,
                  "U_kappa" = - U.kappa              ,
                  "U_nu"    = - U.nu                 ))
    }

    else if(raw){
      return(-c(U.beta, U.delta, U.kappa, U.nu))
    }

    else{
      return(list("U_beta"  = - U.beta               ,
                  "U_delta" = - U.delta              ,
                  "U_kappa" = - U.kappa              ,
                  "U_nu"    = - U.nu                 ))
    }
  }

  Surv_opt_1 = function( data , init_params = NA , opt_params = 'default' , TT = 1000, raw = F){

    if(is.na(init_params)){
      params = param_func(data)
    }
    else{
      params = init_params
    }

    if(opt_params == 'default'){
      b1 = 0.9
      b2 = 0.999
      t_0 = 1
      eps = 1e-8
      alpha = 1
      decay_1 = 0.75
      decay_2 = 0.95
      decay_steps = 100
    }
    else{
      b1 = opt_params$b1
      b2 = opt_params$b2
      t_0 = opt_params$t_0
      eps = opt_params$eps
      alpha = opt_params$alpha
      decay_1 = opt_params$decay_1
      decay_2 = opt_params$decay_2
      decay_steps = opt_params$decay_steps
    }

    L = Flex_derivative_1( data,params , likelihood = T )
    ell = L$ell
    Utot = c(L$U_beta,L$U_delta,L$U_kappa,L$U_nu)
    m.0 = (1-b1)* (Utot)
    v.0 = (1-b2)* (Utot^2)
    mhat = m.0/(1-(b1^t_0))
    vhat = v.0/(1-(b2^t_0))
    dir = -mhat / (eps + sqrt(vhat))
    alpha_0 = alpha

    num = 0
    for(tt in 1:TT){

      params_new = params + dir * alpha

      ell_new = Flex_lkl_1(data , params_new)

      if( ell_new <= ( ell - abs ( 1e-4 * alpha * ( dir %*% Utot ) ) ) ) {  #1st Wolfe condition

        params = params_new
        L = Flex_derivative_1( data , params , likelihood = T )
        ell = L$ell
        Utot_new  = c(L$U_beta,L$U_delta,L$U_kappa,L$U_nu)
        num = num + 1
        if(num %% decay_steps == 0){
          alpha_0 = alpha_0 * decay_2
        }
        alpha = alpha_0
        m.0 = b1*m.0 + (1-b1)*Utot_new
        v.0 = b2*v.0 + (1-b2)*(Utot_new^2)
        t_0 = t_0+1
        mhat = m.0/(1-(b1^t_0))
        vhat = v.0/(1-(b2^t_0))
        dir = -mhat / (eps + sqrt(vhat))
        Utot = Utot_new
      }

      else{
        alpha = alpha*decay_1
      }
    }

    if(raw){
      return( params )
    }
    return( list( 'beta'  =      params[ 1: ncol(data$pred) ]                             ,
                  'delta' = exp( params[ (ncol(data$pred) + 1) : (length(params) - 2) ] ) ,
                  'nu'    = exp( params[ length(params) - 1 ]                           ) ,
                  'kappa' = exp( params[ length(params) ]                               ) ))
  }

  Surv_optim_1 = function( data , init_params = NA , n.iter = 1e4, grad = F, raw = T){

    if(is.na(init_params[1])){
      params = param_func(data)
    }
    else{
      params = init_params
    }

    surv_fn = function(params){
      return(Flex_lkl_1(data, params))
    }

    surv_gr = function(params){
      return(Flex_derivative_1(data, params, raw = T))
    }

    lower = rep(-Inf, length(params))
    upper = rep(Inf, length(params))

    if(grad){
      optimized_params = optim(par = params, fn = surv_fn, gr = surv_gr, method = 'L-BFGS',lower = lower, upper = upper, control = list(maxit = n.iter))
    }
    else{
      optimized_params = optim(par = params, fn = surv_fn, method = 'L-BFGS',lower = lower, upper = upper, control = list(maxit = n.iter))
    }

    if(raw){
      return(optimized_params$par)
    }
    else{
      params = optimized_params$par
      return( list( 'beta'  =      params[ 1: ncol(data$pred) ]                             ,
                    'delta' = exp( params[ (ncol(data$pred) + 1) : (length(params) - 2) ] ) ,
                    'nu'    = exp( params[ length(params) - 1 ]                           ) ,
                    'kappa' = exp( params[ length(params) ]   ) ))
    }
  }

  Flex_lkl_2 = function(data, params){

    pred = data$pred
    resp = data$resp
    cens = data$cens

    length_beta = ncol(pred)

    beta = params[1:length_beta]
    exp_delta = exp( params[(length_beta + 1) : (length(params) - 2)] )
    nu_0 = exp( params[ length(params) - 1 ] )
    kappa_0 = exp( params[ length(params) ] )

    n.person = nrow(resp)
    n.incid  = ncol(resp)

    cum_baseline_hazard = cumsum(exp_delta)

    exp_pred_beta = exp(aperm(sapply(1:dim(pred)[3], function(i) pred[,,i] %*% beta), c(2,1)))

    cum_hazard = array( cum_baseline_hazard[ resp ] , dim= c(n.person , n.incid) )

    cum_hazard_ahead = array( cum_baseline_hazard[ resp+1 ] , dim= c(n.person , n.incid) )

    nu_kappa = function(xi_ij , xi_plus_ij , nu){

      eps_ij = (xi_plus_ij - xi_ij) / xi_plus_ij
      eps_ij_rec = 1 - eps_ij

      eps_rec_nu    = eps_ij_rec   ^ nu
      eps_rec_nu_f  = eps_rec_nu   * eps_ij_rec
      eps_rec_nu_ff = eps_rec_nu_f * eps_ij_rec

      rec_rec       = 1 - eps_rec_nu
      rec_rec_f     = 1 - eps_rec_nu_f
      rec_rec_ff    = 1 - eps_rec_nu_ff

      denominator = ( rec_rec * rec_rec_ff ) - nu * ( eps_rec_nu * ( eps_ij ^ 2 ) )

      nu_new    = nu    * ( rec_rec_f^2         / denominator )
      kappa_new = xi_ij * ( rec_rec * rec_rec_f / denominator )

      return( list ( 'nu_new' = nu_new , 'kappa_new' = kappa_new ) )

    }

    nu    = matrix( rep( nu_0    , n.person*n.incid ) , nrow = n.person )
    kappa = matrix( rep( kappa_0 , n.person*n.incid ) , nrow = n.person )

    for(i in 2:n.incid){

      xi      = kappa[,i-1] + exp_pred_beta[,i-1] * cum_hazard[,i-1]

      xi_plus = kappa[,i-1] + exp_pred_beta[,i-1] * cum_hazard_ahead[,i-1]

      updated_params = nu_kappa(xi, xi_plus, nu[,i-1])

      nu[,i]    =  ifelse(cens[,i-1], nu[,i-1], updated_params$nu_new)
      kappa[,i] =  ifelse(cens[,i-1], xi, updated_params$kappa_new)

      }

    kappa_exp_pred_beta = kappa^(-1) * exp_pred_beta
    xi      = 1 + kappa_exp_pred_beta * cum_hazard
    xi_plus = 1 + kappa_exp_pred_beta * cum_hazard_ahead

    ind_likelihood =  (( xi )^(-nu) - ( xi_plus )^(-nu)*(1-cens))
    approx = (ind_likelihood<=0)
    approx = ifelse(is.na(approx),FALSE,approx)

    if(sum(approx)){

      xi_approx = exp_pred_beta * cum_hazard
      xi_plus_approx = exp_pred_beta * cum_hazard_ahead * (1-cens)

      #Laurent series
      ind_likelihood_approx = (nu * ( xi_plus_approx - xi_approx ) / kappa) + ( nu * (nu+1) * (xi_approx^2 - xi_plus_approx^2) ) / ( 2 * kappa^2 ) - ( nu * (nu+1) * (nu+2) * (xi_approx^3 - xi_plus_approx^3) ) / ( 6 * kappa^3 )

      ind_likelihood[approx] = ind_likelihood_approx[approx]

    }

    ind_log_likelihood = log(ind_likelihood)

    joint_log_likelihood  =  sum(  ind_log_likelihood ,na.rm = T )
    joint_log_likelihood = ifelse ( joint_log_likelihood==0, sum(ind_log_likelihood), joint_log_likelihood)

    return(-joint_log_likelihood)
  }

  Surv_optim_2 = function( data , init_params = NA , n.iter = 1e4, raw = T){

    if(is.na(init_params[1])){
      params = param_func(data)
    }
    else{
      params = init_params
    }

    surv_fn = function(params){
      return(Flex_lkl_2(data, params))
    }

    lower = rep(-Inf, length(params))
    upper = rep(Inf, length(params))

    optimized_params = optim(par = params, fn = surv_fn, method = 'L-BFGS',lower = lower, upper = upper,hessian = TRUE, control = list(maxit = n.iter))

    if(raw){
      return(optimized_params$par)
    }
    else{
      return(optimized_params)
    }

  }

  bootstrap_function = function(X, start_params, n_bootstrap = 1e3, bootstrap_n_iter = 5e3) {

    num_params = length(start_params)
    bootstrap_estimates = matrix(NA, n_bootstrap, num_params)

    for(i in 1:n_bootstrap) {
      resampled_X = X[sample(nrow(X), replace = TRUE), ]
      resampled_data = Array_func(X,intercept = intercept, psi = psi)

      bootstrap_estimates[i, ] = Surv_optim_2(resampled_data, start_params, n.iter = bootstrap_n_iter, raw = T)
    }
    bootstrap_se = apply(bootstrap_estimates, 2, sd)

    return(bootstrap_estimates)
  }

  data = Array_func(X,intercept = intercept , psi = psi)

  params_1 = Surv_optim_1(data, init_params = init_params, n.iter = as.integer(round(n_iter/10)))

  surv_est   = Surv_optim_2(data, init_params = params_1, n.iter = n_iter, raw = F)

  table_builder = function(est, positions, exponentiate = F){
    par_est  = est$par [ positions ]
    se       = sqrt ( diag ( solve ( est$hessian [ positions , positions ] ) ) )
    if(bootstrap){
      se = bootstrap_function(data, start_params = est$par, n_bootstrap = n_bootstrap, bootstrap_n_iter = bootstrap_n_iter)
    }
    CI_lower = par_est - qnorm(.975) * se
    CI_upper = par_est + qnorm(.975) * se
    wald_sts = par_est / se
    p_values = 2 * pnorm(-abs (wald_sts) )

    if(exponentiate){
      par_est = exp(par_est)
      se = NA
      CI_lower = exp(CI_lower)
      CI_upper = exp(CI_upper)
    }

    return(list('est' = par_est ,
                'se'  = se      ,
                'low' = CI_lower,
                'high'= CI_upper,
                'pval'= p_values))


  }

  table_builder = function(est, positions, exponentiate = F, init_lambda = 1e-8, lambda_decr = 0.5, max_tries = 1e3){
    par_est  = est$par [ positions ]
    se = NA
    current_lambda = init_lambda
    attempts = 0
    singular = TRUE

    while(singular && attempts < max_tries){
      attempts = attempts + 1
      hessian = est$hessian [ positions , positions ]
      diag(hessian) = diag(hessian) + current_lambda
      tryCatch({
        se = sqrt(diag(solve(hessian)))
        singular = FALSE
      },
      error = function(e){
        current_lambda = current_lambda / lambda_decr
      })
    }
    if(any(is.na(se))){
      stop('Unable to compute the standard errors after', max_tries, 'attempts.')
    }
    if(bootstrap){
      se = bootstrap_function(data, start_params = est$par, n_bootstrap = n_bootstrap, bootstrap_n_iter = bootstrap_n_iter)
    }
    CI_lower = par_est - qnorm(.975) * se
    CI_upper = par_est + qnorm(.975) * se
    wald_sts = par_est / se
    p_values = 2 * pnorm(-abs (wald_sts) )

    if(exponentiate){
      par_est = exp(par_est)
      se = NA
      CI_lower = exp(CI_lower)
      CI_upper = exp(CI_upper)
    }

    return(list('est' = par_est ,
                'se'  = se      ,
                'low' = CI_lower,
                'high'= CI_upper,
                'pval'= p_values))


  }

  params_length = length(surv_est$par)

  beta            = table_builder(surv_est, 1: ncol(data$pred))
  delta_positions = head(which(surv_est$par[ (ncol(data$pred) + 1) : (params_length - 2) ] !=-1e4),-1) + ncol(data$pred)
  delta           = table_builder(surv_est, delta_positions, exponentiate = T)
  nu_kappa        = table_builder(surv_est, (params_length - 1):params_length, exponentiate = T)

  delta_tot = exp( surv_est$par[ (ncol(data$pred) + 1) : (params_length - 2) ] )
  delta_positions = delta_positions - ncol(data$pred)
  format_number <- function(x) {
    if (abs(x) < 0.001 || abs(x) > 1000) {
      return(formatC(x, format = "e", digits = 3))
    } else {
      return(formatC(x, format = "f", digits = 3))
    }
  }

  add_significance_indicators <- function(p_value) {
    if (p_value < 0.001) {
      return(paste0(format_number(p_value), " ***"))
    } else if (p_value < 0.01) {
      return(paste0(format_number(p_value), " **"))
    } else if (p_value < 0.05) {
      return(paste0(format_number(p_value), " *"))
    } else if (p_value < 0.1) {
      return(paste0(format_number(p_value), " ."))
    } else {
      return(format_number(p_value))
    }
  }


  beta_table = cbind(format_number(beta$est),
                     format_number(beta$se),
                     paste('(',format_number(beta$low),',',format_number(beta$high),')',sep=''),
                     sapply(beta$pval, add_significance_indicators))
  colnames(beta_table) = c('Est.','Std. Error','95% C.I.', 'P-Value')
  rownames(beta_table) = c(colnames(X[, -match(c("id", "response","censor"), names(X))]))

  nu_table = cbind(format_number(nu_kappa$est),
                   paste('(',format_number(nu_kappa$low),',',format_number(nu_kappa$high),')',sep=''))
  colnames(nu_table) = c('Est', '95% CI')
  rownames(nu_table) = c('alpha: shape parameter', 'kappa: rate parameter')

  delta_est = delta_low = delta_high = delta_tot
  delta_low [delta_positions] = delta$low
  delta_high[delta_positions] = delta$high

  delta_high[is.infinite(delta_high)] = 1e10

  delta_est = cumsum(delta_est)
  delta_low = cumsum(delta_low)
  delta_high = cumsum(delta_high)

  delta_table = cbind(delta_est, delta_low, delta_high)

  times_axis = seq ( 1 , length(delta_est) ) * psi

  invisible ( plot ( times_axis , delta_est ,
              col  = adjustcolor('orange',1.5) , lwd = 2 , type = 'l' ,
              ylim = c(0 , 1.5 * max ( delta_est ) ) , xlim = c(0, times_axis[length(times_axis)]-.5 ) ,
              xlab = 'time' , ylab = 'baseline cumulative hazard' ) )
  invisible( polygon ( c( times_axis, rev(times_axis) ) , c( delta_low, rev(delta_high) ),
                       col = 'moccasin' , border = NA ) )
  invisible(lines  ( times_axis , delta_est , lwd = 2,
                     col = adjustcolor('orange' , 1.5) ) )
  invisible( grid() )

  delta_plot = recordPlot()

  return(list ( 'Coefficients' = beta_table ,
                'Frailty'      = nu_table   ,
                'Baseline Hazard' = delta_est[2:length(delta_est)],
                'Baseline Hazard Plot' = delta_plot ) )

}

