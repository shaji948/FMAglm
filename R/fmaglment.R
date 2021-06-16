fma.glmnet <-  function( x, y, focus, family, nlambda, nfolds, grouped, 
                         type.measure, rule, solnptol ){
  
  #### SJ: Only jackknife is supported in this version.
  #### SJ: Define some constants
  n <- nrow(x) #### Sample size
  nbeta <- ncol(x) + 1 #### Number of coefficients, including intercept
  res <- list() # Initiate returned list
  
  ########## Functions needed for weight estimations ##########
  #### SJ: Only useful if focus = "beta"
  if( focus == "beta" ){
    
    fun_prob = function(w, B, x, y, family){
      mu = matrix(NA, nrow = nrow(x), ncol = 1)
      for(i in 1 : nrow(x)){
        
        beta_bar = B[, , i] %*% w
        
        if( family == "binomial" ){
          mu[i] = plogis( sum( c(beta_bar) * c(1, x[i,]) ) )
        } else if( family == "poisson" ){
          mu[i] = exp( sum( c(beta_bar) * c(1, x[i,]) ) )
        }
        
      }
      return(mean((y - mu)^2))
    }
    
    #### VZ: New function needed for Singleton and hybrid averaging, generally a bit
    #### VZ: messier than the one above, since regressors are not the same in all models
    
    if("Singleton" %in% rule){
      fun_prob_sing <- function(w, B, x, y, family){
        mu <- matrix(0, nrow = nrow(x))
        
        for(i in 1:nrow(x)){   # Looping through folds
          
          if(family == "binomial"){
            mu[i] <- c(w[1]) * sum( c(B[, 1, i]) * c(1, x[i, 1]))
            for(j in 2:(nbeta - 1) ){  # Looping through individual Singletons
              mu[i] <-  mu[i] + c(w[j]) * sum( c(B[, j, i]) * c(1, x[i, j]))
            }
            mu[i] <- plogis(mu[i])
          }
          
          if(family == "poisson"){
            mu[i] <- c(w[1]) * sum( c(B[, 1, i]) * c(1, x[i, 1]))
            for(j in 2:(nbeta - 1) ){  # Looping through individual Singletons
              mu[i] <-  mu[i] + c(w[j]) * sum( c(B[, j, i]) * c(1, x[i, j]))
            }
            mu[i] <- exp(mu[i])
          }
        }
        return(mean((y - mu)^2))
      }
    }
    
    
    if("Hybrid" %in% rule){
      fun_prob_hybrid <- function(w, B, x, y, family, comb){
        mu <- matrix(0, nrow = nrow(x))
        
        for(i in 1:nrow(x)){   # Looping through folds
          
          #### VZ: Since the betas are not all of the same dimension, I could
          #### VZ: not find a better way of storing them than in a list
          
          if(family == "binomial"){
            for(j in 1:ncol(comb) ){  # Looping through individual models
              sel <- which(comb[, j] %in% 1) 
              mu[i] <-  mu[i] + c(w[j] %*% cbind(1, x)[i, sel] %*% B[[i]][[j]])
            }
            mu[i] <- plogis(mu[i])
          }
          
          if(family == "poisson"){
            for(j in 1:ncol(comb) ){  # Looping through individual models
              sel <- which(comb[, j] %in% 1) 
              mu[i] <-  mu[i] + c(w[j] %*% cbind(1, x)[i, sel] %*% B[[i]][[j]])
            }
            mu[i] <- exp(mu[i])
          }
        }
        return(mean((y - mu)^2))
      }
      
      fun_eq_hybrid = function(w, B, x, y, family, comb){
        z = sum(w)
        return(z)
      }
    }
    
    fun_eq = function(w, B, x, y, family){
      z = sum(w)
      return(z)
    }
    
  }
  
  #### ---------------------------------------------------- ####
  ####  SJ: Full model estimates
  if( family == "binomial" ){
    fulcoef <- glm.fit( x = cbind(1, x), y = y, family = binomial() )$coefficients
  } else if( family == "poisson" ){
    fulcoef <- glm.fit( x = cbind(1, x), y = y, family = poisson() )$coefficients
  }
  res$fullcoef <- fulcoef
  #### ---------------------------------------------------- ####
  
  #############################################################
  #### SJ: Run glmnet() to get the lambda sequence.
  #### Note: lambda sequence is decreasing.
  fit_alldata = glmnet(x, y, family = family, nfolds = nfolds,
                       grouped = grouped, nlambda = nlambda,
                       type.measure = type.measure)
  lambda_raw <- fit_alldata$lambda
  nlambda_raw <- length(lambda_raw)
  
  #### ---------------------------------------------------- ####
  #### SJ: We need to create the "correct" lambda sequence from the lambda sequence in glmnet()
  ####     This step depends on the quadrature rule that we want to use.
  ####     If we use "raw" (the discrete weight unit-simplex) or "Riemann" for integral, we do not need to introduce new lambda values.
  ####     At this stage, we do not include the full model (lambda = 0) into the lambda sequence
  if( "raw" %in% rule | "Riemann" %in% rule ){
    lambda <- lambda_raw
  }
  if( "Simpson(1/3)" %in% rule | "Simpson(3/8)" %in% rule ){
    
    #### SJ: We need to augment the raw lambda sequence
    if( "Simpson(1/3)" %in% rule ){
      
      lambda <- lambda_raw
      ####  Add the mid points to the lambda sequence
      l_seq = sort(c(lambda, 0))
      S = 2 : length(l_seq)
      l_mid = c((l_seq[S] + l_seq[S - 1])/2)
      lambda = numeric(length(l_mid) + length(l_seq))
      even = c(seq(from = 2, to = length(lambda) - 2, by = 2), length(lambda) - 1)
      lambda[even] = l_mid
      lambda[-even] = l_seq
      lambda <- lambda[-1] ####  We temporarily remove the full model.
      lambda <- sort(lambda, decreasing = TRUE) ####  glmnet() works with decreasing sequence.
      
    }
    if( "Simpson(3/8)" %in% rule ){
      
      lambda <- c( 0.0, sort(lambda_raw, decreasing = FALSE) )
      ####  Add the mid points to the lambda sequence
      lambdamat <- cbind( lambda[-length(lambda)], NA, NA, lambda[-1] )
      lambdamat[, 2] <- (2.0 * lambdamat[,1] + lambdamat[,4]) / 3.0
      lambdamat[, 3] <- (lambdamat[,1] + 2.0 * lambdamat[,4]) / 3.0
      lambda <- c(lambda_raw, lambdamat[,2], lambdamat[,3]) ####  We do not need the full model for cv.glmnet()
      lambda <- sort(lambda, decreasing = TRUE) ####  glmnet() works with decreasing sequence.
      
    }
    
  }
  nlambda <- length(lambda)
  #### ---------------------------------------------------- ####
  
  
  #### ---------------------------------------------------- ####
  #### SJ: Using the above lambda sequence, we can start the CV step.
  if( focus == "beta" ) {
    
    B <-  array(0, dim = c(ncol(x) + 1, nlambda + 1, nfolds))  # Array for betas
    B_raw <- array(0, dim = c(ncol(x) + 1, length(lambda_raw) + 1, nfolds))
    
    #### VZ: Adding corresponding storage for hybrid and singletons
    if("Hybrid" %in% rule){
      betas_lasso <-  as.matrix(coef(fit_alldata, s = lambda_raw))  # Coefficients for hybrid lambda
      u <-  betas_lasso != 0                    # Logical, TRUE if coefficient is not 0
      hybrid_comb <-  unique(u, MARGIN = 2)     # Filtering unique model forms
      #hybrid_comb <-  cbind(hybrid_comb, 1) 
      
      if(all(hybrid_comb[, ncol(hybrid_comb)] == 1) & 
         all(hybrid_comb[, ncol(hybrid_comb) - 1] == 1)){
        hybrid_comb <- hybrid_comb[, - length(hybrid_comb)]
      }
      
      B_hybrid = list()
    }
    if("Singleton" %in% rule){
      B_sing = array(0, dim = c(2, (nbeta - 1), nfolds))
    }
    
    for(i in 1 : nfolds){
      
      x_i <- x[-i, , drop = FALSE]
      y_i <- y[-i]
      
      ####  glmnet()
      fit <-  glmnet(x = x_i, y = y_i, family = family, lambda = lambda)
      fit_raw <- glmnet(x = x_i, y = y_i, family = family, lambda = lambda_raw)
      ####  full model
      if( family == "binomial" ){
        cvfit <- glm.fit( x = cbind(1, x_i), y = y_i, family = binomial() )
        
        #### VZ: Adding hybrid and singleton CV
        if( "Hybrid" %in% rule ){
          
          B_hybrid[[i]] <- list()
          
          for( j in 1:ncol(hybrid_comb) ){                      # Looping through model forms
            
            sel <- which(hybrid_comb[, j] %in% 1)               # Indices of regressors to use
            sel_fit <- glm.fit(cbind(1, x_i)[, sel], y_i, 
                               family = binomial())             # Fitting model
            B_hybrid[[i]][[j]] <- coef(sel_fit)                 # Extracting coefficients
          }
        }
        
        if( "Singelton" %in% rule ){                                    
          for( l in 1:(nbeta - 1) ){                              # Looping through Singletons
            
            sing_fit <- glm.fit(cbind(1, x_i[ , l]), y_i,         # Fitting Singletons
                                family = binomial())
            B_sing[, l, i] <-  matrix(coef(sing_fit))             # Extracting coefficients
          }
        }
      } else if( family == "poisson" ){
        cvfit <- glm.fit( x = cbind(1, x_i), y = y_i, family = poisson() )
        
        #### VZ: Adding hybrid and singleton CV
        if( "Hybrid" %in% rule ){
          
          B_hybrid[[i]] <- list()
          
          for( j in 1:ncol(hybrid_comb) ){                      # Looping through model forms
            
            sel <- which(hybrid_comb[, j] %in% 1)               # Indices of regressors to use
            sel_fit <- glm.fit(cbind(1, x_i)[, sel], y_i, 
                               family = poisson())             # Fitting model
            B_hybrid[[i]][[j]] <- coef(sel_fit)                 # Extracting coefficients
          }
        }
        
        if( "Singelton" %in% rule ){                                    
          for( l in 1:(nbeta - 1) ){                              # Looping through Singletons
            
            sing_fit <- glm.fit(cbind(1, x_i[ , l]), y_i,         # Fitting Singletons
                                family = poisson())
            B_sing[, l, i] <-  matrix(coef(sing_fit))             # Extracting coefficients
          }
        }
        
      }
      
      #### SJ: We add the full model as the last column, as lambda is a decreasing sequence
      B[, , i] = cbind( matrix(coef(fit), ncol = length(lambda)),
                        cvfit$coefficients )
      B_raw[, , i] = cbind( matrix(coef(fit_raw), ncol = length(lambda_raw)),
                            cvfit$coefficients )
      
    }
    
    #### SJ: In order to estimate the weights, we add 0 to the lambda sequence here.
    ####     Need to make sure that lambda remains a decreasing sequence!
    lambda <- c(lambda, 0.0)
    nlambda <- length(lambda)
    init_w = matrix( 1 / nlambda, nrow = nlambda )
    nlopt <- try( solnp(pars = init_w, fun = fun_prob, eqfun = fun_eq, eqB = c(1),
                        LB = rep(0, nlambda), UB = rep(1, nlambda),
                        B = B, x = x, y = y, family = family,
                        control = list(trace = FALSE, tol = solnptol)),
                  TRUE )
    
    lambda_raw_0 <- c(lambda_raw, 0.0)
    nlambda_raw_0 <- length(lambda_raw_0)
    init_w_raw = matrix( 1 / nlambda_raw_0, nrow = nlambda_raw_0 )
    nlopt_raw <- try( solnp(pars = init_w_raw, fun = fun_prob, eqfun = fun_eq, eqB = c(1),
                            LB = rep(0, nlambda_raw_0), UB = rep(1, nlambda_raw_0),
                            B = B_raw, x = x, y = y, family = family,
                            control = list(trace = FALSE, tol = solnptol)),
                      TRUE )
    
    if("Singleton" %in% rule){
      init_w_sing = matrix( 1 / (nbeta - 1), nrow = (nbeta - 1) )
      nlopt_sing <- try( solnp(pars = init_w_sing, fun = fun_prob_sing, eqfun = fun_eq, 
                               eqB = c(1),
                               LB = rep(0, (nbeta - 1)), UB = rep(1, (nbeta - 1)),
                               B = B_sing, x = x, y = y, family = family,
                               control = list(trace = FALSE, tol = solnptol)),
                         TRUE)
      
      
      res[["Singleton"]] <- list()
      if( inherits(nlopt_sing, "try-error") == FALSE){
        res[["Singleton"]]$w <- nlopt_sing$pars
        res[["Singleton"]]$convergence <- nlopt_sing$convergence
        res[["Singleton"]]$elapsed <- nlopt_sing$elapsed
      }else{
        res[["Singleton"]] <- nlopt_sing
      }
    }
    
    if("Hybrid" %in% rule){
      init_w_hybrid = matrix( 1 / ncol(hybrid_comb), nrow = ncol(hybrid_comb) )
      nlopt_hybrid <- try( solnp(pars = init_w_hybrid, fun = fun_prob_hybrid, 
                                 eqfun = fun_eq_hybrid, eqB = c(1), 
                                 LB = rep(0, ncol(hybrid_comb)), 
                                 UB = rep(1, ncol(hybrid_comb)), 
                                 B = B_hybrid, x = x, y = y,
                                 family = family, comb = hybrid_comb,
                                 control = list(trace = FALSE, tol = solnptol)),
                           TRUE)
      
      res[["Hybrid"]] <- list()
      if( inherits(nlopt_hybrid, "try-error") == FALSE){
        res[["Hybrid"]]$w <- nlopt_hybrid$pars
        res[["Hybrid"]]$convergence <- nlopt_hybrid$convergence
        res[["Hybrid"]]$elapsed <- nlopt_hybrid$elapsed
        res[["Hybrid"]]$models <- hybrid_comb 
      }else{
        res[["Hybrid"]] <- nlopt_hybrid
      }
    }
    
    if( inherits(nlopt_raw, "try-error") == FALSE ){
      
      #### Note: lambda sequence is decreasing.
      fit_alldata = glmnet(x, y, family = family, lambda = lambda[-nlambda] )
      coef_alldata <- as.matrix( coef( object = fit_alldata ) )
      coef_alldata <- cbind(coef_alldata, fulcoef)
      
      ####
      res$lambda <- lambda
      res$cvval <- fun_prob( w = nlopt_raw$pars, B = B_raw, x = x, y = y, family = family )
      weights_path <- as.matrix(nlopt_raw$pars)
      convergence <- nlopt_raw$convergence #### 0 = converged
      #avecoef <- coef_alldata %*% weights_path
      elapsed <- nlopt_raw$elapsed
      #### SJ: Put the estimates in the right slot, depending on the quadrature method
      if( "raw" %in% rule ){
        res[["raw"]] <- list()
        res[["raw"]]$w <- nlopt_raw$pars
        res[["raw"]]$convergence <- convergence
        #res[["raw"]]$avecoef <- avecoef
        res[["raw"]]$elapsed <- elapsed
        res[["raw"]]$lambda <- lambda_raw_0
      }
      if( "Riemann" %in% rule ){
        res[["Riemann"]] <- list()
        res[["Riemann"]]$v <- nlopt$pars
        lambdamat <- rbind( lambda[-1], lambda[-nlambda])
        res[["Riemann"]]$w <- rbind( lambda[-1],
                                     nlopt$pars[-1] * 2.0 / (lambdamat[2,] - lambdamat[1,]) )
        rownames(res[["Riemann"]]$w) <- c("lambda", "w")
        res[["Riemann"]]$convergence <- convergence
        #res[["Riemann"]]$avecoef <- avecoef
        res[["Riemann"]]$elapsed <- elapsed
      }
      if( "Simpson(1/3)" %in% rule ){
        res[["Simpson(1/3)"]] <- list()
        res[["Simpson(1/3)"]]$v <- nlopt$pars
        lambda_raw_zero <- c(lambda_raw, 0.0)
        scale <- (lambda_raw_zero[1] - lambda_raw_zero[2]) / 6
        for( i in 1 : nlambda_raw ){
          scale <- c(scale,
                     2 * (lambda_raw_zero[i] - lambda_raw_zero[i + 1]) / 3,
                     (lambda_raw_zero[i] - lambda_raw_zero[i + 2]) / 6 )
        }
        scale[nlambda] <- (lambda_raw_zero[i] - lambda_raw_zero[i + 1]) / 6
        res[["Simpson(1/3)"]]$w <- rbind( lambda,
                                          nlopt$pars / scale )
        #### SJ: You had the following, which is not correct.
        #  w_res[even] = (2/3)*(l_seq[S] - l_seq[S - 1])
        #  w_res[-even] = (1/6)*c((l_seq[S] - l_seq[S - 1]),(l_seq[length(l_seq)] - l_seq[length(l_seq) - 1]))
        rownames(res[["Simpson(1/3)"]]$w) <- c("lambda", "w")
        res[["Simpson(1/3)"]]$convergence <- convergence
        #res[["Simpson(1/3)"]]$avecoef <- avecoef
        res[["Simpson(1/3)"]]$elapsed <- elapsed
      }
      if( "Simpson(3/8)" %in% rule ){
        res[["Simpson(3/8)"]] <- list()
        res[["Simpson(3/8)"]]$v <- nlopt$pars
        #### SJ: Haven't finished this one.
        res[["Simpson(3/8)"]]$convergence <- convergence
        #res[["Simpson(3/8)"]]$avecoef <- avecoef
        res[["Simpson(3/8)"]]$elapsed <- elapsed
      }
    } else {
      res <- nlopt
    }
    
  } else if( focus == "mu" ){
    
    # Fitting the LASSO model
    cvglm = try(cv.glmnet(x, y, family = family, nfolds = nfolds,
                          grouped = grouped, keep = TRUE, lambda = lambda,
                          type.measure = type.measure ), TRUE)
    
    #### VZ: Finding unique model forms
    betas_lasso = as.matrix(coef(cvglm, s = lambda_raw))  # Coefficients for hybrid lambda
    u = betas_lasso != 0                    # Logical, TRUE if coefficient is not 0
    hybrid_comb = unique(u, MARGIN = 2)     # Filtering unique model forms
    # hybrid_comb = cbind(hybrid_comb, 1)   # Adding full model
    
    # if(all(hybrid_comb[, ncol(hybrid_comb)] == 1) & 
    #   all(hybrid_comb[, ncol(hybrid_comb) - 1] == 1)){
    #  hybrid_comb <- hybrid_comb[, - length(hybrid_comb)]
    # }
    
    ####  SJ: fit full model (returning eta)
    ####  VZ: Fit hybrid and Singletons
    fulmod_eta <- rep( NA, nrow(x) )
    hybrid_eta <- matrix(NA, nrow = nrow(x), ncol = ncol(hybrid_comb) ) 
    singleton_eta <- matrix(0, nrow = nrow(x), ncol = (nbeta - 1) )
    for( i in 1 : nrow(x) ){
      
      x_i <- x[-i,, drop = FALSE]
      y_i <- y[-i]
      if( family == "binomial" ){
        cvfit <- glm.fit( x = cbind(1, x_i), y = y_i, family = binomial() )
        
        #### VZ: Hybrid CV
        if( "Hybrid" %in% rule ){
          for( j in 1:ncol(hybrid_comb) ){                         # Looping through model forms
            sel <-  which(hybrid_comb[, j] %in% 1)               # Indices of regressors to use
            sel_fit <-  glm.fit(cbind(1, x_i)[, sel], y_i, 
                                family = binomial())               # Fitting model
            b_sel <-  matrix(coef(sel_fit))                      # Extracting coefficients
            hybrid_eta[i, j] <- cbind(1, x)[i, sel] %*% b_sel    # Predicting# Computing eta
          }
        }
        
        #### VZ: Singleton CV
        if( "Singelton" %in% rule ){                                    
          for( l in 1:(nbeta - 1) ){                              # Looping throug Singletons
            sing_fit <- glm.fit(cbind(1, x_i[ , l]), y_i,          # Fitting Singletons
                               family = binomial())
            b_sing <-  matrix(coef(sing_fit))                     # Extracting coefficients
            singleton_eta[i, l] <- cbind(1, x[i, l]) %*% b_sing   # Predicting
          }
        }
        
      } else if( family == "poisson" ){
        cvfit <- glm.fit( x = cbind(1, x_i), y = y_i, family = poisson() )
        
        #### VZ: Hybrid CV
        if( "Hybrid" %in% rule ){
          for( j in 1:ncol(hybrid_comb) ){                       # Looping through model forms
            sel <- which(hybrid_comb[, j] %in% 1)                # Indices of regressors to use
            sel_fit <- glm.fit(cbind(1, x_i)[, sel], y_i, 
                               family = poisson())               # Fitting model
            b_sel <-  matrix(coef(sel_fit))                      # Extracting coefficients
            hybrid_eta[i, j] <- cbind(1, x)[i, sel] %*% b_sel    # Predicting# Computing eta
            
          }
        }
        
        #### VZ: Singleton CV
        if( "Singelton" %in% rule ){                                    
          for( l in 1:(nbeta - 1) ){                              # Looping through Singletons
            sing_fit <- glm.fit(cbind(1, x_i[ , l]), y_i,         # Fitting Singletons
                                family = poisson())
            b_sing <-  matrix(coef(sing_fit))                     # Extracting coefficients
            singleton_eta[i, l] <- cbind(1, x[i, l]) %*% b_sing   # Predicting
            
          }
          
        }
      }
      fulmod_eta[i] <- sum(cvfit$coefficients * c(1.0, x[i,]) )
      
    }
    
    
    #### VZ: Adding non-augmented lambda weights
    lambda_index <- match(lambda_raw, lambda)
    cv_eta_raw <- cbind(cvglm$fit.preval[, lambda_index], fulmod_eta)
    lambda_raw_0 <- c(lambda_raw, 0.0)
    nlambda_raw_0 <- length(lambda_raw_0)
    if( family == "binomial" ){
      mu_cv_raw = plogis( cv_eta_raw )
    } else if( family == "poisson" ){
      mu_cv_raw = exp( cv_eta_raw )
    }
    
    #### The symmetric matrix in the quadratic programming
    qmat_raw <- (1/n) * crossprod( matrix(y, nrow = n, ncol = nlambda_raw_0) - mu_cv_raw )
    
    ####  Quadratic programming
    qrsolve_raw <- try( kernlab::ipop( c = rep(0.0, nlambda_raw_0),
                                       H = 2.0 * qmat_raw,
                                       A = matrix(1.0, nrow = 1, ncol = nlambda_raw_0),
                                       b = 1.0,
                                       l = rep(0.0, nlambda_raw_0),
                                       u = rep(1.0, nlambda_raw_0),
                                       r = 0.0),
                        TRUE)
    
    ####  In version 4.0-2, it is the linear predictor scale
    ####  We also add the full model.
    cv_eta <- cbind(cvglm$fit.preval, fulmod_eta)
    lambda <- c(lambda, 0.0)
    nlambda <- length(lambda)
    if( family == "binomial" ){
      mu_cv = plogis( cv_eta )
    } else if( family == "poisson" ){
      mu_cv = exp( cv_eta )
    }
    
    #### The symmetric matrix in the quadratic programming
    qmat <- (1/n)*crossprod( matrix(y, nrow = n, ncol = nlambda) - mu_cv )
    
    ####  Quadratic programming
    qrsolve <- try( kernlab::ipop( c = rep(0.0, nlambda),
                                   H = 2.0 * qmat,
                                   A = matrix(1.0, nrow = 1, ncol = nlambda),
                                   b = 1.0,
                                   l = rep(0.0, nlambda),
                                   u = rep(1.0, nlambda),
                                   r = 0.0),
                    TRUE)
    
    
    #### VZ: Adding QP for hybrid averaging
    if( "Hybrid" %in% rule ){
      
      if( family == "binomial" ){
        mu_cv_hybrid = plogis( hybrid_eta )
      } else if( family == "poisson" ){
        mu_cv_hybrid = exp( hybrid_eta )
      }
      
      # Symmetric matrix for QP
      qmat_hybrid <- (1/n)*crossprod( matrix(y, nrow = n, ncol = ncol(hybrid_eta)) - mu_cv_hybrid) 
      #Quadratic programmint
      qrsolve_hybrid <- try(kernlab::ipop( c = rep(0.0, ncol(hybrid_eta)),
                                           H = 2.0 * qmat_hybrid,
                                           A = matrix(1.0, ncol = ncol(hybrid_eta)),
                                           b = 1.0,
                                           l = rep(0.0, ncol(hybrid_eta)),
                                           u = rep(1.0, ncol(hybrid_eta)),
                                           r = 0.0), 
                            TRUE)
      
      #### VZ: Adding results to function output
      res[["Hybrid"]] <- list()
      if( inherits(qrsolve_hybrid, "try-error") == FALSE){
        res[["Hybrid"]]$w <- qrsolve_hybrid@primal
        res[["Hybrid"]]$convergence <- qrsolve_hybrid@how
        res[["Hybrid"]]$models <- hybrid_comb 
        res[["Hybrid"]]$pos_def <- all(eigen(qmat_hybrid)$values > 0) 
        
      }else{
        res[["Hybrid"]] <- qrsolve_hybrid
      }
    }
    
    #### VZ: Adding QP for Singleton averaging
    if( "Singleton" %in% rule){
      
      if( family == "binomial" ){
        mu_cv_singleton = plogis( singleton_eta )
      } else if( family == "poisson" ){
        mu_cv_singleton = exp( singleton_eta )
      }
      
      # Symmetric matrix for QP
      qmat_singleton <- (1/n)*crossprod( matrix(y, nrow = n, ncol = ncol(singleton_eta)) - mu_cv_singleton) 
      
      #Quadratic programming
      qrsolve_singleton <- try(kernlab::ipop( c = rep(0.0, ncol(singleton_eta)),
                                              H = 2.0 * qmat_singleton,
                                              A = matrix(1.0, ncol = ncol(singleton_eta)),
                                              b = 1.0,
                                              l = rep(0.0, ncol(singleton_eta)),
                                              u = rep(1.0, ncol(singleton_eta)),
                                              r = 0.0), 
                               TRUE)
      
      #### VZ: Adding results to function output
      res[["Singleton"]] <- list()
      if( inherits(qrsolve_singleton, "try-error") == FALSE){
        res[["Singleton"]]$w <- qrsolve_singleton@primal
        res[["Singleton"]]$convergence <- qrsolve_singleton@how
        res[["Singleton"]]$pos_def <- all(eigen(qmat_singleton)$values > 0) 
      }else{
        res[["Singleton"]] <- qrsolve_singleton
      }
      
    }
    
    if( inherits(qrsolve_raw, "try-error") == FALSE ){
      
      ####
      res$lambda <- lambda
      res$cvval <- matrix(qrsolve@primal, nrow = 1 ) %*% qmat %*% matrix(qrsolve@primal, ncol = 1 )
      convergence <- qrsolve@how
      
      #### SJ: Put the estimates in the right slot, depending on the quadrature method
      if( "raw" %in% rule ){
        res[["raw"]] <- list()
        res[["raw"]]$w <- qrsolve_raw@primal
        res[["raw"]]$convergence <- convergence
        res[["raw"]]$lambda <- lambda_raw_0
        res[["raw"]]$pos_def <- all(eigen(qmat_raw)$values > 0) 
      }
      if( "Riemann" %in% rule ){
        res[["Riemann"]] <- list()
        res[["Riemann"]]$v <- qrsolve@primal
        lambdamat <- rbind( lambda[-1], lambda[-nlambda])
        res[["Riemann"]]$w <- rbind( lambda[-1],
                                     qrsolve@primal[-1] * 2.0 / (lambdamat[2,] - lambdamat[1,]) )
        rownames(res[["Riemann"]]$w) <- c("lambda", "w")
        res[["Riemann"]]$convergence <- convergence
      }
      if( "Simpson(1/3)" %in% rule ){
        res[["Simpson(1/3)"]] <- list()
        res[["Simpson(1/3)"]]$v <- qrsolve@primal
        lambda_raw_zero <- c(lambda_raw, 0.0)
        scale <- (lambda_raw_zero[1] - lambda_raw_zero[2]) / 6
        for( i in 1 : nlambda_raw ){
          scale <- c(scale,
                     2 * (lambda_raw_zero[i] - lambda_raw_zero[i + 1]) / 3,
                     (lambda_raw_zero[i] - lambda_raw_zero[i + 2]) / 6 )
        }
        scale[nlambda] <- (lambda_raw_zero[i] - lambda_raw_zero[i + 1]) / 6
        res[["Simpson(1/3)"]]$w <- rbind( lambda,
                                          qrsolve@primal / scale )
        rownames(res[["Simpson(1/3)"]]$w) <- c("lambda", "w")
        res[["Simpson(1/3)"]]$convergence <- convergence
        res[["Simpson(1/3)"]]$pos_def <- all(eigen(qmat)$values > 0) 
      }
      if( "Simpson(3/8)" %in% rule ){
        res[["Simpson(3/8)"]] <- list()
        res[["Simpson(3/8)"]]$v <- qrsolve@primal
        #### SJ: Haven't finished this one.
        res[["Simpson(3/8)"]]$convergence <- convergence
      }
      res$cvglmnet <- cvglm
      
    } else {
      res <- qrsolve
    }
    
  }
  #### ---------------------------------------------------- ####
  
  res
  
}


