fma.glm <- function(x, y, focus, family, method, glmnetfit = NULL, nfolds, grouped, nlambda = 100, solnptol = 1e-08 ){

    n <- nrow(x) # Sample size
    k <- ncol(x) # Number of regressors, excluding intercept
    nbeta <- k + 1

    ####  If we do not have the glment() object, we need to compute it here.
    ####  We do not need a cv.glmnet() here, because we only need some model forms, not the estimates
    if( is.null(glmnetfit) == TRUE ){
        glmnetfit <- glmnet( x, y, family = family, nlambda = nlambda )
    }
    lambda <- glmnetfit$lambda

    ####  We only ues the glment package to suggest possible candidate model forms.
    res <- list()
    if( "JMA" %in% method ){

        res[["JMA"]] <- list()

        ####  Extract coefficients from cvglmnet
        beta_glmnet <- as.matrix(coef(glmnetfit, s = lambda, exact = TRUE)) ####  Same if remove exact = TRUE.
        u = (beta_glmnet != 0)              # Logical, TRUE if coefficient is not 0
        u_comb = unique(u, MARGIN = 2)     # Filtering unique model forms
        ####  If the full model is not included, we need to add the full model.
        if( all(colSums(u_comb) != nbeta) ){
            u_comb = cbind(u_comb, 1)
        }
        M_JMA = ncol(u_comb)  # Number of unique model forms

        ####  Initiate beta matrix or CV mu matrix
        if( focus == "beta" ){
            B_JMA = array(0, dim = c(nbeta, M_JMA, nfolds))
        } else if( focus == "mu" ){
            cv_mu_JMA <- matrix( NA, n, M_JMA )
        }

        ####  Fit the model using all observations
        est_allobs <- matrix(0, nbeta, M_JMA)
        for(m in 1 : M_JMA){

            sel = which(u_comb[, m] %in% 1)                       # Indices of regressors to use
            ####  Fit the model
            if( family == "binomial" ){
                sel_fit = glm.fit( x = cbind(1, x)[, sel],
                                   y = y,
                                   family = binomial() )
            } else if( family == "poisson" ){
                sel_fit = glm.fit( x = cbind(1, x),
                                   y = y,
                                   family = poisson() )
            }
            est_allobs[sel, m] <- coef(sel_fit)

        }
        res[["JMA"]]$beta <- est_allobs

    }
    if( "singleton" %in% method ){

        res[["singleton"]] <- list()
        M_sing <- k + 1 # We also allow the intercept only model to be included.

        ####  Initiate beta matrix or CV mu matrix
        if( focus == "beta" ){
            B_sing = array(0, dim = c(nbeta, M_sing, nfolds))
        } else if( focus == "mu" ){
            cv_mu_sing <- matrix( NA, n, M_sing )
        }

        ####  Fit the model using all observations
        est_allobs <- matrix(0, nbeta, M_sing)
        for(m in 0 : k){

            ####  We first add the intercept only model
            if( m == 0 ){

                if( family == "binomial" ){
                    sing_fit = glm.fit( x = matrix(1, nrow(x), 1),
                                        y = y,
                                        family = binomial() )
                } else if( family == "poisson" ){
                    sing_fit = glm.fit( x = matrix(1, nrow(x), 1),
                                        y = y,
                                        family = poisson() )
                }
                est_allobs[1, 1] <- coef(sing_fit)

            }
            ####  We then add other singleton models
            if( m != 0 ){

                if( family == "binomial" ){
                    sing_fit = glm.fit( x = cbind(1, x[, m]),
                                        y = y,
                                        family = binomial() )
                } else if( family == "poisson" ){
                    sing_fit = glm.fit( x = cbind(1, x[, m]),
                                        y = y,
                                        family = poisson() )
                }
                est_allobs[ c(1, m + 1), m + 1 ] <- coef(sing_fit)

            }

        }
        res[["singleton"]]$beta <- est_allobs

    }

    ####  Loop through the folds in CV
    for(i in 1 : n){

        x_i <- x[-i, , drop = F]
        y_i <- y[-i]
        # Looping through model forms
        if( focus == "beta" ){
            if( "JMA" %in% method ){
                beta_JMA <- matrix( 0.0, nbeta, M_JMA )
            }
            if( "singleton" %in% method ){
                beta_sing <- matrix( 0.0, nbeta, M_sing )
            }
        }

        ####  Start with JMA, if applicable
        if( "JMA" %in% method ){

            for(m in 1 : M_JMA){

                sel = which(u_comb[, m] %in% 1)                       # Indices of regressors to use
                ####  Fit the model
                if( family == "binomial" ){
                    sel_fit = glm.fit( x = cbind(1, x_i)[, sel],
                                       y = y_i,
                                       family = binomial() )
                    b_sel = coef(sel_fit)
                    if( focus == "beta" ){
                        beta_JMA[sel, m] <- b_sel
                    } else if( focus == "mu" ){
                        cv_mu_JMA[i, m] <- plogis( sum( c(1, x[i,])[sel] * b_sel ) )
                    }
                } else if( family == "poisson" ){
                    sel_fit = glm.fit( x = cbind(1, x_i)[, sel],
                                       y = y_i,
                                       family = poisson() )
                    b_sel = coef(sel_fit)
                    if( focus == "beta" ){
                        beta_JMA[sel, m] <- b_sel
                    } else if( focus == "mu" ){
                        cv_mu_JMA[i, m] <- exp( sum( c(1, x[i,])[sel] * b_sel ) )
                    }
                }

            }

            if( focus == "beta" ){
                B_JMA[ , , i] <- beta_JMA
            }

        }

        ####  Then, singleton models, if applicable
        if( "singleton" %in% method ){

            for(m in 0 : k){

                ####  We first add the intercept only model
                if( m == 0 ){

                    if( family == "binomial" ){
                        sing_fit = glm.fit( x = matrix(1, nrow(x_i), 1),
                                            y = y_i,
                                            family = binomial() )
                        b_sing = coef(sing_fit)
                        if( focus == "beta" ){
                            beta_sing[1, 1] <- b_sing
                        } else if( focus == "mu" ){
                            cv_mu_sing[i, 1] <- plogis( sum( 1.0 * b_sing ) )
                        }
                    } else if( family == "poisson" ){
                        sing_fit = glm.fit( x = matrix(1, nrow(x_i), 1),
                                            y = y_i,
                                            family = poisson() )
                        b_sing = coef(sing_fit)
                        if( focus == "beta" ){
                            beta_sing[1, 1] <- b_sing
                        } else if( focus == "mu" ){
                            cv_mu_sing[i, 1] <- exp( sum( 1.0 * b_sing ) )
                        }
                    }

                }
                ####  We then add other singleton models
                if( m != 0 ){

                    if( family == "binomial" ){
                        sing_fit = glm.fit( x = cbind(1, x_i[, m]),
                                            y = y_i,
                                            family = binomial() )
                        b_sing = coef(sing_fit)
                        if( focus == "beta" ){
                            beta_sing[ c(1, m + 1), m + 1 ] <- b_sing
                        } else if( focus == "mu" ){
                            cv_mu_sing[i, m + 1] <- plogis( sum( c(1, x[i, m]) * b_sing ) )
                        }
                    } else if( family == "poisson" ){
                        sing_fit = glm.fit( x = cbind(1, x_i[, m]),
                                            y = y_i,
                                            family = poisson() )
                        b_sing = matrix(coef(sing_fit))
                        if( focus == "beta" ){
                            beta_sing[ c(1, m + 1), m + 1 ] <- b_sing
                        } else if( focus == "mu" ){
                            cv_mu_sing[i, m + 1] <- exp( sum( c(1, x[i, m]) * b_sing ) )
                        }
                    }

                }

            }

            if( focus == "beta" ){
                B_sing[ , , i] <- beta_sing
            }

        }

    }

    ####  Depending on the focus, we choose between nonlinear optimization and quadratic programming
    if( focus == "beta" ){

        fun_prob = function(w, B, x, y, family){
            mu = matrix(NA, nrow = nrow(x), ncol = 1)
            for(i in 1 : nfolds){

                beta_bar = B[, , i] %*% w

                if( family == "binomial" ){
                    mu[i] = plogis( sum( c(beta_bar) * c(1, x[i,]) ) )
                } else if( family == "poisson" ){
                    mu[i] = exp( sum( c(beta_bar) * c(1, x[i,]) ) )
                }

            }
            return(sum((y - mu)^2))
        }

        fun_eq = function(w, B, x, y, family){
            z = sum(w)
            return(z)
        }

        if( "JMA" %in% method ){

            init_w = matrix( 1 / M_JMA, nrow = M_JMA )
            nlopt <- try( solnp(pars = init_w, fun = fun_prob, eqfun = fun_eq, eqB = c(1),
                                LB = rep(0, M_JMA), UB = rep(1, M_JMA),
                                B = B_JMA, x = x, y = y, family = family,
                                control = list(trace = FALSE, tol = solnptol)),
                          TRUE )

            if( inherits(nlopt, "try-error") == FALSE ){

                res[["JMA"]]$w <- nlopt$pars
                res[["JMA"]]$cvval <- fun_prob( w = nlopt$pars, B = B_JMA, x = x, y = y, family = family )
                res[["JMA"]]$convergence <- nlopt$convergence #### 0 = converged
                res[["JMA"]]$avecoef <- res[["JMA"]]$beta %*% matrix(res[["JMA"]]$w, ncol = 1)

            }

        }
        if( "singleton" %in% method ){

            fun_prob = function(w, B, x, y, family){
                mu = matrix(NA, nrow = nrow(x), ncol = 1)
                for(i in 1 : nfolds){

                    beta_bar = B[, , i] %*% w

                    if( family == "binomial" ){
                        mu[i] = plogis( sum( c(beta_bar) * c(1, x[i,]) ) )
                    } else if( family == "poisson" ){
                        mu[i] = exp( sum( c(beta_bar) * c(1, x[i,]) ) )
                    }

                }
                return(sum((y - mu)^2))
            }

            fun_eq = function(w, B, x, y, family){
                z = sum(w)
                return(z)
            }

            init_w = matrix( 1 / M_sing, nrow = M_sing )
            nlopt <- try( solnp(pars = init_w, fun = fun_prob, eqfun = fun_eq, eqB = c(1),
                                LB = rep(0, M_sing), UB = rep(1, M_sing),
                                B = B_sing, x = x, y = y, family = family,
                                control = list(trace = FALSE, tol = solnptol)),
                          TRUE )

            if( inherits(nlopt, "try-error") == FALSE ){

                res[["singleton"]]$w <- nlopt$pars
                res[["singleton"]]$cvval <- fun_prob( w = nlopt$pars, B = B_sing, x = x, y = y, family = family )
                res[["singleton"]]$convergence <- nlopt$convergence #### 0 = converged
                res[["singleton"]]$avecoef <- res[["singleton"]]$beta %*% matrix(res[["singleton"]]$w, ncol = 1)

            }

        }

    } else if( focus == "mu" ){

        if( "JMA" %in% method ){

            #### The symmetric matrix in the quadratic programming
            qmat <- crossprod( matrix(y, nrow = n, ncol = M_JMA) - cv_mu_JMA )

            ####  Quadratic programming
            qrsolve <- try( kernlab::ipop( c = rep(0.0, M_JMA),
                                           H = 2.0 * qmat,
                                           A = matrix(1.0, nrow = 1, ncol = M_JMA),
                                           b = 1.0,
                                           l = rep(0.0, M_JMA),
                                           u = rep(1.0, M_JMA),
                                           r = 0.0),
                            TRUE)

            if( inherits(qrsolve, "try-error") == FALSE ){

                res[["JMA"]]$w <- qrsolve@primal
                res[["JMA"]]$cvval <- c( matrix(qrsolve@primal, nrow = 1 ) %*% qmat %*% matrix(qrsolve@primal, ncol = 1 ) )
                res[["JMA"]]$convergence <- qrsolve@how

            }

        }
        if( "singleton" %in% method ){

            #### The symmetric matrix in the quadratic programming
            qmat <- crossprod( matrix(y, nrow = n, ncol = M_sing) - cv_mu_sing )

            ####  Quadratic programming
            qrsolve <- try( kernlab::ipop( c = rep(0.0, M_sing),
                                           H = 2.0 * qmat,
                                           A = matrix(1.0, nrow = 1, ncol = M_sing),
                                           b = 1.0,
                                           l = rep(0.0, M_sing),
                                           u = rep(1.0, M_sing),
                                           r = 0.0),
                            TRUE)

            if( inherits(qrsolve, "try-error") == FALSE ){

                res[["singleton"]]$w <- qrsolve@primal
                res[["singleton"]]$cvval <- c( matrix(qrsolve@primal, nrow = 1 ) %*% qmat %*% matrix(qrsolve@primal, ncol = 1 ) )
                res[["singleton"]]$convergence <- qrsolve@how

            }

        }



    }

    res

}

