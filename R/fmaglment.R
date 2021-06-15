fma.glmnet = function( x, y, focus, family, nlambda, nfolds, grouped, type.measure, rule, solnptol ){

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
            return(sum((y - mu)^2))
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

        B = array(0, dim = c(ncol(x) + 1, nlambda + 1, nfolds))  # Array for betas
        for(i in 1 : nfolds){

            x_i <- x[-i,, drop = FALSE]
            y_i <- y[-i]

            ####  glmnet()
            fit = glmnet(x = x_i, y = y_i, family = family, lambda = lambda)
            ####  full model
            if( family == "binomial" ){
                cvfit <- glm.fit( x = cbind(1, x_i), y = y_i, family = binomial() )
            } else if( family == "poisson" ){
                cvfit <- glm.fit( x = cbind(1, x_i), y = y_i, family = poisson() )
            }

            #### SJ: We add the full model as the last column, as lambda is a decreasing sequence
            B[, , i] = cbind( matrix(coef(fit), ncol = length(lambda)),
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

        if( inherits(nlopt, "try-error") == FALSE ){

            #### Note: lambda sequence is decreasing.
            fit_alldata = glmnet(x, y, family = family, lambda = lambda[-nlambda] )
            coef_alldata <- as.matrix( coef( object = fit_alldata ) )
            coef_alldata <- cbind(coef_alldata, fulcoef)

            ####
            res$lambda <- lambda
            res$cvval <- fun_prob( w = nlopt$pars, B = B, x = x, y = y, family = family )
            weights_path = as.matrix(nlopt$pars)
            convergence <- nlopt$convergence #### 0 = converged
            avecoef <- coef_alldata %*% weights_path
            elapsed <- nlopt$elapsed
            #### SJ: Put the estimates in the right slot, depending on the quadrature method
            if( "raw" %in% rule ){
                res[["raw"]] <- list()
                res[["raw"]]$w <- nlopt$pars
                res[["raw"]]$convergence <- convergence
                res[["raw"]]$avecoef <- avecoef
                res[["raw"]]$elapsed <- elapsed
            }
            if( "Riemann" %in% rule ){
                res[["Riemann"]] <- list()
                res[["Riemann"]]$v <- nlopt$pars
                lambdamat <- rbind( lambda[-1], lambda[-nlambda])
                res[["Riemann"]]$w <- rbind( lambda[-1],
                                             nlopt$pars[-1] * 2.0 / (lambdamat[2,] - lambdamat[1,]) )
                rownames(res[["Riemann"]]$w) <- c("lambda", "w")
                res[["Riemann"]]$convergence <- convergence
                res[["Riemann"]]$avecoef <- avecoef
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
                res[["Simpson(1/3)"]]$avecoef <- avecoef
                res[["Simpson(1/3)"]]$elapsed <- elapsed
            }
            if( "Simpson(3/8)" %in% rule ){
                res[["Simpson(3/8)"]] <- list()
                res[["Simpson(3/8)"]]$v <- nlopt$pars
                #### SJ: Haven't finished this one.
                res[["Simpson(3/8)"]]$convergence <- convergence
                res[["Simpson(3/8)"]]$avecoef <- avecoef
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

        ####  SJ: fit full model (returning eta)
        fulmod_eta <- rep( NA, nrow(x) )
        for( i in 1 : nrow(x) ){

            x_i <- x[-i,, drop = FALSE]
            y_i <- y[-i]
            if( family == "binomial" ){
                cvfit <- glm.fit( x = cbind(1, x_i), y = y_i, family = binomial() )
            } else if( family == "poisson" ){
                cvfit <- glm.fit( x = cbind(1, x_i), y = y_i, family = poisson() )
            }
            fulmod_eta[i] <- sum(cvfit$coefficients * c(1.0, x[i,]) )

        }

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
        qmat <- crossprod( matrix(y, nrow = n, ncol = nlambda) - mu_cv )

        ####  Quadratic programming
        qrsolve <- try( kernlab::ipop( c = rep(0.0, nlambda),
                                       H = 2.0 * qmat,
                                       A = matrix(1.0, nrow = 1, ncol = nlambda),
                                       b = 1.0,
                                       l = rep(0.0, nlambda),
                                       u = rep(1.0, nlambda),
                                       r = 0.0),
                        TRUE)

        if( inherits(qrsolve, "try-error") == FALSE ){

            ####
            res$lambda <- lambda
            res$cvval <- matrix(qrsolve@primal, nrow = 1 ) %*% qmat %*% matrix(qrsolve@primal, ncol = 1 )
            convergence <- qrsolve@how

            #### SJ: Put the estimates in the right slot, depending on the quadrature method
            if( "raw" %in% rule ){
                res[["raw"]] <- list()
                res[["raw"]]$w <- qrsolve@primal
                res[["raw"]]$convergence <- convergence
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

