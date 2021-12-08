fma.glmnet <-  function( x, y, focus, family, nlambda = 100, nfolds = NULL, grouped,
                         penalty.factor = NULL, force.nlambda, singleton.intercept.only = FALSE,
                         type.measure, rule, solnptol, qp.scale ){

    #### Note: Only jackknife is supported in this version.
    ####       x must not contain the intercept column.
    #### SJ: Define some constants
    n <- nrow(x) #### Sample size
    nbeta <- ncol(x) + 1 #### Number of coefficients, including intercept
    res <- list() # Initiate returned list
    if(is.null(penalty.factor) == TRUE){
        penalty.factor <- rep(1.0, ncol(x))
    }
    if(is.null(nfolds) == TRUE){
        nfolds = n
    }

    ########## Functions needed for weight estimations ##########
    #### SJ: Only useful if focus = "beta"
    if( focus == "eta" ){

        ####  SJ: fun_prob computes the cross-validation MSE
        fun_prob = function(w, Eta, y, family){

            ave.eta <- Eta %*% matrix(w, ncol = 1)
            if( family == "binomial" ){
                mu = plogis(ave.eta)
            } else if( family == "poisson" ){
                mu = exp(ave.eta)
            }
            return(mean((y - mu) ^ 2))

        }

        ####  Sum of the weights
        fun_eq = function(w, Eta, y, family){
            return(sum(w))
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
    #### SJ: Run glmnet() to get the lambda sequence. (Any better way?)
    #### Note: lambda sequence is decreasing.
    fit_alldata = glmnet(x, y, family = family, nfolds = nfolds,
                         grouped = grouped, nlambda = nlambda,
                         type.measure = type.measure,
                         penalty.factor = penalty.factor)
    lambda_raw <- fit_alldata$lambda
    nlambda_raw <- length(lambda_raw)
    if(force.nlambda == TRUE & nlambda_raw != nlambda){
        #### I think glmnet use equal space on the log scale.
        fit_alldata = glmnet(x, y, family = family, nfolds = nfolds,
                             grouped = grouped,
                             lambda = exp(seq(max(log(lambda_raw)), min(log(lambda_raw)), length.out = nlambda)),
                             type.measure = type.measure,
                             penalty.factor = penalty.factor)
        lambda_raw <- fit_alldata$lambda
        nlambda_raw <- length(lambda_raw)
    }
    rm(nlambda)

    #### ---------------------------------------------------- ####
    #### SJ: We need to create the "correct" lambda sequence from the lambda sequence in glmnet()
    ####     This step depends on the quadrature rule that we want to use.
    ####     If we use "raw" (the discrete weight unit-simplex) or "Riemann" for integral, we do not need to introduce new lambda values.
    ####     At this stage, we do not include the full model (lambda = 0) into the lambda sequence
    lambda <- list()
    if("raw" %in% rule){
        lambda[["raw"]] <- lambda_raw
    }
    if("Riemann" %in% rule){
        lambda[["Riemann"]] <- lambda_raw
    }
    if( "Simpson(1/3)" %in% rule ){

        lambda_temp <- lambda_raw
        ####  Add the mid points to the lambda sequence
        l_seq = sort(c(lambda_temp, 0))
        S = 2 : length(l_seq)
        l_mid = c((l_seq[S] + l_seq[S - 1])/2)
        lambda_temp = numeric(length(l_mid) + length(l_seq))
        even = c(seq(from = 2, to = length(lambda_temp) - 2, by = 2), length(lambda_temp) - 1)
        lambda_temp[even] = l_mid
        lambda_temp[-even] = l_seq
        lambda_temp <- lambda_temp[-1] ####  We temporarily remove the full model.
        lambda[["Simpson(1/3)"]] <- sort(lambda_temp, decreasing = TRUE) ####  glmnet() works with decreasing sequence.

    }
    if( "Simpson(3/8)" %in% rule ){

        #### SJ: Not finished yet!!!!
        lambda_temp <- c( 0.0, sort(lambda_raw, decreasing = FALSE) )
        ####  Add the mid points to the lambda sequence
        lambdamat <- cbind( lambda_temp[-length(lambda_temp)], NA, NA, lambda_temp[-1] )
        lambdamat[, 2] <- (2.0 * lambdamat[,1] + lambdamat[,4]) / 3.0
        lambdamat[, 3] <- (lambdamat[,1] + 2.0 * lambdamat[,4]) / 3.0
        lambda_temp <- c(lambda_raw, lambdamat[,2], lambdamat[,3]) ####  We do not need the full model for cv.glmnet()
        lambda[["Simpson(3/8)"]] <- sort(lambda_temp, decreasing = TRUE) ####  glmnet() works with decreasing sequence.

    }

    #### SJ: If Simpson is considered, its lambda sequence always nests "raw".
    ####     So we only need to figure out what the "longest" lambda sequence is.
    if("Simpson(3/8)" %in% rule){
        lambda_long <- lambda[["Simpson(3/8)"]]
    } else if("Simpson(1/3)" %in% rule){
        lambda_long <- lambda[["Simpson(1/3)"]]
    } else if("raw" %in% rule){
        lambda_long <- lambda[["raw"]]
    } else if("Riemann" %in% rule){
        lambda_long <- lambda[["Riemann"]]
    }
    #### ---------------------------------------------------- ####


    #### ---------------------------------------------------- ####
    #### SJ: Using the above lambda sequence, we can start the CV step.
    if( focus == "eta" ) {

        #### VZ: Adding corresponding storage for hybrid and singletons
        if("Hybrid" %in% rule){

            betas_lasso <-  as.matrix(coef(fit_alldata, s = lambda_raw))  # Coefficients for hybrid lambda
            u <-  betas_lasso != 0                    # Logical, TRUE if coefficient is not 0
            hybrid_comb <- unique(u, MARGIN = 2)     # Filtering unique model forms
            #### SJ: Check the full model. If it is already added, we do nothing.
            ####     If not added, we add the full model.
            if(all(colSums(hybrid_comb) != nbeta)){
                #### SJ: Full model is absent, we add the full model
                hybrid_comb <- cbind(hybrid_comb, 1)
            } else {
                #### SJ: If the full model is present, we make sure it is not added more than once
                full.index <- which(colSums(hybrid_comb) == nbeta)
                if(length(full.index) >= 2){
                    hybrid_comb <- hybrid_comb[, -full.index[-1]]
                }
            }
            eta_hybrid = matrix(0, nfolds, ncol(hybrid_comb))

        }
        if("Singleton" %in% rule){
            if(singleton.intercept.only == TRUE){
                eta_sing = matrix(0, nfolds, nbeta)
            } else {
                eta_sing = matrix(0, nfolds, nbeta - 1)
            }
        }

        ####  CV part. This part takes a long time, especially Hybrid and singleton.
        fuleta <- rep(NA, nfolds)
        for(i in 1 : nfolds){

            x_i <- x[-i, , drop = FALSE]
            y_i <- y[-i]

            ####  full model and others.
            if( family == "binomial" ){

                cvfit <- glm.fit( x = cbind(1, x_i), y = y_i, family = binomial() )

                #### VZ: Adding hybrid and singleton CV
                if( "Hybrid" %in% rule ){

                    for( j in 1:ncol(hybrid_comb) ){                      # Looping through model forms

                        sel <- which(hybrid_comb[, j] %in% 1)               # Indices of regressors to use
                        sel_fit <- glm.fit(cbind(1, x_i)[, sel], y_i,
                                           family = binomial())             # Fitting model
                        eta_hybrid[i, j] <- sum(c(1.0, x[i,])[sel] * coef(sel_fit))                     # Extracting coefficients
                    }

                }

                if( "Singleton" %in% rule ){

                    for( l in 1:(nbeta - 1) ){                              # Looping through Singletons

                        sing_fit <- glm.fit(cbind(1, x_i[ , l]), y_i,         # Fitting Singletons
                                            family = binomial())
                        eta_sing[i, l] <- sum(c(1.0, x[i, l]) * coef(sing_fit))                     # Extracting coefficients

                    }
                    if(singleton.intercept.only == TRUE){
                        eta_sing[i, nbeta] = qlogis(mean(y_i))
                    }

                }

            } else if( family == "poisson" ){

                ####  Full model.
                cvfit <- glm.fit( x = cbind(1, x_i), y = y_i, family = poisson() )

                #### VZ: Adding hybrid and singleton CV
                if( "Hybrid" %in% rule ){

                    for( j in 1 : ncol(hybrid_comb) ){                      # Looping through model forms

                        sel <- which(hybrid_comb[, j] %in% 1)               # Indices of regressors to use
                        sel_fit <- glm.fit(cbind(1, x_i)[, sel], y_i,
                                           family = poisson())             # Fitting model
                        eta_hybrid[i, j] <- sum(c(1.0, x[i,])[sel] * coef(sel_fit))                 # Extracting coefficients

                    }

                }

                if( "Singleton" %in% rule ){

                    for( l in 1:(nbeta - 1) ){                              # Looping through Singletons

                        sing_fit <- glm.fit(cbind(1, x_i[ , l]), y_i,         # Fitting Singletons
                                            family = poisson())
                        eta_sing[i, l] <- sum(c(1.0, x[i, l]) * coef(sing_fit))             # Extracting coefficients

                    }
                    if(singleton.intercept.only == TRUE){
                        eta_sing[i, nbeta] = log(mean(y_i))
                    }

                }

            }

            fuleta[i] <- sum(c(1.0, x[i,]) * cvfit$coefficients)

        }
        ####  SJ: The essence of averaging beta is just to average the linear predictor.
        ####      Hence, we can use cv.glmnet to help us to get the CV linear predictor.
        cvglm = try(cv.glmnet(x, y, family = family, nfolds = nfolds,
                              grouped = grouped, keep = TRUE, lambda = lambda_long,
                              type.measure = type.measure,
                              penalty.factor = penalty.factor), TRUE)
        glmnet.eta <- cbind(cvglm$fit.preval, fuleta)

        #### VZ: Adds cvglmnet object to output
        res[["cvglmnet"]] <- cvglm

        #### SJ: In order to estimate the weights, we add 0 to the lambda sequence here.
        ####     Need to make sure that lambda remains a decreasing sequence!
        if("raw" %in% rule){

            lambda_index <- match(lambda[["raw"]], lambda_long)
            eta.raw <- cbind(glmnet.eta[, lambda_index], fuleta)
            coef_raw = cbind(as.matrix(coef(cvglm, s = lambda[["raw"]])),
                                       fulcoef)
            lambda[["raw"]] <- c(lambda[["raw"]], 0)
            init_w = rep( 1 / length(lambda[["raw"]]), length(lambda[["raw"]]) )
            nlopt_raw <- try( Rsolnp::solnp(pars = init_w, fun = fun_prob, eqfun = fun_eq, eqB = c(1),
                                            LB = rep(0, length(lambda[["raw"]])),
                                            UB = rep(1, length(lambda[["raw"]])),
                                            Eta = eta.raw, y = y, family = family,
                                            control = list(trace = FALSE, tol = solnptol)),
                              TRUE )
            res[["raw"]] <- list()
            if( inherits(nlopt_raw, "try-error") == FALSE ){

                res[["raw"]]$w <- nlopt_raw$pars
                res[["raw"]]$convergence <- nlopt_raw$convergence #### 0 = converged
                res[["raw"]]$avecoef <- coef_raw %*% matrix(nlopt_raw$pars, ncol = 1)
                res[["raw"]]$elapsed <- nlopt_raw$elapsed
                res[["raw"]]$lambda <- lambda[["raw"]]

            }

        }
        if("Riemann" %in% rule){
            ####  SJ: Not finished yet
            lambda[["Riemann"]] <- c(lambda[["Riemann"]], 0)
        }
        if("Simpson(1/3)" %in% rule){

            lambda_index <- match(lambda[["Simpson(1/3)"]], lambda_long)
            eta.simp <- cbind(glmnet.eta[, lambda_index], fuleta)
            coef_simp = cbind(as.matrix(coef(cvglm, s = lambda[["Simpson(1/3)"]])),
                             fulcoef)
            lambda[["Simpson(1/3)"]] <- c(lambda[["Simpson(1/3)"]], 0)
            init_w = rep( 1 / length(lambda[["Simpson(1/3)"]]), length(lambda[["Simpson(1/3)"]]) )
            nlopt_simp <- try( Rsolnp::solnp(pars = init_w, fun = fun_prob, eqfun = fun_eq, eqB = c(1),
                                             LB = rep(0, length(lambda[["Simpson(1/3)"]])),
                                             UB = rep(1, length(lambda[["Simpson(1/3)"]])),
                                             Eta = eta.simp, y = y, family = family,
                                             control = list(trace = FALSE, tol = solnptol)),
                               TRUE )
            res[["Simpson(1/3)"]] <- list()
            if( inherits(nlopt_simp, "try-error") == FALSE ){

                res[["Simpson(1/3)"]]$v <- nlopt_simp$pars
                res[["Simpson(1/3)"]]$convergence <- nlopt_simp$convergence #### 0 = converged
                res[["Simpson(1/3)"]]$avecoef <- coef_simp %*% matrix(nlopt_simp$pars, ncol = 1)
                res[["Simpson(1/3)"]]$elapsed <- nlopt_simp$elapsed
                lambda_raw_zero <- c(lambda_raw, 0.0)
                scale <- (lambda_raw_zero[1] - lambda_raw_zero[2]) / 6
                for( i in 1 : nlambda_raw ){
                    scale <- c(scale,
                               2 * (lambda_raw_zero[i] - lambda_raw_zero[i + 1]) / 3,
                               (lambda_raw_zero[i] - lambda_raw_zero[i + 2]) / 6 )
                }
                scale[length(scale)] <- (lambda_raw_zero[i] - lambda_raw_zero[i + 1]) / 6
                res[["Simpson(1/3)"]]$w <- nlopt_simp$pars / scale
                res[["Simpson(1/3)"]]$lambda <- lambda[["Simpson(1/3)"]]

            }

        }
        if("Simpson(3/8)" %in% rule){
            ####  SJ: Not finished yet
            lambda[["Simpson(3/8)"]] <- c(lambda[["Simpson(3/8)"]], 0)
        }
        if("Singleton" %in% rule){

            if(singleton.intercept.only == TRUE){
                init_w_sing = matrix( 1 / nbeta, nrow = nbeta )
                nlopt_sing <- try( Rsolnp::solnp(pars = init_w_sing, fun = fun_prob, eqfun = fun_eq,
                                                 eqB = c(1),
                                                 LB = rep(0, nbeta), UB = rep(1, nbeta),
                                                 Eta = eta_sing, y = y, family = family,
                                                 control = list(trace = FALSE, tol = solnptol)),
                                   TRUE)
            } else {
                init_w_sing = matrix( 1 / (nbeta - 1), nrow = (nbeta - 1) )
                nlopt_sing <- try( Rsolnp::solnp(pars = init_w_sing, fun = fun_prob, eqfun = fun_eq,
                                                 eqB = c(1),
                                                 LB = rep(0, (nbeta - 1)), UB = rep(1, (nbeta - 1)),
                                                 Eta = eta_sing, y = y, family = family,
                                                 control = list(trace = FALSE, tol = solnptol)),
                                   TRUE)
            }

            res[["Singleton"]] <- list()
            if( inherits(nlopt_sing, "try-error") == FALSE){

                res[["Singleton"]]$w <- nlopt_sing$pars
                res[["Singleton"]]$convergence <- nlopt_sing$convergence
                res[["Singleton"]]$elapsed <- nlopt_sing$elapsed
                coef_sing <- matrix(0, nbeta, nbeta - 1)
                for( l in 1 : (nbeta - 1) ){
                    if(family == "binomial"){
                        sing_fit <- glm.fit(cbind(1, x[ , l]), y, family = binomial())
                    } else if(family == "poisson"){
                        sing_fit <- glm.fit(cbind(1, x[ , l]), y, family = poisson())
                    }
                    coef_sing[c(1, l + 1), l] <-  coef(sing_fit)
                }
                if(singleton.intercept.only == TRUE){
                    coef_sing = cbind(coef_sing, rep(0, nbeta))
                    coef_sing[1, nbeta] = log(mean(y))
                }
                res[["Singleton"]]$avecoef <- coef_sing %*% matrix(nlopt_sing$pars, ncol = 1)
            }

        }
        if("Hybrid" %in% rule){

            init_w_hybrid = matrix( 1 / ncol(hybrid_comb), nrow = ncol(hybrid_comb) )
            nlopt_hybrid <- try( Rsolnp::solnp(pars = init_w_hybrid, fun = fun_prob,
                                               eqfun = fun_eq, eqB = c(1),
                                               LB = rep(0, ncol(hybrid_comb)),
                                               UB = rep(1, ncol(hybrid_comb)),
                                               Eta = eta_hybrid, y = y, family = family,
                                               control = list(trace = FALSE, tol = solnptol)),
                                 TRUE)

            res[["Hybrid"]] <- list()
            if( inherits(nlopt_hybrid, "try-error") == FALSE){

                res[["Hybrid"]]$w <- nlopt_hybrid$pars
                res[["Hybrid"]]$convergence <- nlopt_hybrid$convergence
                res[["Hybrid"]]$elapsed <- nlopt_hybrid$elapsed
                res[["Hybrid"]]$models <- hybrid_comb
                coef_hybrid <- matrix(0, nbeta, ncol(hybrid_comb))
                for(j in 1 : ncol(hybrid_comb)){

                    sel <- which(hybrid_comb[, j] %in% 1)
                    if(family == "binomial"){
                        sel_fit <- glm.fit(cbind(1, x)[, sel], y, family = binomial())
                    } else if(family == "poisson"){
                        sel_fit <- glm.fit(cbind(1, x)[, sel], y, family = poisson())
                    }
                    coef_hybrid[sel, j] <- coef(sel_fit)

                }
                res[["Hybrid"]]$avecoef <- coef_hybrid %*% matrix(nlopt_hybrid$pars, ncol = 1)

            }

        }

    } else if( focus == "mu" ){

        # Fitting the LASSO model
        cvglm = try(cv.glmnet(x, y, family = family, nfolds = nfolds,
                              grouped = grouped, keep = TRUE, lambda = lambda_long,
                              type.measure = type.measure,
                              penalty.factor = penalty.factor), TRUE)
        #### VZ: Adds cvglmnet object to output
        res[["cvglmnet"]] <- cvglm

        #### VZ: Finding unique model forms
        betas_lasso = as.matrix(coef(cvglm, s = lambda_raw))  # Coefficients for hybrid lambda
        u = betas_lasso != 0                    # Logical, TRUE if coefficient is not 0
        hybrid_comb = unique(u, MARGIN = 2)     # Filtering unique model forms
        #### SJ: Check the full model. If it is already added, we do nothing.
        ####     If not added, we add the full model.
        if(all(colSums(hybrid_comb) != nbeta)){
            #### SJ: Full model is absent, we add the full model
            hybrid_comb <- cbind(hybrid_comb, 1)
        } else {
            #### SJ: If the full model is present, we make sure it is not added more than once
            full.index <- which(colSums(hybrid_comb) == nbeta)
            if(length(full.index) >= 2){
                hybrid_comb <- hybrid_comb[, -full.index[-1]]
            }
        }

        ####  SJ: fit full model (returning eta)
        ####  VZ: Fit hybrid and Singletons
        fulmod_eta <- rep( NA, nrow(x) )
        hybrid_eta <- matrix(NA, nrow = nrow(x), ncol = ncol(hybrid_comb) )
        if(singleton.intercept.only == TRUE){
            singleton_eta <- matrix(0, nrow = nrow(x), ncol = nbeta )
        } else {
            singleton_eta <- matrix(0, nrow = nrow(x), ncol = (nbeta - 1) )
        }

        ####  SJ: It is possible to write the following loop in Rcpp to gain speed.
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
                if( "Singleton" %in% rule ){
                    for( l in 1:(nbeta - 1) ){                              # Looping through Singletons
                        sing_fit <- glm.fit(cbind(1, x_i[ , l]), y_i,          # Fitting Singletons
                                            family = binomial())
                        b_sing <-  matrix(coef(sing_fit))                     # Extracting coefficients
                        singleton_eta[i, l] <- cbind(1, x[i, l]) %*% b_sing   # Predicting
                    }
                    if(singleton.intercept.only == TRUE){
                        singleton_eta[i, nbeta] = qlogis(mean(y_i))
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
                if( "Singleton" %in% rule ){

                    for( l in 1:(nbeta - 1) ){                              # Looping through Singletons
                        sing_fit <- glm.fit(cbind(1, x_i[ , l]), y_i,         # Fitting Singletons
                                            family = poisson())
                        b_sing <-  matrix(coef(sing_fit))                     # Extracting coefficients
                        singleton_eta[i, l] <- cbind(1, x[i, l]) %*% b_sing   # Predicting

                    }
                    if(singleton.intercept.only == TRUE){
                        singleton_eta[i, nbeta] = log(mean(y_i))
                    }

                }
            }
            fulmod_eta[i] <- sum(cvfit$coefficients * c(1.0, x[i,]) )

        }

        #### SJ: need to make sure fit.preval is the linear predictor scale
        ####     which is true in version 4.0-2
        all.coef <- as.matrix(coef(cvglm, s = lambda_long))
        if("raw" %in% rule){

            lambda_index <- match(lambda[["raw"]], lambda_long)
            lambda[["raw"]] <- c(lambda[["raw"]], 0)
            nlambda_raw_0 <- length(lambda[["raw"]])
            cv_eta_raw <- cbind(cvglm$fit.preval[, lambda_index],
                                fulmod_eta)
            if( family == "binomial" ){
                mu_cv_raw = plogis( cv_eta_raw )
            } else if( family == "poisson" ){
                mu_cv_raw = exp( cv_eta_raw )
            }

            #### The symmetric matrix in the quadratic programming
            qmat_raw <- crossprod( matrix(y, nrow = n, ncol = nlambda_raw_0) - mu_cv_raw )
            if(qp.scale == "max"){
                qmat_raw <- qmat_raw / max(qmat_raw)
            } else if(qp.scale == "1/n"){
                qmat_raw <- qmat_raw / n
            }

            ####  Quadratic programming
            qrsolve_raw <- try( kernlab::ipop( c = rep(0.0, nlambda_raw_0),
                                               H = qmat_raw,
                                               A = matrix(1.0, nrow = 1, ncol = nlambda_raw_0),
                                               b = 1.0,
                                               l = rep(0.0, nlambda_raw_0),
                                               u = rep(1.0, nlambda_raw_0),
                                               r = 0.0),
                                TRUE)

            if(inherits(qrsolve_raw, "try-error") == FALSE){
                res[["raw"]] <- list()
                res[["raw"]]$w <- qrsolve_raw@primal
                res[["raw"]]$convergence <- qrsolve_raw@how
                res[["raw"]]$pos_def <- all(eigen(qmat_raw)$values > 0)
                res[["raw"]]$lambda <- lambda[["raw"]]
                res[["raw"]]$coef <- all.coef[, lambda_index]
            }

        }
        if("Riemann" %in% rule){

        }
        if("Simpson(1/3)" %in% rule){

            lambda_index <- match(lambda[["Simpson(1/3)"]], lambda_long)
            lambda[["Simpson(1/3)"]] <- c(lambda[["Simpson(1/3)"]], 0)
            nlambda_simp_0 <- length(lambda[["Simpson(1/3)"]])
            cv_eta_simp <- cbind(cvglm$fit.preval[, lambda_index],
                                fulmod_eta)
            if( family == "binomial" ){
                mu_cv_simp = plogis( cv_eta_simp )
            } else if( family == "poisson" ){
                mu_cv_simp = exp( cv_eta_simp )
            }

            #### The symmetric matrix in the quadratic programming
            qmat_simp <- crossprod( matrix(y, nrow = n, ncol = nlambda_simp_0) - mu_cv_simp )
            if(qp.scale == "max"){
                qmat_simp <- qmat_simp / max(qmat_simp)
            } else if(qp.scale == "1/n"){
                qmat_simp <- qmat_simp / n
            }

            ####  Quadratic programming
            qrsolve_simp <- try( kernlab::ipop( c = rep(0.0, nlambda_simp_0),
                                               H = qmat_simp,
                                               A = matrix(1.0, nrow = 1, ncol = nlambda_simp_0),
                                               b = 1.0,
                                               l = rep(0.0, nlambda_simp_0),
                                               u = rep(1.0, nlambda_simp_0),
                                               r = 0.0),
                                TRUE)

            res[["Simpson(1/3)"]] <- list()
            if(inherits(qrsolve_simp, "try-error") == FALSE){

                res[["Simpson(1/3)"]]$v <- qrsolve_simp@primal
                lambda_raw_zero <- c(lambda_raw, 0.0)
                scale <- (lambda_raw_zero[1] - lambda_raw_zero[2]) / 6
                for( i in 1 : nlambda_raw ){
                    scale <- c(scale,
                               2 * (lambda_raw_zero[i] - lambda_raw_zero[i + 1]) / 3,
                               (lambda_raw_zero[i] - lambda_raw_zero[i + 2]) / 6 )
                }
                scale[length(scale)] <- (lambda_raw_zero[i] - lambda_raw_zero[i + 1]) / 6
                res[["Simpson(1/3)"]]$w <- qrsolve_simp@primal / scale
                res[["Simpson(1/3)"]]$convergence <- qrsolve_simp@how
                res[["Simpson(1/3)"]]$pos_def <- all(eigen(qmat_simp)$values > 0)
                res[["Simpson(1/3)"]]$lambda <- lambda[["Simpson(1/3)"]]
                res[["Simpson(1/3)"]]$coef <- all.coef[, lambda_index]

            }

        }
        #### VZ: Adding QP for hybrid averaging
        if( "Hybrid" %in% rule ){

            if( family == "binomial" ){
                mu_cv_hybrid = plogis( hybrid_eta )
            } else if( family == "poisson" ){
                mu_cv_hybrid = exp( hybrid_eta )
            }

            # Symmetric matrix for QP
            qmat_hybrid <- crossprod( matrix(y, nrow = n, ncol = ncol(hybrid_eta)) - mu_cv_hybrid)
            if(qp.scale == "max"){
                qmat_hybrid <- qmat_hybrid / max(qmat_hybrid)
            } else if(qp.scale == "1/n"){
                qmat_hybrid <- qmat_hybrid / n
            }

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
                res[["Hybrid"]]$coef <- matrix(0, nbeta, ncol(hybrid_comb))
                for( j in 1 : ncol(hybrid_comb) ){
                    sel <-  which(hybrid_comb[, j] %in% 1)
                    if(family == "binomial"){
                        sel_fit <-  glm.fit(cbind(1, x)[, sel], y, family = binomial())
                    } else if(family == "poisson"){
                        sel_fit <-  glm.fit(cbind(1, x)[, sel], y, family = poisson())
                    }
                    res[["Hybrid"]]$coef[sel, j] <- coef(sel_fit)
                }

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
            qmat_singleton <- crossprod( matrix(y, nrow = n, ncol = ncol(singleton_eta)) - mu_cv_singleton)
            if(qp.scale == "max"){
                qmat_singleton <- qmat_singleton / max(qmat_singleton)
            } else if(qp.scale == "1/n"){
                qmat_singleton <- qmat_singleton / n
            }

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
                res[["Singleton"]]$coef <- matrix(0, nbeta, nbeta - 1)
                for( l in 1 : (nbeta - 1) ){
                    if(family == "binomial"){
                        sing_fit <- glm.fit(cbind(1, x[ , l]), y, family = binomial())
                    } else if(family == "poisson"){
                        sing_fit <- glm.fit(cbind(1, x[ , l]), y, family = poisson())
                    }
                    res[["Singleton"]]$coef[c(1, l + 1), l] <- coef(sing_fit) ##  add +1 for the intercept.
                }
                if(singleton.intercept.only == TRUE){

                    res[["Singleton"]]$coef <- cbind(res[["Singleton"]]$coef,
                                                     matrix(0, nbeta, 1))
                    if(family == "binomial"){
                        res[["Singleton"]]$coef[1, nbeta] = qlogis(mean(y))
                    } else if(family == "poisson"){
                        res[["Singleton"]]$coef[1, nbeta] = log(mean(y))
                    }

                }

            }

        }

    }
    #### ---------------------------------------------------- ####

    res

}


