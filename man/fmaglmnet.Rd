\name{fma.glmnet}
\alias{fma.glmnet}
\title{Frequentist model averaging of adaptive LASSO in GLM}
\usage{
fma.glmnet( x, y, focus, family, nlambda = 100, nfolds = NULL, grouped,
            penalty.factor = NULL, force.nlambda,
            singleton.intercept.only = FALSE,
            type.measure = "mse", rule, solnptol, qp.scale,
            glmnetfit = NULL)
}

\arguments{
\item{x}{the matrix of regressors, excluding the intercept.}

\item{y}{the vector of response variable.}

\item{focus}{"mu" or "eta".}

\item{family}{"binomial" or "poisson".}

\item{nlambda}{the number of lambda. Default is 100.}

\item{nfolds}{the number of folds in cross validation. Default equals the sample size.}

\item{grouped}{same as grouped in the glmnet package.}

\item{penalty.factor}{It will be passed to penalty.factor in the glmnet package. Default is a vector of 1, indicating same scaling. If penalty.factor is a numerical vector, it will be directly passed to cv.glment(). If penalty.factor="adaptiveLasso", the weight will be the reciporacal of the absolute value of the full model estimates.}

\item{force.nlambda}{it is common that the user specify a value of nlambda, but the actual number of lambda is lower than nlambda. If force.nlambda = TRUE, the number of lambda values is forced to be nlambda by forming a better lambda sequence.}

\item{singleton.intercept.only}{if TRUE, then the interpcet only model is also included in the singleton models. FALSE means that the intercept only model is left out.}

\item{type.measure}{same as type.measure in the glmnet package.}

\item{rule}{a vector that contains entries from "IMSE", "Singleton", and "Hybrid". "IMSE" averages the lambda sequence produced from the glmnet package, which aims to approximate the integrated MSE.}

\item{solnptol}{tolerance used in solnp().}

\item{qp.scale}{"max" or "1/n". "max" means that the raw symmetric matrix in QP is divided by the max number of the matrix. "1/n" means that it is divided by the sample size.}

\item{glmnetfit}{The user can fit a glmnet object themselves and pass it to the function to save time. It is not provided, a glmnet object will be fitted in the function.}
}

\value{
It returns a list of the following objects
\itemize{
\item{fullcoef:} estimated beta from the full model.
\item{penalty.factor:} the weights for each covariate.
\item{IMSE:} if "IMSE" belongs to rule, then it contains a list of the following objects.
\itemize{
  \item{w:} estimated weights.
  \item{convergence:} 0 or converged = converged.
  \item{avecoef:} averaged beta coefficients, if focus = "eta".
  \item{elapsed:} elapsed time in seconds, if focus = "eta".
  \item{lambda:} the lambda sequence.
  \item{pos_def:} whether QP is positive definite, if focus = "mu".
  \item{coef:} the matrix containing all estimated beta (excluding the full model), if focus = "mu".
}
\item{Hybrid:} if "Hybrid" belongs to rule, then it contains a list of the following objects.
\itemize{
  \item{w:} estimated weights.
  \item{convergence:} 0 or converged = converged.
  \item{avecoef:} averaged beta coefficients, if focus = "eta".
  \item{elapsed:} elapsed time in seconds, if focus = "eta".
  \item{models:} unique model structures.
  \item{pos_def:} whether QP is positive definite, if focus = "mu".
  \item{coef:} the matrix containing all estimated beta, if focus = "mu".
}
\item{Singleton:} if "Singleton" belongs to rule, then it contains a list of the following objects.
\itemize{
  \item{w:} estimated weights.
  \item{convergence:} 0 or converged = converged.
  \item{avecoef:} averaged beta coefficients, if focus = "eta".
  \item{elapsed:} elapsed time in seconds, if focus = "eta".
  \item{pos_def:} whether QP is positive definite, if focus = "mu".
  \item{coef:} the matrix containing all estimated beta, if focus = "mu".
}
}
}

\description{
Cross-validation frequentist model averaging for generalized linear models. The candidate models are selected by adaptive LASSO.
}

\examples{
data("ToyData_Binomial")
fma.glmnet(x = ToyBinomial[1 : 100, paste("x", 1 : 8, sep = "")],
           y = ToyBinomial[1 : 100, "y"],
           focus = "mu",
           family = "binomial",
           nlambda = 50,
           nfolds = 10,
           penalty.factor = "adaptiveLasso",
           grouped = FALSE,
           type.measure = "class",
           rule = c("IMSE", "Hybrid", "Singleton"),
           solnptol = 1e-08,
           force.nlambda = TRUE,
           qp.scale = "max",
           singleton.intercept.only = TRUE)

fma.glmnet(x = ToyBinomial[1 : 100, paste("x", 1 : 8, sep = "")],
           y = ToyBinomial[1 : 100, "y"],
           focus = "eta",
           family = "binomial",
           nlambda = 50,
           nfolds = 10,
           penalty.factor = "adaptiveLasso",
           grouped = FALSE,
           type.measure = "class",
           rule = c("IMSE", "Hybrid", "Singleton"),
           solnptol = 1e-08,
           force.nlambda = TRUE,
           qp.scale = "max",
           singleton.intercept.only = TRUE)

data("ToyData_Poisson")
fma.glmnet(x = ToyPoisson[1 : 100, paste("x", 1 : 8, sep = "")],
           y = ToyPoisson[1 : 100, "y"],
           focus = "mu",
           family = "poisson",
           nlambda = 50,
           nfolds = 10,
           penalty.factor = "adaptiveLasso",
           grouped = FALSE,
           type.measure = "mse",
           rule = c("IMSE", "Hybrid", "Singleton"),
           solnptol = 1e-08,
           force.nlambda = TRUE,
           qp.scale = "max",
           singleton.intercept.only = TRUE)

fma.glmnet(x = ToyPoisson[1 : 100, paste("x", 1 : 8, sep = "")],
           y = ToyPoisson[1 : 100, "y"],
           focus = "eta",
           family = "poisson",
           nlambda = 50,
           nfolds = 10,
           penalty.factor = "adaptiveLasso",
           grouped = FALSE,
           type.measure = "mse",
           rule = c("IMSE", "Hybrid", "Singleton"),
           solnptol = 1e-08,
           force.nlambda = TRUE,
           qp.scale = "max",
           singleton.intercept.only = TRUE)

}
