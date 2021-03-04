
#' Fit semi-parametric GLM
#'
#'@rdname spglm
#'@param formula a one-sided formula of the form \code{cbind(r, s) ~ predictors} where \code{r} and \code{s} give the number of responses and non-responses within each cluster, respectively (so the cluster size is \code{r+s}), and \code{predictors} describes the covariates.
#'@param data  an optional matrix or data frame containing the variables in the formula \code{formula}. By default the variables are taken from \code{environment(formula).}
#'@param subset an optional vector specifying a subset of observations to be used.
#'@param weight an optional vector specifying observation weights.
#'@param link     a link function for the mean.
#'@param start    an optional list with elements named \code{beta}, \code{q}, and/or \code{mu1} giving starting values for estimation. If none or only some are specified, starting values are selected automatically.
#'@param control        a list with parameters controlling the algorithm.
#'@return an object of class \code{spglm} with the fitted model.
#'@export
#' @importFrom stats terms model.matrix

spglm <- function(formula, data, subset, weights, link="logit", start=NULL, control=list(eps=0.001, maxit=100), ...){
    fam <- binomial(link=link)
    
    
       if (missing(formula) || (length(formula) != 3L))
            stop("'formula' missing or incorrect")
       if (missing(data))
            data <- environment(formula)
        mc <- match.call(expand.dots = FALSE)
        m <- match(c("formula", "data", "subset", "weights"), names(mc), 0L)
        m <- mc[c(1L,m)]
        if (is.matrix(eval(m$data, parent.frame())))
            m$data <- as.data.frame(data)
        m[[1L]] <- quote(stats::model.frame)
        mf <- eval(m, parent.frame())

        # extract and check response variable
        Y <- model.response(mf)
        if (ncol(Y) != 2) stop("The response variable should be a 2-column matrix")

        # create model matrix
        mm <- model.matrix(formula, data=mf)

        # extract weights
        weights <- as.vector(model.weights(mf))
        if (!is.null(weights) && !is.numeric(weights))
            stop("'weights' must be a numeric vector")
        if (!is.null(weights) && any(weights < 0))
            stop("negative weights not allowed")
    
    
     while (iter < control$maxit & difference > control$eps) {
        iter <- iter + 1
        referencef0Pre <- referencef0
        tiltingWsPre <- tiltingWs
        betasPre <- betas
        # passing parameter updates to the next iteration
        if (iter == 1) {
          startUpdates <- list(iter = iter,
                               f0StartValue = NULL,
                               thetaStartValue = NULL,
                               betaStartValue = NULL)
        } else {
          startUpdates <- list(iter = iter,
                               f0StartValue = referencef0Pre,
                               thetaStartValue = tiltingWsPre,
                               betaStartValue = betasPre)
        }
        # EM NR algorithm
        res <- gldrMaximize(DesignMatrix.CONST, CBData, referencef0Pre, tiltingWsPre, startUpdates)
        referencef0 <- res$density.ref 
        # notice that this tilting parameter estimates are w.r.t the intermediate dataset, original dataset extended by max(cluster size)+1 times
        tiltingWs <- res$tiltParam[(0:(nrow(CBData)-1))*length(referencef0)+1]
        betas <- res$regressionEst
        # monitor parameter estimates for baseline distribution, tilting parameters, regression coefficients
        difference <- sum(abs(referencef0Pre - referencef0)) + sum(abs(betasPre - betas))
      }
      
      list(referencef0 = referencef0,
           tiltingWs = tiltingWs,
           betas = betas,
           EMgldrmFitRes = res)
    

    mt <- attr(mf, "terms")
    res <- list(coefficients = beta_new, q=q_new, niter = iter, loglik=logl,
                link = link, call = mc, terms = mt,
                xlevels = .getXlevels(mt, mf),
                data_object=list(model_matrix=mm, resp=Y, n=rowSums(Y), weights=weights),
                model_fun=model_fun)
    class(res) <- "spglm"
    res

}
