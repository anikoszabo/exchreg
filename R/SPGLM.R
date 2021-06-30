
#' Fit semi-parametric GLM
#'
#'@rdname spglm
#'@param formula a one-sided formula of the form \code{cbind(r, s) ~ predictors} where \code{r} and \code{s} give the number of responses and non-responses within each cluster, respectively (so the cluster size is \code{r+s}), and \code{predictors} describes the covariates.
#'@param data  an optional matrix or data frame containing the variables in the formula \code{formula}. By default the variables are taken from \code{environment(formula).}
#'@param subset an optional vector specifying a subset of observations to be used.
#'@param weight an optional vector specifying observation weights.
#'@param link     a link function for the mean.
#'@param mu0      an optional numeric value constraining the mean of the baseline distribution
#'@param control        a list with parameters controlling the algorithm.
#'@return an object of class \code{spglm} with the fitted model.
#'@export
#' @importFrom stats terms model.matrix model.frame model.offset model.response model.weights

spglm <- function(formula, data, subset, weights, offset, link="logit", mu0=NULL, 
                  control=list(eps=0.001, maxit=100), ...){

    
       if (missing(formula) || (length(formula) != 3L))
            stop("'formula' missing or incorrect")
       if (missing(data))
            data <- environment(formula)
        mc <- match.call(expand.dots = FALSE)
        m <- match(c("formula", "data", "subset", "weights", "offset"), names(mc), 0L)
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
        
        if (is.character(link)) {
          link <- stats::make.link(link)
        }
        else if (!is.list(link) || !(all(c("linkfun", "linkinv", 
                                         "mu.eta") %in% names(link)))) {
          stop(paste0("link should be a character string or a list containing ", 
                    "functions named linkfun, linkinv, and mu.eta"))
        }

        # extract offset
        offset <- as.vector(model.offset(mf))
        if (is.null(offset)) 
          offset <- rep(0, nrow(mm))
        
        # extract weights
        weights <- as.vector(model.weights(mf))
        if (!is.null(weights) && !is.numeric(weights))
            stop("'weights' must be a numeric vector")
        if (!is.null(weights) && any(weights < 0))
            stop("negative weights not allowed")
        if (is.null(weights))
            weights <- rep(1, nrow(mm))
       
       # define dimensions     
        N <- max(rowSums(Y))
        nobs <- nrow(mm)
        p <- ncol(mm)
      
    
    
    data_object <- list(model_matrix=mm, resp=Y, n=rowSums(Y), weights=weights, offset=offset, maxN=N)
    
    
     
       betas <- NULL
       pooled <- CorrBin::CBData(data.frame(Trt = "All", NResp = Y[,1], ClusterSize = rowSums(Y), Freq=ceiling(weights)), 
                         trt="Trt", clustersize="ClusterSize", nresp="NResp", freq="Freq")
       est <- CorrBin::mc.est(pooled)
       referencef0 <- est$Prob[est$ClusterSize == N]
       # ensure all positive values
       referencef0 <- (referencef0 + 1e-6)/(1+(N+1)*1e-6)
       
       fTiltMatrix <- matrix(rep(referencef0, times=nobs), byrow=TRUE, nrow=nobs, ncol=N+1)
       spt <- 0:N

       if (is.null(mu0))
         mu0 <- weighted.mean(Y[,1]/rowSums(Y), weights)
       
       # hypergeometric terms for log-likelihood calculation
       hp <- sapply(0:N, function(t)dhyper(x=Y[,1], m=rowSums(Y), n=N-rowSums(Y), k=t))   
       
       # initial log-likelihood
       llik <- log(rowSums(fTiltMatrix * hp)) %*% weights

     
     
       # replace each cluster with N+1 clusters of size N
       new <- cbind(y=0:N)
       rep_idx <- rep(1:nrow(Y), each=nrow(new))
       Y2 <- cbind(1:nrow(Y), Y)[rep_idx,]
       colnames(Y2) <- c("i", "resp","nonresp")

       rep_idx2 <- rep(1:nrow(new), times=nrow(Y))
       Ycomb <- cbind(Y2, new[rep_idx2, ,drop=FALSE])
       # select possible combinations
       possible <- (Ycomb[,"y"] >= Ycomb[,"resp"]) & (N-Ycomb[,"y"] >= Ycomb[,"nonresp"])
       Ycomb <- Ycomb[possible,]
       
       obs_start <- match(1:nobs, Ycomb[,"i"])
      
       mm2 <- mm[Ycomb[,"i"], ,drop=FALSE]
       weights2 <- weights[Ycomb[,"i"]]
       offset2 <- offset[Ycomb[,"i"]]
     
     
     iter <- 0
     difference <- 100
     while (iter < control$maxit & difference > control$eps) {
        iter <- iter + 1
        referencef0Pre <- referencef0
        betasPre <- betas
        llikPre <- llik
        
        
          # convert fTiltMatrix to long vector, watching out for potentially different support
          y_idx <- match(Ycomb[,"y"], spt)
          fTilt2 <- ifelse(!is.na(y_idx), fTiltMatrix[cbind(Ycomb[,"i"], y_idx)], 0)
          
          # numerator of weights
          pp_num <- choose(n=Ycomb[,"y"], k=Ycomb[,"resp"]) * 
                   choose(n=N-Ycomb[,"y"], k=Ycomb[,"nonresp"]) * 
                   fTilt2
          # denominator 
          pp_denom <- tapply(pp_num, list(i=Ycomb[,"i"]), sum, simplify=TRUE)
          pp <- c(pp_num / pp_denom[Ycomb[,"i"]])
          if (!is.null(weights)) pp <- weights2 * pp
        
        
          gldrmControl0 <- gldrm.control(returnfTiltMatrix = TRUE, returnf0ScoreInfo = FALSE, print=FALSE,
                                        betaStart = betasPre, f0Start = referencef0Pre)

          mod <- gldrmFit(x = mm2, y=Ycomb[,"y"]/N, linkfun=link$linkfun, linkinv = link$linkinv,
                          mu.eta = link$mu.eta, mu0 = mu0, offset = offset2, weights = pp, 
                          gldrmControl = gldrmControl0,  thetaControl=theta.control(),
                          betaControl=beta.control(), f0Control=f0.control())
                          
          fTiltMatrix <- mod$fTiltMatrix[obs_start,]
          betas <- mod$beta
          referencef0 <- mod$f0
          spt <- round(mod$spt * N)
          llik <- c(log(rowSums(fTiltMatrix * hp)) %*% weights)
        
        
        difference <- abs(llik - llikPre)
      }
      
      
        ll <- function(x){
          spglm_loglik(beta=x[1:p], f0 = exp(x[-(1:p)]), data_object=data_object, link=link)
        }

        hess <- numDeriv::hessian(func=ll,  x=c(betas, log(referencef0)))
        
        # revert to unlogged f0
        grad <- c(rep(1, p), 1/referencef0)
        hess1 <- diag(grad) %*% hess %*% diag(grad)
        
        # create bordered hessian
        border1 <- c(rep(0, p), rep(1, N+1))  # gradient of sum-to-one constraint
        border2 <- c(rep(0, p), 0:N)          # gradient of fixed-mean constraint
        bhess <- rbind(cbind(hess1, border1, border2), c(border1,0,0), c(border2,0,0))
        
        # calculate variance-covariance matrix
        vc <- solve(-bhess)
        SEbeta <- sqrt(diag(vc)[1:p])
        SEf0 <- sqrt(diag(vc)[p+1+(0:N)])
      
      
    

    mt <- attr(mf, "terms")
    names(betas) <- colnames(mm)
    names(referencef0) <- 0:N
    res <- list(coefficients = betas, SE = SEbeta, f0=referencef0, SEf0 = SEf0, mu0=mu0, niter = iter, 
                loglik=llik, link = link, call = mc, terms = mt,
                xlevels = .getXlevels(mt, mf),
                data_object=data_object)
    class(res) <- "spglm"
    res

}

#' @rdname spglm
#' @export
#' @importFrom stats printCoefmat
print.spglm <- function(x, digits = max(3L, getOption("digits") - 3L),...){
  # based on print.lm, summary.lm, and print.summary.lm
  cat("\nA semi-parametric generalized linear regression model fit\n")
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
  if (length(x$coefficients)) {
    beta <- x$coefficients
    SEbeta <- x$SE
    z <- beta / SEbeta
    pval <- 2 * stats::pnorm(abs(z), lower.tail=FALSE)
    coefmat <-  cbind(Estimate = beta, `Std. Error` = SEbeta, 
                      `z value` = z, `Pr(>|z|)` = pval)
    cat("Coefficients:\n")
    printCoefmat(coefmat, digits = digits, na.print = "NA", ...)
    }
    else cat("No coefficients\n")

    if (length(x$f0)){
        cat("\n Baseline probabilities (q0):\n")        
        q0 <- x$f0
         print.default(format(q0, digits = digits),
            print.gap = 2, quote = FALSE)
    }
    else cat("No baseline joint probabilities\n")

   cat("\n Log-likelihood: ", format(x$loglik, digits=digits), "\n")
   invisible(x)
}

#' Predict methods for SPGLM fits
#' Obtains predictions from a fitted semi-parametric generalized linear model. Note that \code{offset}
#' and \code{weight} terms are not implemented for predicting from new data.
#' @param object fitted model of class \code{spglm}
#' @param newdata optionally, a data frame in which to look for covariate values for prediction. 
#'  If NULL, the original data set is used
#' @param type the type of prediction requested. The default is "mean", the mean event probability;
#' "prob" requests the probability of the observation (given number of responses with given cluster size), 
#' "tilt" returns the tilting parameter which achieves the modeled mean from the baseline distribution;
#' and "lp" requests the linear predictor.
#' @param newn if \code{newdata} is provided and \code{type="prob"}, an integer or integer vector specifying the 
#' @param newevents if \code{newdata} is provided and \code{type="prob"}, an integer or integer vector specifying the 
#'  number of events for the predictions
#' @export
#' @importFrom stats  .checkMFClasses .getXlevels delete.response

predict.spglm <- function(object, newdata=NULL,
                              type=c("mean", "prob", "tilt", "lp"),
                              newn=NULL, newevents=NULL, ...){
  type <- match.arg(type)
  tt <- terms(object)
  if (!missing(newdata)){    
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses")))
            .checkMFClasses(cl, m)
    mm <- model.matrix(Terms, m)
    data_object <- list(model_matrix=mm, offset = rep(0, nrow(mm)), weights=rep(1, nrow(mm)),
                        maxN = object$data_object$maxN)
    if (!missing(newn)){
       if (!(length(newn) == 1L || length(newn) == nrow(newdata)))
          stop("'newn' should have length 1 or equal to the number of rows of 'newdata'")
        data_object$n <- rep(newn, length=nrow(newdata))
        if (!all(data_object$n %in% 0:data_object$maxN))
          stop("Values in 'newn' should be integers between 0 and the maximum cluster size in the original data")
       }
    if (!missing(newevents)){
       if (!(length(newevents) == 1L || length(newevents) == nrow(newdata)))
          stop("'newevents' should have length 1 or equal to the number of rows of 'newdata'")
        data_object$resp <- cbind(rep(newevents, length=nrow(newdata)))
       }
  } else {
    data_object <- object$data_object
  }
  
  if (type=="mean"){
    pred <- spglm_pred_mean(beta=object$coefficients, data_object, link=object$link)
  } else
  if (type=="prob"){
    if (!missing(newdata) && (missing(newn) || missing(newevents)))
       stop("For prediction of probability vectors with new data, cluster sizes should be specified in 'newn' and number of events in 'newevents'.")     
     pred <- spglm_probs(beta=object$coefficients, f0=object$f0, 
                              data_object=data_object, link=object$link)
  } else
  if (type == "lp"){
     pred <- spglm_lp(beta=object$coefficients, data_object=data_object)
  } else
  if (type == "tilt"){
     pred <- spglm_tilt(beta=object$coefficients, f0=object$f0, 
                       data_object=data_object, link=object$link)
  }
  return(pred)    
}

#' @keywords internal
#' @importFrom stats dhyper
spglm_lp <- function(beta, data_object){
  eta <- c(data_object$model_matrix %*% beta + data_object$offset)
  eta
}

spglm_pred_mean <- function(beta, data_object, link){
  eta <- c(data_object$model_matrix %*% beta + data_object$offset)
  mu <- link$linkinv(eta)
  mu
}

spglm_tilt <- function(beta, f0, data_object, link){

    mu <- spglm_pred_mean(beta=beta, data_object=data_object, link=link)
    N <- data_object$maxN
    nobs <- nrow(data_object$model_matrix)
    
    # ySptIndex is only used to calculate the log-likelihood, which we will not be using
    th <- getTheta(spt=(0:N)/N, f0=f0, mu=mu, weights=data_object$weights, ySptIndex=rep(1, nobs),
                   thetaStart=NULL, thetaControl=theta.control())
    th$theta
}

spglm_probs <- function(beta, f0, data_object, link){

    mu <- spglm_pred_mean(beta=beta, data_object=data_object, link=link)
    N <- data_object$maxN
    nobs <- nrow(data_object$model_matrix)
    
    # ySptIndex is only used to calculate the log-likelihood, which we will not be using
    th <- getTheta(spt=(0:N)/N, f0=f0, mu=mu, weights=data_object$weights, ySptIndex=rep(1, nobs),
                   thetaStart=NULL, thetaControl=theta.control())
                   
    hp <- sapply(0:N, function(t)dhyper(x=data_object$resp[,1], m=data_object$n, n=N-data_object$n, k=t))   
    probs <- rowSums(t(th$fTilt) * hp)
    probs
}

spglm_loglik <- function(beta, f0, data_object, link){
  
    probs <- spglm_probs(beta, f0, data_object, link)
    llik <- log(probs) %*% data_object$weights
    c(llik)
}
