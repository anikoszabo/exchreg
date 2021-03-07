
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
#' @importFrom stats terms model.matrix

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
    
    
     
       N <- max(rowSums(Y))
       nobs <- nrow(mm)
       p <- ncol(mm)
       
       betas <- rep(0, times= p)

       pooled <- CBData(data.frame(Trt = "All", NResp = Y[,1], ClusterSize = rowSums(Y), Freq=ceiling(weights)), 
                         trt="Trt", clustersize="ClusterSize", nresp="NResp")
       est <- mc.est(pooled)
       referencef0 <- est$Prob[est$ClusterSize == N]
       # ensure all positive values
       referencef0 <- (referencef0 + 1e-6)/(1+(N+1)*1e-6)
       
       fTiltMatrix <- matrix(rep(referencef0, times=nobs), byrow=TRUE, nrow=nobs, ncol=N+1)
       spt <- 0:N

       if (is.null(mu0))
         mu0 <- weighted.mean(Y[,1]/rowSums(Y), weights)
       
     
     
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
      
       mm2 <- mm[Ycomb[,"i"], ,drop=FALSE]
       weights2 <- weights[Ycomb[,"i"]]
       offset2 <- offset[Ycomb[,"i"]]
     
     
     iter <- 0
     difference <- 100
     while (iter < control$maxit & difference > control$eps) {
        iter <- iter + 1
        referencef0Pre <- referencef0
        betasPre <- betas
        
        
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
                          
          fTiltMatrix <- mod$fTiltMatrix
          betas <- mod$beta
          referencef0 <- mod$f0
          spt <- round(mod$spt * N)
        
        
        difference <- sum(abs(referencef0Pre - referencef0)) + sum(abs(betasPre - betas))
      }
      
    

    mt <- attr(mf, "terms")
    res <- list(coefficients = betas, f0=referencef0, mu0=mu0, niter = iter, loglik=mod$llik,
                link = link, call = mc, terms = mt,
                xlevels = .getXlevels(mt, mf),
                data_object=list(model_matrix=mm, resp=Y, n=rowSums(Y), weights=weights, offset=offset))
    class(res) <- "spglm"
    res

}

spglm_pred_mean <- function(beta, data_object, link){
  eta <- c(data_object$model_matrix %*% beta + data_object$offset)
  mu <- link$linkinv(eta)
  mu
}

spglm_loglik <- function(beta, f0, data_object, link){
  
    mu <- spglm_pred_mean(beta=beta, data_object=data_object, link=link)
    N <- max(data_object$n)
    nobs <- nrow(data_object$model_matrix)
    
    # ySptIndex is only used to calculate the log-likelihood, which we will not be using
    th <- getTheta(spt=(0:N)/N, f0=f0, mu=mu, weights=data_object$weights, ySptIndex=rep(1, nobs),
                   thetaStart=NULL, thetaControl=theta.control())
    
    hp <- sapply(0:N, function(t)dhyper(x=data_object$resp[,1], m=data_object$n, n=N-data_object$n, k=t))   
    llik_term <- rowSums(t(th$fTilt) * hp)
    llik <- llik_term %*% data_object$weights
    c(llik)
}
