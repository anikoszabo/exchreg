
require(CorrBin)

# lambda-to-p weight matrix, rows:r, cols:s
# (-1)^(s-r)*choose(n,r)*choose(n-r,s-r)

weight_mat <- function(n){
  s <- r <- 0:n
  res <- outer(r, s, function(x,y)(-1)^(y-x)*choose(n,x)*choose(n-x,y-x))
  res
}

p_from_lambda <- function(lambda.vec, n){
  H <- weight_mat(n)
  c(H %*% lambda.vec[1:(n+1)])
}

# p-to-lambda weight matrix, rows:k, cols:r
# (choose(n-k, n-r)/choose(n,r)

weight_mat2 <- function(n){
  k <- r <- 0:n
  res <- outer(k, r, function(x,y)choose(n-x, n-y)/choose(n,y))
  res
}

lambda_from_p <- function(p.vec, n=length(p.vec)-1){
  H <- weight_mat2(n)
  c(H %*% p.vec)
}

#'@rdname sprr
#'@param formula a one-sided formula of the form \code{cbind(r, s) ~ predictors} where \code{r} and \code{s} give the number of responses and non-responses within each cluster, respectively (so the cluster size is \code{r+s}), and \code{predictors} describes the covariates.
#'@param data  an optional matrix or data frame containing the variables in the formula \code{formula}. By default the variables are taken from \code{environment(formula).}
#'@param subset an optional vector specifying a subset of observations to be used.
#'@param weight an optional vector specifying observation weights.
#'@param link     a link function for the binomial distribution, the inverse of which models the covariate effects.
#'@param mu1      an optional value between 0 and 1 giving the maximal predicted marginal probability. If set to NULL (default), the algorithm will try to estimate it from the data.
#'@param start    an optional list with elements named \code{beta}, \code{q}, and/or \code{mu1} giving starting values for estimation. If none or only some are specified, starting values are selected automatically.
#'@param control        a list with parameters controlling the algorithm.
#'@return an object of class \code{sprr} with the fitted model.
#'@export
#' @importFrom stats terms model.matrix

#' Semi-parametric relative risk model
sprr <- function(formula, data, subset, weights, link="cloglog", mu1=NULL, start=NULL, control=list(eps=0.001, maxit=100), ...){
    fam <- binomial(link=link)
    model_fun <- fam$linkinv

    
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
    
    
        
        mean_constrained_probs <- function(cc, m=NULL){
            if (is.null(m)) return(cc/sum(cc))          
            N <- length(cc) - 1
            r_eq <- Vectorize(function(rho){i <- 0:N; sum(cc * (i - m)/(1+(exp(rho)-1/N)*i))})
            rho <- uniroot(r_eq, interval=c(-10, 10), extendInt="yes")
            r <- exp(rho$root) - 1/N
            q <- cc / (1 + r * (0:N))
            q <- q/sum(q)
            q
        }
        
        
         N <- max(rowSums(Y))

         if (is.null(start$beta)){   
            p0 <- (Y[,1] + 0.5)/(Y[,1] + Y[,2]+ 1)
            if (is.null(mu1) && is.null(start$mu1))
               lp0 <- fam$linkfun(pmin(1-.Machine$double.eps, p0))
            else if (!is.null(mu1))
               lp0 <- fam$linkfun(pmin(1-.Machine$double.eps, p0/mu1))
            else if (!is.null(start$mu1))
               lp0 <- fam$linkfun(pmin(1-.Machine$double.eps, p0/start$mu1))
            if (is.null(weights))       
                lm0 <- lm.fit(x=mm, y=lp0)
            else         
                lm0 <- lm.wfit(x=mm, y=lp0, w=weights)
            beta_new <- coef(lm0)
         } else {
            beta_new <- start$beta
         }
         
         if (is.null(start$q)){
            pooled <- CBData(data.frame(Trt = "All", NResp = Y[,1], ClusterSize = rowSums(Y)), trt="Trt",
                             clustersize="ClusterSize", nresp="NResp")
            est <- mc.est(pooled)
            q0 <- est$Prob[est$ClusterSize == N]
            if (is.null(mu1) && is.null(start$mu1))
                q_new <- mean_constrained_probs(pmax(q0, 0.01))
             else if (!is.null(mu1))
                q_new <- mean_constrained_probs(pmax(q0, 0.01), N*mu1)
             else if (!is.null(start$mu1))
                q_new <- mean_constrained_probs(pmax(q0, 0.01), N*start$mu1)
         } else {
           q_new <- start$q
           if (!is.null(start$mu1))
              q_new <- mean_constrained_probs(q_new, N*start$mu1)
          }
        

        
            # replace each cluster with N(N-1)/2 clusters of size N
            new <- cbind(s=rep(0:N, each=N+1), y=rep(0:N, times=N+1))
            new <- new[new[,"s"] <= new[,"y"],]

            rep_idx <- rep(1:nrow(Y), each=nrow(new))
            Y2 <- cbind(1:nrow(Y), Y)[rep_idx,]
            colnames(Y2) <- c("i", "resp","nonresp")
            mm2 <- mm[rep_idx,,drop=FALSE]

            rep_idx2 <- rep(1:nrow(new), times=nrow(Y))
            Ycomb <- cbind(Y2, new[rep_idx2,])
        
        iter <- 0
        diff <- 100
        while ((diff > control$eps) & (iter < control$maxit)){
            iter <- iter + 1
            beta_old <- beta_new
            q_old <- q_new

            
                # calculate a_{iys}^k
                lp <- mm2 %*% beta_old
                theta <- model_fun(lp)
                a <- dbinom(x=Ycomb[,"s"], size=Ycomb[,"y"], prob=theta) * q_old[Ycomb[,"y"]+1]

                # calculate e_{is}^k w_{iys}^k
                ew_num <- dhyper(x=Ycomb[,"resp"], m=Ycomb[,"s"], n=N - Ycomb[,"s"],
                                 k=Ycomb[,"resp"] + Ycomb[,"nonresp"]) * a
                ew_denom <- tapply(ew_num, list(i=Ycomb[,"i"]), sum, simplify=TRUE)
                ew <-  ew_num / ew_denom[Ycomb[,"i"]]
                if (!is.null(weights)) ew <- weights * ew
            
            
                mod <- glm(cbind(Ycomb[,"s"], Ycomb[,"y"]-Ycomb[,"s"]) ~ mm2+0, 
                           family=fam, weights=ew, start=beta_old)
                beta_new <- coef(mod)
            
            
                c_vec <- tapply(ew, list(y=Ycomb[,"y"]), sum, simplify=TRUE)
                if (is.null(mu1))
                   q_new <- mean_constrained_probs(c_vec)
                else
                   q_new <- mean_constrained_probs(c_vec, N*mu1)
            

            diff <- sum(abs(beta_old - beta_new)) + sum(abs(q_old - q_new))
        }
        logl <- loglik(beta_new, q_new, list(model_matrix=mm, resp=Y, n=rowSums(Y), weights=weights), model_fun=model_fun)
        names(beta_new) <- colnames(mm)

    

    mt <- attr(mf, "terms")
    res <- list(coefficients = beta_new, q=q_new, niter = iter, loglik=logl,
                link = link, call = mc, terms = mt,
                xlevels = .getXlevels(mt, mf),
                data_object=list(model_matrix=mm, resp=Y, n=rowSums(Y), weights=weights),
                model_fun=model_fun)
    class(res) <- "sprr"
    res

}

pred_lp <- function(beta, data_object){
    lp <- data_object$model_matrix %*% beta
    c(lp)
}
pred_theta <- function(beta, data_object, model_fun){
    lp <- pred_lp(beta, data_object)
    theta <- model_fun(lp)
    theta
}

pred_lambda <- function(beta, q, data_object, model_fun){
    lp <- data_object$model_matrix %*% beta
    theta <- model_fun(lp)
    N <- length(q)-1
    lambda_N <- lambda_from_p(q)
    th <- sapply(0:N, function(k)theta^k)
    lambda <- apply(th, 1, function(t)t * lambda_N)
    rownames(lambda) <- 0:N
    t(lambda)
}
pred_pvec <- function(beta, q, data_object, model_fun){
   ll <- pred_lambda(beta, q, data_object, model_fun)
   cs <- data_object$n
   pp <- lapply(seq_along(cs), 
                function(i){res <- p_from_lambda(ll[i,], n=cs[i]);
                            names(res) <- 0:cs[i];
                            res})
   pp
}
pred_p <- function(beta, q, data_object, model_fun){
   pvec <- pred_pvec(beta, q, data_object, model_fun)
   rvec <- data_object$resp[,1]
   pp <- sapply(seq_along(rvec), function(i)pvec[[i]][rvec[i]+1])
   pp
}
loglik <- function(beta, q, data_object, model_fun){
    p <- pred_p(beta=beta, q=q, data_object = data_object, model_fun=model_fun)
    w <- data_object$weights
    if (is.null(w)) ll <- sum(log(p))
    else ll <- sum(ifelse(w==0, 0, w*log(p)))
    ll
}

# based on print.lm
print.sprr <- function(x, digits = max(3L, getOption("digits") - 3L),...){
  cat("\nA semi-parametric relative risk regression model fit\n")
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(x$coefficients, digits = digits),
            print.gap = 2, quote = FALSE)
    }
    else cat("No coefficients\n")

    if (length(x$q)){
        cat("Baseline joint probabilities (mu):\n")     
        N <- length(x$q)-1
        ll <- c(lambda_from_p(x$q))
        names(ll) <- 0:N
        print.default(format(ll, digits = digits),
            print.gap = 2, quote = FALSE)
    }
    else cat("No baseline joint probabilities\n")

   cat("Log-likelihood: ", format(x$loglik, digits=digits), "\n")
   invisible(x)
}

predict.sprr <- function(object, newdata=NULL,
                              type=c("mean", "relrisk", "likelihood", "probvec", "lvec", "lp"),
                              newn=NULL, ...){
  type <- match.arg(type)
  tt <- terms(object)
  if (!missing(newdata)){    
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses")))
            .checkMFClasses(cl, m)
    mm <- model.matrix(Terms, m)
    data_object <- list(model_matrix=mm)
    if (!missing(newn)){
       if (!(length(newn) == 1L || length(newn) == nrow(newdata)))
          stop("'newn' should have length 1 or equal to the number of rows of 'newdata'")
        data_object$n <- rep(newn, length=nrow(newdata))
       }
  } else {
    data_object <- object$data_object
  }
  

  if (type=="likelihood"){
     if (!missing(newdata))
        stop("Type = 'likelihood' is not available for new data. Consider using type='probvec' to get vector of likelihood values.")
     pred <- pred_p(beta=object$coefficients, q=object$q,
                    data_object=data_object,
                    model_fun=object$model_fun)
  } else
  if (type %in% c("mean", "lvec")){
     ll <- pred_lambda(beta=object$coefficients, q=object$q,
                    data_object=data_object,
                    model_fun=object$model_fun)
     pred <- if (type=="mean") ll[,2] else ll           
  } else
  if (type == "lp"){
     pred <- pred_lp(beta=object$coefficients, data_object=data_object)
  } else
  if (type == "relrisk"){
     pred <- pred_theta(beta=object$coefficients, data_object=data_object,
                        model_fun=object$model_fun)
  }
  if (type == "probvec"){
    if (!missing(newdata) && missing(newn))
       stop("For prediction of probability vectors with new data, cluster sizes should be specified in 'newn'")
    pred <- pred_pvec(beta=object$coefficients, q=object$q, data_object=data_object,
                        model_fun=object$model_fun) 
  }
  return(pred)    
}
