\documentclass[reqno]{amsart}
\usepackage[margin=1in]{geometry}
\usepackage[colorlinks=true,linkcolor=blue]{hyperref}
\renewcommand{\NWtarget}[2]{\hypertarget{#1}{#2}}
\renewcommand{\NWlink}[2]{\hyperlink{#1}{#2}}
\newcommand{\bv}{\mathbf{v}}
\newcommand{\bq}{\mathbf{q}}
\newcommand{\bpi}{\text{\boldmath $\pi$}}
\newcommand{\leqst}{\mathrel{\preceq^{st}}}
\newcommand{\geqst}{\mathrel{\succeq^{st}}}
\newcommand{\qz}{\boldsymbol{q}^{(0)}}
\DeclareMathOperator{\E}{\mathbb{E}}

\usepackage{dsfont}



\title{Semi-parametric generalized linear model for correlated binary data}
\author{Aniko Szabo}
\date{\today}
\begin{document}
\maketitle

Note: this setup assumes that parts defined in \texttt{SPregress.w} are available.

\section{Introduction}

We extend the Rathouz-Gao semi-parametric generalized linear model to clustered binary outcomes with varying cluster sizes. Treating the number of events $Y_i$ from a cluster scaled by corresponding cluster size $n_i$ as the response variable, we use a parametric regression model for the marginal response probability $Y_/n_i$ and assume that for the maximum cluster size $N$ the probability density distribution of $Y_i$ is an exponentially tilted version of the reference density distribution.

The conditional mean model is as follows:
\begin{equation}\label{M:conditionalMeanSPGLM}
     E(\dfrac{Y_i}{n_i} | \boldsymbol{Z}_i; \boldsymbol{\beta})=\mu(\boldsymbol{Z}_i,\boldsymbol{\beta}) \equiv \mu_i = h^{-1}( \boldsymbol{Z}_i^{T} \boldsymbol{\beta})
\end{equation}
where $\boldsymbol{Z}_i$ is a matrix of cluster-level covariates, the function $h(\cdot)$ is a known strictly increasing link function.

The assumed probability density function for a cluster of maximal size $N$ is:
\begin{equation}\label{M:densitySPGLM}
    q_{y,N} (\boldsymbol{Z}, \boldsymbol{\beta}) \propto q_{y,N}^{(0)} \times \exp [ - \omega(\boldsymbol{Z}, \boldsymbol{\beta}) \times y]
\end{equation}
where $q_{y,N}^{(0)}$ is a non-parametric reference distribution, and $\omega(\boldsymbol{Z}, \boldsymbol{\beta})$ is the tilting parameter chosen so that the expectation of $Y/N$ equals to the value defined by the mean model \eqref{M:conditionalMeanSPGLM}, i.e.\  $\sum_{y=0}^{N} \dfrac{y}{N} q_{y,N}(\boldsymbol{Z}, \boldsymbol{\beta}) = \mu = h^{-1}( \boldsymbol{Z}^{T} \boldsymbol{\beta}) $. 
%We define the tilting parameter for the reference density to be $0$, i.e., $\omega_0 \equiv 0$.

For the purpose of identifiability, the mean of marginal probability $\{ q_{y,N}^{(0)} \}$ is set at an arbitrary fixed value, which can be the mean of the reference group, or the overall marginal mean of the data. We extend the conditional mean and probability density model for  cluster sizes smaller than the maximal size $N$ by assuming marginal compatibility.

\section{Main function}

@O ../R/SPGLM.R @{
#' Fit semi-parametric GLM
#'
#'@@rdname spglm
#'@@param formula a one-sided formula of the form \code{cbind(r, s) ~ predictors} where \code{r} and \code{s} give the number of responses and non-responses within each cluster, respectively (so the cluster size is \code{r+s}), and \code{predictors} describes the covariates.
#'@@param data  an optional matrix or data frame containing the variables in the formula \code{formula}. By default the variables are taken from \code{environment(formula).}
#'@@param subset	an optional vector specifying a subset of observations to be used.
#'@@param weights	an optional vector specifying observation weights.
#'@@param offset	an optional vector specifying an offset term (predictor with a fixed coefficient of 1).
#'@@param link	  a link function for the mean.
#'@@param mu0	  an optional numeric value constraining the mean of the baseline distribution
#'@@param control	a list with parameters controlling the algorithm.
#'@@return an object of class \code{spglm} with the fitted model.
#'@@export
#' @@importFrom stats terms model.matrix model.frame model.offset model.response model.weights

spglm <- function(formula, data, subset, weights, offset, link="logit", mu0=NULL, 
                  control=list(eps=0.001, maxit=100)){

    @< Create model matrix from formula and data@>
    
    data_object <- list(model_matrix=mm, resp=Y, n=rowSums(Y), weights=weights, offset=offset, maxN=N, spt=0:N)
    
    @< Fit model@>

    mt <- attr(mf, "terms")
    names(betas) <- colnames(mm)
    names(referencef0) <- 0:N
    rownames(vc) <- colnames(vc) <- c(names(betas), names(referencef0))
    res <- list(coefficients = betas, f0=referencef0, vcov = vc, mu0=mu0, niter = iter, 
                loglik=llik, link = link, call = mc, terms = mt,
                xlevels = .getXlevels(mt, mf),
                data_object=data_object)
    class(res) <- "spglm"
    res

}
@}



@D Create model matrix from formula and data @{
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
  
@}

\subsection{Methods for the \texttt{spglm} class}

First, define a printing method which does not show the saved data and model matrix

@O ../R/SPGLM.R @{
#' @@export
#' @@importFrom stats printCoefmat
print.spglm <- function(x, digits = max(3L, getOption("digits") - 3L),...){
  # based on print.lm, summary.lm, and print.summary.lm
  cat("\nA semi-parametric generalized linear regression model fit\n")
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
  if (length(x$coefficients)) {
    beta <- x$coefficients
    p <- length(beta)
    SEbeta <- sqrt(diag(x$vcov[1:p, 1:p]))
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
@}

Custom extractors for the coefficients and variance-covariance matrix, with option to get only beta's (default), only the backbone pdf, or both.

@O ../R/SPGLM.R @{
#' Extracting components of fitted `spglm` object
#'
#' @@param object fitted model of class \code{spglm}
#' @@param beta logical, whether the regression parameter estimates / variances should be returned
#' @@param f0 logical, whether the backbone density estimates / variances should be returned
#' @@param ... not used, present fo consistency with generics
#' @@importFrom stats coef
#' @@export

coef.spglm <- function(object, beta=TRUE, f0=FALSE, ...){
  res <- NULL
  if (beta) res <- c(res, object$coefficients)
  if (f0) res <- c(res, object$f0)
  res
}

#' @@rdname coef.spglm
#'
#' @@importFrom stats vcov
#' @@export
vcov.spglm <- function(object, beta=TRUE, f0=FALSE, ...){
  p <- length(object$coefficients)
  m <- nrow(object$vcov)
  
  idx <- NULL
  if (beta) idx <- c(idx, 1:p)
  if (f0) idx <- c(idx, (p+1):m)
  
  vc <- object$vcov[idx,idx]
  vc
}

#' @@rdname coef.spglm
#'
#' @@importFrom stats logLik
#' @@export
logLik.spglm <- function(object, ...){
  object$loglik
}
                              
@}

The prediction method will predict for a variety of scenarios:
\begin{itemize}
\item Input data set:
  \begin{itemize}
    \item the data used in the fitting;
    \item new data.
  \end{itemize}
\item Results:
  \begin{itemize}
    \item $\mu(z,\beta)$: the mean event probability, given $z$;
    \item $p_{r,n}(z)$: the probability of observing the given $r$ responses with cluster size $n$, given $z$;
    \item $\{p_{\cdot,n}(z)\}$: the entire vector of response probabilities for cluster size $n$, given $z$ (will be a list due to varying lengths);
    \item $\beta'z$: the linear predictor value $z$;
    \item $\omega(\beta'z)$: the tilting parameter at predictor value $z$.
  \end{itemize}
\end{itemize}

@O ../R/SPGLM.R @{
#' Predict methods for SPGLM fits
#'
#' Obtains predictions from a fitted semi-parametric generalized linear model. Note that \code{offset}
#' and \code{weight} terms are not implemented for predicting from new data.
#' @@param object fitted model of class \code{spglm}
#' @@param newdata optionally, a data frame in which to look for covariate values for prediction. 
#'  If NULL, the original data set is used
#' @@param type the type of prediction requested. The default is "mean", the mean event probability;
#' "prob" requests the probability of the observation (given number of responses with given cluster size), 
#' "tilt" returns the tilting parameter which achieves the modeled mean from the baseline distribution;
#' and "lp" requests the linear predictor.
#' @@param newn if \code{newdata} is provided and \code{type="prob"}, an integer or integer vector specifying the clustersize for the predictors
#' @@param newevents if \code{newdata} is provided and \code{type="prob"}, an integer or integer vector specifying the 
#'  number of events for the predictions
#' @@param ... not used, present for consistency with generic
#' @@export
#' @@importFrom stats  .checkMFClasses .getXlevels delete.response

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
@}

\section{EM algorithm for model fitting}
We will combine an Expectation Maximization with the Newton-Raphson algorithm proposed by Rathouz and Gao (2009) to estimate the model parameters. Recall that the observed data consists of the number of events, cluster size, and cluster-level covariates: $(y_i, N_i, \boldsymbol{Z}_i)$, $i=1,\cdots,I$. 

\subsection{Observed data log-likelihood}
Under the marginal compatibility assumption, the observed log-likelihood function is:
\begin{equation}\label{E:loglikelihood}
\begin{split}
    \ell (q_{y,N}, \boldsymbol{\beta}) &= \sum_{i=1}^I \log \{ q_{y_i,N_i} (\boldsymbol{Z}_i, \boldsymbol{\beta}) \} \\
    &= \sum_{i=1}^I \log \bigg[ \binom{N_i}{y_i} \sum_{t=y}^{N-N_i+y_i} \dfrac{\binom{N-N_i}{t-y_i}}{\binom{N}{t}} \times \dfrac{q_{t,N}^{(0)} \times \exp \{ \omega(\boldsymbol{Z}_i; \boldsymbol{\beta}) t \}} {\sum_{t'=0}^{N} q_{t',N}^{(0)} \times \exp \{ \omega(\boldsymbol{Z}_i; \boldsymbol{\beta}) t' \}}  \bigg].
\end{split}
\end{equation}

\subsection{Model prediction and likelihood}
We define internal function to calculate predicted values and the log-likelihood.

@O ../R/SPGLM.R @{
#' @@keywords internal
#' @@importFrom stats dhyper
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
    spt <- data_object$spt
    
    # ySptIndex is only used to calculate the log-likelihood, which we will not be using
    th <- getTheta(spt=spt/N, f0=f0[spt+1], mu=mu, weights=data_object$weights, ySptIndex=rep(1, nobs),
                   thetaStart=NULL, thetaControl=theta.control())
    th$theta
}

spglm_probs <- function(beta, f0, data_object, link){

    mu <- spglm_pred_mean(beta=beta, data_object=data_object, link=link)
    N <- data_object$maxN
    nobs <- nrow(data_object$model_matrix)
    spt <- data_object$spt
   
    # ySptIndex is only used to calculate the log-likelihood, which we will not be using
    th <- getTheta(spt=spt/N, f0=f0[spt+1], mu=mu, weights=data_object$weights, ySptIndex=rep(1, nobs),
                   thetaStart=NULL, thetaControl=theta.control())
                   
    hp <- sapply(0:N, function(t)dhyper(x=data_object$resp[,1], m=data_object$n, n=N-data_object$n, k=t))  
    
    # fill in 0's for values with no support
    fTilt <- matrix(0, ncol = nobs, nrow = N+1)
    fTilt[spt+1,] <- th$fTilt
    #probs <- rowSums(t(th$fTilt) * hp)
    probs <- rowSums(t(fTilt) * hp)
    probs
}

spglm_loglik <- function(beta, f0, data_object, link){
  
    probs <- spglm_probs(beta, f0, data_object, link)
    llik <- log(probs) %*% data_object$weights
    c(llik)
}
@}

It is difficult to utilize the observed log-likelihood directly for parameter estimation. Its score function is complex because there is a logarithm of a sum in \eqref{E:loglikelihood}. Therefore, we implement an EM algorithm for parameter estimation.


@D Fit model @{
 @< Set starting values@>
 @< Set up for E-step @>
 @< Identify zero support@>

 iter <- 0
 difference <- 100
 while (iter < control$maxit & difference > control$eps) {
    iter <- iter + 1
    referencef0Pre <- referencef0
    betasPre <- betas
    llikPre <- llik
    
    @< E-step @>
    @< M-step @>
    
    difference <- abs(llik - llikPre)
  }
  
  @< Calculate standard errors @>
  
@}


\subsection{Missing data setup}
Stefanescu et al (2003) have proved that the marginal compatibility assumption is equivalent to assuming that each cluster of size $N'$ arises from a cluster of maximal cluster size $N$ with some binary observations missing completely missing at random (MCAR). The missing data setup and the expectation step from EM algorithm are based on this MCAR interpretation. 

For the setup of the EM algorithm, the complete data are clusters of size $N$ with a corresponding number of events $s_i$, where $s_i \in \{y_i, \cdots, N \}$. The missing data correspond to the unobserved outcomes in of the $(N - N_i)$ cluster components. The missing data can be summarized as $(s_i - y_i)$ events out of $(N - N_i)$ elements.

The complete data log-likelihood is:
\begin{equation}\label{E:logLikelihoodCompleteMapFromObs}
\begin{split}
    \ell_{\text{complete}} (q_{y,N}, \boldsymbol{\beta} \mid N_i,s_i) &= \sum_{i=0}^I \log \{ q_{s_i,N}(\boldsymbol{Z}_i; \boldsymbol{\beta}) \} \\
    &= \sum_{i=0}^I \sum_{s=y_i}^N {\mathds{1}}(s_i=s) \times \log \{ q_{s,N}(\boldsymbol{Z}_i; \boldsymbol{\beta}) \}, 
\end{split}
\end{equation}
where $q_{s_i,N}$ is the probability of achieving $s_i$ events in a cluster of size $N$ in the complete data set. The $\mathds{1}$ denotes the indicator function.


Equation \ref{E:logLikelihoodCompleteMapFromObs} can be interpreted as the log-likelihood of a model with fixed cluster size $N$ for an expanded data set: each observation with cluster size $N_i$ and number of responses $y_i$ is expanded to $(N+1)$ replicates with $0, 1, \ldots, N$, responses, respectively.

@D Set up for E-step @{
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

@}

\subsection{Expectation step}

We take the expectation of complete data log-likelihood \eqref{E:logLikelihoodCompleteMapFromObs} given parameters from the previous $k^{th}$ iteration, which is shown as follows:
\begin{equation}\label{E:ExpectedlogLikelihoodCompleteMapFromObs}
\begin{split}
    Q \bigg( q_{y,N}^{(k)}, \boldsymbol{\beta}^{(k)} \bigg) &= 
     \E \bigg\{ \ell_{\text{complete}} (q_{y,N}, \boldsymbol{\beta} \mid N_i, s_i) \mid q_{y,N}^{(k)}, \boldsymbol{\beta}^{(k)}, N_i, y_i \bigg\} \\
    &= \E \Bigg[ \sum_{i=0}^I \sum_{s=0}^N {\mathds{1}}(s=s_i) \times \log \{ q_{s,N}(\boldsymbol{Z}_i; \boldsymbol{\beta}) \} \mid q_{y,N}^{(k)}, y_i,N_i \Bigg] \\
    &= \sum_{i=0}^I \sum_{s=0}^N \underbrace{ \E \big\{ {\mathds{1}}(s=s_i) \mid q_{y,N}^{(k)}, \boldsymbol{\beta}^{(k)}, y_i, N_i \big\}}_{p_{N_iy_is_i}^{(k)}} \times \log \{q_{s,N}(\boldsymbol{Z}_i; \boldsymbol{\beta}) \},
\end{split}
\end{equation}
where $p_{N'y'y} = \textbf{Pr}_N(Y=y \mid Y'=y',N')$ denotes the conditional probability of achieving $y$ events in the complete cluster of size $N$, given that $y'$ events have been observed in the original cluster of size $N'$. 
The expression for $p_{N'y'y}$ can be derived from $q_{y,N}$ using Bayes theorem and the MCAR interpretation as follows:
\begin{equation}\label{D:prst}
\begin{split}
p_{N'y'y} &= \textbf{Pr}_N(Y=y \mid Y'=y',N') \\ 
          &= \frac{\textbf{Pr}_{N'}(Y'=y' \mid Y=y) \textbf{Pr}_N(Y=y)}{\sum_{t=0}^{N} \textbf{Pr}_{y'}(Y'=y' \mid Y=t) \textbf{Pr}_N(Y=t)} \\ 
          &= \frac{\binom{y}{y'} \binom{N-y}{N'-y'} q_{y,N}}{\sum_{t=y'}^{N-N'+y'} \dbinom{t}{y'} \binom{N-t}{N'-y'} q_{t,N}},
\end{split}
\end{equation}
and $p_{N_iy_is_i}^{(k)}$ at the $k$th step of the iteration can be obtained by using $q_{y,N}^{(k)}$ on the right-hand side of the expression. 

The value of $q_{y,N}^{(k)}(\boldsymbol{Z}_i; \boldsymbol{\beta}^{(k)})$ can be obtained from ${q_{y,N}^{(0)}}^{(k)}$ using \ref{M:densitySPGLM}. But it is also returned by \texttt{gldrmFit} as \texttt{fTiltMatrix}.


@D E-step @{
  # convert fTiltMatrix to long vector, watching out for potentially different support
  y_idx <- match(Ycomb[,"y"], 0:N)
  fTilt2 <- ifelse(!is.na(y_idx), fTiltMatrix[cbind(Ycomb[,"i"], y_idx)], 0)
  
  # numerator of weights
  pp_num <- choose(n=Ycomb[,"y"], k=Ycomb[,"resp"]) * 
           choose(n=N-Ycomb[,"y"], k=Ycomb[,"nonresp"]) * 
           fTilt2
  # denominator 
  pp_denom <- tapply(pp_num, list(i=Ycomb[,"i"]), sum, simplify=TRUE)
  pp <- c(pp_num / pp_denom[Ycomb[,"i"]])
  if (!is.null(weights)) pp <- weights2 * pp
@}

\subsection{Maximization step}

Equation \eqref{E:ExpectedlogLikelihoodCompleteMapFromObs} can be interpreted as the log-likelihood for the SPGLM with fixed cluster size $N$ for an expanded data set: each observation with cluster size $N_i$ and number of responses $y_i$ is expanded to $(N+1)$ replicates with $0, 1, \ldots, N$, responses, respectively. The probabilities $p_{N_iy_is_i}^{(k)}$ serve as cluster weights in the expanded data set. From here we can use a slight generalization of the Newton-Raphson algorithm developed by Wurm and Rathouz (2018) to estimate the baseline density distribution $ \boldsymbol{q}_N^{(0)} $ and the regression coefficients $\boldsymbol{\beta}$. 

Default settings are used for the algorithm, except the starting values is updated based on the previous iteration

@D M-step @{
  gldrmControl0 <- gldrm.control(returnfTiltMatrix = TRUE, returnf0ScoreInfo = FALSE, print=FALSE,
                                betaStart = betasPre, f0Start = referencef0Pre[spt+1])

  mod <- gldrmFit(x = mm2, y=Ycomb[,"y"]/N, linkfun=link$linkfun, linkinv = link$linkinv,
                  mu.eta = link$mu.eta, mu0 = mu0, offset = offset2, weights = pp, 
                  gldrmControl = gldrmControl0,  thetaControl=theta.control(),
                  betaControl=beta.control(), f0Control=f0.control())
                  
  fTiltMatrix[,spt+1] <- mod$fTiltMatrix[obs_start,]
  betas <- mod$beta
  referencef0[spt+1] <- mod$f0
#  spt <- round(mod$spt * N)
  llik <- c(log(rowSums(fTiltMatrix * hp)) %*% weights)
@}

\subsection{Starting values}
We start the regression coefficients are set to 0, so $q^{(0)}_N$ would be the overall estimate under marginally compatibility. The default value of \texttt{mu0} is set to the mean of $y_i/n_i$.

@D Set starting values @{
  betas <- NULL
  pooled <- CorrBin::CBData(data.frame(Trt = "All", NResp = Y[,1], ClusterSize = rowSums(Y), Freq=ceiling(weights)), 
                    trt="Trt", clustersize="ClusterSize", nresp="NResp", freq="Freq")
  est <- CorrBin::mc.est(pooled)
  referencef0 <- est$Prob[est$ClusterSize == N]
  
  
@}

Values of $y$ that cannot be a possible number of responses in a cluster of size $N$ corresponding to an observed cluster will always have $q_{y,N}=0$, so they need to be removed from the support during the `gldrm` estimation process (which will use $\log q_{y,N}$)

@D Identify zero support @{
  
  # identify possible y-values
  spt <- sort(unique(Ycomb[,"y"]))
  data_object$spt <- spt
  
  # ensure all positive values within support
  referencef0[spt+1] <- referencef0[spt+1] + 1e-6
  referencef0 <- referencef0 / sum(referencef0)

  
  fTiltMatrix <- matrix(rep(referencef0, times=nobs), byrow=TRUE, 
      nrow=nobs, ncol=N+1)

  if (is.null(mu0))
    mu0 <- weighted.mean(Y[,1]/rowSums(Y), weights)
  
  # hypergeometric terms for log-likelihood calculation
  hp <- sapply(0:N, function(t)dhyper(x=Y[,1], m=rowSums(Y), n=N-rowSums(Y), k=t))   
  
  # initial log-likelihood
  llik <- log(rowSums(fTiltMatrix * hp)) %*% weights

@}

\subsection{Standard errors}

We will use numeric derivation to obtain the hessian of the observed data likelihood. The variance-covariance matrix can be estimated as the inverse of the negative hessian. The sum-to-one and fixed-mean constraints of the reference distribution are included by expanding the hessian to a bordered hessian before inversion by including the derivatives of the constraints (this can be derived based on the Lagrange multiplier method). The reference distribution is included in the hessian on the log-scale, and is multiplied by the gradient before inversion to revert to the original scale.

@D Calculate standard errors @{
  ll <- function(x){
    spglm_loglik(beta=x[1:p], f0 = exp(x[-(1:p)]), data_object=data_object, link=link)
  }

  hess <- numDeriv::hessian(func=ll,  x=c(betas, log(referencef0)))
  hess_idx <- c(1:p, p+spt+1) # betas and f0 elements with non-zero support
  
  # revert to unlogged f0
  grad <- c(rep(1, p), 1/referencef0[spt+1])
  hess1 <- diag(grad) %*% hess[hess_idx, hess_idx] %*% diag(grad)
  
  # create bordered hessian
  border1 <- c(rep(0, p), rep(1, length(spt)))  # gradient of sum-to-one constraint
  border2 <- c(rep(0, p), spt)          # gradient of fixed-mean constraint
  bhess <- rbind(cbind(hess1, border1, border2), c(border1,0,0), c(border2,0,0))
  
  # calculate variance-covariance matrix
  bvc <- solve(-bhess)
  # remove border and set to 0 for zero-support values
  vc <- matrix(0, nrow=p+N+1, ncol=p+N+1)
  vc[hess_idx, hess_idx] <- bvc[1:(p+length(spt)), 1:(p+length(spt))]     
@}


\end{document}