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

\title{Semi-parametric generalized linear model for correlated binary data}
\author{Aniko Szabo}
\date{\today}


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

@O R/SPGLM.R

#'@@rdname spglm
#'@@param formula a one-sided formula of the form \code{cbind(r, s) ~ predictors} where \code{r} and \code{s} give the number of responses and non-responses within each cluster, respectively (so the cluster size is \code{r+s}), and \code{predictors} describes the covariates.
#'@@param data  an optional matrix or data frame containing the variables in the formula \code{formula}. By default the variables are taken from \code{environment(formula).}
#'@@param subset	an optional vector specifying a subset of observations to be used.
#'@@param weight	an optional vector specifying observation weights.
#'@@param link	  a link function for the mean.
#'@@param start	  an optional list with elements named \code{beta}, \code{q}, and/or \code{mu1} giving starting values for estimation. If none or only some are specified, starting values are selected automatically.
#'@@param control	a list with parameters controlling the algorithm.
#'@@return an object of class \code{spglm} with the fitted model.
#'@@export
#' @@importFrom stats terms model.matrix

#' Semi-parametric generalized linear model
spglm <- function(formula, data, subset, weights, link="logit", start=NULL, control=list(eps=0.001, maxit=100), ...){
    fam <- binomial(link=link)
    
    @< Create model matrix from formula and data@>
    @< Fit model@>

    mt <- attr(mf, "terms")
    res <- list(coefficients = beta_new, q=q_new, niter = iter, loglik=logl,
                link = link, call = mc, terms = mt,
                xlevels = .getXlevels(mt, mf),
                data_object=list(model_matrix=mm, resp=Y, n=rowSums(Y), weights=weights),
                model_fun=model_fun)
    class(res) <- "spglm"
    res

}
@}

\end{document}

@D Create model matrix from formula and data @{
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
@}
@}

\begin{document}
\maketitle
\end{document}