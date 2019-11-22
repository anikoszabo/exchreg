\documentclass[reqno]{amsart}
\usepackage[margin=1in]{geometry}
\usepackage[colorlinks=true,linkcolor=blue]{hyperref}
\renewcommand{\NWtarget}[2]{\hypertarget{#1}{#2}}
\renewcommand{\NWlink}[2]{\hyperlink{#1}{#2}}
\newcommand{\vv}[1]{\boldsymbol{#1}}
\newcommand{\ff}[1]{\dot{#1}}
\newcommand{\Z}[2][K]{\mathbb{Z}_{#1}^{#2}}

\title{Gaussian mixture model for clustered observations}
\author{Aniko Szabo}
\date{\today}


\begin{document}
\maketitle

\section{Mixture model}
Let $Y_{ij}$, $i=1,\ldots,I$, $j=1,\ldots,n_i$ denote a continuous random variable 
measured in subject $j$ of cluster $i$. We want to model the distribution of $Y$
as a mixture of $K$ latent normal distributions:
\begin{equation}
[Y_{ij} \mid \xi_{ij}=k] \sim N(\mu_k, \sigma_k),
\end{equation}
where $\xi_{ij}$ is the latent class indicator for subject $j$ of cluster $i$
taking values from 1 to $K$.

Due to the clustering we cannot assume that $\xi_{ij}$ are independent, however we 
will assume that they follow and exchangeable multinomial distribution.

Denoting $\phi(y; \mu, \sigma)$  the normal pdf with mean $\mu$ and standard deviation 
$\sigma$, the likelihood is
\begin{align}
L = & \prod_{i=1}^I \sum_{\vv{x}\in \Z{n_i}} Pr(\vv{\xi}_{i}=\vv{x})
\times \prod_{j=1}^{n_i} \phi(y_{ij}; \mu_{x_j}, \sigma_{x_j})\\
 = &  \prod_{i=1}^I \sum_{\vv{x}\in \Z{n_i}} q_{\ff{\vv{x}}|n_i} 
\times \prod_{j=1}^{n_i} \phi(y_{ij}; \mu_{x_j}, \sigma_{x_j}),
\end{align}
where $\vv{x} = (x_1,\ldots,x_{n_i})$, $\vv{\xi}_{i}=(\xi_{i1},\ldots,\xi_{in_i})$, 
$\ff{\vv{x}} = (\sum_j I(x_j=1),\ldots,\sum_j I(x_j=K))$ is the vector of latent category 
frequencies, and
\begin{equation}
q_{\vv{r}|n_i} = \sum_{\vv{t}\in\Z{N}} h(\vv{r}, \vv{t}, N, n_i) q_{\vv{t}|L}
\end{equation}
are the exchangeable and marginally reproducible outcome probabilities of the 
exchangeable multinomial distribution with maximal cluster size $N$ with 
$h(\vv{r}, \vv{t}, N, n_i)$ defined as multivariate hypergeometric probabilities.



\section{EM algorithm}

We can consider $\xi_{ij}$ missing data. Then the complete data likelihood is
\begin{equation}
L_c(\mu, \sigma, q | \xi_{ij}, y_{ij}) = 
 \sum_{i=1}^I  q_{\ff{\vv{\xi}}|n_i} \prod_{j=1}^{n_i} \log \phi(y_{ij}; \mu_{\xi_{ij}}, \sigma_{\xi_{ij}}) =
 \sum_{i=1}^I \prod_{\vv{z}\in\Z{n_i}} \big[q_{\ff{\vv{z}}|n_i}
    \prod_{j=1}^{n_i} \log \phi(y_{ij}; \mu_{z_j}, \sigma_{z_j})\big]^{I(\vv{\xi}_i = \vv{z})}.
\end{equation}

Its logarithm is

\begin{equation}
\log L_c(\mu, \sigma, q | \xi_{ij}, y_{ij}) = 
  \sum_{i=1}^I \sum_{\vv{z}\in\Z{n_i}} I(\vv{\xi}_i = \vv{z}) \big[\log q_{\ff{\vv{z}}|n_i} +  
        \sum_{j=1}^{n_i} \log \phi(y_{ij}; \mu_{z_j}, \sigma_{z_j})\big]
\end{equation}

The expected complete data log-likelihood, given previous estimates $\mu^{(m)}, \sigma^{(m)}, q^{(m)}$ is
\begin{multline}
Q(\mu, \sigma, q) = E\big[\log L_c(\mu, \sigma, q \mid \mu^{(m)}, \sigma^{(m)}, q^{(m)}, y_{ij})\big] = \\
 \sum_{i=1}^I \sum_{\vv{z}\in\Z{n_i}}  Pr(\vv{\xi}_{i}=\vv{z} \mid \mu^{(m)}, \sigma^{(m)}, q^{(m)}, y_{ij} ) 
 [\log q_{\ff{\vv{z}}|n_i} +\sum_{j=1}^{n_i}\log \phi(y_{ij}; \mu_{x_j}, \sigma_{x_j})],
\end{multline}

\subsection{E-step}
Using the Bayes theorem
\begin{align}
 e_{i\vv{z}}^{(m)} = & Pr(\vv{\xi}_{i}=\vv{z} \mid \mu^{(m)}, \sigma^{(m)}, q^{(m)}, y_{ij} )=\\
  & \frac{Pr(\vv{y}_{i} \mid \vv{\xi}_i=\vv{z},  \mu^{(m)}, \sigma^{(m)}, q^{(m)}) Pr(\vv{\xi}_{i}=\vv{z} \mid \mu^{(m)}, \sigma^{(m)}, q^{(m)})}%
  {\sum_{\vv{w}} Pr(\vv{y}_{i} \mid \vv{\xi}_i=\vv{w}, \mu^{(m)}, \sigma^{(m)}, q^{(m)}) Pr(\vv{\xi}_{i}=\vv{w} \mid \mu^{(m)}, \sigma^{(m)}, q^{(m)}) } = 
  \frac{\prod_{j=1}^{n_i} \phi(y_{ij}; \mu_{z_j}^{(m)}, \sigma^{(m)}_{z_j})  q^{(m)}_{\ff{\vv{z}}|n_i}}%
  {\sum_{\vv{w}}\prod_{j=1}^{n_i} \phi(y_{ij}; \mu_{w_j}^{(m)}, \sigma^{(m)}_{w_j})  q^{(m)}_{\ff{\vv{w}}|n_i}  }.  
\end{align}

\subsection{M-step}

\begin{align}
Q(\mu, \sigma, q) = & \sum_{i=1}^I\sum_{\vv{z}\in\Z{n_i}} e_{i\vv{z}}^{(m)}   
 [\log q_{\ff{\vv{z}}|n_i} +\sum_{j=1}^{n_i}\log \phi(y_{ij}; \mu_{z_j}, \sigma_{z_j})] =\\
 &  \sum_{i=1}^I\sum_{\vv{z}\in\Z{n_i}} e_{i\vv{z}}^{(m)} \log q_{\ff{\vv{z}}|n_i} +
  \sum_{i=1}^I\sum_{\vv{z}\in\Z{n_i}} e_{i\vv{z}}^{(m)}  \sum_{j=1}^{n_i}\log \phi(y_{ij}; \mu_{z_j}, \sigma_{z_j})  =\\
 &  \sum_{i=1}^I\sum_{\vv{z}\in\Z{n_i}} e_{i\vv{z}}^{(m)} \log q_{\ff{\vv{z}}|n_i} +
  \sum_{k=1}^K \sum_{i=1}^I \sum_{j=1}^{n_i}\sum_{\vv{z}: z_j=k} e_{i\vv{z}}^{(m)}  \log \phi(y_{ij}; \mu_k, \sigma_k),
\end{align}
where the two components can be maximized separately, and within the second component the terms corresponding to different $k$'s
can be maximized separately..


%For now, let's assume $K=2$ -- a two-component mixture.

\end{document}

\section{Initial setup}
@O R/MixMod.R @{
require(CorrBin)
@}
