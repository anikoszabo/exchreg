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

\title{Gaussian mixture model for clustered observations}
\author{Aniko Szabo}
\date{\today}


\begin{document}
\maketitle

Let $Y_ij$, $i=1,\ldots,I$, $j=1,\ldots,n_i$ denote a continuous random variable 
measured in subject $j$ of cluster $i$. We want to model the distribution of $Y$
as a mixture of $K$ latent normal distributions:
\begin{equation}
Y_{ij} \mid \xi_{ij}=k \sim N(\mu, \sigma_k^2),
\end{equation}
where $\xi_{ij}$ is the latent class indicato for subject $j$ of cluster $i$.

Due to the clustering we cannot assume that $\xi_{ij}$ are independent, however we 
will assume that they are \textit{exchangeable} with \textit{marginal compatibility}.

For now, let's assume $K=2$ -- a two-component mixture.

The likelihood is
\begin{equation}
L = \prod_{i=1}^I \sum_{x_1,\ldots,x_{n_i}} Pr(\xi_{i1}=x_1, \ldots, \xi_{in_i}=x_{n_i})
\times \prod_{j=1}^{n_i} \phi(y_ij; \mu_{x_j}, \sigma^2_{x_j}),
\end{equation}
where $\phi(y; \mu, \sigma^2)$ is the normal pdf with mean $\mu$ and variance
$\sigma^2$.


\section{Initial setup}
@O R/MixMod.R @{
require(CorrBin)
@}

\end{document}
