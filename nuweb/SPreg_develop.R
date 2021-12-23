library(devtools)
source('nuweb/Nuweb.R')

#  setting up the package and its infrastructure
create("../exchreg")
use_package("CorrBin")
use_package("numDeriv")
use_build_ignore("nuweb")
use_gpl3_license()
desc::desc_add_author(given = "Xinran", family = "Qi", email = "xinqi@mcw.edu", role = "aut",
               comment = NULL, file = ".", normalize = FALSE)
use_readme_rmd()
use_rcpp_armadillo()

# add data
boric_acid <- CorrBin::read.CBData("z:/EOGeorge/Data/Binary/BoricAcidMousedata_processed.csv",
                  sep=",", skip=1)
use_data(boric_acid)

#############
ex <- as.package("../exchreg")
nuweb(ex)

#shell("cd c:/exchreg/nuweb/ && texify --pdf --quiet --run-viewer  SPGLM.tex")
shell("cd c:/exchreg/nuweb/ && texify --pdf --quiet SPGLM.tex")
shell("cd c:/exchreg/nuweb/ && texify --pdf --quiet wgldrm.tex")


document(ex)
load_all(ex)

build_readme()

check(ex, check_dir = "c:/Temp", cran = TRUE, manual=TRUE)


## test SPGLM
data(boric_acid)
m <- spglm(cbind(NResp, ClusterSize - NResp) ~ Trt, data=boric_acid, link="logit",
           weights=boric_acid$Freq, control = list(eps=1e-5, maxit=100))


ba2 <- boric_acid[rep(1:nrow(boric_acid), boric_acid$Freq),]
m2 <- spglm(cbind(NResp, ClusterSize - NResp) ~ Trt, data=ba2, link="logit", control = list(eps=1e-5, maxit=1000))

all.equal(coef(m), coef(m2))
all.equal(m$loglik, m2$loglik)

# values with 0 support for in-between values
m12 <- spglm(cbind(NResp, ClusterSize - NResp) ~ Trt, 
             data=subset(boric_acid, ClusterSize < 13),
             link="logit", weights=Freq, control = list(eps=1e-5, maxit=100))

# values with 0 support for 0 responses
mzero <- spglm(cbind(NResp, ClusterSize - NResp) ~ Trt, 
             data=subset(boric_acid, NResp > 0),
             link="logit", weights=Freq, control = list(eps=1e-5, maxit=100))

# values with 0 support for 0 & N responses 
mzero2 <- spglm(cbind(NResp, ClusterSize - NResp) ~ Trt, 
               data=subset(boric_acid, (NResp > 0) & (NResp < ClusterSize)),
               link="logit", weights=Freq, control = list(eps=1e-5, maxit=100))


# offset
moff <- spglm(cbind(NResp, ClusterSize - NResp) ~ Trt, offset=Freq, data=boric_acid)
# prediction
a <- spglm_pred_mean(m$coefficients, m$data_object, m$link)
b <- spglm_loglik(beta=m$coefficients, f0=m$f0, m$data_object, m$link)
b2 <- spglm_loglik(beta=m2$coefficients, f0=m2$f0, m$data_object, m$link)

nd <- data.frame(Trt = unique(boric_acid$Trt))
nd$p <- predict(m, newdata=nd, type="mean")
nd$p2 <- predict(moff, newdata=nd, type="mean")
nd

b3 <- predict(m, type="prob")
predict(m, newdata=nd, newn=5, newevents=0, type="prob")
nd_long <- merge(nd, data.frame(Resp=0:5))
nd_long$p3 <- predict(m, newdata=nd_long, newevents = nd_long$Resp, newn=5, type="prob")
p3mat <- matrix(nd_long$p3, nrow=nrow(nd), 
                dimnames = list(Trt = nd$Trt, Resp=0:5))
matplot(0:5, t(p3mat), type="l")

mn <- spglm(cbind(NResp, ClusterSize - NResp) ~ Trt, data=boric_acid, link="log",
              weights=boric_acid$Freq, control = list(eps=1e-5, maxit=1000))

theta <- predict(m, newdata=nd, type="tilt")
predict(moff, newdata=nd, type="lp")

# laplace link
laplace.link <- function(mu) ifelse(mu>0.5, -log(2-2*mu), log(2*mu))
laplace.invlink <- function(eta) 0.5 + 0.5*sign(eta)*(1-exp(-abs(eta)))
laplace.mu.eta <- function(eta) pmax(0.5 * exp(-abs(eta)), .Machine$double.eps)

mll <- spglm(cbind(NResp, ClusterSize - NResp) ~ Trt, data=boric_acid, 
             link=list(linkfun=laplace.link, linkinv=laplace.invlink, mu.eta=laplace.mu.eta),
            weights=boric_acid$Freq, control = list(eps=1e-5, maxit=1000))

# random data
dd <- ran.spglm(1:5, means=seq(0.2, 0.6, length.out = 5), q0 = dbinom(0:5, size=5, prob=0.3))
dd

# profiling -> getTheta takes most of the time originally
# after rewriting it in C, the dhyper calls within spglm_probs stands out
profvis::profvis({m <- spglm(cbind(NResp, ClusterSize - NResp) ~ Trt, data=boric_acid, link="logit",
                   weights=boric_acid$Freq, control = list(eps=1e-5, maxit=100))})

######### Testing wGLDRM
run_gldrm <- 
function (formula, data = NULL, link = "identity", mu0 = NULL,  weights=NULL,
          offset = NULL, gldrmControl = gldrm.control(), thetaControl = theta.control(), 
          betaControl = beta.control(), f0Control = f0.control()) {
  if (is.null(weights)) weights <- rep(1, nrow(data))
  mf <- model.frame(formula, data)
  x <- stats::model.matrix(attr(mf, "terms"), mf)
  attributes(x)[c("assign", "contrasts")] <- NULL
  y <- stats::model.response(mf, type = "numeric")
  if (is.null(offset)) 
    offset <- rep(0, nrow(x))
  if (length(offset) != nrow(x)) 
    stop("offset should be NULL or a vector with length equal to the number of observations.")
  if (is.character(link)) {
    link <- stats::make.link(link)
  }
  else if (!is.list(link) || !(all(c("linkfun", "linkinv", 
                                     "mu.eta") %in% names(link)))) {
    stop(paste0("link should be a character string or a list containing ", 
                "functions named linkfun, linkinv, and mu.eta"))
  }
  yMin <- min(y)
  yMax <- max(y)
  yMed <- (yMin + yMax)/2
  Z <- function(y) (y - yMed) * 2/(yMax - yMin)
  Y <- function(z) z * (yMax - yMin)/2 + yMed
  z <- Z(y)
  linkfunZ <- function(muZ) link$linkfun(Y(muZ))
  linkinvZ <- function(eta) Z(link$linkinv(eta))
  mu.etaZ <- function(eta) 2/(yMax - yMin) * link$mu.eta(eta)
  if (is.null(mu0)) {
    mu0Z <- NULL
  }
  else {
    mu0Z <- Z(mu0)
  }
  modZ <- gldrmFit(x = x, y = z, linkfun = linkfunZ, linkinv = linkinvZ, 
                   mu.eta = mu.etaZ, mu0 = mu0Z, offset = offset, weights = weights, 
                   gldrmControl = gldrmControl, thetaControl = thetaControl, 
                   betaControl = betaControl, f0Control = f0Control)
  modZ$mu <- Y(modZ$mu)
  muetaAdj <- link$mu.eta(modZ$eta)/mu.etaZ(modZ$eta)
  modZ$seMu <- modZ$seMu * muetaAdj
  modZ$mu0 <- Y(modZ$mu0)
  modZ$spt <- Y(modZ$spt)
  modZ$theta <- modZ$theta * 2/(yMax - yMin)
  modZ$bPrime <- Y(modZ$bPrime)
  modZ$bPrime2 <- modZ$bPrime2 * ((yMax - yMin)/2)^2
  names(modZ$beta) <- colnames(x)
  modZ$formula <- formula
  modZ$data <- data.frame(mf)
  modZ$link <- link
  modZ$offset <- offset
  modZ
}

data(shelltox)
sh2 <- shelltox[rep(1:nrow(shelltox), shelltox$Freq), ]

# original code
res0 <- gldrm::gldrm(I(NResp/ClusterSize) ~ Trt, data=sh2, link="logit")
# full dataset without frequency weights
res1 <- run_gldrm(I(NResp/ClusterSize) ~ Trt, data=sh2, link="logit")
# reduced dataset with frequency weights
res2 <- run_gldrm(I(NResp/ClusterSize) ~ Trt, data=shelltox, link="logit", weights=shelltox$Freq)
# rescale weights to a sum of 1
res3 <- run_gldrm(I(NResp/ClusterSize) ~ Trt, data=shelltox, link="logit", weights=shelltox$Freq/sum(shelltox$Freq))

all.equal(res0[c("conv", "iter", "beta", "f0", "mu0", "llik")],
          res1[c("conv", "iter", "beta", "f0", "mu0", "llik")])


all.equal(res1[c("conv", "iter", "beta", "f0", "mu0", "llik")],
          res2[c("conv", "iter", "beta", "f0", "mu0", "llik")])

all.equal(res1[c("conv", "iter", "beta", "f0", "mu0", "llik")],
          c(res3[c("conv", "iter", "beta", "f0", "mu0")], llik = res3$llik * sum(shelltox$Freq)))


## Testing getTheta C-version
getThetaR <- function(spt, f0, mu, thetaStart=NULL, thetaControl=theta.control())
{
  ## Extract control arguments
  if (class(thetaControl) != "thetaControl")
    stop("thetaControl must be an object of class \'thetaControl\' returned by
             thetaControl() function.")
  logit <- thetaControl$logit
  eps <- thetaControl$eps
  maxiter <- thetaControl$maxiter
  maxhalf <- thetaControl$maxhalf
  maxtheta <- thetaControl$maxtheta

  ## Define value from inputs
  sptN <- length(spt)
  m <- min(spt)
  M <- max(spt)
  n <- length(mu)
  
  ## Format arguments
  spt <- as.vector(spt)
  f0 <- as.vector(f0)
  mu <- as.vector(mu)
  thetaStart <- as.vector(thetaStart)
  
  ## Value does not change
  gMu <- exchreg:::g(mu, m, M)
  
  ## Initialize values
  theta <- thetaStart  # initial values required
  thetaOld <- bPrimeErrOld <- rep(NA, n)
  conv <- rep(FALSE, n)
  maxedOut <- rep(FALSE, n)
  

    fUnstd <- f0 * exp(tcrossprod(spt, theta))  # |spt| x n matrix of tilted f0 values
    b <- colSums(fUnstd)
    fTilt <- fUnstd / rep(b, each=sptN)  # normalized
  
  bPrime <- colSums(spt*fTilt)  # mean as a function of theta
  bPrime2 <- colSums(outer(spt, bPrime, "-")^2 * fTilt)  # variance as a function of theta
  bPrimeErr <- bPrime - mu  # used to assess convergence
  
  
  ## Update theta until convergence
  conv <- (abs(bPrimeErr) < eps) | (theta==maxtheta & bPrimeErr<0) |
  (theta==-maxtheta & bPrimeErr>0)
  s <- which(!conv)
  iter <- 0
  while(length(s)>0 && iter<maxiter) {
    iter <- iter + 1
    bPrimeErrOld[s] <- bPrimeErr[s]  # used to assess convergence
     
    ## 1) Update theta
    thetaOld[s] <- theta[s]
    thetaS <- theta[s] - bPrimeErr[s] / bPrime2[s]
     thetaS[thetaS > maxtheta] <- maxtheta
     thetaS[thetaS < -maxtheta] <- -maxtheta
     theta[s] <- thetaS
  #   
  ## 2) Update fTilt, bPrime, and bPrime2 and take half steps if bPrimeErr not improved
   ss <- s
   nhalf <- 0
   while(length(ss)>0 && nhalf<maxhalf) {
     ## 2a) Update fTilt, bPrime, and bPrime2
      fUnstd[, ss] <- f0*exp(tcrossprod(spt, theta[ss]))  # |spt| x n matrix of tilted f0 values
      b[ss] <- colSums(fUnstd[, ss, drop=FALSE])
      fTilt[, ss] <- fUnstd[, ss, drop=FALSE] / rep(b[ss], each=sptN)  # normalized
   
     bPrime[ss] <- colSums(spt*fTilt[, ss, drop=FALSE])  # mean as a function of theta
     bPrime2[ss] <- colSums(outer(spt, bPrime[ss], "-")^2 * fTilt[, ss, drop=FALSE])  # variance as a function of theta
     bPrimeErr[ss] <- bPrime[ss] - mu[ss]  # used to assess convergence
       
     ## 2b) Take half steps if necessary
    ss <- ss[abs(bPrimeErr[ss]) > abs(bPrimeErrOld[ss])]
    if (length(ss) > 0) nhalf <- nhalf + 1
    theta[ss] <- (theta[ss] + thetaOld[ss]) / 2
   }
   ## If maximum half steps are exceeded, set theta to previous value
   maxedOut[ss] <- TRUE
   theta[ss] <- thetaOld[ss]
     
   ## 3) Check convergence
   conv[s] <- (abs(bPrimeErr[s]) < eps) |
     (theta[s]==maxtheta & bPrimeErr[s]<0) |
     (theta[s]==-maxtheta & bPrimeErr[s]>0)
   s <- s[!conv[s] & !maxedOut[s]]
  }
 
  list(theta=theta, fTilt=fTilt, bPrime=bPrime, bPrime2=bPrime2,
       conv=conv, iter=iter)
}

.spt <- 0:5
.f0 <- dbinom(0:5, size=5, prob=0.2)
.mu <- seq(2, 4, by=0.5)
.thetaStart <- rep(0, length(.mu))

  
Rcpp::sourceCpp("src/getTheta.cpp")
(a <- getThetaC(.spt, .f0, .mu, .thetaStart, exchreg:::theta.control()))
(b <- getThetaR(.spt, .f0, .mu, .thetaStart, exchreg:::theta.control()))

all.equal(lapply(a,c), lapply(b,c))
microbenchmark::microbenchmark(
  getThetaC(.spt, .f0, .mu, .thetaStart, exchreg:::theta.control()),
  getThetaR(.spt, .f0, .mu, .thetaStart, exchreg:::theta.control())
)



ab <- sprr(cbind(NResp, ClusterSize-NResp) ~ Trt, data=shelltox, weights=shelltox$Freq,
           control=list(eps=0.01,maxit=100), start=list(beta=c(0,-1,1,1)))

aa <- sprr(cbind(NResp, ClusterSize-NResp) ~ Trt, data=shelltox, 
          control=list(eps=0.01,maxit=100))
a0 <- sprr(cbind(NResp, ClusterSize-NResp) ~ 1, 
                 data=shelltox, weights=shelltox$Freq,
                 control=list(eps=0.01,maxit=100))
a1 <- sprr(cbind(NResp, ClusterSize-NResp) ~ as.numeric(Trt), 
                 data=shelltox, weights=shelltox$Freq,
                 control=list(eps=0.01,maxit=100))

data(egde)
b <- sprr(cbind(NResp, ClusterSize-NResp) ~ Trt+log(ClusterSize),
                data=egde, weights=egde$Freq,
                control=list(eps=0.01,maxit=100))

#
data("boric_acid")
ba.mod <- sprr(cbind(NResp, ClusterSize-NResp) ~ Trt,
          data=boric_acid, weights=boric_acid$Freq,link="log",
          control=list(eps=0.01,maxit=100))


ba.mod1 <- sprr(cbind(NResp, ClusterSize-NResp) ~ Trt + log(ClusterSize),
               data=boric_acid, weights=boric_acid$Freq, link="log",
               control=list(eps=0.01,maxit=100),
               start=list(q=ba.mod$q, beta=c(coef(ba.mod), 0),
                          mu1=lambda_from_p(ba.mod$q)[2]))
ba.mod2 <- sprr(cbind(ClusterSize-NResp, NResp) ~ Trt ,
                data=ba, weights=ba$Freq, link="log",
                control=list(eps=0.01,maxit=100))

zz <- mc.est(ba)
ba_zz <- merge(ba, zz, by=c("Trt","ClusterSize","NResp"), all.x = TRUE,
               all.y=FALSE)
with(ba_zz, sum(Freq * log(Prob)))

nd <- data.frame(Trt = unique(ba$Trt))
nd$p <- predict(ba.mod, newdata=nd)
nd

l <- predict(ba.mod, newdata=nd, newn=15, type="lvec")
matplot(t(l), type="l")

#########


# create native routine registration file
#tools::package_native_routine_registration_skeleton("z://RForge/corrbin/pkg/CorrBin",
#                                                    con="z://RForge/corrbin/pkg/CorrBin/src/init.c")

#document(cb)
#run_examples(cb)  # or dev_example("ran.CMData")
#load_all(cb)

#build(cb)
#check(cb,  check_dir = "z:/RForge", check_version = TRUE, cran = TRUE)


