library(devtools)
<<<<<<< HEAD:nuweb/SPreg_develop.R
source('nuweb/Nuweb.R')
=======
# create("../exchreg")
source('z:/RForge/Nuweb.R')
#source('/home/aszabo/RForge/Nuweb_linux.R')
>>>>>>> e91f19d9246eeff8aac80ff853ad9aaf61539ae9:SPreg_develop.R

#  setting up the package and its infrastructure
create("../exchreg")
use_package("CorrBin")
use_build_ignore("nuweb")


ex <- as.package("../exchreg")
nuweb(ex)

<<<<<<< HEAD:nuweb/SPreg_develop.R
shell("cd c:/exchreg/nuweb/ && texify --pdf --quiet --run-viewer  SPregress.tex")
shell("cd c:/exchreg/nuweb/ && texify --pdf --quiet --run-viewer  SPGLM.tex")
shell("cd c:/exchreg/nuweb/ && texify --pdf --quiet --run-viewer  wgldrm.tex")


document(ex)
load_all(ex)
=======
shell("cd c:/exchreg/ && texify --pdf --quiet --run-viewer --clean SPregress.tex")
shell("cd c:/exchreg/ && texify --pdf --quiet --run-viewer --clean SPGLM.tex")
>>>>>>> e91f19d9246eeff8aac80ff853ad9aaf61539ae9:SPreg_develop.R

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

######### Testing SPreg


a <- sprr(cbind(NResp, ClusterSize-NResp) ~ Trt, data=shelltox, weights=shelltox$Freq,
                control=list(eps=0.001, maxit=1000), 
          start = list(mu1=0.4))

nd <- data.frame(Trt = unique(shelltox$Trt))
nd$p <- predict(a, newdata=nd)
nd

a2 <- sprr(cbind(NResp, ClusterSize-NResp) ~ Trt + log(ClusterSize), 
           data=shelltox, weights=shelltox$Freq,
           control=list(eps=0.001,maxit=100),
           start = list(beta=c(coef(a),0), q=a$q))


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
ba <- read.CBData("z:/EOGeorge/Data/Binary/BoricAcidMousedata_processed.csv",
                  sep=",", skip=1)
ba.mod <- sprr(cbind(NResp, ClusterSize-NResp) ~ Trt,
          data=ba, weights=ba$Freq,link="log",
          control=list(eps=0.01,maxit=100))


ba.mod1 <- sprr(cbind(NResp, ClusterSize-NResp) ~ Trt + log(ClusterSize),
               data=ba, weights=ba$Freq, link="log",
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


