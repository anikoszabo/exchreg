library(devtools)
source('z:/RForge/Nuweb.R')
#source('/home/aszabo/RForge/Nuweb_linux.R')

run_nuweb(file="SPregress.w", path=getwd())

source('R/SPreg.R', echo=FALSE)

data(shelltox)
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
tools::package_native_routine_registration_skeleton("z://RForge/corrbin/pkg/CorrBin",
                                                    con="z://RForge/corrbin/pkg/CorrBin/src/init.c")

document(cb)
run_examples(cb)  # or dev_example("ran.CMData")
load_all(cb)

build(cb)
check(cb,  check_dir = "z:/RForge", check_version = TRUE, cran = TRUE)


