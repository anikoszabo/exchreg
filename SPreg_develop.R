library(devtools)
source('z:/RForge/Nuweb.R')
source('/home/aszabo/RForge/Nuweb_linux.R')

run_nuweb(file="SPregress.w", path=getwd())

source('R/SPreg.R', echo=FALSE)

data(shelltox)
a <- sprr(cbind(NResp, ClusterSize-NResp) ~ Trt, data=shelltox, weights=shelltox$Freq,
                control=list(eps=0.01,maxit=100))

nd <- data.frame(Trt = unique(shelltox$Trt))
p <- predict(a, newdata=nd)

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



#########


# create native routine registration file
tools::package_native_routine_registration_skeleton("z://RForge/corrbin/pkg/CorrBin",
                                                    con="z://RForge/corrbin/pkg/CorrBin/src/init.c")

document(cb)
run_examples(cb)  # or dev_example("ran.CMData")
load_all(cb)

build(cb)
check(cb,  check_dir = "z:/RForge", check_version = TRUE, cran = TRUE)


