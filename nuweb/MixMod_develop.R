library(devtools)
source('z:/RForge/Nuweb.R')
#source('/home/aszabo/RForge/Nuweb_linux.R')

run_nuweb(file="Clustering.w", path=getwd())

source('R/MixMod.R', echo=FALSE)


### initial setup
create_package("z:/exchreg")
use_package("CorrBin")
