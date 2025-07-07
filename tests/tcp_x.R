## this will render the output independent from the version of the package
suppressPackageStartupMessages(library(rrcov3way))

## IGNORE_RDIFF_BEGIN
##
##  ...
##
## IGNORE_RDIFF_END

##  cp_gen() =================================================
try(xdat <- cp_gen(nsim=1, eps=0, type="bl"))           # wrong eps-type combination
try(xdat <- cp_gen(nsim=1, eps=0.2, type="none"))       # wrong eps-type combination
try(xdat <- cp_gen(nsim=1, eps=0.2, type="xx"))         # wrong type

xdat <- cp_gen(I=50, J=100, K=10, nsim=1, nf=2,
     noise=0.15, noise1=0.10, Acol=TRUE, Bcol=TRUE, Ccol=TRUE,
     congA=0.5, congB=0.5, congC=0.5,
     eps=0.2, type="bl")
names(xdat)

xdat <- cp_gen(nsim=1, eps=0, type="none")
xdat <- cp_gen(nsim=1, eps=0.2, type="gl")
xdat <- cp_gen(nsim=1, eps=0.2, type="og")
xdat <- cp_gen(nsim=1, eps=0.2, type="og", c1=2, c2=0.1)
xdat <- cp_gen(nsim=1, eps=0.2, type="gl", c1=1, c2=0.2)




