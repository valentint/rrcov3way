## this will render the output independent from the version of the package
suppressPackageStartupMessages(library(rrcov3way))

## VT::14.01.2020 - added ## IGNORE_RDIFF_BEGIN/END
##  when printing the A, B and C matrices - because this could cause
##  differences on some platform, e.g. see for example sign indeterminacy
##

## Example with the UNIDO Manufacturing value added data
data(va3way)
dim(va3way)

## Treat quickly and dirty the zeros in the data set (if any)
va3way[va3way==0] <- 0.001

(res <- Parafac(va3way))
print(res$fit)
## IGNORE_RDIFF_BEGIN
print(res)
print(res$A)
print(res$B)
print(res$C)
## IGNORE_RDIFF_END
print(res$rd)
print(res$cutoff.rd)

## IGNORE_RDIFF_BEGIN
(res.r <- Parafac(va3way, robust=TRUE, center=TRUE, scale=TRUE, trace=TRUE))
## IGNORE_RDIFF_END

print(res.r$fit)
## IGNORE_RDIFF_BEGIN
print(res.r$A)
print(res.r$B)
print(res.r$C)
## IGNORE_RDIFF_END
print(res.r$rd)
print(res$cutoff.rd)
print(res.r$Hset)
print(res.r$flag)

(res.c <- Parafac(va3way, coda.transform="ilr"))
print(res.c$fit)
## IGNORE_RDIFF_BEGIN
print(res.c$A)
print(res.c$B)
print(res.c$Bclr)
print(res.c$C)
## IGNORE_RDIFF_END
print(res.c$rd)
print(res$cutoff.rd)

(res.rc <- Parafac(va3way, robust=TRUE, coda.transform="ilr", center=TRUE,
    scale=TRUE))
print(res.rc$fit)
## IGNORE_RDIFF_BEGIN
print(res.rc$A)
print(res.rc$B)
print(res.rc$Bclr)
print(res.rc$C)
## IGNORE_RDIFF_END
print(res.rc$rd)
print(res$cutoff.rd)
print(res.rc$Hset)
print(res.rc$flag)

set.seed(1)
(res.int2 <- Parafac(va3way, center=TRUE, scale=TRUE, optim="int2"))

set.seed(1)
(res.atld <- Parafac(va3way, center=TRUE, scale=TRUE, optim="atld"))

##================================================================
##
##  Example with the TV data from ThreeWay: ratings of 15 American television
##  shows on 16 bipolar scales made by 30 students: 15x16x30 array
##  in A-mode, i.e. 15x450.
data(TV, package="ThreeWay")

## Transform to 3-way array, to pass to Parafac() and set the dimnames
tv <- toArray(TV[[1]], 30, 16, 15, mode="B")
dimnames(tv)[[1]] <- TV[[4]]
dimnames(tv)[[2]] <- TV[[2]]
dimnames(tv)[[3]] <- TV[[3]]

(tvcp <- Parafac(tv, 2))

## IGNORE_RDIFF_BEGIN
tvcp$A; tvcp$B; tvcp$C; tvcp$GA
## IGNORE_RDIFF_END

tvcp$cutoff.rd
sort(tvcp$rd)
tvcp$cutoff.sd
sort(tvcp$sd)

(rtvcp <- Parafac(tv, 2, robust=TRUE))
## IGNORE_RDIFF_BEGIN
rtvcp$A; rtvcp$B; rtvcp$C; rtvcp$GA
## IGNORE_RDIFF_END
rtvcp$cutoff.rd
sort(rtvcp$rd)
rtvcp$cutoff.sd
sort(rtvcp$sd)

## ===============================================================
##
##  Compositional data and robustness

data(ulabor)

(res0 <- Parafac(ulabor))

## IGNORE_RDIFF_BEGIN
res0$A; res0$B; res0$C; res0$GA
## IGNORE_RDIFF_END

res0$cutoff.rd
sort(res0$rd)
res0$cutoff.sd
sort(res0$sd)

(res <- Parafac(ulabor, robust=TRUE, coda.transform="ilr"))

## IGNORE_RDIFF_BEGIN
res$A; res$B; res$C; res$GA
## IGNORE_RDIFF_END

res$cutoff.rd
sort(res$rd)
res$cutoff.sd
sort(res$sd)

(res1 <- Parafac(ulabor, coda.transform="clr"))

## IGNORE_RDIFF_BEGIN
res1$A; res1$B; res1$C; res1$GA
## IGNORE_RDIFF_END

res1$cutoff.rd
sort(res1$rd)
res1$cutoff.sd
sort(res1$sd)
