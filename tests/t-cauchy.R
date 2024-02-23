require(OneStep)

n <- 1e3
set.seed(1)
x <- rcauchy(n)
onestep(x, "cauchy")
