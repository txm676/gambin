## ------------------------------------------------------------------------
library("gambin")
data(moths, package="gambin")

##unimodal model
fit = fit_abundances(moths)
barplot(fit)
points(fit)
AIC(fit)

##unimodal model (fit to a subsample of 1000 individuals)
fit2 = fit_abundances(moths, subsample = 1000)
barplot(fit2)
points(fit2)
AIC(fit2)

##bimodal model (using 3 cores)
fit3 = fit_abundances(moths, no_of_components = 2, cores = 3)
barplot(fit3)
points(fit3)
AIC(fit3)

