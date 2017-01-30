library("gambin")

ll_w = function(alpha, x, maxoct) {
  if(any(alpha < 0)|| alpha[3] > 1) return(Inf)
  dgamb <- alpha[3]*dgambin(x$Var1, alpha[1], maxoct, log = FALSE) +
    (1-alpha[3])*dgambin(x$Var1, alpha[2], maxoct, log = FALSE)
  -sum(x$Freq * log(dgamb))
}

set.seed(1)
n = 1e4
x1 = rgambin(n, 4, 10)
x2 = rgambin(n, 1, 10)
x = as.data.frame(table(c(x1, x2)), stringsAsFactors = FALSE)
x$Var1 = as.numeric(x$Var1)
optim(c(1, 1, 0.1), ll_w, x=x,  maxoct=10) # Set away from the truth
# Parameter estimates are
# 0.8341 3.6372 0.4108
# Looks OK as these parameter estimates have a lower 
# minumum than the true values
ll_w(c(0.8341,3.6372,0.4108), x, 10)
ll_w(c(1,4,0.5), x, 10)


### Let's reduce the sample size: n=200 per component
set.seed(2)
n = 500
x1 = rgambin(n, 4, 10)
x2 = rgambin(n, 1, 10)
x = as.data.frame(table(c(x1, x2)), stringsAsFactors = FALSE)
x$Var1 = as.numeric(x$Var1)
(est = optim(c(1, 1, 0.1), ll_w, x=x,  maxoct=10))
# Parameter estimates are
#2.386e-05 2.914e+00 1.454e-01
# Compare to the truth - again true minimum
ll_w(est$par, x, 10)
ll_w(c(1,4,0.5), x, 10)

## Let's simulate the from the new parameter values
x1_sim = rgambin(ceiling(2*n*est$par[3]), est$par[1], 10)
x2_sim = rgambin(2*n * (1-est$par[3]), est$par[2], 10)
par(mfrow=c(1, 2))
hist(c(x2, x1))
hist(c(x2_sim, x1_sim))

