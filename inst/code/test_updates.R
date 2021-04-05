---
title: Test ED updates with general V
author: Yunqi Yang
date: 3/15/2021
output: html_document
---

```{r simulate-data}
set.seed(1)
n <- 500
V <- rbind(c(0.8,0.2),
           c(0.2,1.5))
U <- list(none   = rbind(c(0,0),
                         c(0,0)),
          shared = rbind(c(1.0,0.9),
                         c(0.9,1.0)),
          only1  = rbind(c(1,0),
                         c(0,0)),
          only2  = rbind(c(0,0),
                         c(0,1)))
w <- c(0.8,0.1,0.075,0.025)
X <- simulate_ud_data(n,w,U,V)
```

## Compare unconstrained updates for common V vs. general V

Here, I specify a common V and extend it to a list of Vs. fit1 and
fit2 both use `unconstrained.update = "ed"`. In this case, they are
doing the same update with the same Vs. So results are expected to be
exactly the same.
 
```{r ed-update-1}
set.seed(1)
fit1 <- ud_init(X,V = rep(list(V),n))
fit1 <- ud_fit(fit1,
               control = list(unconstrained.update = "ed",maxiter = 10,
                              resid.update = "none",scaled.update = "none",
                              rank1.update = "none"))
```

```{r ed-update-2}
set.seed(1)
fit2 <- ud_init(X,V = V)
fit2 <- ud_fit(fit2,
               control = list(unconstrained.update = "ed",maxiter = 10,
                              resid.update = "none",scaled.update = "none",
                              rank1.update = "none"))
```

```{r check-ed-updates}
plot(fit1$progress$loglik,type = "l",col = "darkblue",lwd = 2)
lines(fit2$progress$loglik,col = "darkorange",lwd = 2,lty = "dashed")
range(simplify2array(fit1$U) - simplify2array(fit2$U))
```

## Compare rank1 updates for common V and general V

If `rank1.update = "ed"`, we are using David's ultimate deconvolution
algorithm.  If `rank1.update = 'teem'`, we are using TEEM rank1
update.

Those results are not expected to be the same.

```{r}
set.seed(1)
fit0 <- ud_init(X,V = rep(list(V),n))
fit1 <- ud_fit(fit0,control = list(unconstrained.update = "none", scaled.update = 'none', rank1.update = 'teem', minval = 0, maxiter = 5000 ))
```

```{r}
set.seed(1)
fit2 <- ud_init(X, V = V)
fit3 <- ud_fit(fit2, control = list(unconstrained.update = "none", scaled.update = 'none', rank1.update = 'teem', resid.update = 'none', minval = 0, maxiter = 5000 ))

```

```{r  }
plot(fit1$progress$loglik)
points(fit3$progress$loglik, col = 'blue')
```
