## -----------------------------------------------------------------------------
knitr::opts_chunk$set(comment = "#",
                      results = "hold",
                      collapse = TRUE,
                      fig.align = "center",
                      warning = FALSE,
                      message = FALSE)

## -----------------------------------------------------------------------------
library(udr)

## -----------------------------------------------------------------------------

set.seed(1)
n <- 400
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


