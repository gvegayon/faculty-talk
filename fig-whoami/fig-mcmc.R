library(coda)
# dat <- readRDS("fig-mcmc/mcmc_partially_annotated_no_prior.rds")

library(fmcmc)
data(logit, package = "mcmc")
out <- glm(y ~ x1 + x2 + x3 + x4, data = logit, family = binomial, x = TRUE)
beta.init <- as.numeric(coefficients(out))
lupost_factory <- function(x, y) function(beta) {
  eta <- as.numeric(x %*% beta)
  logp <- ifelse(eta < 0, eta - log1p(exp(eta)), - log1p(exp(- eta)))
  logq <- ifelse(eta < 0, - log1p(exp(eta)), - eta - log1p(exp(- eta)))
  logl <- sum(logp[y == 1]) + sum(logq[y == 0])
  return(logl - sum(beta^2) / 8)
}

lupost <- lupost_factory(out$x, out$y)
khaario <- kernel_adapt(freq = 1, warmup = 2000)

set.seed(12)

out_harrio_1 <- MCMC(
  initial   = rbind(beta.init, beta.init + rnorm(5))[1,],
  fun       = lupost,
  nsteps    = 6000,    # We will only run the chain for 100 steps                    
  kernel    = khaario, # We passed the predefined kernel
  thin      = 1,       # No thining here
  nchains   = 1L,      # A single chain
  multicore = FALSE    # Running in serial
)

graphics.off()
svg("fig-mcmc/fig-mcmc.svg", bg = "transparent")
traceplot(out_harrio_1[,5], main = "Trace of Adaptive Transition Kernel")
abline(v = 2000, col = "red", lwd = 2, lty=2)
dev.off()