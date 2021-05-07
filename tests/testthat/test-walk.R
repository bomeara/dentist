test_that("multivariate walking works", {
 sims <- rlnorm(100, meanlog=1, sdlog=3)
 
 dlnorm_to_run <- function(par, sims) {
   return(-sum(dlnorm(sims, meanlog=par[1], sdlog=par[2], log=TRUE)))
 }
 
 optimized_results <- optim(c(meanlog=.9, sdlog=2.9), dlnorm_to_run, sims=sims)
 best_par <- optimized_results$par
 best_neglnL <- optimized_results$value
 
 dented_results <- dent_walk(par=best_par, fn=dlnorm_to_run, best_neglnL=best_neglnL,  nsteps=1000, sims=sims)
 expect_equal(nrow(dented_results$results), 1001)
 expect_true(any(dented_results$acceptances))
})

test_that("univariate walking works", {
 sims <- rnorm(100, mean=17)
 possible_means <- seq(from=16, to=18, length.out=100) # for optimize
 
  dnorm_to_run <- function(par, sims) {
   return(-sum(dnorm(x=sims, mean=par, log=TRUE)))
 }
 
 optimized_results <- optimize(dnorm_to_run,interval=range(possible_means), sims=sims, maximum=FALSE)
 best_par <- optimized_results$minimum
 names(best_par) <- "mean"
 best_neglnL <- optimized_results$objective
 
 dented_results <- dent_walk(par=best_par, fn=dnorm_to_run, best_neglnL=best_neglnL,  nsteps=1000, sims=sims)

 expect_equal(nrow(dented_results$results), 1001)
 expect_true(any(dented_results$acceptances))
})

