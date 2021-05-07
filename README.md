
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dentist

<!-- badges: start -->

[![R-CMD-check](https://github.com/bomeara/dentist/workflows/R-CMD-check/badge.svg)](https://github.com/bomeara/dentist/actions)
<!-- badges: end -->

`dentist` is an R package to sample points around a specified distance
from the maximum likelihood estimates. This should be a better way to
estimate uncertainty than using the Hessian of the likelihood equation.
It works by “denting” the likelihood surface to make a ridge at your
desired ∆lnL and then “walks” around this dented surface, sampling
points.

<https://bomeara.github.io/dentist/> for a website

<https://github.com/bomeara/dentist> for the source code

## Installation

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("bomeara/dentist")
```

## Example

Imagine we had empirical data from some distribution:

``` r
sims <- stats::rlnorm(100, meanlog=1, sdlog=3)
```

We could write a function for the likelihood of the data and optimize
this function:

``` r
# Define the likelihood function

 dlnorm_to_run <- function(par, sims) {
   return(-sum(stats::dlnorm(sims, meanlog=par[1], sdlog=par[2], log=TRUE)))
 }
 
 # Optimize the model given the empirical data. We guess at the starting values
 optimized_results <- stats::optim(c(meanlog=.5, sdlog=1), dlnorm_to_run, sims=sims)
 best_par <- optimized_results$par
 best_neglnL <- optimized_results$value
```

That gives us a point estimate of the best values:

``` r
print(best_par)
#>   meanlog     sdlog 
#> 0.8496578 2.6874693
```

But how confident should we be? For familiar distributions like a
binomial distribution we can compute the confidence interval for a
parameter estimate. For less familiar ones, we can approximate it using
the curvature at the peak. Another approach is to vary one parameter at
a time while holding others at their maximum likelihood estimate until
the likelihood gets much worse (typically a cutoff of ∆2 log likelihood
units is used). That can underestimate the uncertainty, since there
could be a ridge in parameter space. A better approach would be to try
many values and return all those within some specified likelihood
bounds. This could be done with latin hypercubes, but it can be a large
part of parameter space to explore. A different approach would be to
focus on sampling points right at the boundary of “good enough.” This
package does this. If a likelihood surface is a peak, this package wants
to sample points around a specified height below the peak to fully
sample the uncertainty. It does it by “walking” with a
Metropolis-Hastings algorithm around a dented likelihood surface, so
that values better than the ideal threshold are reflected back. First,
to show how the original surface (left) is dented (right):

<img src="man/figures/README-awesomeness-1.png" width="100%" />

Note that for plotting, it’s shown with the regular log likelihood:
typically values &lt;&lt; 0, and higher is better. In the package, it
wants you to use negative log likelihood, so it would be the mirror
image of these plots (lower better).

And now to sample around that ring on the right:

``` r
library(dentist)
 dented_results <- dent_walk(par=best_par, fn=dlnorm_to_run, best_neglnL=best_neglnL,  nsteps=1000, print_freq=250, sims=sims)
#> [1] "Done replicate 250"
#> [1] "CI of values"
#>            X1       X2       X3
#> [1,] 325.7501 0.358857 2.367714
#> [2,] 327.7347 1.286969 3.106623
#> [1] "Done replicate 500"
#> [1] "CI of values"
#>            X1       X2       X3
#> [1,] 325.7501 0.358857 2.358058
#> [2,] 327.7347 1.286969 3.106623
#> [1] "Done replicate 750"
#> [1] "CI of values"
#>            X1       X2       X3
#> [1,] 325.7501 0.358857 2.358058
#> [2,] 327.7347 1.376533 3.106623
#> [1] "Done replicate 1000"
#> [1] "CI of values"
#>            X1        X2       X3
#> [1,] 325.7501 0.3272925 2.358058
#> [2,] 327.7365 1.3765331 3.106623
```

This generates information about the confidence:

``` r
print(dented_results)
#> This ran 1000 steps looking for all points within 2 negative log likelihood units of the best parameter values.
#> 
#> Parameters: 
#>                    meanlog    sdlog
#> best             0.8496578 2.687469
#> lower.CI         0.3272925 2.358058
#> upper.CI         1.3765331 3.106623
#> lowest.examined  0.1642158 2.267184
#> highest.examined 1.4858718 3.325918
```

And also has a way to visualize the results:

``` r
plot(dented_results)
```

<img src="man/figures/README-summarized_results-1.png" width="100%" />

You do not need to use this package for simple distributions, but in
phylogenetics, where programs like OUwie, corHMM, hisse, OUCH, and more
give a point estimate, this lets you get confidence for the parameters
if there’s a way to have a function where you pass in a set of parameter
values and get back a negative log likelihood.
