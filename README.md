
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dentist

<!-- badges: start -->

[![R-CMD-check](https://github.com/bomeara/dentist/workflows/R-CMD-check/badge.svg)](https://github.com/bomeara/dentist/actions)
[![R-CMD-check](https://github.com/bomeara/dentist/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bomeara/dentist/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`dentist` is an R package to sample points around a specified distance
from the maximum likelihood estimates. This should be a better way to
estimate uncertainty than using the Hessian of the likelihood equation.
It works by “denting” the likelihood surface to make a ridge at your
desired delta lnL and then “walks” around this dented surface, sampling
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
#> 0.9035365 2.9155160
```

But how confident should we be? For familiar distributions like a
binomial distribution we can compute the confidence interval for a
parameter estimate. For less familiar ones, we can approximate it using
the curvature at the peak. Another approach is to vary one parameter at
a time while holding others at their maximum likelihood estimate until
the likelihood gets much worse (typically a cutoff of delta 2 log
likelihood units is used). That can underestimate the uncertainty, since
there could be a ridge in parameter space. A better approach would be to
try many values and return all those within some specified likelihood
bounds. This could be done with latin hypercubes, but it can be a large
part of parameter space to explore. A different approach would be to
focus on sampling points right at the boundary of “good enough.” This
package does this. If a likelihood surface is a peak, this package wants
to sample points around a specified height below the peak to fully
sample the uncertainty. It does it by “walking” with a
Metropolis-Hastings algorithm around a dented likelihood surface, so
that values better than the ideal threshold are reflected back. First,
to show how the original surface (left) is dented (right):

<img src="man/figures/README-awesomeness-1.png" width="100%" /><img src="man/figures/README-awesomeness-2.png" width="100%" />

Note that for plotting, it’s shown with the regular log likelihood:
typically values \<\< 0, and higher is better. In the package, it wants
you to use negative log likelihood, so it would be the mirror image of
these plots (lower better).

And now to sample around that ring on the right:

``` r
library(dentist)
 dented_results <- dent_walk(par=best_par, fn=dlnorm_to_run, best_neglnL=best_neglnL,  nsteps=1000, print_freq=250, sims=sims)
#> [1] "Calculating intervals at a confidence level of 95%"
#> [1] "Done replicate 250"
#> [1] "CI of values (the 105 replicates within 2.99573227355399 neglnL of the optimum)"
#>        neglnL   meanlog    sdlog
#> [1,] 339.3150 0.3174435 2.511672
#> [2,] 342.2822 1.5966083 3.461012
#> [1] "Rough volume of good region is 1.21436251821118"
#> [1] "Done replicate 500"
#> [1] "CI of values (the 180 replicates within 2.99573227355399 neglnL of the optimum)"
#>        neglnL   meanlog    sdlog
#> [1,] 339.3150 0.2776119 2.499837
#> [2,] 342.2822 1.6138289 3.495837
#> [1] "Rough volume of good region is 1.33087214415585"
#> [1] "Done replicate 750"
#> [1] "CI of values (the 255 replicates within 2.99573227355399 neglnL of the optimum)"
#>        neglnL  meanlog    sdlog
#> [1,] 339.3150 0.192684 2.499837
#> [2,] 342.2822 1.613829 3.495837
#> [1] "Rough volume of good region is 1.41546033655253"
#> [1] "Done replicate 1000"
#> [1] "CI of values (the 318 replicates within 2.99573227355399 neglnL of the optimum)"
#>        neglnL  meanlog    sdlog
#> [1,] 339.3150 0.192684 2.482288
#> [2,] 342.2822 1.613829 3.495837
#> [1] "Rough volume of good region is 1.44040005320574"
```

This generates information about the confidence:

``` r
print(dented_results)
#> This ran 1000 steps looking for all points within 2.99573227355399 negative log likelihood units of the best parameter values.
#> 
#> Parameters: 
#>                      meanlog     sdlog
#> best             0.903536538 2.9155160
#> lower.CI         0.192684020 2.4822875
#> upper.CI         1.613828879 3.4958366
#> lowest.examined  0.006333159 0.3362002
#> highest.examined 2.603547188 5.1997580
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
