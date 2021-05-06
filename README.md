# dentist
R package to sample points around a specified distance from the maximum likelihood estimates. This should be a better way to estimate uncertainty than using the Hessian of the likelihood equation. It works by "denting" the likelihood surface to make a ridge at your desired ∆lnL and then "walks" around this dented surface, sampling points.
