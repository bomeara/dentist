#' Sample points from along a ridge
#' This "dents" the likelihood surface by reflecting points better than a threshold back across the threshold (think of taking a hollow plastic model of a mountain and punching the top so it's a volcano). It then uses essentially a Metropolis-Hastings walk to wander around the new rim. It adjusts the proposal width so that it samples points around the desired likelihood. 
#' This is better than using the curvature at the maximum likelihood estimate since it can actually sample points in case the assumptions of the curvature method do not hold. It is better than varying one parameter at a time while holding others constant because that could miss ridges: if I am fitting 5=x+y, and get a point estimate of (3,2), the reality is that there are an infinite range of values of x and y that will sum to 5, but if I hold x constant it looks like y is estimated very precisely. Of course, one could just fully embrace the Metropolis-Hastings lifestyle and use a full Bayesian approach.
#' 
#' While running, it will display the current range of likelihoods in the desired range (by default, the best negative log likelihood + 2 negative log likelihood units) and the parameter values falling in that range. If things are working well, the range of values will stabilize during a search.
#' 
#' The algorithm tunes: if it is moving too far away from the desired likelihoods, it will decrease the proposal width; if it staying in areas better than the desired likelihood, it will increase the proposal width. It will also expand the proposal width for parameters where the extreme values still appear good enough to try to find out the full range for these values. 
#' 
#' In general, the idea of this is not to give you a pleasingly narrow range of possible values -- it is to try to find the actual uncertainty, including finding any ridges that would not be seen in univariate space.
#' @param par Starting parameter vector. Ideally giving the desired likelihood, but starting at the optimum is likely ok. If named, the vector names are used to label output parameters.
#' @param fn The likelihood function, assumed to return negative log likelihoods
#' @param best_neglnL The negative log likelihood at the optimum; other values will be greater than this.
#' @param delta How far from the optimal negative log likelihood to focus samples
#' @param nsteps How many steps to take in the analysis
#' @param print_freq Output progress every print_freq steps.
#' @param lower_bound Minimum parameter values to try. One for all or a vector of the length of par.
#' @param upper_bound Maximum parameter values to try. One for all or a vector of the length of par.
#' @param ... Other arguments to fn. 
#' @return A dentist object containing results, the data.frame of negative log likelihoods and the parameters associated with them, and acceptances, the vector of whether a proposed move was accepted each step.
#' @export
#' @examples
#' sims <- rnorm(100, mean=17)
#' possible_means <- seq(from=16, to=18, length.out=100)
#' 
#' # Make sure we have a function that takes in a parameters vector, other arguments if needed,
#' # and returns the negative log likelihood
#' dnorm_to_run <- function(par, sims) {
#'   return(-sum(dnorm(x=sims, mean=par, log=TRUE)))
#' }
#' 
#' optimized_results <- optimize(dnorm_to_run,interval=range(possible_means), sims=sims, maximum=FALSE)
#' best_par <- optimized_results$minimum
#' best_neglnL <- optimized_results$objective
#' 
#' dented_results <- dent_walk(par=best_par, fn=dnorm_to_run, best_neglnL=best_neglnL,  nsteps=1000, sims=sims)
#' plot(dented_results$results$mean, dented_results$results$neglnL)
dent_walk <- function(par, fn, best_neglnL, delta=2, nsteps=1000, print_freq=50, lower_bound=0, upper_bound=Inf, ...) {
  results <- data.frame(matrix(NA, nrow=nsteps+1, ncol=length(par)+1))
  results[1,] <- c(best_neglnL, par)
  sd_vector <- 0.1*par
  acceptances <- rep(NA, nsteps)
  for (rep_index in sequence(nsteps)) {
    old_params <- as.numeric(results[rep_index,-1])
    new_params <- dent_propose(old_params, lower_bound=lower_bound, upper_bound=upper_bound, sd=sd_vector) 
    new_neglnL <- fn(par=new_params, ...)
    
    new_dented_neglnL <- dent_likelihood(new_neglnL, best_neglnL, delta)
    old_dented_neglnL <- dent_likelihood(results[rep_index,][1], best_neglnL, delta)
    
    
    if(new_dented_neglnL<=old_dented_neglnL) {
      results[rep_index+1,] <- c(new_neglnL, new_params)
      acceptances[rep_index]<-TRUE
    } else {
      if(new_dented_neglnL-old_dented_neglnL < runif(1)) {
        results[rep_index+1,] <- c(new_neglnL, new_params)
        acceptances[rep_index]<-TRUE
      } else {
        results[rep_index+1,] <- results[rep_index,]
        acceptances[rep_index]<-FALSE
      }
    }
    if(rep_index%%100==0) { # adaptively change proposal width for all params at once
      acceptances_run <- tail(acceptances[!is.na(acceptances)],100)
      if(sum(acceptances_run)/length(acceptances_run) > 0.3) {
        sd_vector <- sd_vector * 1.5
      }
      if(sum(acceptances_run)/length(acceptances_run) < 0.1) {
        sd_vector <- sd_vector * 0.8
      }     
    }
    if(rep_index%%200==0) { # if we haven't found values of some parameters that are outside the CI, widen the search for just those
      good_enough_results_range <- apply(results[which(results[1:(rep_index+1),1]-min(results[1:(rep_index+1),1])<=2),], 2, range)[,-1]
      all_results_range <- apply(results[1:(rep_index+1),], 2, range)[,-1]
      not_past_bounds <- apply(all_results_range==good_enough_results_range, 2, any)
      if(any(not_past_bounds)) {
        sd_vector[not_past_bounds] <- sd_vector[not_past_bounds]*1.5
      }
    }
    if(rep_index%%print_freq==0) {
      print(paste("Done replicate",rep_index))
      print("CI of values")
      intermediate_results <- apply(results[which(results[1:(rep_index+1),1]-min(results[1:(rep_index+1),1])<=2),], 2, range)
      #colnames(intermediate_results) <- c("neglnL", names(par))
      print(intermediate_results)
    }
  }
  colnames(results) <- c("neglnL", names(par))
  final_results <- list(results=results, acceptances=acceptances)
  class(final_results) <- c("dentist", "list")
  return(final_results)
}

#' Propose new values
#' This proposes new values using a normal distribution centered on the original parameter values, with desired standard deviation. If any proposed values are outside the bounds, it will propose again.
#' @param old_params The original parameter values
#' @param lower_bound Minimum parameter values to try. One for all or a vector of the length of par.
#' @param upper_bound Maximum parameter values to try. One for all or a vector of the length of par.
#' @param sd Standard deviation to use for the proposals. One for all or a vector of the length of par.
#' @return A vector of the new parameter values
dent_propose <- function(old_params, lower_bound=0, upper_bound=Inf, sd=1) {
  new_params <- rnorm(length(old_params), old_params, sd)
  while(any(new_params<lower_bound) | any(new_params>upper_bound)) {
    new_params <- rnorm(length(old_params), old_params, sd)
  }
  return(new_params)
}

#' Dents the likelihood surface
#' This takes any values that are better (lower) than the desired negative log likelihood and reflects them across the best_neglnL + delta line, "denting" the likelihood surface.
#' @param neglnL The original negative log likelihood
#' @param best_neglnL The negative log likelihood at the optimum; other values will be greater than this.
#' @param delta How far from the optimal negative log likelihood to focus samples
#' @return The transformed negative log likelihood
#' @export
#' @examples
#' sims <- rnorm(100, mean=17)
#' possible_means <- seq(from=16, to=18, length.out=100)
#' results_normal <- rep(NA, length(possible_means))
#' results_dented <- results_normal
#' dnorm_to_run <- function(par, sims) {
#'   return(-sum(dnorm(x=sims, mean=par, log=TRUE)))
#' }
#' best_neglnL <- optimize(dnorm_to_run,interval=range(possible_means), sims=sims, maximum=FALSE)$objective
#' for (i in seq_along(possible_means)) {
#'   results_normal[i] <- dnorm_to_run(possible_means[i], sims=sims)
#'   results_dented[i] <- dent_likelihood(par=possible_means[i], fn=dnorm_to_run, best_neglnL=best_neglnL, delta=2, sims=sims)
#' }
#' 
#' plot(possible_means, results_normal)
#' lines(possible_means, results_dented)
dent_likelihood <- function(neglnL, best_neglnL, delta=2) {
  difference <- neglnL - (best_neglnL+delta)
  if(difference<0) {
    neglnL <- (best_neglnL+delta) - difference
  }
  return(neglnL)
}