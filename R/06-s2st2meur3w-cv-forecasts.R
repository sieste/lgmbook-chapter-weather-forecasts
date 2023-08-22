library(tidyverse)
library(Matrix)
load('../data/hindcast-eur3w.Rdata')
load('../data/s2st2meur3w-crossvalidation-results.Rdata')
source('00-functions.R')

# calculate crossvalidation predictive samples
eur3w = hindcast_eur3w
eur3w_loo = s2st2meur3w_cv_result
years = as.numeric(names(eur3w_loo))
s2st2meur3w_cv_forecasts = list()

# loop over left-out years
for (yr in years) {

  # calculate covariate (forecast in yr minus mean over forecasts w/o yr) in
  # correct order
  fcst_mean = hindcast_eur3w %>%
              filter(year != yr) %>%
              group_by(long, lat) %>%
              summarise(fcst_mean = mean(fcst), .groups='drop') %>%
              arrange(long, desc(lat)) %>%
              pull(fcst_mean)
              
  fcst_yr = hindcast_eur3w %>% 
         filter(year == yr) %>% 
         arrange(long, desc(lat)) %>%
         pull(fcst)

  z_yr = fcst_yr - fcst_mean

  # extract left-out observations
  y_yr = hindcast_eur3w %>%
         filter(year == yr) %>%
         arrange(long, desc(lat)) %>%
         pull(obs)

  # extract posterior parameter samples
  x_sampls_post = eur3w_loo[[paste(yr)]][['xstar_sampls']]

  # dimension constants
  n_s = length(y_yr)
  n_x = 3 * n_s
  n_sampls = ncol(x_sampls_post)
  
  # generate parameter samples without spatial prior
  x_hat_yr = eur3w_loo[[paste(yr)]][['x_hat']]
  x_hat_sd_yr = eur3w_loo[[paste(yr)]][['x_hat_sd']]
  set.seed(yr)
  x_sampls_mle = replicate(n_sampls, rnorm(n_x, x_hat_yr, x_hat_sd_yr))

  # generate forecasts from spatial model for this year
  alpha_star = x_sampls_post[1:n_s, ]
  beta_star = x_sampls_post[1:n_s + n_s, ]
  tau_star = x_sampls_post[1:n_s + 2*n_s, ]
  y_star = vapply(X = 1:n_sampls, 
                  FUN.VALUE = numeric(n_s),
                  FUN = function(ii) {
                    mu = alpha_star[, ii] + beta_star[, ii] * z_yr
                    sigma = exp(0.5 * tau_star[, ii])
                    rnorm(n = n_s, mean = mu, sd = sigma)
                  })

  # generate forecasts with mles for this year
  alpha_hat = x_sampls_mle[1:n_s, ]
  beta_hat = x_sampls_mle[1:n_s + n_s, ]
  tau_hat = x_sampls_mle[1:n_s + 2*n_s, ]
  y_hat = vapply(X = 1:n_sampls, 
                  FUN.VALUE = numeric(n_s),
                  FUN = function(ii) {
                    mu = alpha_hat[, ii] + beta_hat[, ii] * z_yr
                    sigma = exp(0.5 * tau_hat[, ii])
                    rnorm(n = n_s, mean = mu, sd = sigma)
                  })

  # save forecasts and observations

  s2st2meur3w_cv_forecasts[[paste(yr)]][['y_star']] = y_star
  s2st2meur3w_cv_forecasts[[paste(yr)]][['y_hat']] = y_hat
  s2st2meur3w_cv_forecasts[[paste(yr)]][['y']] = y_yr
 
}

# save coordinate grid
grid = hindcast_eur3w %>% 
       filter(year == years[1]) %>% 
       arrange(long, desc(lat)) %>%
       select(long, lat)

s2st2meur3w_cv_forecasts[['grid']] = grid

save(file = '../data/s2st2meur3w-crossvalidation-forecasts.Rdata', 
     list = 's2st2meur3w_cv_forecasts')

