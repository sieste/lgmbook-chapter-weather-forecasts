library(tidyverse)
library(Matrix)
library(INLA)
load('../data/hindcast-eur.Rdata')
source('00-functions.R')

hindcast = hindcast_eur %>%
  mutate(year = as.integer(substr(date, 1, 4))) %>%
  select(-ens_sd) %>%
  rename(fcst = ens_mean)

N_i = hindcast %>% pull(lat) %>% unique %>% length
N_j = hindcast %>% pull(long) %>% unique %>% length
years = hindcast %>% pull(year) %>% unique %>% sort
lead_times = hindcast %>% pull(lead_time) %>% unique %>% sort

cv_res = NULL

for (yr_loo in years) {

  for (lt in lead_times) {
    
    cat('year: ', yr_loo, '   lead time: ', lt, '\n')
  
    hindcast_train = hindcast %>% 
      filter(year != yr_loo, lead_time == lt)


    # calculate mles and asymptotic variances per grid point and lead time
    df_mle = hindcast_train %>%
      group_by(long, lat) %>%
      mutate(obs_anom = obs - mean(obs)) %>%
      mutate(fcst_anom = fcst - mean(fcst)) %>%
      summarise(.groups='drop',
                vyy = mean(obs_anom^2),
                vzz = mean(fcst_anom^2),
                vyz = mean(obs_anom * fcst_anom),
                ybar = mean(obs),
                alpha_hat = ybar,
                beta_hat = vyz / vzz,
                tau_hat = 0.5*log(mean((obs - alpha_hat - beta_hat * fcst_anom)^2)),
                n = n(),
                var_alpha_hat = exp(2 * tau_hat) / n,
                var_beta_hat = exp(2 * tau_hat) / (n * vzz),
                var_tau_hat = 1 / (2 * n)) %>%
      select(-vyy, -vzz, -vyz, -ybar) %>%
      mutate(i = as.numeric(factor(lat)),
             j = as.numeric(factor(long)),
             s = (j - 1) * N_i + i) %>%
      arrange(s)

    inla_alpha = inla(
      formula = alpha_hat ~ -1 + f(s, model='rw2d', nrow=N_i, ncol=N_j, constr=FALSE),
      data = df_mle, family='gaussian', 
      control.predictor=list(compute=TRUE),
      control.family=list(hyper=list(prec=list(initial=0, fixed=TRUE))),
      scale=1/df_mle$var_alpha_hat)
  
    inla_beta = inla(
      formula = beta_hat ~ -1 + f(s, model='rw2d', nrow=N_i, ncol=N_j, constr=FALSE),
      data = df_mle, family='gaussian', 
      control.predictor=list(compute=TRUE),
      control.family=list(hyper=list(prec=list(initial=0, fixed=TRUE))),
      scale=1/df_mle$var_beta_hat)
  
    inla_tau = inla(
      formula = tau_hat ~ -1 + f(s, model='rw2d', nrow=N_i, ncol=N_j, constr=FALSE),
      data = df_mle, family='gaussian', 
      control.predictor=list(compute=TRUE),
      control.family=list(hyper=list(prec=list(initial=0, fixed=TRUE))),
      scale=1/df_mle$var_tau_hat)
    
    df_mle = df_mle %>% 
      select(-i, -j, -s) %>%
      mutate(year = yr_loo, 
             lead_time = lt) %>%
      mutate(alpha_smooth = inla_alpha$summary.fitted.values$mean,
             alpha_smooth_var = inla_alpha$summary.fitted.values$sd^2,
             beta_smooth = inla_beta$summary.fitted.values$mean,
             beta_smooth_var = inla_beta$summary.fitted.values$sd^2,
             tau_smooth = inla_tau$summary.fitted.values$mean,
             tau_smooth_var = inla_tau$summary.fitted.values$sd^2)
  
    if (is.null(cv_res)) {
      cv_res = df_mle
    } else {
      cv_res = bind_rows(cv_res, df_mle)
    }
  }
}

s2st2meur_cv_result = cv_res

save(file='../data/s2st2meur-crossvalidation-results.Rdata', list='s2st2meur_cv_result')


