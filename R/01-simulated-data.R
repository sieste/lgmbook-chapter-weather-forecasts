#!/usr/bin/env R

library(tidyverse)
library(Matrix)
library(mvtnorm)

source('00-functions.R')

# specify grid dimension
nlon = 80
nlat = 100
nt   = 20


set.seed(12345)

# simulate true values of regression parameters (latent variables) 
nmod_lon   = 4
nmod_lat   = 4
alpha_true = simulate_spatial_noise(nlat, nlon, nmod_lat, nmod_lon, c(10, 20))
beta_true  = simulate_spatial_noise(nlat, nlon, nmod_lat, nmod_lon, c(-0.5, 2))
tau_true   = simulate_spatial_noise(nlat, nlon, nmod_lat, nmod_lon, c(-1, 1))

# stacked into one vector
x_true = c(alpha_true, beta_true, tau_true)
n_x  = length(x_true)

# simulate forecast 2d fields
fcst = array(dim=c(nlat, nlon, nt), dimnames=list(x=1:nlat, y=1:nlon, t=1:nt))
for (tt in seq_len(nt)) {
  fcst[ , , tt] = simulate_spatial_noise(nlat, nlon, nmod_lon, nmod_lat, c(10, 20))
}

# format as lat/lon/time data frame and calculate local forecast means
data_true = expand.grid(lat = 1:nlat, lon=1:nlon, t=1:nt) %>% 
  dplyr::mutate(alpha_true = alpha_true[cbind(lat, lon)],
                beta_true  = beta_true[cbind(lat, lon)],
                tau_true   = tau_true[cbind(lat, lon)],
                fcst       = fcst[cbind(lat, lon, t)]) %>%
  dplyr::group_by(lon, lat) %>%
  dplyr::mutate(fcst_mean = mean(fcst)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(fcst_anom = fcst - fcst_mean) %>%
  # ... simulate observations
  dplyr::mutate(sigma_true = exp(tau_true),
                epsilon    = rnorm(n(), 0, sigma_true),
                obs        = alpha_true + beta_true * fcst_anom + epsilon)
                        
# separate out the data that is used for inference
data = data_true %>% 
  dplyr::select(lon, lat, fcst_anom, obs)

# summary measures of data set required for inference
data_summary = data %>%
  group_by(lon, lat) %>%
  mutate(obs_anom = obs - mean(obs)) %>%
  summarise(.groups='drop',
            vyy = mean(obs_anom^2),
            vzz = mean(fcst_anom^2),
            vyz = mean(obs_anom * fcst_anom),
            ybar = mean(obs))



################################################################################
# calculate linear regression mles per grid point

lm_data = data %>% 
  dplyr::group_by(lon, lat) %>%
  tidyr::nest(hindcast=c(fcst_anom, obs)) %>%
  dplyr::mutate(lm = purrr::map(hindcast, ~ lm(obs ~ fcst_anom, data = .)))

# extract the MLEs, and also sxx and n to calculate the Hessian
mle_data = lm_data %>%
  dplyr::mutate(alpha_hat = purrr::map_dbl(lm, ~ coef(.)[1]),
                beta_hat = purrr::map_dbl(lm, ~ coef(.)[2]),
                sigma_hat = purrr::map_dbl(lm, ~ sqrt(mean(resid(.)^2))),
                tau_hat = log(sigma_hat),
                sxx = purrr::map_dbl(lm, ~ mean(.$model$fcst_anom^2)),
                n = purrr::map_int(lm, ~ length(.$residuals)),
                var_alpha_hat = exp(2 * tau_hat) / n,
                var_beta_hat = exp(2 * tau_hat) / (n * sxx),
                var_tau_hat = 1 / (2 * n))

# vector of MLEs and negative Hessian 
x_hat = with(mle_data, c(alpha_hat, beta_hat, tau_hat))
nH_diag = with(mle_data, c(n * exp(-2 * tau_hat),
                           n * sxx * exp(-2 * tau_hat),
                           2 * n))
nH_y = sparseMatrix(i=1:length(nH_diag), j=1:length(nH_diag), x=nH_diag)


################################################################################
# Spatial inference 
################################################################################

# 2d random walk precision matrices with fixed hyperparameters
RR = make_rw2d_structure_matrix(nlon, nlat)
kappa = c(alpha = 8, beta = 150, tau = 150)
Q_x = bdiag(list(kappa['alpha'] * RR, 
                 kappa['beta'] * RR, 
                 kappa['tau'] * RR))


# posterior mean 
cQ_xy = Cholesky(Q_x + nH_y, LDL=FALSE) # Cholesky of posterior precision
x_star = drop(solve(cQ_xy, nH_y %*% x_hat))

# conditional posterior samples
n_sampls = 500
zz = solve(cQ_xy, matrix(rnorm(n_x * n_sampls), n_x, n_sampls), system='Lt')
postsampls_x = drop(x_star + solve(cQ_xy, zz, system='Pt'))
postvar_x = apply(postsampls_x, 1, 'var')


################################################################################
# infer kappa
################################################################################

params = list(SS=nlon * nlat, TT=nt, 
              vyy=data_summary$vyy, vzz=data_summary$vzz, 
              vyz=data_summary$vyz, ybar=data_summary$ybar, 
              RR=RR, nH_y=nH_y, x_hat=x_hat)

# find posterior mode of kappa with optim
opt_lkappa = optim(par = c(alpha=0, beta=0, tau=0),
                   fn = lp_kappa_I_y, params=params,
                   control=list(fnscale=-1))

# find hessian of the log posterior at the mode for laplace approximation
hess_lkappa = optimHess(par=opt_lkappa[['par']], fn = lp_kappa_I_y, params=params)


################################################################################
# approximate marginal posterior of latents
################################################################################

# draw log kappas from the normal approximation of the posterior
n_sampls = 100
mu_lkappa = opt_lkappa$par
sigma_lkappa = solve(-hess_lkappa)

# Draw an x samples from the conditional p(x|y,kappa) for each kappa sample
# (NOTE: Using makeCluster(7), clusterExport(...), parLapply(...) can reduce
# the time from 2 minutes to about 1 minute. Might be useful when drawing more
# samples.)
x_sampls = vapply(
  X=1:n_sampls, 
  FUN.VALUE = numeric(n_x),
  FUN = function(k) {
    lkappa = mvtnorm::rmvnorm(1, mean=mu_lkappa, sigma=sigma_lkappa)
    kappa = exp(drop(lkappa))
    Q_x = Matrix::bdiag(list(kappa['alpha'] * RR, 
                             kappa['beta'] * RR, 
                             kappa['tau'] * RR))
    cQ_xy = Matrix::Cholesky(Q_x + nH_y, LDL=FALSE) 
    x_star = drop(Matrix::solve(cQ_xy, nH_y %*% x_hat))
    zk = Matrix::solve(cQ_xy, rnorm(n_x), system='Lt')
    xk = drop(x_star + Matrix::solve(cQ_xy, zk, system='Pt'))
    return(xk)
})

postmean_x = rowMeans(x_sampls)
postsd_x = apply(x_sampls, 1, sd)
post_data = mle_data %>% ungroup %>% select(lon, lat) %>%
  mutate(alpha_star = postmean_x[1:SS],
         beta_star = postmean_x[1:SS + SS],
         tau_star = postmean_x[1:SS + 2*SS])


################################################################################
# generate testing data (forecasts and observations
################################################################################
# simulate forecast 2d field
fcst_test = simulate_spatial_noise(nlat, nlon, nmod_lon, nmod_lat, c(10, 20))

# get forecast means from training data to be able to calculate forecast
# anomalies in test data
fcst_mean_train = data_true %>% 
  group_by(lat, lon) %>% 
  summarise(.groups='drop', fcst_mean = fcst_mean[1])

# simulate test data: assemble true regression parameter values and test
# forecasts, and simulate observations
data_test = 
  expand.grid(lat = 1:nlat, lon=1:nlon) %>% 
  as_tibble() %>%
  dplyr::mutate(alpha_true = alpha_true[cbind(lat, lon)],
                beta_true  = beta_true[cbind(lat, lon)],
                tau_true   = tau_true[cbind(lat, lon)],
                fcst       = fcst_test[cbind(lat, lon)]) %>%
  dplyr::right_join(fcst_mean_train, by=c('lat', 'lon')) %>%
  dplyr::mutate(fcst_anom = fcst - fcst_mean) %>%
  # ... simulate observations
  dplyr::mutate(sigma_true = exp(tau_true),
                epsilon    = rnorm(n(), 0, sigma_true),
                obs        = alpha_true + beta_true * fcst_anom + epsilon)


# add post-processed forecasts using MLEs of regression parameters
data_test = 
  data_test %>% 
  right_join(select(mle_data, lon, lat, alpha_hat, beta_hat), by=c('lat', 'lon')) %>%
  mutate(fcst_hat = alpha_hat + beta_hat * fcst_anom)

# add post-processed forecasts using posterior means of regression parameters
data_test =
  data_test %>%
  right_join(post_data, by=c('lon', 'lat')) %>%
  mutate(fcst_star = alpha_star + beta_star * fcst_anom)

# add "perfectly post-processed" forecasts, using the true parameter values
data_test =
  data_test %>%
  mutate(fcst_perf = alpha_true + beta_true * fcst_anom)

# add the climatological forecast
data_test = 
  data_test %>%
  right_join(select(data_summary, lon, lat, ybar), by=c('lon', 'lat')) %>%
  rename(fcst_clim = ybar)


# calculate "global" mse 
mse_test = data_test %>% summarise(
  mse_hat = mean((obs-fcst_hat)^2),
  mse_star = mean((obs-fcst_star)^2),
  mse_perf = mean((obs - fcst_perf)^2),
  mse_clim = mean((obs - fcst_clim)^2))

# Conclusion: postprocessed forecasts much better for smoothed than for
# non-smoothed parameters, MSE almost same as when using the true parameters


# save
save(file = '01-simulated-data.Rdata', 
     list=c('data_true', 'data', 'data_summary', 'nlon', 'nlat', 'nt', 'n_x'))


# bring data into convenient format for plotting
plt_data_true = 
  tidyr::crossing(variable = c('alpha', 'beta', 'tau'),
                  lon = 1:nlon, 
                  lat = 1:nlat) %>%
  dplyr::mutate(true = x_true)

plt_data_lik = mle_data %>%
  dplyr::select(-hindcast, -lm, -sxx, -n, -sigma_hat) %>%
  tidyr::gather('variable', 'value', -lon, -lat) %>%
  dplyr::mutate(what = ifelse(stringr::str_detect(variable, '^var_'), 
                              'variance', 'estimate')) %>%
  dplyr::mutate(variable = stringr::str_replace_all(variable, 
                                                    '^var_|_hat$', '')) %>%
  tidyr::spread(what, value) %>%
  dplyr::mutate(method = 'likelihood')

plt_data_post = 
  tidyr::crossing(variable = c('alpha', 'beta', 'tau'),
                  lon = 1:nlon, 
                  lat = 1:nlat) %>%
  dplyr::mutate(estimate = x_star,
                variance = postvar_x,
                method = 'posterior')

plt_data_all = bind_rows(plt_data_lik, plt_data_post) %>%
               right_join(plt_data_true, by=c('variable', 'lat', 'lon'))


################################################################################
# some summary plots (figure files produced at the end) 
################################################################################

# scatter plot of estimate +/- 2sd vs true value
plt_scatter = list()
for (var in c('alpha', 'beta', 'tau')) {
  plt_scatter[[var]] = 
    ggplot(plt_data_all %>% filter(variable == var)) + 
      geom_point(aes(x=true, y=estimate)) + 
      geom_segment(aes(x=true, xend = true, 
                       y = estimate - 2*sqrt(variance), 
                       yend = estimate + 2*sqrt(variance))) +
      facet_grid(~method, scales='free') +
      geom_abline(intercept=0, slope=1, col='white', lty=2)
}

# spatial plots of mean
plt_spatial = list()
for (var in c('alpha', 'beta', 'tau')) {
  plt_spatial[[var]] = 
    ggplot(plt_data_all %>% filter(variable == var)) + 
      geom_raster(aes(x=lon, y=lat, fill=estimate)) + 
      scale_fill_gradientn(colours = hcl.colors(10)) +
      facet_grid(~method) +
      labs(x=NULL, y=NULL, fill=var)
}

# posterior standard deviation with vs without spatial prior
plt_stdev = list()
for (var in c('alpha', 'beta', 'tau')) {

  df_ = plt_data_all %>% 
    filter(variable == var) %>%
    select(lon, lat, method, variance) %>%
    spread(method, variance)

  plt_stdev[[var]] =
    ggplot(df_) +
    geom_point(aes(x=sqrt(likelihood), y=sqrt(posterior))) +
    geom_abline(intercept = 0, slope = 1)

}

outdir = 'fig/01-simulated-data'
for (var in c('alpha', 'beta', 'tau')) {

  ggsave(file.path(outdir, paste('compare-spatial-est-', var, '.pdf', sep='')), 
       plt_spatial[[var]],
       width=7, height=3)

  ggsave(file.path(outdir, paste('compare-scatter-', var, '.pdf', sep='')),
         plt_scatter[[var]],
         width=7, height=3)
}


# mse with and without spatial prior
mses = plt_data_all %>% 
  ungroup %>%
  group_by(method, variable) %>%
  summarise(.groups='drop',
            mse = mean((estimate - true)^2)) %>%
  arrange(variable)





