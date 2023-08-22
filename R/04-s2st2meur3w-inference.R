library(tidyverse)
library(Matrix)
library(ggpubr)
load('../data/hindcast-eur3w.Rdata')
source('00-functions.R')

eur3w = hindcast_eur3w

world = rnaturalearth::ne_coastline() %>% fortify
lat_range = range(eur3w$lat)
long_range = range(eur3w$long)




# calculate forecast/obs summary measures per grid point
eur3w_summary = eur3w %>%
  group_by(long, lat) %>%
  mutate(obs_anom = obs - mean(obs)) %>%
  summarise(.groups='drop',
            vyy = mean(obs_anom^2),
            vzz = mean(fcst_anom^2),
            vyz = mean(obs_anom * fcst_anom),
            ybar = mean(obs))
            
# fit linear regression models obs ~ fcst_anom per grid point
eur3w_lm_data = 
  eur3w %>% 
  dplyr::group_by(long, lat) %>%
  dplyr::select(-year, -fcst) %>%
  tidyr::nest(hindcast=c(fcst_anom, obs)) %>%
  dplyr::mutate(lm = purrr::map(hindcast, ~ lm(obs ~ fcst_anom, data = .)))

# extract the MLEs, and also sxx and n to calculate the Hessian
eur3w_mle_data = eur3w_lm_data %>%
  dplyr::mutate(alpha_hat = purrr::map_dbl(lm, ~ coef(.)[1]),
                beta_hat = purrr::map_dbl(lm, ~ coef(.)[2]),
                sigma_hat = purrr::map_dbl(lm, ~ sqrt(mean(resid(.)^2))),
                tau_hat = log(sigma_hat),
                sxx = purrr::map_dbl(lm, ~ mean(.$model$fcst_anom^2)),
                n = purrr::map_int(lm, ~ length(.$residuals)),
                var_alpha_hat = exp(2 * tau_hat) / n,
                var_beta_hat = exp(2 * tau_hat) / (n * sxx),
                var_tau_hat = 1 / (2 * n))

# arrange by long and lat such that the estimators are arranged as a vectorised
# 2d field (stacked columns) from upper left to lower right
eur3w_mle_data = eur3w_mle_data %>% arrange(long, desc(lat))

# extract vector of MLEs and negative Hessian 
x_hat = with(eur3w_mle_data, c(alpha_hat, beta_hat, tau_hat))
nH_diag = with(eur3w_mle_data, c(n * exp(-2 * tau_hat),
                                 n * sxx * exp(-2 * tau_hat),
                                 2 * n))
nH_y = sparseMatrix(i=1:length(nH_diag), j=1:length(nH_diag), x=nH_diag)

# dimension paramaters
n_lon = eur3w_mle_data %>% pull(long) %>% n_distinct
n_lat = eur3w_mle_data %>% pull(lat) %>% n_distinct
n_t = nrow(eur3w_lm_data$hindcast[[1]])
n_s = n_lon * n_lat
n_x = 3 * n_s

# 2d random walk structure matrix (Q = kappa * R)
RR = make_rw2d_structure_matrix(n_lon, n_lat)

# # log prior for precision kappa is exponential for sigma = 1/sqrt(kappa) with
# # lambda chosen such that p(sigma > sigma_total) = 0.01, where sigma_total is
# # the total standard deviation of the corresponding MLE
# lambda = c(alpha = -log(0.01) / sd(eur3w_mle_data$alpha_hat),
#            beta = -log(0.01) / sd(eur3w_mle_data$beta_hat),
#            tau = -log(0.01) / sd(eur3w_mle_data$tau_hat))
# 
# lp_kappa = function(lkappa) {
#   kappa = exp(lkappa)
#   lp = sum(-lambda / sqrt(kappa) - 1.5 * kappa)
#   return(lp)
# }

# another approach to motivating the prior for kappa: our "null model" is to
# use the mles, i.e. no smoothing, which corresponds to kappa = 0, so if we
# want to penalise complexity, this should be our prior mode. A kappa that is
# too large oversmoothes, which we can define as the corresponding standard
# deviation being smaller than 1% of the total standard deviation. So we should
# use an exponential prior for kappa with lambda chosen as before. 
# lambda = c(alpha = -log(0.01) * var(eur3w_mle_data$alpha_hat) / .1,
#            beta = -log(0.01) * var(eur3w_mle_data$beta_hat) / .1,
#            tau = -log(0.01) * var(eur3w_mle_data$tau_hat) / .1)
# 
# # joint log prior as function of log kappa (exponential distribution)
# lp_lkappa = function(lkappa, lambda) {
#   lp = sum(-lambda * exp(lkappa))
#   return(lp)
# }
# NOTE: the exponential prior leads to too much smoothing. Despite assigning only 1% pprior probability to kappa > kappa_total, the prior allows for the posterior to extend far into positive values. Therefore, I will use a prior with a shorter tail, namely a half-Gaussian. Lambda here is the precision 

lambda = c(alpha = (var(eur3w_mle_data$alpha_hat) * qnorm(0.005))^2,
           beta = (var(eur3w_mle_data$beta_hat) * qnorm(0.005))^2,
           tau = (var(eur3w_mle_data$tau_hat) * qnorm(0.005))^2)

# joint log prior as function of log kappa (half normal distribution with
# precision lambda, such that p(kappa > kappa_total) = 0.01)
lp_lkappa = function(lkappa, lambda) {
  lp = sum(-lambda/2 * exp(lkappa)^2)
  return(lp)
}






# save parameters as list for the fitting functions
params = list(SS=n_lon * n_lat, TT=n_t, 
              vyy=eur3w_summary$vyy, vzz=eur3w_summary$vzz, 
              vyz=eur3w_summary$vyz, ybar=eur3w_summary$ybar, 
              RR=RR, nH_y=nH_y, x_hat=x_hat, lambda=lambda)



# approximate samples from p(x|y) by sampling x from the conditional
# p(x|y,kappa') for kappa' values sampled from p(kappa|y)
sample_lkappa_I_y = makefunction_sample_lkappa_I_y(params, lp_lkappa)
sample_x_I_lkappa_y = makefunction_sample_x_I_lkappa_y(params)
n_sampls = 500

x_sampls = replicate(n_sampls, {
             lk = sample_lkappa_I_y()
             x = sample_x_I_lkappa_y(lk)
           })


# calculate marginal posterior expectation and stdev of the latents
postexp_x = rowMeans(x_sampls)
postsd_x = apply(x_sampls, 1, sd)


eur3w_spat_data = eur3w_mle_data %>%
  ungroup() %>%
  mutate(alpha_star = postexp_x[1:n_s],
         beta_star = postexp_x[1:n_s + n_s],
         tau_star = postexp_x[1:n_s + 2 * n_s]) %>%
  mutate(sd_alpha_star = postsd_x[1:n_s],
         sd_beta_star = postsd_x[1:n_s + n_s],
         sd_tau_star = postsd_x[1:n_s + 2 * n_s])
  
plt_alpha = 
  eur3w_spat_data %>% select(long, lat, alpha_hat, alpha_star) %>%
  gather('alpha', 'value', alpha_hat, alpha_star) %>%
  ggplot() + 
    geom_raster(aes(x=long, y=lat, fill=value)) +
    facet_wrap(~alpha) + 
    scale_fill_gradientn(colors=hcl.colors(10)) +
    geom_path(data=world, aes(x=long, y=lat, group=group)) +
    coord_cartesian(xlim=long_range, ylim=lat_range)

plt_beta = 
eur3w_spat_data %>% select(long, lat, beta_hat, beta_star) %>%
  gather('beta', 'value', beta_hat, beta_star) %>%
  ggplot() + 
    geom_raster(aes(x=long, y=lat, fill=value)) +
    facet_wrap(~beta) + 
    scale_fill_gradient2(limits=c(-1.5,1.5)) +
    geom_path(data=world, aes(x=long, y=lat, group=group)) +
    coord_cartesian(xlim=long_range, ylim=lat_range)

plt_sigma = 
eur3w_spat_data %>% select(long, lat, tau_hat, tau_star) %>%
  gather('tau', 'value', tau_hat, tau_star) %>%
  rename(sigma = tau) %>% 
  mutate(value = exp(0.5*value)) %>%
  ggplot() + 
    geom_raster(aes(x=long, y=lat, fill=value)) +
    facet_wrap(~sigma) + 
    scale_fill_gradientn(colors=hcl.colors(10)) +
    geom_path(data=world, aes(x=long, y=lat, group=group)) +
    coord_cartesian(xlim=long_range, ylim=lat_range)

show(ggarrange(plt_alpha, plt_beta, plt_sigma, ncol=1))


