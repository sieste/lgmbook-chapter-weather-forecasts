library(tidyverse)
suppressPackageStartupMessages(library(INLA))
library(Matrix)
library(ggpubr)
load('../data/hindcast-eur.Rdata')
source('00-functions.R')
outdir = 'fig/04-inference'

hindcast = hindcast_eur %>%
  select(-ens_sd) %>%
  rename(fcst = ens_mean)

world = rnaturalearth::ne_coastline() %>% fortify
lat_range = range(hindcast$lat)
long_range = range(hindcast$long)


# calculate mles and asymptotic variances per grid point and lead time
hindcast_mle_data = hindcast %>%
  group_by(long, lat, lead_time) %>%
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
  select(-vyy, -vzz, -vyz, -ybar)


# arrange by long and lat such that the estimators are arranged as a vectorised
# 2d field (stacked columns) from upper left to lower right
# hindcast_mle_data = hindcast_mle_data %>% arrange(lead_time, long, desc(lat))


# smooth alpha and beta for lead time 1

N_i = hindcast %>% pull(lat) %>% unique %>% length
N_j = hindcast %>% pull(long) %>% unique %>% length
df_mle = hindcast_mle_data %>% 
  filter(lead_time == 1) %>% 
  mutate(i = as.numeric(factor(lat)),
         j = as.numeric(factor(long)),
         s = (j - 1) * N_i + i) %>%
  #select(s, long, lat, alpha_hat, var_alpha_hat) %>%
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
  mutate(alpha_smooth = inla_alpha$summary.fitted.values$mean,
         alpha_smooth_var = inla_alpha$summary.fitted.values$sd^2,
         beta_smooth = inla_beta$summary.fitted.values$mean,
         beta_smooth_var = inla_beta$summary.fitted.values$sd^2,
         tau_smooth = inla_tau$summary.fitted.values$mean,
         tau_smooth_var = inla_tau$summary.fitted.values$sd^2)

ggplot(df_mle) + 
  geom_abline() + 
  geom_point(aes(x=alpha_hat, y=alpha_smooth)) +
  geom_segment(aes(x=alpha_hat, xend=alpha_hat, y=alpha_smooth - 2 * sqrt(alpha_smooth_var), yend=alpha_smooth + 2 * sqrt(alpha_smooth_var)))


ggplot(df_mle) + 
  geom_abline() + 
  geom_point(aes(x=beta_hat, y=beta_smooth)) +
  geom_segment(aes(x=beta_hat, xend=beta_hat, y=beta_smooth - 2 * sqrt(beta_smooth_var), yend=beta_smooth + 2 * sqrt(beta_smooth_var)))

ggplot(df_mle) + 
  geom_abline() + 
  geom_point(aes(x=tau_hat, y=tau_smooth)) +
  geom_segment(aes(x=tau_hat, xend=tau_hat, y=tau_smooth - 2 * sqrt(tau_smooth_var), yend=tau_smooth + 2 * sqrt(tau_smooth_var)))


post_df = df_mle %>%
  filter(lead_time == 1) %>%
  select(long, lat, 
         alpha_hat, alpha_smooth, 
         beta_hat, beta_smooth, 
         tau_hat, tau_smooth) %>%
  pivot_longer(col=c(-long, -lat), names_to='par', values_to='value') %>%
  separate(par, into=c('par', 'type'), sep='_') %>%
  mutate(type=stringr::str_replace_all(type, c('hat'='maximum likelihood', 'smooth'='posterior mean')))

post_plt = ggarrange(
  post_df %>% filter(par == 'alpha') %>%
  ggplot() +
    geom_raster(aes(x=long, y=lat, fill=value)) +
    scale_fill_viridis_c() +
    geom_path(data=world, aes(x=long, y=lat, group=group), colour='white') +
    xlim(long_range) + ylim(lat_range) +
    facet_wrap(~type) +
    labs(x='Longitude', y='Latitude', fill='alpha'),
  post_df %>% filter(par == 'beta') %>%
  ggplot() +
    geom_raster(aes(x=long, y=lat, fill=value)) +
    scale_fill_viridis_c() +
    geom_path(data=world, aes(x=long, y=lat, group=group), colour='white') +
    xlim(long_range) + ylim(lat_range) +
    facet_wrap(~type) +
    labs(x='Longitude', y='Latitude', fill='beta'),
  post_df %>% filter(par == 'tau') %>%
  ggplot() +
    geom_raster(aes(x=long, y=lat, fill=value)) +
    scale_fill_viridis_c() +
    geom_path(data=world, aes(x=long, y=lat, group=group), colour='white') +
    xlim(long_range) + ylim(lat_range) +
    facet_wrap(~type) + 
    labs(x='Longitude', y='Latitude', fill='tau'),
  ncol=1
)
ggsave(paste(outdir, 'mle-vs-postmean.png', sep='/'), post_plt, width=7, height=9)


hyper_plt = 
bind_rows(
  inla_alpha$marginals.hyperpar$`Precision for s` %>% 
    as_tibble() %>%
    mutate(par = 'alpha'),
  inla_beta$marginals.hyperpar$`Precision for s` %>% 
    as_tibble() %>%
    mutate(par = 'beta'),
  inla_tau$marginals.hyperpar$`Precision for s` %>% 
    as_tibble() %>%
    mutate(par = 'tau')
) %>%
ggplot() +
geom_line(aes(x=x, y=y, colour=par), show.legend=FALSE) + facet_wrap(~par, scale='free', ncol=1) +
labs(x='Precision hyperparameter', y='Posterior density')
ggsave(paste(outdir, 'hyper-posterior.png', sep='/'), hyper_plt, width=5, height=6)



# # extract vector of MLEs and negative Hessian 
# x_hat = with(eur3w_mle_data, c(alpha_hat, beta_hat, tau_hat))
# nH_diag = with(eur3w_mle_data, c(n * exp(-2 * tau_hat),
#                                  n * sxx * exp(-2 * tau_hat),
#                                  2 * n))
# nH_y = sparseMatrix(i=1:length(nH_diag), j=1:length(nH_diag), x=nH_diag)
# 
# # dimension paramaters
# n_lon = eur3w_mle_data %>% pull(long) %>% n_distinct
# n_lat = eur3w_mle_data %>% pull(lat) %>% n_distinct
# n_t = nrow(eur3w_lm_data$hindcast[[1]])
# n_s = n_lon * n_lat
# n_x = 3 * n_s
# 
# # 2d random walk structure matrix (Q = kappa * R)
# RR = make_rw2d_structure_matrix(n_lon, n_lat)

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
# NOTE: the exponential prior leads to too much smoothing. Despite assigning
# only 1% pprior probability to kappa > kappa_total, the prior allows for the
# posterior to extend far into positive values. Therefore, I will use a prior
# with a shorter tail, namely a half-Gaussian. Lambda here is the precision 
# 
# lambda = c(alpha = (var(eur3w_mle_data$alpha_hat) * qnorm(0.005))^2,
#            beta = (var(eur3w_mle_data$beta_hat) * qnorm(0.005))^2,
#            tau = (var(eur3w_mle_data$tau_hat) * qnorm(0.005))^2)
# 
# # joint log prior as function of log kappa (half normal distribution with
# # precision lambda, such that p(kappa > kappa_total) = 0.01)
# lp_lkappa = function(lkappa, lambda) {
#   lp = sum(-lambda/2 * exp(lkappa)^2)
#   return(lp)
# }
# 
# 
# # save parameters as list for the fitting functions
# params = list(SS=n_lon * n_lat, TT=n_t, 
#               vyy=eur3w_summary$vyy, vzz=eur3w_summary$vzz, 
#               vyz=eur3w_summary$vyz, ybar=eur3w_summary$ybar, 
#               RR=RR, nH_y=nH_y, x_hat=x_hat, lambda=lambda)
# 
# 
# # approximate samples from p(x|y) by sampling x from the conditional
# # p(x|y,kappa') for kappa' values sampled from p(kappa|y)
# sample_lkappa_I_y = makefunction_sample_lkappa_I_y(params, lp_lkappa)
# sample_x_I_lkappa_y = makefunction_sample_x_I_lkappa_y(params)
# n_sampls = 500
# 
# x_sampls = replicate(n_sampls, {
#              lk = sample_lkappa_I_y()
#              x = sample_x_I_lkappa_y(lk)
#            })
# 
# 
# # calculate marginal posterior expectation and stdev of the latents
# postexp_x = rowMeans(x_sampls)
# postsd_x = apply(x_sampls, 1, sd)
# 
# 
# eur3w_spat_data = eur3w_mle_data %>%
#   ungroup() %>%
#   mutate(alpha_star = postexp_x[1:n_s],
#          beta_star = postexp_x[1:n_s + n_s],
#          tau_star = postexp_x[1:n_s + 2 * n_s]) %>%
#   mutate(sd_alpha_star = postsd_x[1:n_s],
#          sd_beta_star = postsd_x[1:n_s + n_s],
#          sd_tau_star = postsd_x[1:n_s + 2 * n_s])
#   
# plt_alpha = 
#   eur3w_spat_data %>% select(long, lat, alpha_hat, alpha_star) %>%
#   gather('alpha', 'value', alpha_hat, alpha_star) %>%
#   ggplot() + 
#     geom_raster(aes(x=long, y=lat, fill=value)) +
#     facet_wrap(~alpha) + 
#     scale_fill_gradientn(colors=hcl.colors(10)) +
#     geom_path(data=world, aes(x=long, y=lat, group=group)) +
#     coord_cartesian(xlim=long_range, ylim=lat_range)
# 
# plt_beta = 
# eur3w_spat_data %>% select(long, lat, beta_hat, beta_star) %>%
#   gather('beta', 'value', beta_hat, beta_star) %>%
#   ggplot() + 
#     geom_raster(aes(x=long, y=lat, fill=value)) +
#     facet_wrap(~beta) + 
#     scale_fill_gradient2(limits=c(-1.5,1.5)) +
#     geom_path(data=world, aes(x=long, y=lat, group=group)) +
#     coord_cartesian(xlim=long_range, ylim=lat_range)
# 
# plt_sigma = 
# eur3w_spat_data %>% select(long, lat, tau_hat, tau_star) %>%
#   gather('tau', 'value', tau_hat, tau_star) %>%
#   rename(sigma = tau) %>% 
#   mutate(value = exp(0.5*value)) %>%
#   ggplot() + 
#     geom_raster(aes(x=long, y=lat, fill=value)) +
#     facet_wrap(~sigma) + 
#     scale_fill_gradientn(colors=hcl.colors(10)) +
#     geom_path(data=world, aes(x=long, y=lat, group=group)) +
#     coord_cartesian(xlim=long_range, ylim=lat_range)
# 
# show(ggarrange(plt_alpha, plt_beta, plt_sigma, ncol=1))





