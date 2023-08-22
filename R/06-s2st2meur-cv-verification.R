library(tidyverse)
library(Matrix)
load(here::here('data', 'hindcast-eur.Rdata'))
load(here::here('data', 's2st2meur-crossvalidation-results.Rdata'))
source(here::here('R', '00-functions.R'))
outdir = here::here('R', 'fig', '06-s2st2m-cv-verification')

world = rnaturalearth::ne_coastline() %>% fortify

weighted.sd = function(x, w) {
  sqrt(weighted.mean(x^2, w) - weighted.mean(x, w)^2)
}

# calculate leave-one-out forecast anomalies, and leave-one-out predictions
# with mle, smoothed, and climatology
hindcast = hindcast_eur %>% 
  mutate(year = as.integer(substr(date, 1, 4))) %>%
  rename(fcst = ens_mean) %>%
  group_by(long, lat, lead_time) %>%
  # leave-one-out anomaly x(i) - xbar(-i) = n/(n-1) (x(i) - xbar), and
  # leave-one-out climatology ybar(-i)
  mutate(fcst_anom = n() / (n() - 1) * (fcst - mean(fcst)),
         ybar = n() / (n() - 1) * (mean(obs) - obs / n())) %>%
  select(long, lat, year, lead_time, fcst_anom, obs, ybar)

# join with verifying observations
hindcast_cv = s2st2meur_cv_result %>%
  inner_join(hindcast, by=c('long', 'lat', 'year', 'lead_time')) %>%
  mutate(yhat = alpha_hat + beta_hat * fcst_anom,
         ystar = alpha_smooth + beta_smooth * fcst_anom,
         yhat_var = var_alpha_hat + var_beta_hat * fcst_anom^2 + exp(2*tau_hat),
         ystar_var = alpha_smooth_var + beta_smooth_var * fcst_anom^2 + exp(2*tau_smooth))

# summarise MSE
mse_cv = 
hindcast_cv %>% 
  group_by(lead_time) %>%
  mutate(ww = cos(lat * pi / 180)) %>%
  summarise(mse_hat = weighted.mean((obs - yhat)^2, ww),
            mse_hat_se = weighted.sd((obs - yhat)^2, ww) / sqrt(n()),
            mse_smooth = weighted.mean((obs - ystar)^2, ww),
            mse_smooth_se = weighted.sd((obs - yhat)^2, ww) / sqrt(n()),
            mse_clim = weighted.mean((obs - ybar)^2, ww),
            mse_clim_se = weighted.sd((obs - ybar)^2, ww) / sqrt(n()))

# summarise relative improvement at high lead times (no plot)
mse_cv %>% 
  filter(between(lead_time, 10, 14)) %>% 
  select(-ends_with('_se')) %>% 
  mutate(skill = (mse_hat - mse_smooth) / (mse_clim - mse_hat)) %>%
  round(digits=2) %>%
  knitr::kable()


# plot MSE vs lead time
mse_plt =
ggplot(mse_cv) + 
geom_line(aes(x=lead_time, y=mse_hat, colour='unsmoothed')) +
geom_segment(aes(x=lead_time, xend=lead_time, 
                 y=mse_hat - 2*mse_hat_se, yend=mse_hat + 2 * mse_hat_se, 
                 colour='unsmoothed')) +
geom_line(aes(x=lead_time, y=mse_smooth, colour='smoothed')) +
geom_segment(aes(x=lead_time, xend=lead_time, 
                 y=mse_smooth - 2*mse_smooth_se, yend=mse_smooth + 2 * mse_smooth_se, 
                 colour='smoothed')) +
geom_line(aes(x=lead_time, y=mse_clim, colour='climatology')) +
geom_segment(aes(x=lead_time, xend=lead_time, 
                 y=mse_clim - 2*mse_clim_se, yend=mse_clim + 2 * mse_clim_se, 
                 colour='climatology')) +
labs(x='lead time [days]', y='Mean squared error', colour=NULL)
ggsave(paste(outdir, 'mse.png', sep='/'), mse_plt, width=5, height=3)


# plot MSE difference vs lead time
msediff_plt =
hindcast_cv %>% 
  group_by(lead_time) %>%
  mutate(ww = cos(lat * pi / 180)) %>%
  summarise(mse_diff = weighted.mean((obs - yhat)^2 - (obs-ystar)^2, ww),
            mse_diff_se = weighted.sd((obs - yhat)^2 - (obs - ystar)^2, ww) / sqrt(n())) %>%
ggplot() + 
geom_ribbon(aes(x=lead_time, ymin=mse_diff - 2*mse_diff_se, ymax=mse_diff + 2*mse_diff_se), fill='lightblue') +
geom_line(aes(x=lead_time, y=mse_diff)) +
labs(x='lead time [days]', y='Mean squared error difference', colour='NULL') +
geom_abline(slope=0, lty='dashed')
ggsave(paste(outdir, 'mse-diff.png', sep='/'), msediff_plt, width=5, height=3)


# relative improvement of smoothing over MLE measured in units of improvement
# of MLE over clim
ggplot(mse_cv) +
geom_line(aes(x=lead_time, y=(mse_smooth - mse_hat) / (mse_hat - mse_clim))) +
coord_cartesian(xlim=c(1,25), ylim=c(0,1))


# MSE improvement map (red = smoothing improves)
msemap_plt = 
hindcast_cv %>% 
  mutate(ww = cos(lat * pi / 180)) %>%
  group_by(lead_time, long, lat) %>%
  summarise(mse_diff = weighted.mean((yhat - obs)^2 - (ystar - obs)^2, w=ww)) %>%
  ggplot() + 
    geom_raster(aes(x=long, y=lat, fill=mse_diff)) + 
    facet_wrap(~lead_time) + 
    scale_fill_gradient2(low='blue', mid='white', high='red', limits=c(-2,2)) +
    geom_path(data=world, aes(x=long, y=lat, group=group), colour='black') +
    xlim(range(hindcast_cv$long)) + ylim(range(hindcast_cv$lat)) +
    labs(x='Longitude', y='Latitude', fill='MSE\ndifference')
ggsave(paste(outdir, 'msediff-map.png', sep='/'), msemap_plt, width=13, height=10)


# coverage of the approximate 95% prediction interval
qq = qt(0.975, df=19)
covg_plt = 
hindcast_cv %>% 
mutate(ww = cos(lat * pi / 180)) %>%
mutate(lwr_hat = yhat - qq * sqrt(yhat_var), 
       upr_hat = yhat + qq * sqrt(yhat_var),
       lwr_star = ystar - qq * sqrt(ystar_var), 
       upr_star = ystar + qq * sqrt(ystar_var)) %>%
group_by(lead_time) %>%
summarise(cov95_hat = weighted.mean(obs > lwr_hat & obs < upr_hat, ww),
          cov95_star = weighted.mean(obs > lwr_star & obs < upr_star, ww)) %>%
ggplot() + 
  geom_line(aes(x=lead_time, y=cov95_hat, colour='unsmoothed')) +
  geom_line(aes(x=lead_time, y=cov95_star, colour='smoothed')) +
  labs(x='lead time [days]', y='Coverage of the 95% CI', colour=NULL)
ggsave(paste(outdir, 'coverage.png', sep='/'), covg_plt, width=5, height=3)


# log score difference of the predictive t distribution with df=19
logs_plt = 
hindcast_cv %>% 
mutate(ww = cos(lat * pi / 180)) %>%
mutate(ign_hat = -tau_hat + dt((obs - yhat)*exp(-tau_hat), df=19, log=TRUE),
       ign_star = -tau_smooth + dt((obs - ystar)*exp(-tau_smooth), df=19, log=TRUE)) %>%
group_by(lead_time) %>%
summarise(ign_diff = weighted.mean(ign_star - ign_hat, ww),
          ign_diff_sd = weighted.sd(ign_star - ign_hat, ww) / sqrt(n())) %>%
ggplot() + 
  geom_ribbon(aes(x=lead_time, ymin=ign_diff - 2*ign_diff_sd, ymax=ign_diff + 2*ign_diff_sd), fill='lightblue') +
  geom_line(aes(x=lead_time, y=ign_diff)) +
  labs(x='lead time [days]', y='log score difference') +
  geom_abline(slope=0, lty='dashed')
ggsave(paste(outdir, 'logscore-diff.png', sep='/'), logs_plt, width=5, height=3)



