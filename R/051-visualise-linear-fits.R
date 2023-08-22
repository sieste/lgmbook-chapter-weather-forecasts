library(tidyverse)
library(lubridate)
load('../data/s2st2meur-crossvalidation-results.Rdata')
load('../data/hindcast-eur.Rdata')
outdir = here::here('R', 'fig', '05-fit-plots')

long0 = 0
dlong = 3
lat0 = 51
dlat = 3
year0 = 2019
lead0 = 5

params = s2st2meur_cv_result %>%
  filter(between(long, long0, long0+dlong), 
         between(lat, lat0, lat0+dlat),
         year == year0,
         lead_time == lead0) %>%
  select(long, lat, alpha_hat, beta_hat, tau_hat, alpha_smooth, beta_smooth, tau_smooth)

dat = hindcast_eur %>%
  filter(between(long, long0, long0+dlong), 
         between(lat, lat0, lat0+dlat),
         #year(as.Date(date)) != year0,
         lead_time == lead0) 

xbar = dat %>% group_by(long, lat) %>% summarise(xbar = mean(ens_mean))
params = params %>% inner_join(xbar, by=c('long','lat'))


plt = 
ggplot(dat) +
  geom_point(aes(x=ens_mean, y=obs, pch='hindcast data')) + 
  facet_grid(lat~long, labeller=labeller(lat=\(x) paste(x, 'deg. N'), long=\(x) paste(x, 'deg. W'))) +
  geom_abline(data=params, aes(intercept=alpha_hat-beta_hat*xbar, slope=beta_hat, colour='MLE')) +
  geom_abline(data=params, aes(intercept=alpha_smooth-beta_smooth*xbar, slope=beta_smooth, colour='smoothed')) +
  #geom_abline(aes(intercept=0, slope=1, colour='0-1 diagonal'), lty=2) +
  xlim(10,25) + ylim(10,25) +
  labs(colour=NULL, pch=NULL, x='forecast temperature [C]', y='observed temperature [C]')
  

ggsave(filename=here::here(outdir, 'linear-fits.png'), plot=plt, width=6, height=5)


