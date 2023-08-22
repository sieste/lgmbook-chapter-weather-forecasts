library(tidyverse)
library(Matrix)
load('../data/hindcast-s2secmwfens-eraint.Rdata')
source('00-functions.R')

hindcast = hindcast_s2secmwfens_eraint

world = rnaturalearth::ne_coastline() %>% fortify

# fit linear model y_t = a + b * (x_t - xm) + sigma * e_t
# and plot b and sigma on a map

# regression
mos_coefs = 
hindcast %>% 
  select(-ens_sd) %>%
  rename(fcst = ens_mean) %>%
  filter(lead_time %in% c(3,7,14,28)) %>%
  group_by(long, lat, lead_time) %>%
  summarise(ybar = mean(obs), 
            xbar = mean(fcst), 
            sxx = sum((fcst-xbar)^2), 
            sxy = sum((obs - ybar) * (fcst-xbar)), 
            beta0_hat = ybar, 
            beta1_hat = sxy/sxx, 
            sigma_hat = sqrt(mean((obs - beta0_hat - beta1_hat * (fcst - xbar))^2)),
            .groups='drop')


outdir = 'fig/02-s2st2m-eda'

# intercept estimates on a map
###########################################################

n_colrs = 11
colrs = RColorBrewer::brewer.pal(n=n_colrs, 'OrRd')

beta0_plt = 
  mos_coefs %>% 
  ggplot() + 
    geom_raster(aes(x=long, y=lat, fill=beta0_hat)) + 
    facet_wrap(~lead_time, ncol=1, 
               labeller=labeller(lead_time = function(x) paste(x,'days'))) + 
    scale_fill_stepsn(colours=colrs,
                      breaks=seq(-10,45,len=n_colrs+1),
                      limits=c(-10,45)) +
#    scale_fill_gradientn(colours=colrs, limits=c(-2,2)) +
    geom_path(data=world, aes(x=long, y=lat, group=group)) +
    ylim(-70,80) +
    theme(legend.position='bottom') +
    guides(fill=guide_colourbar(barwidth=unit(0.7, 'npc'), 
           title.position='bottom', 
           title.hjust=.5, 
           barheight=unit(.01,'npc'))) +
    labs(fill=expression(hat(beta)[0]), x='Longitude', y='Latitude')
ggsave(paste(outdir, 'beta0hatmap.png', sep='/'), 
       beta0_plt, width=4, height=8)


# slope estimates on a map
###########################################################

n_colrs = 11
colrs = rev(RColorBrewer::brewer.pal(n=n_colrs, 'PuOr'))
beta1_plt = 
  mos_coefs %>% 
  ggplot() + 
    geom_raster(aes(x=long, y=lat, fill=beta1_hat)) + 
    facet_wrap(~lead_time, ncol=1, 
               labeller=labeller(lead_time = function(x) paste(x,'days'))) + 
    scale_fill_stepsn(colours=colrs,
                      breaks=pretty(seq(-1,1,len=n_colrs+1),n=n_colrs+1),
                      limits=c(-1,1)) +
#    scale_fill_gradientn(colours=colrs, limits=c(-2,2)) +
    geom_path(data=world, aes(x=long, y=lat, group=group)) +
    ylim(-70,80) +
    theme(legend.position='bottom') +
    guides(fill=guide_colourbar(barwidth=unit(0.7, 'npc'), 
           title.position='bottom', 
           title.hjust=.5, 
           barheight=unit(.01,'npc'))) +
    labs(fill=expression(hat(beta)[1]), x='Longitude', y='Latitude')

ggsave(paste(outdir, 'beta1hatmap.png', sep='/'), 
       beta1_plt, width=4, height=8)

# residual stdev estimates on a map
###########################################################

n_colrs=10
colrs = RColorBrewer::brewer.pal(n=n_colrs, 'OrRd')
sigma_plt = 
  mos_coefs %>% 
  ggplot() + 
    geom_raster(aes(x=long, y=lat, fill=sigma_hat)) + 
    facet_wrap(~lead_time, ncol=1, 
               labeller=labeller(lead_time = function(x) paste(x,'days'))) + 
    scale_fill_stepsn(colours=colrs,
                      breaks=seq(0,5,len=n_colrs+1),
                      limits=c(0,5)) +
#    scale_fill_gradientn(colours=colrs, limits=c(-2,2)) +
    geom_path(data=world, aes(x=long, y=lat, group=group)) +
    ylim(-70,80) +
    theme(legend.position='bottom') +
    guides(fill=guide_colourbar(barwidth=unit(0.7, 'npc'), 
           title.position='bottom', 
           title.hjust=.5, 
           barheight=unit(.01,'npc'))) +
    labs(fill=expression(hat(sigma)), x='Longitude', y='Latitude')

ggsave(paste(outdir, 'sigmahatmap.png', sep='/'), 
       sigma_plt, width=4, height=8)



