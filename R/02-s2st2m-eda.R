library(tidyverse)
library(Matrix)
load('../data/hindcast-s2secmwfens-eraint.Rdata')
source('00-functions.R')

hindcast = hindcast_s2secmwfens_eraint %>%
  mutate(year = as.integer(substr(date, 1, 4)))


world = rnaturalearth::ne_coastline() %>% fortify

# correlation and bias of daily forecasts 
summary_daily = hindcast %>% 
  select(-ens_sd) %>%
  rename(fcst = ens_mean) %>%
  filter(lead_time %in% c(3,7,21,42)) %>%
  group_by(long, lat, lead_time) %>%
  summarise(cor = cor(fcst, obs), bias=mean(fcst - obs), .groups='drop')

outdir = 'fig/02-s2st2m-eda'

# worldwide correlation map on a few days 
n_colrs = 5
colrs = rev(RColorBrewer::brewer.pal(n=n_colrs, 'PuOr'))

cor_daily_plt = 

summary_daily %>% 
  filter(lead_time %in% c(3,7,21,42)) %>%
  ggplot() + 
    geom_raster(aes(x=long, y=lat, fill=cor)) + 
    facet_wrap(~lead_time, ncol=1, 
               labeller=labeller(lead_time = function(x) paste(x,'days'))) + 
    scale_fill_stepsn(colours=colrs,
                      breaks=seq(-1,1,len=n_colrs+1),
                      limits=c(-1,1)) +
    geom_path(data=world, aes(x=long, y=lat, group=group)) +
    ylim(-70,80) +
    theme(legend.position='bottom') +
    guides(fill=guide_colourbar(barwidth=unit(0.7, 'npc'), 
           title.position='bottom', 
           title.hjust=.5, 
           barheight=unit(.01,'npc'))) +
    labs(fill='correlation', x='Longitude', y='Latitude')

ggsave(paste(outdir, 'cormap-daily.png', sep='/'), 
       cor_daily_plt, width=4, height=8)


# worldwide bias map on a few days 
n_colrs = 10
colrs = rev(RColorBrewer::brewer.pal(n=n_colrs, 'PuOr'))

bias_daily_plt = 
summary_daily %>% 
  #filter(lead_time %in% c(3,7,14,28)) %>%
  ggplot() + 
    geom_raster(aes(x=long, y=lat, fill=bias)) +
    facet_wrap(~lead_time, ncol=1,
               labeller=labeller(lead_time = function(x) paste(x,'days'))) + 
    scale_fill_stepsn(colours=colrs, breaks=seq(-5,5,len=n_colrs+1),
                      limits=c(-5,5),
                      labels = function(x) ifelse(abs(x) <= 4, x, '')
                      ) +
    geom_path(data=world, aes(x=long, y=lat, group=group)) +
    ylim(-70,80) +
    theme(legend.position='bottom') +
    guides(fill=guide_colourbar(barwidth=unit(0.7, 'npc'), 
           title.position='bottom', 
           title.hjust=.5, 
           barheight=unit(.01,'npc'))) +
    labs(fill='forecast bias [K]', x='Longitude', y='Latitude')

ggsave(paste(outdir, 'biasmap-daily.png', sep='/'), 
       bias_daily_plt, width=4, height=8)


# # calculate weekly averages of observations and forecasts 
# hindcast_weekly =
#   hindcast %>% 
#   mutate(lead_week = ceiling(lead_time / 7)) %>%
#   group_by(long, lat, year, lead_week) %>%
#   summarise(obs = mean(obs), 
#             fcst = mean(fcst), 
#             .groups = 'drop')
#
# # correlation and bias of weekly averages
# summary_weekly =
#   hindcast_weekly %>%
#   group_by(long, lat, lead_week) %>%
#   summarise(cor = cor(fcst, obs), bias=mean(fcst - obs), .groups='drop') 

# # europe correlation
# ggplot(summary_weekly) + 
#   geom_raster(aes(x=long, y=lat, fill=cor)) + 
#   facet_wrap(~lead_week) + 
#   scale_fill_gradient2(low='blue', mid='white', high='red', limits=c(-1,1)) + 
#   geom_path(data=world, aes(x=long, y=lat, group=group)) + 
#   xlim(-15,35) + ylim(30,70)
# 
# # africa correlation
# ggplot(summary_weekly) + 
#   geom_raster(aes(x=long, y=lat, fill=cor)) + 
#   facet_wrap(~lead_week) + 
#   scale_fill_gradient2(low='blue', mid='white', high='red', limits=c(-1,1)) + 
#   geom_path(data=world, aes(x=long, y=lat, group=group)) + 
#   xlim(-20,55) + ylim(-40,40)
# 
# 
# # china/se-asia correlation
# ggplot(summary_weekly) + 
#   geom_raster(aes(x=long, y=lat, fill=cor)) + 
#   facet_wrap(~lead_week) + 
#   scale_fill_gradient2(low='blue', mid='white', high='red', limits=c(-1,1)) + 
#   geom_path(data=world, aes(x=long, y=lat, group=group)) + 
#   xlim(90,150) + ylim(-15,40)
# 
# 
# # subsets:
# # EUR3: 15W-30E, 35N-70N, week 3
# # AFR3: 20W-55E, 40SN-40N, week 3
# # SEA3: 90E-150E, 15S-40N, week 3
# 
# 
# # TODO: 
# # cv error of obs ~ fcst vs obs ~ 1
# # any further generic facts about weather forecasting that I want to illustrate



