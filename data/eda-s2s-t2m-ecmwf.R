load('hindcast-s2secmwf-eraint.Rdata')

library(tidyverse)
library(rnaturalearth)

# load coastline data from rnaturalearth, fortify() so the spatialdataframe is
# converted into a lat/long data frame that can be plotted with geom_path
world = rnaturalearth::ne_coastline() %>% fortify() 

# calculate summary statistics (e.g. regression parameters) for each grid point
# and lead time
hindcast = hindcast_s2secmwf_eraint %>% group_by(lat, long, lead_time)
smry = hindcast %>%
  mutate(fcst_anom = fcst - mean(fcst)) %>%
  summarise(alpha_hat = mean(obs),
            beta_hat = cov(fcst, obs) / var(fcst),
            tau_hat = 0.5 * log(mean((obs - alpha_hat - beta_hat * fcst_anom)^2)))


# plot
smry %>% filter(lead_time %% 7 == 1) %>%
ggplot() +
  geom_raster(aes(x=long, y=lat, fill=beta_hat)) + 
  geom_path(data=world, aes(x=long, y=lat, group=group)) + 
  facet_wrap(~lead_time) + 
  scale_fill_gradient2(low='blue', mid='white', high='red')


