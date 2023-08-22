library(tidyverse)
load('../data/s2st2meur3w-crossvalidation-forecasts.Rdata')
load('../data/hindcast-eur3w.Rdata')

# from the saved data, extract the different forecast types for each year, the
# grid, and the observations
forecasts = s2st2meur3w_cv_forecasts
forecast_grid = forecasts[['grid']]
forecasts = forecasts[ names(forecasts) != 'grid' ]

y_hat = forecasts %>% map('y_hat')   # mle based
y_star = forecasts %>% map('y_star') # spatially smoothed

obs = forecasts %>% 
      map('y') %>% 
      map( ~ tibble(obs = .) ) %>%
      map( ~ bind_cols(., forecast_grid) ) %>%
      map_dfr(~ ., .id = 'year')

# combine posterior predictive means and observations
y_hat_bar = y_hat %>% 
  map( ~ tibble(y_hat_bar = rowMeans(.)) ) %>%
  map_dfr( ~ bind_cols(., forecast_grid), .id = 'year' ) %>%
  right_join(obs, by=c('year', 'long', 'lat')) %>%
  mutate(sqerr_hat = (y_hat_bar - obs)^2)
  
y_star_bar = y_star %>% 
  map( ~ tibble(y_star_bar = rowMeans(.)) ) %>%
  map_dfr( ~ bind_cols(., forecast_grid), .id = 'year' ) %>%
  right_join(obs, by=c('year', 'long', 'lat')) %>%
  mutate(sqerr_star = (y_star_bar - obs)^2)

# calculate mse
mse_hat = y_hat_bar %>% 
  group_by(long, lat) %>%
  summarise(mse_hat = mean(sqerr_hat), .groups='drop')

mse_star = y_star_bar %>% 
  group_by(long, lat) %>%
  summarise(mse_star = mean(sqerr_star), .groups='drop')

mse = mse_hat %>% right_join(mse_star, by=c('long', 'lat'))

cat('### coslat weighted mse ### \n')
mse %>% summarise(mse_hat = weighted.mean(mse_hat, w=cos(lat/180*pi)),
                  mse_star = weighted.mean(mse_star, w=cos(lat/180*pi))) %>%
print

cat('### fraction of forecasts improved ### \n')

right_join(
  y_hat_bar %>% select(year, long, lat, sqerr_hat),
  y_star_bar %>% select(year, long, lat, sqerr_star),
  by = c('year', 'long', 'lat')
) %>%
summarise(frac_impr = weighted.mean(sqerr_star < sqerr_hat, w=cos(lat/180*pi))) %>% 
print

cat('### coverage of the 90% CI ###\n')


