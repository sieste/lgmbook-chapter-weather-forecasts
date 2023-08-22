library(tidyverse)
load('../data/hindcast-s2secmwf-eraint.Rdata')

# africa: 20W-55E, 40SN-40N
filter_africa = function(x) {
  filter(x, between(long, -20, 55), between(lat, -40, 40))
}
# southeast asia: 90E-150E, 15S-40N
filter_se_asia = function(x) {
  filter(x, between(long, 90, 150), between(lat, -15, 40))
}
# europe 15W-30E, 30N-70N :
filter_europe = function(x) {
  filter(x, between(long, -15, 30), between(lat, 30, 70))
}
# british islands
filter_uk = function(x) {
  filter(x, between(long, -10.5, 2.5), between(lat, 48, 60))
}

hindcast_eur3w = hindcast_s2secmwf_eraint %>%
  mutate(lead_week = ceiling(lead_time / 7)) %>%
  filter_europe() %>%
  filter(lead_week == 3) %>%
  select(-lead_week, -lead_time) %>%
  mutate(year = as.integer(substr(date, 1, 4))) %>%
  #filter(year >= 2015) %>%
  group_by(long, lat, year) %>%
  summarise(obs = mean(obs), 
            fcst = mean(fcst), 
            .groups = 'drop') %>%
  group_by(long, lat) %>%
  mutate(fcst_anom = fcst - mean(fcst)) %>%
  ungroup

save(file = '../data/hindcast-eur3w.Rdata', list='hindcast_eur3w')

