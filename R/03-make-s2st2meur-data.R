library(tidyverse)
load('../data/hindcast-s2secmwfens-eraint.Rdata')

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
  filter(x, between(long, -24, 60), between(lat, 34.5, 69))
}
# british islands
filter_uk = function(x) {
  filter(x, between(long, -10.5, 2.5), between(lat, 48, 60))
}

hindcast_eur = hindcast_s2secmwfens_eraint %>%
  filter_europe()

save(file = '../data/hindcast-eur.Rdata', list='hindcast_eur')

