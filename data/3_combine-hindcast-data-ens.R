library(tidyverse)
library(data.table)
library(lubridate)

# combine all forecasts and observations into a hindcast data set with columns
# lat, long, start_date, lead_time, fcst, obs
load('eraint-t2m-for-s2s.Rdata')
load('s2s_ecmwf_t2m_ens.Rdata')

meta = s2s_ecmwf_t2m_ens_meta
meta[, date := .(paste(start_date, veri_date, sep='::'))]

obs = eraint_t2m_for_s2s
ens = s2s_ecmwf_t2m_ens
ens = setnames(ens, old=meta[, band], new=meta[, date])

obs = melt(obs, id.vars=c('long', 'lat'), 
           variable.name='date', value.name='obs')
ens = melt(ens, id.vars=c('long', 'lat'),
            variable.name='date', value.name='fcst')

# summarise ensemble (mean and sd) 
fcst = ens[, .(ens_mean = mean(fcst), ens_sd = sd(fcst)), by=.(long, lat, date)]


fcst[, start_date := .(substr(date, 1, 10))]
fcst[, date := .(substr(date, 13,22))]
fcst[, lead_time := .(lubridate::yday(date) - lubridate::yday(start_date))]
fcst[, start_date := NULL]

hindcast = merge(obs, fcst, by=c('long', 'lat', 'date'))
hindcast_s2secmwfens_eraint = as_tibble(hindcast)
save(file='hindcast-s2secmwfens-eraint.Rdata', 
     list='hindcast_s2secmwfens_eraint')








