library(tidyverse)
library(data.table)

# combine all forecasts and observations into a hindcast data set with columns
# lat, long, start_date, lead_time, fcst, obs
load('eraint-t2m-for-s2s.Rdata')
load('s2s-ecmwf-t2m.Rdata')

meta = s2s_ecmwf_t2m_meta
meta[, fcst_dates := .(paste(start_date, veri_date, sep='::'))]

obs = eraint_t2m_for_s2s
fcst = s2s_ecmwf_t2m
fcst = setnames(fcst, old=meta[, band], new=meta[, fcst_dates])

obs = melt(obs, id.vars=c('long', 'lat'), 
           variable.name='date', value.name='obs')
fcst = melt(fcst, id.vars=c('long', 'lat'),
            variable.name='date', value.name='fcst')

fcst[, start_date := .(substr(date, 1, 10))][
     , date := .(substr(date, 13,22))][
     , lead_time := .(as.numeric(as.Date(date)) - 
                      as.numeric(as.Date(start_date)))][
     , start_date := NULL]

hindcast = merge(obs, fcst, by=c('long', 'lat', 'date'))
hindcast_s2secmwf_eraint = as_tibble(hindcast)
save(file='hindcast-s2secmwf-eraint.Rdata', 
     list='hindcast_s2secmwf_eraint')








