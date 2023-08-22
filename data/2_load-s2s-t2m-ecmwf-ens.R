# This R script load the t2m forecast data from the grib file with
# rgdal::readGDAL and extracts meta information (forecast date and lead time)
# from the grib header with gdalUtils::gdalinfo. 

library('tidyverse')
library('data.table')
library('rgdal')
library('gdalUtils')

grib_file = 's2s_ecmwf_t2m_ens.grib'

# load the grib file and transform to data table
out = rgdal::readGDAL(grib_file)
out = as.data.table(out)

# `out` has 9200 "bands": 10 bands per day times 46 forecast days per year times
# 20 hindcast years. There are also coordinates columns x and y
head(names(out)) # band1, band2, ...
tail(names(out)) # ..., band9199, band9200, x, y

# The resolution is 1.5x1.5 degree, so there are 240 grid points around a
# latitude circle, and 121 grid points from pole to pole. Each band is a vector
# of length 240*121=29040, holding the data of one spatial output field. 
nrow(out) # 29040

# To make the output data frame play nicer with plotting functions later, we
# rename x to long and y to lat, and wrap longitudes around so they range from
# -180 to 180 instead of 0 to 360.
setnames(out, old=c('x', 'y'), new=c('long', 'lat'))
out[, long := ifelse(long > 180, long - 360, long)]


# It is not clear from any of the data returned by rgdal::readGDAL which year
# and forecast day corresponds to which band. But in the original grib file
# this metadata is stored. So we have to parse the output of
# gdalUtils::gdalinfo to get this information.
outinfo = gdalUtils::gdalinfo(grib_file)

# outinfo is a vector of about 140000 lines of text. Some lines start on "Band "
# thus indicating that a new band is documented below that line. Let's get the
# indices of these beginning-of-block lines.
band_lines = grep('^Band ', outinfo)
length(band_lines) # 9200 bands, as expected
diff(band_lines) # each descriptor block is 15 lines long

# in line 8 of each block there is a string similar to
# "REF_TIME=2000-06-29T00:00:00Z" which indicates the forecast start date of
# the corresponding band. Let's extract everything between GRIB_REF_TIME= and T
start_dates = grep('GRIB_REF_TIME=', outinfo, value=TRUE)
start_dates = strsplit(start_dates, ' +') %>% map_int(~ as.integer(.[3]))
start_dates = as.POSIXct(start_dates, origin='1970-01-01') %>% format('%Y-%m-%d')


# In line 10 of each block, there is a long string of numbers separated by
# spaces. Positions 25,26,27 correspond to the year/month/day on which the
# forecast verifies.
veri_dates = grep('GRIB_PDS_TEMPLATE_ASSEMBLED_VALUES', outinfo, value=TRUE)
veri_dates = strsplit(veri_dates, ' +') %>% 
             map_chr(~ paste(.[26:28], collapse='-')) %>%
             as.Date %>% format('%Y-%m-%d')


# next we match the initialisation dates and verification dates with their
# band
nbands = ncol(out) - 2 
meta = data.table(band = paste('band', 1:nbands, sep=''),
                  start_date = start_dates,
                  veri_date = veri_dates)

s2s_ecmwf_t2m_ens_meta = meta
s2s_ecmwf_t2m_ens = out


save(file='s2s_ecmwf_t2m_ens.Rdata', list=c('s2s_ecmwf_t2m_ens', 's2s_ecmwf_t2m_ens_meta'))


