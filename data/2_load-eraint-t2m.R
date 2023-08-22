# This R script loads the (large) ERA-Interim temperature reanalysis grib data
# file (2000-2019, 1.5*1.5deg resolution, 4 values per day) using the
# rgdal::readGDAL function.  It extracts the dates from the grib metadata with
# gdalUtils::gdalinfo. The final result is a data table with two colums
# corresponding to long and lat and the other columns corresponding to dates
# (YYYY-MM-DD), and each column is a time instance of daily average temperature
# data indexed by longitude and latitude in the long and lat column. The
# generated file is a 500Mb heavily bzip-compressed Rdata file.

library('tidyverse')
library('rgdal')
library('gdalUtils')
library('data.table')

grib_file = 'eraint-t2m.grib'

# First we parse the output of gdalUtils::gdalinfo to get the dates for each
# band.
outinfo = gdalUtils::gdalinfo(grib_file)

# outinfo is a vector of almost 300k lines of text. Some lines start on "Band "
# followed by a band id, thus indicating that a new band is documented below
# that line. Let's extract the band ids from these beginning-of-block lines.
band_lines = grep('^Band ', outinfo, value=TRUE)
band_id = strsplit(band_lines, ' ') %>% map_int(~ as.integer(.[2]))

# in line 6 of each block there is a string similar to "GRIB_REF_TIME=
# 946706400 sec UTC" which indicates the date of the corresponding band in
# seconds since 1970-01-01. Let's extract everything between GRIB_REF_TIME= and
# sec. There are 4 values per day, but we only save the YMD date so it is
# easier later to use date as a grouping variable and calculate daily averages.
reftime_lines = grep('GRIB_REF_TIME', outinfo, value=TRUE)
reftime_sec = strsplit(reftime_lines, ' +') %>% map_int(~ as.integer(.[3]))
reftime_ymd = as.POSIXct(reftime_sec, origin='1970-01-01') %>% format('%Y-%m-%d')

# clean up
rm(outinfo)

# next we match the dates and bands
meta = data.table(date = reftime_ymd, 
                  iband = band_id,
                  band = paste('band', band_id, sep=''))

# we will extract data from the grib file year by year (R crashed when I tried
# to load everything at once)

meta[, year := as.integer(format(as.Date(date), '%Y'))]
setkey(meta, band)
result = list()
years = sort(unique(meta[, year]))

for (the_year in years) {

  # get bands corresponding to that year
  meta_sub = meta[year == the_year, ]
  meta_sub = meta_sub[order(iband)]
  ibands = meta_sub[, iband]
  bands = meta_sub[, band]

  # load the grib file and transform to data frame
  out = rgdal::readGDAL(grib_file, band = ibands)
  names(out) = bands
  out = as.data.table(out)
  
  # replace column names band -> date
  out_bands = grep('^band.*$', colnames(out), value=TRUE)
  out_dates = meta[out_bands, date]
  setnames(out, old=out_bands, new=out_dates)

  # melt wide into long data.table
  out = melt(out, id.vars = c('x', 'y'), variable.name='date',
             value.name = 't2m')

  # now can group by coordinate and date and calculate daily average
  # temperature
  out = out[, .(t2m = mean(t2m)), list(x, y, date)]

  # cast back into wide table format
  out = dcast(out, x + y ~ date, value.var='t2m')

  # append result and clean up
  result[[paste(the_year)]] = out
  rm(meta_sub, out)
  gc()
}

# merge all data tables in result into one
eraint_t2m = Reduce(merge, result)

# some cosmetics: rename x and y to long and lat, wrap longitude so it ranges
# [-180, 180] instead of [0,360], and convert temperature from K to C
setnames(eraint_t2m, old=c('x', 'y'), new=c('long', 'lat'))
eraint_t2m[, long := ifelse(long > 180, long - 360, long)]
date_cols = grep("long|lat", names(eraint_t2m), value=TRUE, invert=TRUE)
eraint_t2m[, (date_cols) := .SD - 273.15, .SDcols = date_cols]

# save Rdata with decent compression
save(file='eraint-t2m.Rdata', 
     list='eraint_t2m', 
     compress='bzip2', compression_level=9)

# clean up
rm(result, meta, eraint_t2m)
gc()


# also save smaller eraint dataset with only dates for which there are S2S
# forecasts
load('s2s-ecmwf-t2m.Rdata')
veri_dates = unique(s2s_ecmwf_t2m_meta[, veri_date])
eraint_t2m_for_s2s = eraint_t2m[, c('long', 'lat', veri_dates), with=FALSE]

# save Rdata with decent compression
save(file='eraint-t2m-for-s2s.Rdata', 
     list='eraint_t2m_for_s2s', 
     compress='bzip2', compression_level=9)







