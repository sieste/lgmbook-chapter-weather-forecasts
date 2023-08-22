#!/usr/bin/env python
# This script downloads ERA-Interim 2m temperature data from the ECMWF data
# portal. Four values daily from 01-01-2000 to the end of the ERA-Interim
# project period are downloaded on a 1.5x1.5 degree lat-lon grid.  The
# downloaded data file is about 1.6Gb.

from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()

server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "2000-01-01/to/2019-08-31",
    "expver": "1",
    "grid": "1.5/1.5",
    "levtype": "sfc",
    "param": "167.128",
    "step": "0",
    "stream": "oper",
    "time": "00:00:00/06:00:00/12:00:00/18:00:00",
    "type": "an",
    "target": "eraint-t2m.grib",
})


