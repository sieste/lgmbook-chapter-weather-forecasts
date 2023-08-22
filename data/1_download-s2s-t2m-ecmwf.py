#!/usr/bin/env python

# This script downloads hindcast data from the S2S database for daily averaged
# 2m temperature produced by the ECMWF model. The forecast start dates are in
# variable "hdate" and the forecast steps in hours are in variable "step".

from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "s2",
    "dataset": "s2s",
    "date": "2020-06-29",
    "expver": "prod",
    "hdate": "2000-06-29/2001-06-29/2002-06-29/2003-06-29/2004-06-29/2005-06-29/2006-06-29/2007-06-29/2008-06-29/2009-06-29/2010-06-29/2011-06-29/2012-06-29/2013-06-29/2014-06-29/2015-06-29/2016-06-29/2017-06-29/2018-06-29/2019-06-29",
    "levtype": "sfc",
    "model": "glob",
    "origin": "ecmf",
    "param": "167",
    "step": "0-24/24-48/48-72/72-96/96-120/120-144/144-168/168-192/192-216/216-240/240-264/264-288/288-312/312-336/336-360/360-384/384-408/408-432/432-456/456-480/480-504/504-528/528-552/552-576/576-600/600-624/624-648/648-672/672-696/696-720/720-744/744-768/768-792/792-816/816-840/840-864/864-888/888-912/912-936/936-960/960-984/984-1008/1008-1032/1032-1056/1056-1080/1080-1104",
    "stream": "enfh",
    "time": "00:00:00",
    "type": "cf",
    "target": "s2s-t2m-ecmwf.grib"
})

