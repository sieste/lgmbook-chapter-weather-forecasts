# Downloading S2S data

The S2S forecast data is available at [apps.ecmwf.int/datasets](apps.ecmwf.int/datasets)

To be able to download data, you have to register for an account with your email address.

Then login and follow the link to the S2S data. 

Use the sidebar to navigate to the forecast product you want to download. 

For example, for the first data set I downloaded I used the sidebar options 

- Reforecasts
- Daily averaged
- ECMWF
- Surface
- Control forecast 

and in the user interface I selected 

- the latest model version
- all hindcast dates (2000-2019)
- all forecast steps
- and the 2-metre temperature forecast variable.

Then I clicked "View data retreival request" and copied the python code into a
text file `download-s2s-t2m-ecmwf.py`. Running python script downloads the
selected data from the data base.

First of all, to be able to run the python script I had to install the "ECMWF
API client" python module with

```
pip install --user ecmwf-api-client
```

There is lots of useful information on the api client at
[pypi.org/project/ecmwf-api-client](pypi.org/project/ecmwf-api-client). Most
importantly, there is a link back to the ECMWF site to look up your API key
[api.ecmwf.int/v1/key/](https://api.ecmwf.int/v1/key/). Follow the instructions
and save the api key information in a file called `.ecmwfapirc` in your home
directory.

Then the python script can be executed with 

```
python s2s-t2m-ecmwf.py
```

The downloaded reforecast data is 50Mb.


# Loading the data in R

The file `load-s2s-t2m-ecmwf.R` is heavily commented and shows how to load the
grib file and rearrange the data into a tidy data frame format in R.





