#!/usr/bin/env python
import os
import sys
import datetime
import numpy
import math
import cdsapi
from joblib import Parallel, delayed
from time import sleep
from random import randint

#disable unnecessary warnings
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


#declare function to get data for each day
def getECMWF(Year,Day):
  
  Year         = int(Year)
  Day          = int(Day)


  FilePath = format(Year, '04d') +'/'
  FileName = 'era5_' + format(Year, '04d') + 'd' + format(Day,  '03d')
  FileName += '.nc'

  DateString = format(Year, '04d')+'-'+format(1, '02d')+'-'+format(Day, '02d')
  
  DateString = datetime.datetime(Year,1,1)+datetime.timedelta(days=Day-1)
  DateString = DateString.strftime('%Y-%m-%d')

  print('================================================')
  print(FilePath+FileName)
  print('================================================')


  if not os.path.exists(format(Year, '04d') +'/'):
    os.makedirs(format(Year, '04d') +'/')
    
  if not os.path.exists(FilePath+FileName):

      #space out requests
      sleep(randint(10,180))

      #try:
        #get the data!
      c = cdsapi.Client()
      c.retrieve('reanalysis-era5-complete', {    # do not change this!
                "class": "ea",
                "dataset": "era5",
                "date": DateString,
                "expver": "1",
                "levelist": "1/to/137",
                "levtype": "ml",
                "param": "130/131/132",
                "stream": "oper",
                "time": "00:00:00/03:00:00/06:00:00/09:00:00/12:00:00/15:00:00/18:00:00/21:00:00",
                "type": "an",
                "format": "netcdf",
                "grid": "1.5/1.5",
                "area": "90/-180/-90/180",
                }, FilePath+FileName)
        
      #except:
        #print('Error downloading file '+FileName)
      
  return



##call function to download data
#Years = [2018]
#Days  = range(213,243)
#Steps = [0,3,6,12,24,48];
#Times  = [0,12]

#call function to download data
Years = range(2019,2018,-1)
Days  = range(244,350,1)




DaysList, YearList = numpy.meshgrid(Days,Years)
DaysList = numpy.ndarray.flatten(DaysList)
YearList = numpy.ndarray.flatten(YearList)




results = Parallel(n_jobs=60, verbose=11, backend="threading")(map(delayed(getECMWF),YearList,DaysList))
