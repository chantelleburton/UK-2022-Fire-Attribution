# Testing 1960-2013 HadGEM3 vs ERA5 toolbox FWI data, in maps and in histograms

#module load scitools/default-current
#python3
#-*- coding: iso-8859-1 -*-


import numpy as np
import iris
import datetime
import matplotlib
#matplotlib.use('Agg')
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import iris.analysis
import iris.plot as iplt
import iris.coord_categorisation
import os
import cartopy.io.shapereader as shpreader
from ascend import shape
import iris.analysis.stats
import scipy.stats
from scipy import stats
import ascend
from ascend import shape
import cf_units
import seaborn as sns
import pandas as pd
import glob


############# User inputs here #############
Country = 'United Kingdom'
# Options: 'United Kingdom', 'England', 'Wales', 'Scotland', 'Northern Ireland'
############# User inputs end here #############

### Functions
#Set up the 2022 files and months automatically
def SetCountry(Country):
    if Country == 'United Kingdom' or Country == 'England' or Country == 'Scotland' or Country == 'Wales' or Country == 'Northern Ireland':
        print('Running '+Country)
        daterange = iris.Constraint(time=lambda cell: 6<= cell.point.month <=8)
        month = 'JJA'
        percentile = 90
        ERA5_2022 = iris.load_cube('/scratch/cburton/impactstoolbox/Data/era5/Fire-Weather/FWI-2-day/FWI-2-day_ERA5_std_reanalysis_2022-06-01-2022-08-31_-9.0.50.0.2.0.59.0_day_initialise-from-copernicus:True-and-use-numpy=False.nc')
        return month,percentile,daterange,ERA5_2022

def CountryConstrain(cube, Country):
    if Country == 'United Kingdom':
       shpfilename = str(shpreader.natural_earth(resolution='110m', category='cultural', name='admin_0_countries'))
       natural_earth_file = shape.load_shp(str(shpfilename))
       CountryMask = shape.load_shp(shpfilename, Name=Country)
       Country_shape = CountryMask.unary_union()
       Country1 = Country_shape.mask_cube(cube)
    elif Country == 'England' or Country == 'Scotland' or Country == 'Wales':
       shpfilename = str(shpreader.natural_earth(resolution='110m', category='cultural', name='admin_0_map_units'))
       natural_earth_file = shape.load_shp(str(shpfilename))
       CountryMask = shape.load_shp(shpfilename, Name=Country)
       Country_shape = CountryMask.unary_union()
       Country1 = Country_shape.mask_cube(cube)
    elif Country == 'Northern Ireland':
       shpfilename = str(shpreader.natural_earth(resolution='110m', category='cultural', name='admin_0_map_units'))
       natural_earth_file = shape.load_shp(str(shpfilename))
       CountryMask = shape.load_shp(shpfilename, GEOUNIT=Country)
       Country_shape = CountryMask.unary_union()
       Country1 = Country_shape.mask_cube(cube)
    return Country1 

def CountryMean(cube):
    coords = ('longitude', 'latitude')
    for coord in coords:
        if not cube.coord(coord).has_bounds():
            cube.coord(coord).guess_bounds()
    grid_weights = iris.analysis.cartography.area_weights(cube)
    cube = cube.collapsed(coords, iris.analysis.MEAN, weights = grid_weights)
    return cube 

def TimeMean(cube):
    cube = cube.collapsed('time', iris.analysis.MEAN)
    return cube 

def CountryMax(cube):
    coords = ('longitude', 'latitude')
    cube = cube.collapsed(coords, iris.analysis.MAX)
    return cube 

def CountryPercentile(cube, percentile):
    coords = ('longitude', 'latitude')
    cube = cube.collapsed(coords, iris.analysis.PERCENTILE, percent=percentile)
    return cube 

def TimeMean(cube):
    cube = cube.collapsed('time', iris.analysis.MEAN)
    return cube 

def TimeMax(cube):
    cube = cube.collapsed('time', iris.analysis.MAX)
    return cube 

def TimePercentile(cube, percentile):
    cube = cube.collapsed(['time'], iris.analysis.PERCENTILE, percent=percentile)
    return cube 




def compute_average(*args):
    num_arrays = len(args)
    lengths = set(arr.shape[0] for arr in args)
    #if len(lengths) > 1:
    #    raise ValueError("Arrays have different lengths")
    averages = np.mean(np.stack(args), axis=0)
    return averages


def GetERA5(ERA5_2022,Country):
#Get the ERA5 2022 data for the threshold line
    if Country != 'SAM': #(Already constrained to box before making the data for SAM)
        ERA5_2022 = CountryConstrain(ERA5_2022, Country)
    ERA5_2022 = CountryPercentile(ERA5_2022, percentile)
    ERA5_2022 = TimePercentile(ERA5_2022, percentile)
    ERA5_2022 = np.array(ERA5_2022.data)
    return ERA5_2022


###Set Country###
month,percentile,daterange,ERA5_2022 = SetCountry(Country)
ERA5_2022 = GetERA5(ERA5_2022,Country)
###Set Country###

'''
#Get vars first
vars = ('10u','10v', 'pr', 'rh', 'tasmax')
for var in vars:
    print(var)
    ERA5_Arr = []
    for year in np.arange(1960, 2014):
        ERA5_ImpactsToolBox_Arr = []
        print('ERA5',year)
        ERA5_06 = iris.load_cube('/project/applied/Data/OBS_ERA5/daily/'+var+'/*'+str(year)+'06.nc')
        ERA5_07 = iris.load_cube('/project/applied/Data/OBS_ERA5/daily/'+var+'/*'+str(year)+'07.nc')
        ERA5_08 = iris.load_cube('/project/applied/Data/OBS_ERA5/daily/'+var+'/*'+str(year)+'08.nc')
        cubes = iris.cube.CubeList([ERA5_06, ERA5_07,ERA5_08])
        ERA5 = cubes.concatenate_cube()     
        ERA5 = TimePercentile(ERA5, percentile)
        ERA5 = CountryConstrain(ERA5, Country)
        ERA5 = CountryPercentile(ERA5, percentile)
        ERA5 = np.ravel(ERA5.data)
        ERA5_Arr.append(ERA5)
        f = open('/scratch/cburton/scratch/FWI/2022/BiasCorrection/ERA5_'+var+'_'+Country+'1960-2013_90%.dat','a')
        np.savetxt(f,(ERA5))
        f.close() 
    #Save text out to a file
    f = open('/scratch/cburton/scratch/FWI/2022/BiasCorrection/ERA5_'+var+'_'+Country+'1960-2013_90%_ALL.dat','a')
    np.savetxt(f,(ERA5_Arr))
    f.close() 

print('Done ERA5')


members = ('aojaa', 'aojab', 'aojac', 'aojad', 'aojae', 'aojaf', 'aojag', 'aojah', 'aojai', 'aojaj','dlrja', 'dlrjb', 'dlrjc', 'dlrjd', 'dlrje')
for member in members:
    #HadGEM3_Arr = []
    for year in np.arange(1961, 2014):
        print('HadGEM',member,year)
        HadGEM3 = iris.load('/scratch/cburton/scratch/FWI/2022/historical/1960-2013/'+member+'_apa.pp')
        print("cubelist",HadGEM3)
        daterange = iris.Constraint(time=lambda cell: cell.point.year == year)
        HadGEM3 = HadGEM3.extract(daterange)
        print("One Year",HadGEM3)
        HadGEM3 = HadGEM3.extract(iris.Constraint(name="y_wind"))[0] 
        print("One Var",HadGEM3)
        HadGEM3 = HadGEM3.extract(daterange)
        print("One Month",HadGEM3)
        HadGEM3 = TimePercentile(HadGEM3, percentile)
        print("Time percentile",HadGEM3)
        HadGEM3 = CountryConstrain(HadGEM3, Country)
        HadGEM3 = CountryPercentile(HadGEM3, percentile)
        print("Country percentile",HadGEM3.data)
        HadGEM3 = np.ravel(HadGEM3.data)
        print(HadGEM3)
        #HadGEM3_Arr.append(HadGEM3.data)

        #Save text out to a file
        f = open('/scratch/cburton/scratch/FWI/2022/BiasCorrection/HadGEM_'+Country+'1960-2013_'+member+'ywind_90%.dat','a')
        np.savetxt(f,(HadGEM3))
        f.close()  

exit()

#    pr = data2process.extract(iris.Constraint(name="precipitation_flux"))[0] #Mean Precip
#    tas = data2process.extract(iris.Constraint(name="air_temperature"))[1] #Maximum temp (as per Perry et al 2022 https://nhess.copernicus.org/articles/22/559/2022/ )
#    hurs = data2process.extract(iris.Constraint(name="relative_humidity"))[0] #Mean RH
#    uas = data2process.extract(iris.Constraint(name="x_wind"))[0] #eastward stash #3225 
#    vas = data2process.extract(iris.Constraint(name="y_wind"))[0] #northward stash #3226


'''

def main():
    #Vars first
    Vars = ('Temp', 'Precip', 'RH', 'xwind', 'ywind')
    n=1
    for Var in Vars:
        print(Var)
        ERA5_ImpactsToolBox_File = ('/scratch/cburton/scratch/FWI/2022/BiasCorrection/ERA5_'+Var+'_United Kingdom1960-2013_90%_ALL.dat')
        data = []
        with open(ERA5_ImpactsToolBox_File, 'r') as f:
            d = f.readlines()
            for i in d:
                data.append([float(i)]) 
        ERA5_ImpactsToolBox_Arr = np.array(data, dtype='O')

        folder_path = os.path.join('/scratch/cburton/scratch/FWI/2022/BiasCorrection/') 
        filenames = ['aojaa','aojab', 'aojac', 'aojad']#, 'aojae', 'aojaf', 'aojag', 'aojah', 'aojai', 'aojaj','dlrja', 'dlrjb', 'dlrjc', 'dlrjd', 'dlrje']
        Data = []
        for filename in filenames:
            data = []
            files = (folder_path+'HadGEM_United Kingdom1960-2013_'+filename+Var+'_90%.dat')
            with open(files, 'r') as f:
                d = f.readlines()
                for i in d:
                    data.append([float(i)]) 
            array = np.array(data, dtype='O')
            Data.append(array)
        averages = compute_average(*Data)

        if Var == 'Temp':
            averages = averages-273.15
            ylabel = ('deg C')
            title = ('Maximum Daily Temperature')
        if Var == 'Precip':
            averages = averages*86400
            ylabel = ('mm/day')
            title = ('Daily Precip')
        if Var == 'RH':
            title = ('RH')
            ylabel = ('%')
        if Var == 'xwind':
            ylabel = ('m/s')
            title = ('xwind')
        if Var == 'ywind':
            ylabel = ('m/s')
            title = ('ywind')

        plt.subplot(3,2,n)
        plt.plot(averages, color = 'orange', label='HadGEM3')
        plt.plot(ERA5_ImpactsToolBox_Arr, color = 'black', label='ERA5')
        plt.ylabel(ylabel)
        plt.title(title)
        if n == 5:
            years = ('1960','1970','1980','1990','2000','2010')
            x_pos = (0, 10, 20, 30, 40, 50)
            plt.xticks(x_pos, years)
        else:
            years = (' ',' ',' ',' ',' ',' ')
            x_pos = (0, 10, 20, 30, 40, 50)
            plt.xticks(x_pos, years)
        n = n+1

    ##Now for FWI
    print('FWI')
    folder_path = '/scratch/cburton/scratch/FWI/2022/historical/Array/' 
    Data = []
    for filename in filenames:
        data = []
        files = (folder_path+'HadGEM3_Arr_'+Country+'1960-2013_'+filename+'_90%.dat')
        with open(files, 'r') as f:
            d = f.readlines()
            for i in d:
                data.append([float(i)]) 
        array = np.array(data, dtype='O')
        Data.append(array)

    averages = compute_average(*Data)


    ERA5_ImpactsToolBox_File = ('/scratch/cburton/scratch/FWI/2022/ERA5/Array/ERA5_ImpactsToolBox_Arr_'+Country+'1960-2013_90%.dat')
    data = []
    with open(ERA5_ImpactsToolBox_File, 'r') as f:
        d = f.readlines()
        for i in d:
            data.append([float(i)]) 
    ERA5_ImpactsToolBox_Arr = np.array(data, dtype='O')


    plt.subplot(3,2,6)
    plt.plot(averages, color = 'orange', label='HadGEM3')
    plt.plot(ERA5_ImpactsToolBox_Arr, color = 'black', label='ERA5')
    plt.title('FWI')
    years = ('1960','1970','1980','1990','2000','2010')
    x_pos = (0, 10, 20, 30, 40, 50)
    plt.xticks(x_pos, years)
    plt.legend(loc='best')

    plt.show()


if __name__ == "__main__":
    main()









