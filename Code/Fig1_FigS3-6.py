#module load scitools/default-current
#python3
#-*- coding: iso-8859-1 -*-


import numpy as np
import iris
import datetime
import matplotlib
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import iris.analysis
import iris.plot as iplt
import iris.coord_categorisation
import os
import cartopy.io.shapereader as shpreader
import shapefile as shape
from ascend import shape
import iris.analysis.stats
import scipy.stats
from scipy import stats
import ascend
from ascend import shape
import cf_units
import seaborn as sns
import pandas as pd
import statsmodels.api as sm




############# User inputs here #############
Country = 'Wales'
inputs_dir = "/scratch/cburton/scratch/FWI/2022/"
outputs_dat_dir = "/data/dynamic/dkelley/UK-2022-Fire-Attribution/datFiles/"
outputs_dir = '/data/dynamic/dkelley/UK-2022-Fire-Attribution/outputs/'
ERA5_reanalysis_file = '/scratch/cburton/impactstoolbox/Data/era5/Fire-Weather/FWI-2-day/FWI-2-day_ERA5_std_reanalysis_2022-06-01-2022-08-31_-9.0.50.0.2.0.59.0_day_initialise-from-copernicus:True-and-use-numpy=False.nc'
# Options: 'United Kingdom', 'England', 'Wales', 'Scotland', 'Northern Ireland'
############# User inputs end here #############


#Set up the 2022 files and months automatically
if Country == 'United Kingdom' or Country == 'England' or Country == 'Scotland' or Country == 'Wales' or Country == 'Northern Ireland':
    print('Running '+Country)
    daterange = iris.Constraint(time=lambda cell: 6<= cell.point.month <=8)
    month = 'JJA'
    percentile = 90
    ERA5_2022 = iris.load_cube(ERA5_reanalysis_file)


### Functions
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




#############  Subplot (a) - Historical PDF uncorrected ######### 

### Make the .dat files first ###

# ERA5 from toolbox
ERA5_ImpactsToolBox_Arr = []
for year in np.arange(1960, 2014):
    print('ERA5',year)
    ERA5_ImpactsToolBox = iris.load_cube(inputs_dir + 'ERA5/FWI_ccra3_'+str(year)+'.nc')
    ERA5_ImpactsToolBox = ERA5_ImpactsToolBox.extract(daterange)
    ERA5_ImpactsToolBox = TimePercentile(ERA5_ImpactsToolBox, percentile)
    ERA5_ImpactsToolBox = CountryConstrain(ERA5_ImpactsToolBox, Country)
    ERA5_ImpactsToolBox = CountryPercentile(ERA5_ImpactsToolBox, percentile)
    ERA5_ImpactsToolBox = np.ravel(ERA5_ImpactsToolBox.data)
    ERA5_ImpactsToolBox_Arr.append(ERA5_ImpactsToolBox)
    print(ERA5_ImpactsToolBox)

#Save ERA5 text out to a file
f = open(outputs_dat_dir + '/ERA5/Array/ERA5_ImpactsToolBox_Arr_'+Country+'1960-2013_'+str(percentile)+'%.dat','a')
np.savetxt(f,(ERA5_ImpactsToolBox_Arr))
f.close()    
exit()


# HadGEM3
members = ('aojaa', 'aojab', 'aojac', 'aojad', 'aojae', 'aojaf', 'aojag', 'aojah', 'aojai', 'aojaj','dlrja', 'dlrjb', 'dlrjc', 'dlrjd', 'dlrje')   
for member in members:
    HadGEM3_Arr = []
    for year in np.arange(1960, 2014):
        print('HadGEM',member,year)
        HadGEM3 = iris.load_cube(inputs_dir + 'historical/1960-2013/FWI_'+member+'.nc')
        print(HadGEM3)
        if HadGEM3 == None:
            pass
        else:
           daterange = iris.Constraint(time=lambda cell: cell.point.year == year)
        if HadGEM3 == None:
            pass
        else:
            HadGEM3 = HadGEM3.extract(daterange)
        if HadGEM3 == None:
            pass
        else:            
            HadGEM3 = HadGEM3.extract(daterange)
        if HadGEM3 == None:
            pass
        else:
            print('Selected times :\n' + str(HadGEM3.coord('time')))
            HadGEM3 = TimePercentile(HadGEM3, percentile)
            HadGEM3 = CountryConstrain(HadGEM3, Country)
            HadGEM3 = CountryPercentile(HadGEM3, percentile)
            HadGEM3 = np.ravel(HadGEM3.data)
            HadGEM3_Arr.append(HadGEM3)


    #Save HadGEM3 text out to a file
    f = open(outputs_dat_dir + 'historical/Array/HadGEM3_Arr_'+Country+'1960-2013_'+member+'_'+str(percentile)+'%.dat','a')
    np.savetxt(f,(HadGEM3_Arr))
    f.close()  

exit()


### Read the .dat files  ###
ERA5_ImpactsToolBox_File = (outputs_dat_dir + 'ERA5/Array/ERA5_ImpactsToolBox_Arr_'+Country+'1960-2013_'+str(percentile)+'%.dat')
data = []
with open(ERA5_ImpactsToolBox_File, 'r') as f:
    d = f.readlines()
    for i in d:
        data.append([float(i)]) 
ERA5_ImpactsToolBox_Arr = np.array(data, dtype='O')

members = ('aojaa', 'aojab', 'aojac', 'aojad', 'aojae', 'aojaf', 'aojag', 'aojah', 'aojai', 'aojaj', 'dlrja', 'dlrjb', 'dlrjc', 'dlrjd', 'dlrje')
data = []
for member in members:
    HadGEM3_File = (outputs_dat_dir + 'historical/Array/HadGEM3_Arr_'+Country+'1960-2013_'+member+'_'+str(percentile)+'%.dat')
    with open(HadGEM3_File, 'r') as f:
         d = f.readlines()
         for i in d:
            data.append([float(i)]) 
HadGEM3_Arr = np.array(data, dtype='O')

#Get the ERA5 2022 data for the threshold line
if Country != 'SAM': #(Already constrained to box before making the data for SAM)
    ERA5_2022 = CountryConstrain(ERA5_2022, Country)
ERA5_2022 = CountryPercentile(ERA5_2022, percentile)
ERA5_2022 = TimePercentile(ERA5_2022, percentile)
ERA5_2022 = np.array(ERA5_2022.data)
print(ERA5_2022)


###  Make the plot  ###
plt.subplot(2,2,1)
sns.distplot(HadGEM3_Arr, hist=True, kde=True, 
             color = 'yellow', fit_kws={"linewidth":2.5,"color":"orange"}, label='HadGEM3')

sns.distplot(ERA5_ImpactsToolBox_Arr, hist=True, kde=True, 
             color = 'grey', fit_kws={"linewidth":2.5,"color":"black"}, label='ERA5')

plt.axvline(x=ERA5_2022, color='black', linewidth=2.5, label='ERA5 '+month+' 2022')

if Country == 'United Kingdom' or Country == 'England':
    plt.yticks([0.0,0.05,0.1,0.15,0.2])#UK and ENG
elif Country == 'Northern Ireland' or Country == 'Scotland':
    plt.yticks([0.0,0.2,0.4,0.6,0.8])#NI, Scot
elif Country == 'Wales':
    plt.yticks([0.0,0.1,0.2,0.3,0.4])#Wales
plt.xlabel(' ')
plt.title('a) '+month+' 1960-2013 (Uncorrected)')
plt.legend(loc='best')



############ Subplot (b) - Historical PDF bias-correctd and transformed ########### 
BiasCorrDict = {}
FWI_SIM = {}
members = ('aojaa', 'aojab', 'aojac', 'aojad', 'aojae', 'aojaf', 'aojag', 'aojah', 'aojai', 'aojaj','dlrja', 'dlrjb', 'dlrjc', 'dlrjd', 'dlrje')   
for member in members:
    print(member)
    # Step 0; Load fwi data from CSV using pandas
    df_obs = pd.read_csv(outputs_dat_dir + 'ERA5/Array/ERA5_ImpactsToolBox_Arr_'+Country+'1960-2013_'+str(percentile)+'%.dat')
    df_sim = pd.read_csv(outputs_dat_dir + '/historical/Array/HadGEM3_Arr_'+Country+'1960-2013_'+member+'_'+str(percentile)+'%.dat')
    df_obs[np.isnan(df_obs)] = 0.000000000001
    df_sim[np.isnan(df_sim)] = 0.000000000001 

    ####Log transform the data here#### 
    df_obs = np.log(np.exp(df_obs)-1)
    df_sim = np.log(np.exp(df_sim)-1)

    # Extract years and FWI values
    years = np.arange(1960,2013)
    fwi_sim = df_sim.values
    fwi_sim = fwi_sim[:,0]
    fwi_obs = df_obs.values
    fwi_obs = fwi_obs[:,0]

    print(len(fwi_sim))
    print(len(fwi_obs))
    print(len(years))

    # Step 1a: Fit a linear regression model to obs and sim
    t = years - 2022  # shift years to be relative to 2022
    X = sm.add_constant(t)  # add a constant term for intercept
    def find_regression_parameters(fwi):
        model = sm.OLS(fwi, X)
        results = model.fit()

        # Step 1b: Get the coefficients (slope and intercept)
        fwi0, delta = results.params
        return fwi0, delta, np.std(fwi - delta * t) 

    fwi0_sim, delta_sim, std_sim =  find_regression_parameters(fwi_sim)
    fwi0_obs, delta_obs, std_obs =  find_regression_parameters(fwi_obs)

    # Step 2: Detrend the sim and scale to obs
    BiasCorrDict[member] = fwi0_obs + (fwi_sim - delta_sim * t - fwi0_sim) 
    fwi_detrended_sim = BiasCorrDict[member]

fwi_detrended_sim_ENSEMBLE = [BiasCorrDict['aojaa'],BiasCorrDict['aojab'],BiasCorrDict['aojac'],BiasCorrDict['aojad'],BiasCorrDict['aojae'],BiasCorrDict['aojaf'],BiasCorrDict['aojag'],BiasCorrDict['aojah'],BiasCorrDict['aojai'],BiasCorrDict['aojaj'],BiasCorrDict['dlrja'],BiasCorrDict['dlrjb'],BiasCorrDict['dlrjc'],BiasCorrDict['dlrjd'],BiasCorrDict['dlrje']]
fwi_detrended_sim_ENSEMBLE = np.ravel(fwi_detrended_sim_ENSEMBLE)

###Inverse of log transform here###
fwi_detrended_sim_ENSEMBLE = np.log(np.exp(fwi_detrended_sim_ENSEMBLE)+1)
fwi_obs = np.log(np.exp(fwi_obs)+1)
fwi_sim = np.log(np.exp(fwi_sim)+1)

##Make the PDF plot
plt.subplot(2,2,2)
sns.distplot(fwi_detrended_sim_ENSEMBLE, kde_kws={"clip": (0, None)}, color = 'yellow',   label='HadGEM3')
sns.distplot(fwi_obs, kde=True, kde_kws={"clip": (0, None)}, color = 'grey',  label='ERA5')
plt.axvline(x=ERA5_2022, color='black', linewidth=2.5, label='ERA5 '+month+' 2022')
plt.xlabel(' ')
if Country == 'United Kingdom' or Country == 'England':
    plt.yticks([0.0,0.05,0.1,0.15,0.2])#UK and ENG
elif Country == 'Northern Ireland' or Country == 'Scotland':
    plt.yticks([0.0,0.2,0.4,0.6,0.8])#NI, Scot
elif Country == 'Wales':
    plt.yticks([0.0,0.1,0.2,0.3,0.4])#Wales
plt.title('b) '+month+' 1960-2013 (Corrected)')
plt.legend(loc='best')



######### Subplot (c) - Timeseries of bias correction ######### 
plt.subplot(2,2,3)
plt.plot(years, fwi_obs, label='ERA5', color='blue')
plt.plot(years, fwi_sim, label='HadGEM3', color='red')
plt.plot(years, fwi_detrended_sim, label='Detrended & Shifted FWI', color='purple')
plt.xlabel('Year')
plt.ylabel('FWI')
plt.title('c) Time Series Plot of FWI and Detrended & Shifted FWI')
plt.legend()
plt.grid(True)


######### Subplot (d) - Uncorrected 2022 ######### 

### First make the .dat files for hist and histnat (can run each one in paralell to save time) ###

folder = '/scratch/cburton/scratch/FWI/2022/'
index_filestem1 = 'hist/'
index_filestem2 = 'histnat/'
index_name = 'Canadian_FWI'

#HIST
members = np.arange(1,106)
histarray = []
for member in members:
    print ('hist',member)
    for n in np.arange(1,6):
        try:
            if member < 10:
                hist = iris.load_cube(folder+index_filestem1+'/FWI_00'+str(member)+'_'+str(n)+'*.nc', index_name)
            elif member > 9 and member < 100:
                hist = iris.load_cube(folder+index_filestem1+'/FWI_0'+str(member)+'_'+str(n)+'*.nc', index_name)
            else:
                hist = iris.load_cube(folder+index_filestem1+'/FWI_'+str(member)+'_'+str(n)+'*.nc', index_name)
            print(hist)
            hist = CountryConstrain(hist, Country)
            hist = CountryPercentile(hist, percentile)
            hist = TimePercentile(hist, percentile)
            hist = np.ravel(hist.data)
            f = open(outputs_dat_dir + 'hist/Array/'+Country+'_UNCORRECTED_hist'+str(percentile)+'%.dat','a')
            np.savetxt(f,(hist),newline=',',fmt='%s')
            f.write('\n')
            f.close()
            histarray.append(hist)
        except IOError:
             pass 
     
histarray = np.array(histarray)
histarray = np.ravel(histarray)
print(repr(histarray)) 
exit()


folder = '/scratch/cburton/scratch/FWI/2022/'
index_filestem1 = 'hist/'
index_filestem2 = 'histnat/'
index_name = 'Canadian_FWI'

#HISTNAT
histnatarray = []
members = np.arange(1,106)
for member in members:
    print ('histnat',member)
    for n in np.arange(1,6):
        try:
            if member < 10:
                histnat = iris.load_cube(folder+index_filestem2+'/FWI_00'+str(member)+'_'+str(n)+'*.nc', index_name)
            elif member > 9 and member < 100:
                histnat = iris.load_cube(folder+index_filestem2+'/FWI_0'+str(member)+'_'+str(n)+'*.nc', index_name)
            else:
                histnat = iris.load_cube(folder+index_filestem2+'/FWI_'+str(member)+'_'+str(n)+'*.nc', index_name)    
            histnat = CountryConstrain(histnat, Country)
            histnat = CountryPercentile(histnat, percentile)
            histnat = TimePercentile(histnat, percentile)
            histnat = np.ravel(histnat.data)
            f = open(outputs_dat_dir + '/histnat/Array/'+Country+'_UNCORRECTED_histnat'+str(percentile)+'%.dat','a')
            np.savetxt(f,(histnat),newline=',',fmt='  %s')
            f.write('\n')
            f.close()
            histnatarray.append(histnat)
        except IOError:
            pass 
        
histnatarray = np.array(histnatarray)
histnatarray = np.ravel(histnatarray)
print(repr(histnatarray))
exit()

### Then read in the .dat files and make the plot ###

ALL = outputs_dat_dir + 'hist/Array/'+Country+'_UNCORRECTED_hist'+str(percentile)+'%.dat'
NAT = outputs_dat_dir + 'histnat/Array/'+Country+'_UNCORRECTED_histnat'+str(percentile)+'%.dat'

All_array = []
Nat_array = []
with open(ALL) as f:
     all_lines=f.readlines()
All_array += [float(line.rstrip(',\n')) for line in all_lines]

with open(NAT) as f:
     nat_lines=f.readlines()
Nat_array += [float(line.rstrip(',\n')) for line in nat_lines]

NatDict = np.array(Nat_array)
AllDict = np.array(All_array)


plt.subplot(2,2,4)
sns.distplot(AllDict, hist=True, kde=True, 
             color = 'orange', fit_kws={"linewidth":2.5,"color":"orange"}, label='ALL')
sns.distplot(NatDict, hist=True, kde=True, 
             color = 'blue', fit_kws={"linewidth":2.5,"color":"blue"}, label='NAT')
plt.axvline(x=ERA5_2022, color='black', linewidth=2.5, label='ERA5 '+month+' 2022')
plt.xlabel('FWI')
if Country == 'United Kingdom':
    plt.yticks([0.0,0.05,0.1,0.15,0.2])#UK 
elif Country == 'Wales':
    plt.yticks([0.0,0.1,0.2,0.3,0.4])#Wales
plt.title('d) '+month+' 2022 Uncorrected')
plt.legend()
plt.suptitle(str(Country)+' '+str(percentile)+'th percentile FWI')
plt.show()










