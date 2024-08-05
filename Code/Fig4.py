
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
from ascend import shape
import iris.analysis.stats
import scipy.stats
from scipy import stats
import ascend
from ascend import shape
import cf_units
import seaborn as sns
from scipy.stats import genextreme as gev, kstest
import pandas as pd
import statsmodels.api as sm
from pdb import set_trace



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

def CountryMax(cube):
    coords = ('longitude', 'latitude')
    for coord in coords:
        if not cube.coord(coord).has_bounds():
            cube.coord(coord).guess_bounds()
    grid_weights = iris.analysis.cartography.area_weights(cube)
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

def RiskRatio(Alldata,Natdata, Threshold):
    ALL = (np.count_nonzero(Alldata > ERA5_2022)/len(Alldata))
    NAT = (np.count_nonzero(Natdata > ERA5_2022)/len(Natdata))
    RR = ALL/NAT
    return RR

def FAR(Alldata,Natdata, Threshold):
    ALL = (np.count_nonzero(Alldata > ERA5_2022))
    NAT = (np.count_nonzero(Natdata > ERA5_2022))
    FAR = 1 - (NAT/ALL)
    return FAR

def RT(data, Threshold):
    DATA = (np.count_nonzero(data > ERA5_2022)/len(data))
    RT = 1 / DATA
    return RT

def draw_bs_replicates(ALL, NAT, ERA5, func, size):
    """creates a bootstrap sample, computes replicates and returns replicates array"""
    # Create an empty array to store replicates
    RR_replicates = np.empty(size)
    
    # Create bootstrap replicates as much as size
    for i in range(size):
        # Create a bootstrap sample
        ALL_sample = np.random.choice(ALL,size=(int(np.round(len(ALL)-(0.1*len(ALL))))), replace=False)
        ALL_sample = np.random.choice(ALL_sample,size=len(ALL), replace=True)
        NAT_sample = np.random.choice(NAT,size=(int(np.round(len(NAT)-(0.1*len(NAT))))), replace=False)
        NAT_sample = np.random.choice(NAT_sample,size=len(ALL), replace=True)
        # Get bootstrap replicate and append to bs_replicates
        RR_replicates[i] = func(ALL_sample, NAT_sample, ERA5) 
    return RR_replicates

def bootstrap(data, ERA5, func, size):
    """creates a bootstrap sample, computes replicates and returns replicates array"""
    # Create an empty array to store replicates
    RT_replicates = np.empty(size)
    
    # Create bootstrap replicates as much as size
    for i in range(size):
        # Create a bootstrap sample
        DATA_sample = np.random.choice(data,size=(int(np.round(len(data)-(0.1*len(data))))), replace=False)
        DATA_sample = np.random.choice(DATA_sample,size=len(data), replace=True)
        # Get bootstrap replicate and append to bs_replicates
        RT_replicates[i] = func(DATA_sample, ERA5) 
    return RT_replicates



def GetERA5(ERA5_2022,Country):
#Get the ERA5 2022 data for the threshold line
    if Country != 'SAM': #(Already constrained to box before making the data for SAM)
        ERA5_2022 = CountryConstrain(ERA5_2022, Country)
    ERA5_2022 = CountryPercentile(ERA5_2022, percentile)
    ERA5_2022 = TimePercentile(ERA5_2022, percentile)
    ERA5_2022 = np.array(ERA5_2022.data)
    return ERA5_2022



############## 1) Create .dat files and save out to save time in plotting #################
'''
folder = '/scratch/cburton/scratch/FWI/2022/'
index_filestem1 = 'hist/'
index_filestem2 = 'histnat/'
index_name = 'Canadian_FWI'


## For each historical member, bias correct 525 members for 2022 and save out
BiasCorrDict = {}
histmembers = ('aojaa','aojab', 'aojac', 'aojad', 'aojae', 'aojaf','aojag', 'aojah', 'aojai', 'aojaj','dlrja', 'dlrjb', 'dlrjc', 'dlrjd', 'dlrje') 
for histmember in histmembers:
    print(histmember)
    # Step 0; Load fwi data from CSV using pandas
    df_obs = pd.read_csv('/scratch/cburton/scratch/FWI/2022/ERA5/Array/ERA5_ImpactsToolBox_Arr_'+Country+'1960-2013_'+str(percentile)+'%.dat')
    df_sim = pd.read_csv('/scratch/cburton/scratch/FWI/2022/historical/Array/HadGEM3_Arr_'+Country+'1960-2013_'+histmember+'_'+str(percentile)+'%.dat')
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

    
    #### First do for hist array  ####
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

                ####Log transform the data here#### 
                hist = np.log(np.exp(hist)-1)

                # Step 2: Detrend the sim and scale to obs
                Endhist = fwi0_obs + (hist - delta_sim * 0 - fwi0_sim)
                print(Endhist)

                ####inverse Log (exponential) transform here####      
                Endhist = np.log(np.exp(Endhist)+1)
 
                f = open('/scratch/cburton/scratch/FWI/2022/hist/Array/'+Country+'_2022code+BC1'+histmember+'_hist'+str(percentile)+'%_LogTransform_-1+1.dat','a')
                np.savetxt(f,(Endhist),newline=',',fmt='%s')
                f.write('\n')
                f.close()
                histarray.append(hist)
            except IOError:
                 pass 
     
    histarray = np.array(histarray)
    histarray = np.ravel(histarray)
    print(repr(histarray)) 
exit()
    
        
    
    ##### Repeat for histnat array (can run this separately in parallel to save time) ####
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

                ####Log transform the data here#### 
                histnat = np.log(np.exp(histnat)-1)

                # Step 2: Detrend the sim and scale to obs
                Endhist = fwi0_obs + (histnat - delta_sim * 0 - fwi0_sim)

                ####inverse Log (exponential) transform here####      
                Endhist = np.log(np.exp(Endhist)+1)

                f = open('/scratch/cburton/scratch/FWI/2022/histnat/Array/'+Country+'_2022code+BC1'+histmember+'_histnat'+str(percentile)+'%_LogTransform_-1+1.dat','a')
                np.savetxt(f,(Endhist),newline=',',fmt='  %s')
                f.write('\n')
                f.close()
                histnatarray.append(histnat)
            except IOError:
                pass 
        
    histnatarray = np.array(histnatarray)
    histnatarray = np.ravel(histnatarray)

exit()


'''
'''
############## 2a) Create plot for UK #################

month,percentile,daterange,ERA5_2022 = SetCountry(Country)
ERA5_2022 = GetERA5(ERA5_2022,Country)

#Read in data files for each historical member (each one has 525 values for the 2022 data), for ALL and NAT
AllDict = {}
NatDict = {}
AllDict[Country] = []  
NatDict[Country] = []  

members = ('aojaa', 'aojab', 'aojac', 'aojad', 'aojae', 'aojaf', 'aojag', 'aojah', 'aojai', 'aojaj','dlrja', 'dlrjb', 'dlrjc', 'dlrjd', 'dlrje')   
for member in members:
    print(member)
    all_file = '/scratch/cburton/scratch/FWI/2022/hist/Array/United Kingdom_2022code+BC1'+member+'_hist'+str(percentile)+'%_LogTransform_-1+1.dat'
    nat_file = '/scratch/cburton/scratch/FWI/2022/histnat/Array/United Kingdom_2022code+BC1'+member+'_histnat'+str(percentile)+'%_LogTransform_-1+1.dat'
    with open(all_file) as f:
         all_lines=f.readlines()
    AllDict[Country] += [float(line.rstrip(',\n')) for line in all_lines]
    with open(nat_file) as f:
         nat_lines=f.readlines()
    NatDict[Country] += [float(line.rstrip(',\n')) for line in nat_lines]

#Make sure they are arrays, so we can plot them
NatDict[Country] = np.array(NatDict[Country])
AllDict[Country] = np.array(AllDict[Country])
print(len(NatDict[Country]))

### Bootstrap and print the Risk Ratio Results when cycling through each country ###
RR = draw_bs_replicates(AllDict[Country], NatDict[Country], ERA5_2022, RiskRatio, 10000)
print('Risk Ratio')
print(len(RR))
print(np.percentile(RR, 5))
print(np.percentile(RR, 95))

### Bootstrap and print the RT Results when cycling through each country ###
RT_ALL = bootstrap(AllDict[Country], ERA5_2022, RT, 10000)
RT_NAT = bootstrap(NatDict[Country], ERA5_2022, RT, 10000)



#Make box and whisker plot for UK return times, using 5, 50 and 95 percentiles for RT_ALL and RT_NAT
plt.boxplot([RT_ALL, RT_NAT], positions=[1,2])
plt.xticks([1,2], ['ALL', 'NAT'])
plt.ylabel('Return Time (years)')
plt.title('Return Time for UK 2022 JJA')
plt.show()


exit()
'''



############## 2b) Create 5 subplots #################

Countries = ('United Kingdom','England', 'Wales', 'Scotland', 'Northern Ireland')
n = 1
for Country in Countries:
    print(Country)
    month,percentile,daterange,ERA5_2022 = SetCountry(Country)
    ERA5_2022 = GetERA5(ERA5_2022,Country)

    #Read in data files for each historical member (each one has 525 values for the 2022 data), for ALL and NAT
    AllDict = {}
    NatDict = {}
    AllDict[Country] = []  
    NatDict[Country] = []  

    members = ('aojaa', 'aojab', 'aojac', 'aojad', 'aojae', 'aojaf', 'aojag', 'aojah', 'aojai', 'aojaj','dlrja', 'dlrjb', 'dlrjc', 'dlrjd', 'dlrje')   
    for member in members:
        all_file = '/scratch/cburton/scratch/FWI/2022/hist/Array/'+Country+'_2022code+BC1'+member+'_hist'+str(percentile)+'%_LogTransform_-1+1.dat'
        nat_file = '/scratch/cburton/scratch/FWI/2022/histnat/Array/'+Country+'_2022code+BC1'+member+'_histnat'+str(percentile)+'%_LogTransform_-1+1.dat'
        with open(all_file) as f:
             all_lines=f.readlines()
        AllDict[Country] += [float(line.rstrip(',\n')) for line in all_lines]
        with open(nat_file) as f:
             nat_lines=f.readlines()
        NatDict[Country] += [float(line.rstrip(',\n')) for line in nat_lines]

    #Make sure they are arrays, so we can plot them
    NatDict[Country] = np.array(NatDict[Country])
    AllDict[Country] = np.array(AllDict[Country])

    ### Bootstrap and print the RT Results when cycling through each country ###
    RT_ALL = bootstrap(AllDict[Country], ERA5_2022, RT, 10000)
    RT_NAT = bootstrap(NatDict[Country], ERA5_2022, RT, 10000)
    print(np.percentile(RT_ALL, 50))
    print(np.percentile(RT_NAT, 50))

    #Make the plot
    plt.subplot(2,3,n)
    #plt.boxplot([RT_ALL, RT_NAT], positions=[1,2]) #without any changes to whiskers / fliers
    bp1 = plt.boxplot([RT_ALL, RT_NAT], positions=[1,2], widths=0.6, whiskerprops=dict(linewidth=0))
    bp2 = plt.boxplot([RT_ALL, RT_NAT], positions=[1,2], widths=0.6, patch_artist=True, 
                      boxprops=dict(visible=False), medianprops=dict(visible=False), 
                      capprops=dict(visible=False), flierprops=dict(visible=False))
    plt.setp(bp2['whiskers'], linewidth=2)

    plt.xticks([1,2], [' ', ' '])
    plt.yscale('log') # add log scale
    plt.ylabel(' ')
    plt.title(Country)
    if n == 1 or n == 4:
        plt.ylabel('Return Time (years)')
    if n == 4 or n == 5:
        plt.xticks([1,2], ['ALL', 'NAT'])
    n = n+1

plt.show()

exit()




'''
PRINTED RESULTS


- UK  
6.24
11.187802419354835

- Scotland
4.96
11.384615384615385

- NI
3.3419467173193063
4.128583535108959

- Wales
8.867924528301886
15.406525735294114

- England
6.253544776119403
10.487179487179487


'''

