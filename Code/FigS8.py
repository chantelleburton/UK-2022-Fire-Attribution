# Testing 2022 HadGEM3 vs ERA5 in Maps
# Testing ERA5 toolbox vs ERA5 Copernicus FWI data, in maps 


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
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator


daterange = iris.Constraint(time=lambda cell: 6<= cell.point.month <=8)
month = 'JJA'
percentile = 90
ERA5_2022 = iris.load_cube('/scratch/cburton/impactstoolbox/Data/era5/Fire-Weather/FWI-2-day/FWI-2-day_ERA5_std_reanalysis_2022-06-01-2022-08-31_-9.0.50.0.2.0.59.0_day_initialise-from-copernicus:True-and-use-numpy=False.nc')


def CountryConstrain(cube):
    shpfilename = str(shpreader.natural_earth(resolution='110m', category='cultural', name='admin_0_countries'))
    natural_earth_file = shape.load_shp(str(shpfilename))
    Country = shape.load_shp(shpfilename, Name='United Kingdom')
    Country_shape = Country.unary_union()
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

def TimePercentile(cube, percentile):
    cube = cube.collapsed(['time'], iris.analysis.PERCENTILE, percent=percentile)
    return cube 

def TimeMean(cube):
    cube = cube.collapsed('time', iris.analysis.MEAN)
    return cube 


def prepareMap(cube, percentile):
    cube = cube.extract(daterange)
    cube = CountryConstrain(cube)
    cube = cube.collapsed('time', iris.analysis.PERCENTILE, percent=percentile)
    cube = cube.intersection(longitude=(-180, 180))
    return cube





ERA5 = prepareMap(ERA5_2022, percentile=percentile)


folder = '/scratch/cburton/scratch/FWI/2022/'
index_filestem1 = 'hist/'
index_filestem2 = 'histnat/'
index_name = 'Canadian_FWI'

#members = np.arange(1,106)
members = np.arange(1,3)
z=0
for member in members:
    print ('hist',member)
    #for n in np.arange(1,6):
    for n in np.arange(1,2):
        try:
            if member < 10:
                hist = iris.load_cube(folder+index_filestem1+'FWI_00'+str(member)+'_'+str(n)+'-2022.nc', index_name)
            elif member > 9 and member < 100:
                hist = iris.load_cube(folder+index_filestem1+'FWI_0'+str(member)+'_'+str(n)+'-2022.nc', index_name)
            else:
                hist = iris.load_cube(folder+index_filestem1+'FWI_'+str(member)+'_'+str(n)+'-2022.nc', index_name)
            hist=hist.extract(iris.Constraint(latitude=lambda cell: (40.0) < cell < (60.0), longitude=lambda cell: (0) < cell < (360))) #UK
            hist = CountryConstrain(hist)
            hist = hist.extract(daterange)
            hist = TimePercentile(hist, percentile)
            if z == 0:
                allhist = hist.copy()
                z=z+1
            elif z != 0:
                allhist = allhist+hist
                z=z+1
        except IOError:
             pass 
print(z)  
print(allhist)
allhist = allhist/z


z=0
for member in members:
    print ('histnat',member)
    #for n in np.arange(1,6):
    for n in np.arange(1,2):
        try:
            if member < 10:
                histnat = iris.load_cube(folder+index_filestem2+'FWI_00'+str(member)+'_'+str(n)+'-2022.nc', index_name)
            elif member > 9 and member < 100:
                histnat = iris.load_cube(folder+index_filestem2+'FWI_0'+str(member)+'_'+str(n)+'-2022.nc', index_name)
            else:
                histnat = iris.load_cube(folder+index_filestem2+'FWI_'+str(member)+'_'+str(n)+'-2022.nc', index_name)
            histnat=histnat.extract(iris.Constraint(latitude=lambda cell: (40.0) < cell < (60.0), longitude=lambda cell: (0) < cell < (360))) #UK
            histnat = CountryConstrain(histnat)
            histnat = histnat.extract(daterange)
            histnat = TimePercentile(histnat, percentile)
            if z == 0:
                allhistnat = histnat.copy()
                z=z+1
            elif z != 0:
                allhistnat = allhistnat+histnat
                z=z+1
        except IOError:
             pass 
print(z)  
print(allhistnat)
allhistnat = allhistnat/z



admin_boundaries = cfeature.NaturalEarthFeature(
    category='cultural', 
    name='admin_0_map_units', 
    scale='50m', 
    facecolor='none')


fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 8), subplot_kw={'projection': ccrs.PlateCarree()})

ax1 = axes[0]
iplt.pcolormesh(ERA5, vmin = 0, vmax=20.0,axes=ax1)
ax1.set_title('ERA5')
ax1.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=0.5)
ax1.add_feature(admin_boundaries, edgecolor='black', linewidth=0.5)
ax1.set_extent([-9, 2, 50, 59], crs=ccrs.PlateCarree())
ax1.coastlines('50m', linewidth=0.8)

ax2 = axes[1]
mesh=iplt.pcolormesh(allhist, vmin = 0, vmax=20.0,axes=ax2)
ax2.set_title('ALL')
ax2.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=0.5)
ax2.add_feature(admin_boundaries, edgecolor='black', linewidth=0.5)
ax2.set_extent([-9, 2, 50, 59], crs=ccrs.PlateCarree())
ax2.coastlines('50m', linewidth=0.8)

ax3 = axes[2]
iplt.pcolormesh(allhistnat, vmin = 0, vmax=20.0, axes=ax3)
ax3.set_title('NAT')
ax3.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=0.5)
ax3.add_feature(admin_boundaries, edgecolor='black', linewidth=0.5)
ax3.set_extent([-9, 2, 50, 59], crs=ccrs.PlateCarree())
ax3.coastlines('50m', linewidth=0.8)

                         #Position: Left, Bottom, width, height
colorbar_axes = plt.gcf().add_axes([0.3, 0.1, 0.45, 0.04])
colorbar = plt.colorbar(mesh, colorbar_axes, orientation='horizontal', label='FWI')

plt.show()









