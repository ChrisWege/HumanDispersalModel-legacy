# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a Testfile for generating an equidistant grid pattern
"""

import cartopy.crs as ccrs
import numpy as np
import pyproj
from netCDF4 import Dataset
import time

if __name__ == "__main__":
    p_ll = pyproj.Proj('epsg:4326')
    p_mt = pyproj.Proj('EPSG:3857')

    pc = ccrs.PlateCarree()

    stepsize = 50000. #50 km grid scale
    gridpointsx = 100
    gridpointsy = 50

    llcrnrlon = np.array([-15.])
    llcrnrlat = np.array([25.])
    
    llcrnr_mt = pyproj.transform(p_ll, p_mt, llcrnrlon[0], llcrnrlat[0])

    x = np.arange(gridpointsx)*stepsize+llcrnr_mt[1] # [0] and [1] have been switched around, epsg:3857 does something funky here...
    y = np.arange(gridpointsy)*stepsize+llcrnr_mt[0]
    x_mesh,y_mesh = np.meshgrid(x,y)
    
    lats,lons = pyproj.transform(p_mt,p_ll,x_mesh,y_mesh) # lats, lons are switched around here as well...
    
    potential = np.zeros((len(y), len(x)))
    
    path = "/data/sfb806/human_mobility_model/dispersal_model/data/"
    name1 = "Acc_EHEP_test_flat_closed_block_50km_equidis_gog"
    name2 = ".nc"
    
    potential[:,:] = 1.
    
    watermask = np.ones((len(y),len(x)))
    watermask[0:3,:] = 0.
    watermask[-3:,:] = 0.
    watermask[:,0:3] = 0.
    watermask[:,-3:] = 0.
    
    timeunit = "in years"
    
    output = Dataset(path + name1 + name2,mode='w')
    output.description = "EHEP with carrying capacity estimation, test case environment."
    output.history = "Created: "+time.ctime(time.time())
    net_dim_x = output.createDimension("x",len(x))
    net_dim_y = output.createDimension("y",len(y))
    net_lat = output.createVariable("lat","f8",("y","x"))
    net_lon = output.createVariable("lon","f8",("y","x"))
    net_Acc_EHEP = output.createVariable("Acc_EHEP","f8",("y","x"))
    net_watermask = output.createVariable("watermask","f8",("y","x"))
    output.timeslice = str(1000)
    output.timeunit = timeunit
    net_watermask[:,:] = watermask
    net_lat[:,:] = lats
    net_lon[:,:] = lons
    net_Acc_EHEP[:,:] = potential[:,:]
        
        
    # Adding additional documentation    
    net_lat.description = "Latitude"
    net_lat.units = "degrees_north"
    net_lon.description = "Longitude"
    net_lon.units = "degrees_east"
    net_Acc_EHEP.description = "Carrying capacity estimation"
    net_Acc_EHEP.units = "None"
    net_watermask.description = "(Ocean-) Watermask (1=land,0=water)"
    net_watermask.units = "None"
        
    output.close()
