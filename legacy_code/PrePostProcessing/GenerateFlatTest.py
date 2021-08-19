# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 10:51:23 2019

@author: cwegener 

Script for generating several potential fields to be used as testcases in the 
dispersal model. Required to see if the model is working correctly or not.
"""

from netCDF4 import Dataset
import numpy as np
from scipy import interpolate
import time

if __name__ == "__main__":
    path = "/data/sfb806/human_mobility_model/dispersal_model/data/"
    name1 = "Acc_EHEP_test_linear_closed_block_0.52deg"
    name2 = ".nc"
    
    points_x = 81
    points_y = 61
    int_x = 21
    int_y = 21
    lat_start = 30
    lat_stop = 60
    lon_start = -15
    lon_stop = 25
    lat = np.zeros((points_y,points_x))
    lon = np.zeros((points_y,points_x))
    lat_flat = np.linspace(lat_start,lat_stop,points_y)
    lon_flat = np.linspace(lon_start,lon_stop,points_x)
    lon, lat = np.meshgrid(lon_flat,lat_flat)
    carry_capa = np.linspace(0.,1.,points_x)
    potential = np.zeros((points_y,points_x))
    
    rlat = np.zeros((int_y,int_x))
    rlon = np.zeros((int_y,int_x))
    rlat_flat = np.linspace(lat_start,lat_stop,int_y)
    rlon_flat = np.linspace(lon_start,lon_stop,int_x)
    rlon,rlat = np.meshgrid(rlon_flat,rlat_flat)
    
    for i,val in enumerate(lat_flat):
        potential[i,:] = carry_capa
        
    watermask = np.ones((points_y,points_x))
    # Closed boundary
    watermask[0:3,:] = 0.
    watermask[-3:,:] = 0.
    watermask[:,0:3] = 0.
    watermask[:,-3:] = 0.
    # capsulated Island near the center
    #watermask[25:36,50:61] = 0.
    #watermask[26:35,51:60] = 1.
        
    timeunit = "in years"
    
    output = Dataset(path + name1 + name2,mode='w')
    output.description = "EHEP with carrying capacity estimation, test case environment."
    output.history = "Created: "+time.ctime(time.time())
    net_dim_x = output.createDimension("x",points_x)
    net_dim_y = output.createDimension("y",points_y)
    net_lat = output.createVariable("lat","f8",("y","x"))
    net_lon = output.createVariable("lon","f8",("y","x"))
    net_Acc_EHEP = output.createVariable("Acc_EHEP","f8",("y","x"))
    net_watermask = output.createVariable("watermask","f8",("y","x"))
    output.timeslice = str(1000)
    output.timeunit = timeunit
    net_watermask[:,:] = watermask
    net_lat[:,:] = lat
    net_lon[:,:] = lon
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
    
    
