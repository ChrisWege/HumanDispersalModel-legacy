# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 14:42:16 2018

@author: cwegener

Module for the creation of the Accessible Environmental Human Existence
Potential (Accesible EHEP, Phi_S).
Loads in the EHEP out of .nc data (any size) and modifies it with the
topography, natural borders and technology terms.
Calculates a maximum capacity of humans for the grid and applies it on the
Accessible EHEP.
"""

from netCDF4 import Dataset
import numpy as np
from scipy import interpolate
import time
import sys
from geographiclib.geodesic import Geodesic
from itertools import product

### Functions ###
def apply_glacier(Ehep,lat,lon,x=None,y=None):
    """Applys the glacial information on the given potential.
    
    Takes in the information about glaciers and applies it to the EHEP so that
    occurence of them nullifies the potential. Absence does not alter the
    potential."""
    # Loading in the data and preparing interpolation function
    icedata = Dataset(path_ice,mode='r')
    lat_ice = np.array(icedata.variables[var_name_ice_lat][:])
    lon_ice = np.array(icedata.variables[var_name_ice_lon][:])
    fraction = np.array(icedata.variables[var_name_ice_cover][:][:])
    Ice_interpol = interpolate.RegularGridInterpolator((lat_ice,lon_ice),
                   fraction,bounds_error=False,fill_value=0.)
    # Calculating a mask based on the data
    icedata.close()
   
    if gridtype == 'xy':
        combined = list(zip(lat,lon))
    else:
        combined = list(product(lat,lon))
    mask = 1 - Ice_interpol(combined)

    if gridtype == 'xy':
        mask = np.reshape(mask,(y,x))
        alt_Ehep = np.zeros((len(timeslice),y,x))
    else:
        mask = np.reshape(mask,(len(lat),len(lon)))
        alt_Ehep = np.zeros((len(timeslice),len(lat),len(lon)))
    
    # Applying the mask and returning the modified potential
    
    for i,val in enumerate(timeslice): 
        alt_Ehep[i,:,:] = Ehep[i,:,:] * mask
    print("Finished apply_glacier...")
    return alt_Ehep
    
def save_glacier(Ehep,lat,lon,x=None,y=None):
    """Saves the changes made by apply_glacier."""
    # Save file preperations and documentation
    appice = Dataset(path_output + file_out_glac,mode='w',format="NETCDF4")
    appice.description = "Intermediate step of the Accessible EHEP (1/5). "\
"After Glacial modifications."
    appice.history = "Created: " +time.ctime(time.time())
    
    # Creating variables and saving data
    if gridtype == 'xy':
        net_dim_x = appice.createDimension("x",x)
        net_dim_y = appice.createDimension("y",y)
        net_dim_time = appice.createDimension("time",None)
        net_lat = appice.createVariable("lat","f8",("y","x"))
        net_lon = appice.createVariable("lon","f8",("y","x"))
        net_Acc_EHEP = appice.createVariable("Acc_EHEP","f8",("time","y","x"))
        net_lat[:,:] = lat
        net_lon[:,:] = lon
    else:
        net_dim_lon = appice.createDimension("lon",len(lon))
        net_dim_lat = appice.createDimension("lat",len(lat))
        net_dim_time = appice.createDimension("time",None)
        net_lat = appice.createVariable("lat","f8",("lat"))
        net_lon = appice.createVariable("lon","f8",("lon"))
        net_Acc_EHEP = appice.createVariable("Acc_EHEP","f8",("time","lat","lon"))
        net_lat[:] = lat
        net_lon[:] = lon
    
    net_time = appice.createVariable("time","f8",("time"))
    net_Acc_EHEP[:,:,:] = Ehep
    net_time[:] = timeslice
    
    # Additional documentation
    net_lat.description = "Latitude"
    net_lat.units = "degrees_north"
    net_lon.description = "Longitude"
    net_lon.units = "degrees_east"
    net_Acc_EHEP.description = "Accesible EHEP (Step 1/5)"
    net_Acc_EHEP.units = "None"
    net_time.description = "time"
    net_time.unit = timeunit
    appice.close()
    print("Saved glacier changes...")
    
def apply_watermask(Ehep,lat,lon,x=None,y=None):
    """Applies the land mask information on the given potential
    
    Takes in the information about water bodies and masks them in the given
    potential. """
    # Loading the data and preparing interpolation function
    maskdata = Dataset(path_watermask,mode='r')
    lat_mask = np.array(maskdata.variables[var_name_watermask_lat][:])
    lon_mask = np.array(maskdata.variables[var_name_watermask_lon][:])
    fraction = np.array(maskdata.variables[var_name_watermask_mask][:][:])
    mask_interpol = interpolate.RegularGridInterpolator((lat_mask,lon_mask),
                    fraction,method='nearest',bounds_error=False,fill_value=0.)
    maskdata.close()
    # Creating the mask array
    if gridtype == 'xy':
        combined = list(zip(lat,lon))
    else:
        combined = list(product(lat,lon))
    mask = mask_interpol(combined)
    if gridtype == 'xy':
        mask = np.reshape(mask,(y,x))
        alt_Ehep = np.zeros((len(timeslice),y,x))
    else:
        mask = np.reshape(mask,(len(lat),len(lon)))
        alt_Ehep = np.zeros((len(timeslice),len(lat),len(lon)))
   
    # Applying the mask and returning
    for i,val in enumerate(timeslice):
        alt_Ehep[i,:,:] = Ehep[i,:,:] * mask
    print("Finished apply_watermask...")
    return alt_Ehep,mask

def save_watermask(Ehep,lat,lon,x=None,y=None):
    """Saves the changes made by apply_watermask."""
    # Save file preperations and documentation
    appmask = Dataset(path_output + file_out_watermask,mode='w',
                      format="NETCDF4")
    appmask.description = "Intermediate step of the Accessible EHEP (2/5). "\
"After water/land mask modifications."
    appmask.history = "Created: " +time.ctime(time.time())
    
    # Creating variabels and saving data
    if gridtype == 'xy':
        net_dim_x = appmask.createDimension("x",x)
        net_dim_y = appmask.createDimension("y",y)
        net_dim_time = appmask.createDimension("time",None)
        net_lat = appmask.createVariable("lat","f8",("y","x"))
        net_lon = appmask.createVariable("lon","f8",("y","x"))
        net_Acc_EHEP = appmask.createVariable("Acc_EHEP","f8",("time","y","x"))
        net_watermask = appmask.createVariable("watermask","f8",("y","x"))
        net_lat[:,:] = lat
        net_lon[:,:] = lon
    else:
        net_dim_lon = appmask.createDimension("lon",len(lon))
        net_dim_lat = appmask.createDimension("lat",len(lat))
        net_dim_time = appmask.createDimension("time",None)
        net_lat = appmask.createVariable("lat","f8",("lat"))
        net_lon = appmask.createVariable("lon","f8",("lon"))
        net_Acc_EHEP = appmask.createVariable("Acc_EHEP","f8",("time","lat","lon"))
        net_watermask = appmask.createVariable("watermask","f8",("lat","lon"))
        net_lat[:] = lat
        net_lon[:] = lon 
        
    net_time = appmask.createVariable("time","f8",("time"))
    net_watermask[:,:] = mask
    net_Acc_EHEP[:,:] = Ehep
    net_time[:] = timeslice
    
    
    # Additional documentation
    net_lat.description = "Latitude"
    net_lat.units = "degrees_north"
    net_lon.description = "Longitude"
    net_lon.units = "degrees_east"
    net_Acc_EHEP.description = "Accesible EHEP (Step 2/5)"
    net_Acc_EHEP.units = "None"
    net_time.description = "time"
    net_time.unit = timeunit
    appmask.close()
    print("Saved watermask changes...")
    
def apply_forest(Ehep,lat,lon,x=None,y=None):
    """Applies the vegetation forest data to alter the potential.
    
    Loads the information about forests from vegetation dataand uses a function
    to create a grid that alters the incomming potential. High forest densities 
    lower the potential, as humans can not move or hunt through thick forests.
    """
    def forest_func(inx):
        """Mathematical function for the forest fraction."""
        # Currently linear damping
        outy = np.array(inx)
        outy[inx < veg_params[0]] = veg_params[2]
        outy[inx > veg_params[1]] = veg_params[3]
        mask = np.logical_and(inx >= veg_params[0], inx <= veg_params[1])
        outy[mask] = veg_params[2] + (inx[mask]-veg_params[0])*(veg_params[3]
        - veg_params[2])/(veg_params[1] - veg_params[0])
        return outy
    # Loading data
    forestdata = Dataset(path_veg,mode='r')
    lat_forest = np.array(forestdata.variables[var_name_veg_lat][:])
    lon_forest = np.array(forestdata.variables[var_name_veg_lon][:])
    fraction = np.array(forestdata.variables[var_name_veg_data][:][:][:])
    forestdata.close()
    # Sum of different forest types to get one forest fraction
    forestsum = (fraction[0,:,:] + fraction[1,:,:] + fraction[2,:,:] 
                + fraction[3,:,:] + fraction[4,:,:])
    del fraction
    # Creating interpolation function
    forest_interpol = interpolate.RegularGridInterpolator((lat_forest,
                      lon_forest),forestsum,bounds_error=False,fill_value=0.)
    # Calculation of the forest mask
    if gridtype == 'xy':
        combined = list(zip(lat,lon))
    else:
        combined = list(product(lat,lon))
    forest = forest_interpol(combined)
    if gridtype == 'xy':
        forest = np.reshape(forest,(y,x))
        alt_Ehep = np.zeros((len(timeslice),y,x))
        forest = forest_func(forest)
    else:
        forest = np.zeros((len(lat),len(lon)))
        alt_Ehep = np.zeros((len(timeslice),len(lat),len(lon)))
        forest = forest_func(forest)
        
    # Applying the forest alterations on the potential
    for i,val in enumerate(timeslice):
        alt_Ehep[i,:,:] = Ehep[i,:,:] * forest
    print("Finished apply_forest")
    return alt_Ehep
    
    
def save_forest(Ehep,lat,lon,x=None,y=None):
    """Saves the changes made by apply_forest."""
    # Creating save file
    appforest = Dataset(path_output + file_out_forest,mode='w',format="NETCDF4")
    appforest.description = "Intermediate step of the Accessible EHEP (3/5). "\
"After forest modifications."
    appforest.history = "Created: " +time.ctime(time.time())
    
    # Creating variables and saving data
    if gridtype == 'xy':
        net_dim_x = appforest.createDimension("x",x)
        net_dim_y = appforest.createDimension("y",y)
        net_dim_time = appforest.createDimension("time",None)
        net_lat = appforest.createVariable("lat","f8",("y","x"))
        net_lon = appforest.createVariable("lon","f8",("y","x"))
        net_Acc_EHEP = appforest.createVariable("Acc_EHEP","f8",("time","y","x"))
        net_watermask = appforest.createVariable("watermask","f8",("y","x"))
        net_lat[:,:] = lat
        net_lon[:,:] = lon
    else:
        net_dim_lon = appforest.createDimension("lon",len(lon))
        net_dim_lat = appforest.createDimension("lat",len(lat))
        net_dim_time = appforest.createDimension("time",None)
        net_lat = appforest.createVariable("lat","f8",("lat"))
        net_lat = appforest.createVariable("lon","f8",("lon"))
        net_Acc_EHEP = appforest.createVariable("Acc_EHEP","f8",("time","lat","lon"))
        net_watermask = appforest.createVariable("watermask","f8",("lat","lon"))
        net_lat[:] = lat
        net_lon[:] = lon
        
    net_time = appforest.createVariable("time","f8",("time")) 
    net_watermask[:,:] = mask
    net_Acc_EHEP[:,:,:] = Ehep
    net_time[:] = timeslice
    # Additional Documentation
    net_lat.description = "Latitude"
    net_lat.units = "degrees_north"
    net_lon.description = "Longitude"
    net_lon.units = "degrees_east"
    net_Acc_EHEP.description = "Accesible EHEP (Step 3/5)"
    net_Acc_EHEP.units = "None"
    net_time.description = "time"
    net_time.unit = timeunit
    appforest.close()
    print("Saved forest changes...")
    
def apply_technology():
    """Creates a function for the technology term.
    
    Uses the predefined information and creates a multiplier which accounts for
    the technological advancement of humans. This function is mostly a 
    placeholder, as the technological advancement may influence many aspects 
    of the model."""
    
    return
    
def save_technology():
    """Saves the changes made by apply_technology."""
    return

def apply_orography(Ehep,lat,lon,x=None,y=None):
    """Creates a function for the orography term.
    
    Loads the orography information (mean elevation, standard deviation) from
    file and alters the input potential with the two functions height_func and
    std_func based on that data."""
    def height_func(inx):
        """Mathematical function for the height effect on potential."""
        # Currently exponential functions
        outy = np.array(inx)
        outy[inx < elevation_params[0]] = elevation_params[2]
        outy[inx > elevation_params[1]] = elevation_params[3]
        mask = np.logical_and(inx >= elevation_params[0],
                              inx <= elevation_params[1])
        outy[mask] = elevation_params[2] + (inx[mask]-elevation_params[0])*(
        elevation_params[3] - elevation_params[2])/(elevation_params[1]
        - elevation_params[0])
        return outy
        
        
    def std_func(inx):
        """Mathematical function for the roughness effect on potential."""
        outy = np.array(inx)
        # Currently log function
        outy[inx < std_params[0]] = std_params[2]
        outy[inx > std_params[1]] = std_params[3]
        mask = np.logical_and(inx >= std_params[0], inx <= std_params[1])
        outy[mask] = std_params[2] + (inx[mask]-std_params[0])*(std_params[3]
        - std_params[2])/(std_params[1] - std_params[0])
        return outy
    
    # Loading in data and preperations of interpolation functions
    topodata = Dataset(path_topo,mode='r')
    lat_mask = np.array(topodata.variables[var_name_topo_lat][:])
    lon_mask = np.array(topodata.variables[var_name_topo_lon][:])
    height = np.array(topodata.variables[var_name_topo_height][:][:])
    std = np.array(topodata.variables[var_name_topo_std][:][:])
    height_interpol = interpolate.RegularGridInterpolator((lat_mask,
                      lon_mask),height,bounds_error=False,fill_value=0.)
    std_interpol = interpolate.RegularGridInterpolator((lat_mask,
                      lon_mask),std,bounds_error=False,fill_value=0.)
    topodata.close()

    # Reshaping of variables
    if gridtype == 'xy':
        combined = list(zip(lat,lon))
    else:
        combined = list(product(lat,lon))
    height_int = height_interpol(combined)
    std_int = std_interpol(combined)
    if gridtype == 'xy':
        height_int = np.reshape(height_int,(y,x))
        std_int = np.reshape(std_int,(y,x))
        alt_Ehep = np.zeros((len(timeslice),y,x))
        height_int = height_func(height_int)
        std_int = std_func(std_int)
    else:
        height_int = np.reshape(height_int,(len(lat),len(lon)))
        std_int = np.reshape(std_int,(len(lat),len(lon)))
        alt_Ehep = np.zeros((len(timeslice),len(lat),len(lon)))
        height_int = height_func(height_int)
        std_int = std_func(std_int)
        
    
    # Calculations of the damping functions
    height_int = height_func(height_int)
    std_int = std_func(std_int)
    # Applying and returning the modifications
    for i,val in enumerate(timeslice):
        alt_Ehep[i,:,:] = Ehep[i,:,:] * height_int * std_int
    print("Finished apply_orography...")
    return alt_Ehep
    
def save_orography(Ehep,lat,lon,x=None,y=None):
    """Saves the changes made by apply_orography."""
    # Creating save file and documentation
    apporo = Dataset(path_output + file_out_oro,mode='w',format="NETCDF4")
    apporo.description = "Intermediate step of the Accessible EHEP (5/5). "\
"After orography modifications."
    apporo.history = "Created: " +time.ctime(time.time())
    
    # Creating variabels and saving data
    if gridtype == 'xy':
        net_dim_x = apporo.createDimension("x",x)
        net_dim_y = apporo.createDimension("y",y)
        net_dim_time = apporo.createDimension("time",None)
        net_lat = apporo.createVariable("lat","f8",("y","x"))
        net_lon = apporo.createVariable("lon","f8",("y","x"))
        net_Acc_EHEP = apporo.createVariable("Acc_EHEP","f8",("time","y","x"))
        net_watermask = apporo.createVariable("watermask","f8",("y","x"))
        net_lat[:,:] = lat
        net_lon[:,:] = lon
    else:
        net_dim_lon = apporo.createDimension("lon",len(lon))
        net_dim_lat = apporo.createDimension("lat",len(lat))
        net_dim_time = apporo.createDimension("time",None)
        net_lat = apporo.createVariable("lat","f8",("lat"))
        net_lon = apporo.createVariable("lon","f8",("lon"))
        net_Acc_EHEP = apporo.createVariable("Acc_EHEP","f8",("time","lat","lon"))
        net_watermask = apporo.createVariable("watermask","f8",("lat","lon"))
        net_lat[:] = lat
        net_lon[:] = lon
        
    net_time = apporo.createVariable("time","f8",("time"))
    net_watermask[:,:] = mask
    net_Acc_EHEP[:,:,:] = Ehep
    # Additional Documentation
    net_lat.description = "Latitude"
    net_lat.units = "degrees_north"
    net_lon.description = "Longitude"
    net_lon.units = "degrees_east"
    net_Acc_EHEP.description = "Accesible EHEP (Step 5/5)"
    net_Acc_EHEP.units = "None"
    net_time.description = "time"
    net_time.unit = timeunit
    apporo.close()
    print("Saved orography changes...")
    
def estimate_capacity(Ehep,lat,lon,x=None,y=None):
    """Creates an estimation of the maximum capacity for humans.

    Takes in a predefined maximum capacity distribution or creates one on the 
    basis of vegetation data."""
    def calc_capacity(veg,setnum):
        """Subfunction for the carrying capacity."""
        if setnum == 0:
            capacity = (np.array(veg[0,:,:] * 8. + veg[1,:,:] * 16. 
                        + veg[2,:,:] * 8. + veg[3,:,:] * 13. 
                        + veg[4,:,:] * 13. + veg[5,:,:] * 12. 
                        + veg[6,:,:] * 12. + veg[7,:,:] * 18. 
                        + veg[8,:,:] * 0. + veg[9,:,:] * 14. 
                        + veg[10,:,:] * 0.))

        print("Capacity Calculation Subroutine: done")
        return capacity

    # Loading vegetation dataset
    Data = Dataset(path_veg,mode='r')
    vegetation = np.array(Data.variables["landuse"][:,:,:])
    lat_veg = Data.variables["lat"][:]
    lon_veg = Data.variables["lon"][:]
    Data.close()
    if gridtype == "xy":
        capa_int = np.empty([len(vegetation_range),y,x])
        alt_Ehep = np.zeros((len(timeslice),y,x))
    elif gridtype == "latlon":
        capa_int = np.empty([len(vegetation_range),len(lat),len(lon)])
        alt_Ehep = np.zeros((len(timeslice),len(lat),len(lon)))
    
    # Creating an interpolation function for every biome (climate) type and 
    # saving the interpolated biome data in one large array
    for k in vegetation_range:
        veg_interpol = interpolate.RegularGridInterpolator((lat_veg,
                                                        lon_veg),
                                                        vegetation[k,:,:])
        
        combined = list(zip(lat,lon))
        capa_int_step = veg_interpol(combined)
        capa_int_step = np.reshape(capa_int_step,(y,x))
    
        capa_int[k,:,:] = capa_int_step
        print("Step ",k," finished...")
    # Calculation of the carrying capacity and modification of the dimensionless
    # potential into a junction representing maximum human densities
    Capa_calc = calc_capacity(capa_int,setnumber_capacity)

    for i,val in enumerate(timeslice):
        alt_Ehep[i,:,:] = Ehep[i,:,:] * Capa_calc
    print("Finished estimate_capacity")
    return alt_Ehep
    
def save_capa(Ehep,lat,lon,x=None,y=None):
    """Saves the final results with the maximum capacity as value."""
    # Creating the save file and documentation
    appcap = Dataset(path_output + file_out_capa,mode='w',format="NETCDF4")
    appcap.description = "Carrying capacity estimation based on the "\
"Accesible EHEP and the vegetation data."
    appcap.history = "Created: " +time.ctime(time.time())
    
    # creating variables and saving data
    if gridtype == 'xy':
        net_dim_x = appcap.createDimension("x",x)
        net_dim_y = appcap.createDimension("y",y)
        net_dim_time = appice.createDimension("time",None)
        net_lat = appcap.createVariable("lat","f8",("y","x"))
        net_lon = appcap.createVariable("lon","f8",("y","x"))
        net_Acc_EHEP = appcap.createVariable("Acc_EHEP","f8",("time","y","x"))
        net_watermask = appcap.createVariable("watermask","f8",("y","x"))
        net_lat[:,:] = lat
        net_lon[:,:] = lon
    else:
        net_dim_lon = appcap.createDimension("lon",len(lon))
        net_dim_lat = appcap.createDimension("lat",len(lat))
        net_dim_time = appcap.createDimension("time",None)
        net_lat = appcap.createVariable("lat","f8",("lat"))
        net_lon = appcap.createVariable("lon","f8",("lon"))
        net_Acc_EHEP = appcap.createVariable("Acc_EHEP","f8",("time","lat","lon"))
        net_watermask = appcap.createVariable("watermask","f8",("y","x"))
        net_lat[:] = lat
        net_lon[:] = lon
    
    net_time = appcap.createVariable("time","f8",("time"))    
    net_watermask[:,:] = mask
    net_Acc_EHEP[:,:,:] = Ehep
    net_time[:] = timeslice
    
    # Adding additional documentation    
    net_lat.description = "Latitude"
    net_lat.units = "degrees_north"
    net_lon.description = "Longitude"
    net_lon.units = "degrees_east"
    net_Acc_EHEP.description = "Carrying capacity estimation"
    net_Acc_EHEP.units = "None"
    net_time.description = "time"
    net_time.unit = timeunit
    appcap.close()
    print("Saved carrying capacity estimation")
    
def preprocess(EHEP,lat,lon,x=None,y=None):
    """Preprocesses the results to be useable with the dispersal model."""
    geod = Geodesic.WGS84
    if not switch_capa:
        EHEP = EHEP * 10.
    if gridtype == 'xy':
        density = np.zeros((y,x))
        if starting_dist == 'uniform':
            density[:,:] = 2. * mask[:,:]
        elif starting_dist == 'gaus':
            st_x0 = 32.
            st_y0 = 13.
            st_sx0 = 1.
            st_sy0 = 1.
            st_ampli = 20.
            gauß = np.zeros((y,x))
            xarray = np.linspace(0,x,x)
            yarray = np.linspace(0,y,y)
            for i,valx in enumerate(xarray):
                for j,valy in enumerate(yarray):
                    gauß[j,i] = st_ampli* np.exp(-((valx-st_x0)**2/(2.*st_sx0**2)+(valy-st_y0)**2/(2.*st_sy0**2)))
            density[:,:] = gauß * mask
        elif starting_dist == 'firstslice':
            density[:,:] = EHEP[0,:,:]
    else:
        density = np.zeros((len(lat),len(lon)))
        x = len(lon)
        y = len(lat)
        if starting_dist == 'uniform':
            density[:,:] = 2. * mask[:,:]
        elif starting_dist == 'gaus':
            st_x0 = 208.
            st_y0 = 24.
            st_sx0 = 4.
            st_sy0 = 4.
            st_ampli = 20.
            gauß = np.zeros((y,x))
            xarray = np.linspace(0,x,x)
            yarray = np.linspace(0,y,y)
            for i,valx in enumerate(xarray):
                for j,valy in enumerate(yarray):
                    gauß[j,i] = st_ampli* np.exp(-((valx-st_x0)**2/(2.*st_sx0**2)+(valy-st_y0)**2/(2.*st_sy0**2)))
            density[:,:] = gauß * mask
        elif starting_dist == 'firstslice':
            density[:,:] = EHEP[0,:,:]
        
    # create the output file and fill it with structure
    prepro = Dataset(path_output + file_out_mobil,mode='w',format="NETCDF4")
    prepro.description = "Accessible EHEP with additional values for "\
"the mobility model."
    prepro.history = "Created: " +time.ctime(time.time())
    
    net_dim_x = prepro.createDimension("x",x)
    net_dim_y = prepro.createDimension("y",y)
    net_dim_x1 = prepro.createDimension("x_1",x-1)
    net_dim_y1 = prepro.createDimension("y_1",y-1)
    net_dim_x2 = prepro.createDimension("x_2",x+1)
    net_dim_y2 = prepro.createDimension("y_2",y+1)
    net_dim_single = prepro.createDimension("single",1)

    net_lat = prepro.createVariable("lat","f8",("y","x"))
    net_lon = prepro.createVariable("lon","f8",("y","x"))
    net_dens = prepro.createVariable("dens","f8",("y","x"))
    net_hnumber = prepro.createVariable("hnumb","f8",("y","x"))
    net_dx = prepro.createVariable("dx","f8",("single"))
    net_dy = prepro.createVariable("dy","f8",("single"))
    net_dx_grid = prepro.createVariable("dx_grid","f8",("y","x_1"))
    net_dy_grid = prepro.createVariable("dy_grid","f8",("y_1","x"))
    net_posx_lat = prepro.createVariable("posx_lat","f8",("y","x_1"))
    net_posx_lon = prepro.createVariable("posx_lon","f8",("y","x_1"))
    net_posy_lat = prepro.createVariable("posy_lat","f8",("y_1","x"))
    net_posy_lon = prepro.createVariable("posy_lon","f8",("y_1","x"))
    net_area = prepro.createVariable("area","f8",("y","x"))
    net_corner_lat = prepro.createVariable("corner_lat","f8",("y_2","x_2"))
    net_corner_lon = prepro.createVariable("corner_lon","f8",("y_2","x_2"))
    
    
    # Adding additional documentation    
    net_lat.description = "Latitude"
    net_lat.units = "Degrees north"
    net_lon.description = "Longitude"
    net_lon.units = "Degrees east"
    net_dx.description = "Mean distance between two gridpoints, x-direction"
    net_dx.unit = "km"
    net_dy.description = "Mean distance between two gridpoints, y-direction"
    net_dy.unit = "km"
    net_dx_grid.description = "Grid specific distance between two points, x-direction"
    net_dx_grid.uni = "km"
    net_dy_grid.description = "Grid specific distance between two points, y-direction"
    net_dy_grid.uni = "km"
    net_posx_lat.description = "Latitude of border between two points, x-direction"
    net_posx_lon.description = "Longitude of border between two points, x-direction"
    net_posy_lat.description = "Latitude of border between two points, y-direction"
    net_posy_lon.description = "Longitude of border between two points, y-direction"
    net_posx_lat.unit = "Degrees north"
    net_posx_lon.unit = "Degrees east"
    net_posy_lat.unit = "Degrees north"
    net_posy_lon.unit = "Degrees east"
    net_dens.description = "Starting distribution of human density"
    net_dens.unit = "Humans/100km^2"
    net_hnumber.description = "Starting distribution of number of humans"
    net_hnumber.unit = "Humans"
    net_area.description = "Area of the grid cell"
    net_area.unit = "km^2"
    
    # save the first values and initialize others
    if gridtype == 'xy':
        net_lat[:,:] = lat
        net_lon[:,:] = lon
    else:
        latm,lonm = np.meshgrid(lat,lon)
        net_lat[:,:] = latm
        net_lon[:,:] = lonm
    net_dens[:,:] = density
    latxygrid = np.zeros((x,y))
    lonxygrid = np.zeros((x,y))
    #meterslat = np.zeros((x,y))
    #meterslon = np.zeros((x,y))
    
    # define lat/lon coordinates on curvilinear grid (and transpose them)
    if gridtype == 'xy':
        for i in list(range(x)):
            for j in list(range(y)):
                latxygrid[i,j] = lat[j,i]
                lonxygrid[i,j] = lon[j,i]
    else:
        latxygrid = latm
        lonxygrid = lonm
        
    # transform the distances from lowest left grid point into distances between
    # each gridpoint
    metersx = np.zeros((x-1,y))
    metersy = np.zeros((x,y-1))
    posx_lat = np.zeros((x-1,y))
    posx_lon = np.zeros((x-1,y))
    posy_lat = np.zeros((x,y-1))
    posy_lon = np.zeros((x,y-1))
    for i in list(range(x-1)):
        for j in list(range(y-1)):
            #glon = geod.Inverse(latxygrid[i,j],lonxygrid[i,j],latxygrid[i+1,j],lonxygrid[i+1,j])
            mlon = geod.InverseLine(latxygrid[i,j],lonxygrid[i,j],latxygrid[i+1,j],lonxygrid[i+1,j])
            metersx[i,j] = mlon.s13
            mlon_t = mlon.Position(0.5 * mlon.s13)
            posx_lat[i,j] = mlon_t['lat2']
            posx_lon[i,j] = mlon_t['lon2']
            #glat = geod.Inverse(latxygrid[i,j],lonxygrid[i,j],latxygrid[i,j+1],lonxygrid[i,j+1])
            mlat = geod.InverseLine(latxygrid[i,j],lonxygrid[i,j],latxygrid[i,j+1],lonxygrid[i,j+1]) 
            metersy[i,j] = mlat.s13
            mlat_t = mlat.Position(0.5 * mlat.s13)
            posy_lat[i,j] = mlat_t['lat2']
            posy_lon[i,j] = mlat_t['lon2']
            
    for i in list(range(x-1)):
        mlon = geod.InverseLine(latxygrid[i,-1],lonxygrid[i,-1],latxygrid[i+1,-1],lonxygrid[i+1,-1])
        metersx[i,-1] = mlon.s13
        mlon_t = mlon.Position(0.5 * mlon.s13)
        posx_lat[i,-1] = mlon_t['lat2']
        posx_lon[i,-1] = mlon_t['lon2']
    for j in list(range(y-1)):
        mlat = geod.InverseLine(latxygrid[-1,j],lonxygrid[-1,j],latxygrid[-1,j+1],lonxygrid[-1,j+1])
        metersy[-1,j] = mlat.s13
        mlat_t = mlat.Position(0.5 * mlat.s13)
        posy_lat[-1,j] = mlat_t['lat2']
        posy_lon[-1,j] = mlat_t['lon2']
        
    # take the mean of the distances for dx and dy, introduces small error but 
    # neglectiable for most cases. No error if given grid is equidistant.
    dx = metersx.mean()/1000.
    dy = metersy.mean()/1000.
    
    # save the rest of the initial values
    net_dx[:] = dx
    #net_dx[:] = 50.
    net_dy[:] = dy
    #net_dy[:] = 50.
    net_dx_grid[:,:] = metersx.T/1000.
    net_dy_grid[:,:] = metersy.T/1000.
    net_posx_lat[:,:] = posx_lat.T
    net_posx_lon[:,:] = posx_lon.T
    net_posy_lat[:,:] = posy_lat.T
    net_posy_lon[:,:] = posy_lon.T
    
    # calculate the area that every gridpoint represents
    ## if the grid is equidistant, the corner points are exact. Otherwise, an
    ## approximation is used by estimating the corner through lat/lon direction
    ## independently and averaging over both estimations.
    
    corner_points_lat = np.zeros((x+1,y+1))
    corner_points_lon = np.zeros((x+1,y+1))
    
    ### find the positions of the corners of each gridbox:
    # inner points
    for i in list(range(x-1)):
        for j in list(range(y-1)):
            mlat = geod.InverseLine(posy_lat[i,j],posy_lon[i,j],posy_lat[i+1,j],posy_lon[i+1,j])
            mlat_t = mlat.Position(0.5* mlat.s13)
            
            mlon = geod.InverseLine(posx_lat[i,j],posx_lon[i,j],posx_lat[i,j+1],posx_lon[i,j+1])
            mlon_t = mlon.Position(0.5* mlon.s13)
            #print(mlat_t,mlon_t)
            corner_points_lat[i+1,j+1] = np.mean(np.array([mlat_t['lat2'],mlon_t['lat2']]))
            corner_points_lon[i+1,j+1] = np.mean(np.array([mlat_t['lon2'],mlon_t['lon2']]))
            
    # edges
    for i in list(range(x-1)):
        edge_low = geod.InverseLine(posx_lat[i,0],posx_lon[i,0],posx_lat[i,1],posx_lon[i,1])
        edge_high = geod.InverseLine(posx_lat[i,-2],posx_lon[i,-2],posx_lat[i,-1],posx_lon[i,-1])
        edge_low_t = edge_low.Position(-0.5*edge_low.s13)
        edge_high_t = edge_high.Position(1.5*edge_high.s13)
        corner_points_lat[i+1,0] = edge_low_t['lat2']
        corner_points_lon[i+1,0] = edge_low_t['lon2']
        corner_points_lat[i+1,-1] = edge_high_t['lat2']
        corner_points_lon[i+1,-1] = edge_high_t['lon2']
    
    for j in list(range(y-1)):
        edge_left = geod.InverseLine(posy_lat[0,j],posy_lon[0,j],posy_lat[1,j],posy_lon[1,j])
        edge_right = geod.InverseLine(posy_lat[-2,j],posy_lon[-2,j],posy_lat[-1,j],posy_lon[-1,j])
        edge_left_t = edge_left.Position(-0.5*edge_left.s13)
        edge_right_t = edge_right.Position(1.5*edge_right.s13)
        corner_points_lat[0,j+1] = edge_left_t['lat2']
        corner_points_lon[0,j+1] = edge_left_t['lon2']
        corner_points_lat[-1,j+1] = edge_right_t['lat2']
        corner_points_lon[-1,j+1] = edge_right_t['lon2']
        
    # corners
    ## lower left
    corn_x = geod.InverseLine(corner_points_lat[1,0],corner_points_lon[1,0],
                                  corner_points_lat[2,0],corner_points_lon[2,0])
    corn_y = geod.InverseLine(corner_points_lat[0,1],corner_points_lon[0,1],
                                  corner_points_lat[0,2],corner_points_lon[0,2])
    corn_x_t = corn_x.Position(-1*corn_x.s13)
    corn_y_t = corn_y.Position(-1*corn_y.s13)
    corner_points_lat[0,0] = np.mean(np.array([corn_x_t['lat2'],corn_y_t['lat2']]))
    corner_points_lon[0,0] = np.mean(np.array([corn_x_t['lon2'],corn_y_t['lon2']]))
    ## lower right
    corn_x = geod.InverseLine(corner_points_lat[-3,0],corner_points_lon[-3,0],
                                     corner_points_lat[-2,0],corner_points_lon[-2,0])
    corn_y = geod.InverseLine(corner_points_lat[-1,1],corner_points_lon[-1,1],
                                     corner_points_lat[-1,2],corner_points_lon[-1,2])
    corn_x_t = corn_x.Position(2*corn_x.s13)
    corn_y_t = corn_y.Position(-1*corn_y.s13)
    corner_points_lat[-1,0] = np.mean(np.array([corn_x_t['lat2'],corn_y_t['lat2']]))
    corner_points_lon[-1,0] = np.mean(np.array([corn_x_t['lon2'],corn_y_t['lon2']]))
    ## upper left
    corn_x = geod.InverseLine(corner_points_lat[1,-1],corner_points_lon[1,-1],
                                  corner_points_lat[2,-1],corner_points_lon[2,-1])
    corn_y = geod.InverseLine(corner_points_lat[0,-3],corner_points_lon[0,-3],
                                  corner_points_lat[0,-2],corner_points_lon[0,-2])
    corn_x_t = corn_x.Position(-1*corn_x.s13)
    corn_y_t = corn_y.Position(2*corn_y.s13)
    corner_points_lat[0,-1] = np.mean(np.array([corn_x_t['lat2'],corn_y_t['lat2']]))
    corner_points_lon[0,-1] = np.mean(np.array([corn_x_t['lon2'],corn_y_t['lon2']]))
    ## upper right
    corn_x = geod.InverseLine(corner_points_lat[-3,-1],corner_points_lon[-3,-1],
                                  corner_points_lat[-2,-1],corner_points_lon[-2,-1])
    corn_y = geod.InverseLine(corner_points_lat[-1,-3],corner_points_lon[-1,-3],
                                  corner_points_lat[-1,-2],corner_points_lon[-1,-2])
    corn_x_t = corn_x.Position(2*corn_x.s13)
    corn_y_t = corn_y.Position(2*corn_y.s13)
    corner_points_lat[-1,-1] = np.mean(np.array([corn_x_t['lat2'],corn_y_t['lat2']]))
    corner_points_lon[-1,-1] = np.mean(np.array([corn_x_t['lon2'],corn_y_t['lon2']]))
     
    net_corner_lat[:,:] = corner_points_lat.T
    net_corner_lon[:,:] = corner_points_lon.T
    ### actual area calculation
    point_area = np.zeros((x,y))
    shape = geod.Polygon()
    for i in list(range(x)):
        for j in list(range(y)):
            shape.AddPoint(corner_points_lat[i+1,j+1],corner_points_lon[i+1,j+1])
            shape.AddPoint(corner_points_lat[i,j+1],corner_points_lon[i,j+1])
            shape.AddPoint(corner_points_lat[i,j],corner_points_lon[i,j])
            shape.AddPoint(corner_points_lat[i+1,j],corner_points_lon[i+1,j])
            number, perimeter, area = shape.Compute()
            point_area[i,j] = area/1000000. #conversion to square km
            shape.Clear()
        
    net_area[:,:] = point_area.T
    
    net_dim_steps = prepro.createDimension("time_s",len(timeslice))
    net_time_steps = prepro.createVariable("time_steps","i4",("time_s"))
    net_watermask = prepro.createVariable("watermask","f8",("time_s","y","x"))
    net_watermask.description = "(Ocean-) Watermask (1=land,0=water)"
    net_watermask.units = "None"
    net_time_steps.description = "Times of the intermediate steps"
    net_time_steps.unit = "Before Present (BP)"
    net_Acc_EHEP = prepro.createVariable("EHEP","f8",("time_s","y","x"))
    net_Acc_EHEP.description = "Existence Potential for every major and intermediate step"
    net_Acc_EHEP.unit = "Humans / 100 km^2"
    net_Acc_EHEP_num = prepro.createVariable("EHEPnumb","f8",("time_s","y","x"))
    net_Acc_EHEP_num.description = "Existence Potential for every major and intermediate step, total numbers"
    net_Acc_EHEP_num.unit = "Humans"
    
    net_hnumber[:,:] = density[:,:] * net_area[:,:]/100.
    
    for i,val in enumerate(timeslice):
        net_time_steps[i] = val
        net_Acc_EHEP[i,:,:] = EHEP[i,:,:]
        net_Acc_EHEP_num[i,:,:] = EHEP[i,:,:] * net_area/100.
        net_watermask[i,:,:] = mask[:,:]

    prepro.close()
    print("Saved preprocessing for the mobility model")
    
# Main Body of the Program
if __name__ == "__main__":
    ### Initialisation ###

    # File locations (path) required for alteration of the potential
    path_input = "/data/sfb806/human_mobility_model/masoud/dispersal_model/data/ehep/ehep_lategrav_all_norm2.nc"
    path_ice = "/data/sfb806/human_mobility_model/data/climate/pmip3_modified/lgm/ice_30s_exeu.nc"
    path_watermask = "/data/sfb806/human_mobility_model/data/climate/pmip3_modified/lgm/landmask.lgm_snipexeu.nc"
    path_veg = "/data/sfb806/human_mobility_model/data/climate/pmip3_modified/lgm/landuse_lgm_5min_exeu.nc"
    path_topo = "/data/sfb806/human_mobility_model/data/climate/topo_fine.nc"
    
    path_output = "/data/sfb806/human_mobility_model/christian/"
    
    # Output file names for intermediate steps and end result
    file_out_glac = "Acc_EHEP_glac_all_iso.nc"
    file_out_watermask ="Acc_EHEP_watermask_all_iso.nc"
    file_out_forest = "Acc_EHEP_forest_all_iso.nc"
    file_out_oro = "Acchep_lategrav_all_norm2.nc"
    file_out_capa = "Acc_EHEP_capacity_all_iso.nc"
    file_out_mobil = "Acc_HEP_Heinrich_warm-cold_jump_quick_0.5deg.nc"
    
    # Miscellaneous settings for the calculations
    save_interm = False
    gridtype = 'latlon'
    starting_dist = 'uniform'
    setnumber_capacity = 0
    vegetation_range = [0,1,2,3,4,5,6,7,8,9,10,]
    veg_params = [0.5,1.0, # left, right boundary
                 1.0,0.8] # y value for boundary
    elevation_params = [350,2000, # left, right boundary
                        1.0,0.8] # y value for boundary
    std_params = [50.,400., # left, right boundary
                  1.0,0.8] # y value for boundary
    geod = Geodesic.WGS84
    timeunit = "number of run"

    # variable name definitions found in the .nc files loaded above, adjust
    # if using different files
    var_name_lat = "lat"
    var_name_lon = "lon"
    var_name_EHEP = "ehep"
    var_name_x = "x"
    var_name_y = "y"    
    var_name_ice_lat = "lat"
    var_name_ice_lon = "lon"
    var_name_ice_cover = "ice_cover_fractional"
    var_name_watermask_lat = "lat"
    var_name_watermask_lon = "lon"
    var_name_watermask_mask = "land_sea_mask"
    var_name_veg_lat = "lat"
    var_name_veg_lon = "lon"
    var_name_veg_data = "landuse"
    var_name_topo_lat = "lat"
    var_name_topo_lon = "lon"
    var_name_topo_height = "height"
    var_name_topo_std = "roughness"
    #var_name_time= 'runs'
    
    # Loading in Input potential (EHEP)
    print("Starting Accessible EHEP calculations...")
    inputdata = Dataset(path_input,mode='r')
    if gridtype == 'xy':
        x = len(inputdata.dimensions[var_name_x])
        y = len(inputdata.dimensions[var_name_y])
        lat = np.array(inputdata.variables[var_name_lat][:][:])
        lon = np.array(inputdata.variables[var_name_lon][:][:])
        lat_flat = lat.flatten()
        lon_flat = lon.flatten()
    elif gridtype == 'latlon':
        x = None
        y = None
        lat = np.array(inputdata.variables[var_name_lat][:])
        lon = np.array(inputdata.variables[var_name_lon][:])
        lat_flat = lat[:]
        lon_flat = lon[:]
    else:
        sys.exit("Wrong gridtype! "+gridtype+ " not supported.")
    alt_EHEP = np.array(inputdata.variables[var_name_EHEP][:][:][:])
    #timeslice = np.array(inputdata.variables[var_name_time])
    timeslice = np.linspace(1,100,100)
    inputdata.close()
    
    # Starting the glacial mask module
    print("Applying glacier data...")
    alt_EHEP = apply_glacier(alt_EHEP,lat_flat,lon_flat,x,y)
    if save_interm:
        save_glacier(alt_EHEP,lat,lon,x,y)
    
    # Starting the water mask module
    print("Applying water bodies / land mask...")
    alt_EHEP,mask = apply_watermask(alt_EHEP,lat_flat,lon_flat,x,y)
    if save_interm:
        save_watermask(alt_EHEP,lat,lon,x,y)
    
    # Starting the forest damping module
    #print("Applying forest data...")
    #alt_EHEP = apply_forest(alt_EHEP,lat_flat,lon_flat,x,y)
    #if save_interm:
    #    save_forest(alt_EHEP,lat,lon,x,y)
    
    # Starting the orography damping module
    print("Applying orography restrictions...")
    alt_EHEP = apply_orography(alt_EHEP,lat_flat,lon_flat,x,y)
    save_orography(alt_EHEP,lat,lon,x,y)
    
    # Starting the miscellaneous technological advances module
    #print("Applying technological advances...")
    #alt_EHEO = apply_technology(alt_EHEP,lat_flat,lon_flat,x,y)
    #if save_interm:
        #save_technology(alt_EHEP,lat,lon,x,y) 
    
    # Starting the Carying capacity module
    switch_capa = False
    #print("Estimation Maximum Capacity...")
    #alt_EHEP = estimate_capacity(alt_EHEP,lat_flat,lon_flat,x,y) 
    #switch_capa = True
    #save_capa(alt_EHEP,lat,lon,x,y)
    
    # Additional preprocessing that is required for the Fortran script
    #print("Starting the preprocessing for the dispersal model...")
    #preprocess(alt_EHEP,lat,lon,x,y)
    
    print("Finished calculations. Success!")
    
