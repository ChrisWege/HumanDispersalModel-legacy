# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 14:16:22 2020

@author: cwegener
"""

from netCDF4 import Dataset
import numpy as np
#from scipy import interpolate
#import time
#import sys
#from geographiclib.geodesic import Geodesic
#from itertools import product

path = "/data/sfb806/human_mobility_model/christian/Applic_Auri_Grav/ReRunExp2/"
name = "Idealized_Block_Test_HEP.nc"
data = Dataset(path+name,mode='w')

points_x = 81
points_x_1 = points_x - 1
points_y = 61
points_y_1 = points_y - 1

lat = np.zeros((points_y,points_x))
lon = np.zeros((points_y,points_x))
potential = np.zeros((points_y,points_x))
pot_numb = np.zeros((points_y,points_x))
watermask = np.ones((points_y,points_x))
dens = np.zeros((points_y,points_x))
numdens = np.zeros((points_y,points_x))
area = np.zeros((points_y,points_x))
posx_lat = np.zeros((points_y,points_x_1))
posx_lon = np.zeros((points_y,points_x_1))
posy_lat = np.zeros((points_y_1,points_x))
posy_lon = np.zeros((points_y_1,points_x))
dx_grid = np.zeros((points_y,points_x_1))
dy_grid = np.zeros((points_y_1,points_x))

y_flat = np.linspace(0.,1.,points_y)
capa = np.linspace(0.,1.,points_x)

#for i,val in enumerate(y_flat):
#    potential[i,:] = capa
potential[:,:] = 1.

dx = 50.
dy = 50.
area[:,:] = dx*dy

#watermask[0:3,:] = 0.
#watermask[-3:,:] = 0.
#watermask[:,0:3] = 0.
#watermask[:,-3:] = 0.
watermask[0,:] = 0.
watermask[-1:,:] = 0.
watermask[:,0] = 0.
watermask[:,-1:] = 0.
#watermask[26:37,46:57] = 0.
#watermask[27:36,47:56] = 1.


potential[:,:] = potential * 0.05 * watermask
#dens[:,:] = 2. * watermask
yy,xx = np.meshgrid(np.linspace(-2,2,21), np.linspace(-2,2,21))
dg = np.sqrt(yy*yy + xx*xx)
sigmag, mug = 1.00, 0.0
gg = np.exp(-( (dg-mug)**2 / (2.0* sigmag**2)))
dens[40:61,40:61] = gg * 0.1
dens[dens <= 0.001] = 0.
dens[:,:] = dens[:,:] * watermask
#dens[35:46,25:36] = 10.
numdens = area * dens
pot_numb = area * potential


timeunit = "in years"

net_dim_x = data.createDimension("x",points_x)
net_dim_y = data.createDimension("y",points_y)
net_dim_x1 = data.createDimension("x_1",points_x_1)
net_dim_y1 = data.createDimension("y_1",points_y_1)
net_dim_single = data.createDimension("single",1)

net_lat = data.createVariable("lat","f8",("y","x"))
net_lon = data.createVariable("lon","f8",("y","x"))
net_dens = data.createVariable("dens","f8",("y","x"))
net_hnumber = data.createVariable("hnumb","f8",("y","x"))
net_dx = data.createVariable("dx","f8",("single"))
net_dy = data.createVariable("dy","f8",("single"))
net_dx_grid = data.createVariable("dx_grid","f8",("y","x_1"))
net_dy_grid = data.createVariable("dy_grid","f8",("y_1","x"))
net_posx_lat = data.createVariable("posx_lat","f8",("y","x_1"))
net_posx_lon = data.createVariable("posx_lon","f8",("y","x_1"))
net_posy_lat = data.createVariable("posy_lat","f8",("y_1","x"))
net_posy_lon = data.createVariable("posy_lon","f8",("y_1","x"))
net_area = data.createVariable("area","f8",("y","x"))
net_watermask = data.createVariable("watermask","f8",("y","x"))
net_Acc_EHEP = data.createVariable("EHEP","f8",("y","x"))
net_Acc_EHEP_num = data.createVariable("EHEPnumb","f8",("y","x"))

net_lat[:,:] = lat
net_lon[:,:] = lon
net_dens[:,:] = dens
net_hnumber[:,:] = numdens
net_dx[:] = dx
net_dy[:] = dy
net_dx_grid[:,:] = dx_grid
net_dy_grid[:,:] = dy_grid
net_posx_lat[:,:] = posx_lat
net_posy_lat[:,:] = posy_lat
net_posx_lon[:,:] = posx_lon
net_posy_lon[:,:] = posy_lon
net_area[:,:] = area
net_watermask[:,:] = watermask
net_Acc_EHEP[:,:] = potential
net_Acc_EHEP_num[:,:] = pot_numb
    
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
net_dens.unit = "Humans/km^2"
net_hnumber.description = "Starting distribution of number of humans"
net_hnumber.unit = "Humans/gridcell"
net_area.description = "Area of the grid cell"
net_area.unit = "km^2"
net_watermask.description = "(Ocean-) Watermask (1=land,0=water)"
net_watermask.units = "None"
net_Acc_EHEP.description = "Existence Potential for every major and intermediate step"
net_Acc_EHEP.unit = "Humans/km^2"
net_Acc_EHEP_num.description = "Existence Potential for every major and intermediate step, total numbers"
net_Acc_EHEP_num.unit = "Humans/gridcell"

data.close()
