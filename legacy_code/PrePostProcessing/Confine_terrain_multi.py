# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 10:25:58 2018

Script for smoothing the terrain on a given global grid and calculating the
roughness through a variance of the elevation around every gridpoint.

@author: cwegener
"""

import numpy as np
from netCDF4 import Dataset
import time
import sys
from pathos.multiprocessing import ProcessingPool as Pool
import itertools

def Confine_Terrain(path,lat,lon,inc,switch='regional'):
    """Smoothes the terrain on given grid and calculates variance"""
    
    def loop_multi(params):
        """Multiprocessing loop."""
        i = params[0]
        j = params[1]
        
        currentloclon = min(range(len(terralon)),key=lambda x: abs(terralon[x]-j))
        currentloclat = min(range(len(terralat)),key=lambda x: abs(terralat[x]-i))
        if terratopo[currentloclat,currentloclon] == 0.:
            rough = 0.
            smooth = 0.
        else:
            try:
                lonsnip = np.where((terralon>=j-inc)&(terralon<=j+inc))
                latsnip = np.where((terralat>=i-inc)&(terralat<=i+inc))
                rough = np.var(terratopo[latsnip,lonsnip])
                smooth = np.mean(terratopo[latsnip,lonsnip])
            except:
                print("Failure at ",i,j)
                rough = 0.
                smooth = 0.
                
        if (i in response_i) & (j in response_j):
            print("/",i,j,"/",end='', flush=True)
            
                        
        return [rough, smooth]
        
     
    if switch == 'global':
        print("Called Confine_Terrain in global mode")
        terra = Dataset(path,mode='r')
        terralon = np.array(terra.variables["lon"][:])
        terralat = np.array(terra.variables["lat"][:])
        terratopo = np.array(terra.variables["topo"][:][:])
        Smooth = np.empty([len(lat),len(lon)])
        Rough = np.empty([len(lat),len(lon)])
        
        ## Dateline (-180,:)
        lonsnip = np.where(((terralon>=180.0-inc)&(terralon<=180.))
        |((terralon>=-180.)&(terralon<=-180.0+inc)))    
        for indj,j in enumerate(lat[1:-1]):
            latsnip = np.where((terralat>=j-inc)&(terralat<=j+inc))
            try:
                Rough[indj,0] = np.var(terratopo[latsnip,lonsnip])
                Smooth[indj,0] = np.mean(terratopo[latsnip,lonsnip])
            except:
                print("Failure at ",indj,j)
                raise
                sys.exit()

        vali = lat[1:-1]
        valj = lon[1:]
        paramlist = list(itertools.product(vali,valj))
        
        # center points
        print("Starting Loop")
        pool = Pool(3)
        rough = []
        smooth = []
        response_i = set(np.linspace(-90.,90.,19))
        response_j = set(np.linspace(-180.,180.,37))
        result = pool.map(loop_multi,paramlist)
        result2 = np.array(result)
        rough = np.array(result2[:,0])
        smooth = np.array(result2[:,1])
        print(rough.shape,len(rough))
        print(smooth.shape,len(smooth))
        print(Rough.shape,len(Rough))
        print(Smooth.shape,len(Smooth))
        Rough[:,:] = np.reshape(rough,(len(lat),len(lon)))
        Smooth[:,:] = np.reshape(smooth,(len(lat),len(lon)))

                        
        ## Poles (:,-90) and (:,90)
        print("Special Cases")
        Rough[0,:] = np.var(terratopo[0,:])
        Smooth[0,:] = np.mean(terratopo[0,:])
        Rough[-1,:] = np.var(terratopo[-1,:])
        Smooth[-1,:] = np.mean(terratopo[-1,:])
        return Smooth,Rough
    
    if switch == 'regional':
        print("Called Confine_Terrain in regional mode")
        terra = Dataset(path,mode='r')
        terralon = np.array(terra.variables["lon"][:])
        terralat = np.array(terra.variables["lat"][:])
        terratopo = np.array(terra.variables["present_day_topography"][:][:])
        terratopo = terratopo + 90
        Smooth = np.empty([len(lat),len(lon)])
        Rough = np.empty([len(lat),len(lon)])


        paramlist = list(itertools.product(lat,lon))

        # center points
        print("Starting Loop")
        pool = Pool(3)
        rough = []
        smooth = []
        response_i = set(np.linspace(-15.,54.,70))
        response_j = set(np.linspace(32.,55.,24))
        result = pool.map(loop_multi,paramlist)
        result2 = np.array(result)
        rough = np.array(result2[:,0])
        smooth = np.array(result2[:,1])
        print(rough.shape,len(rough))
        print(smooth.shape,len(smooth))
        print(Rough.shape,len(Rough))
        print(Smooth.shape,len(Smooth))
        Rough[:,:] = np.reshape(rough,(len(lat),len(lon)))
        Smooth[:,:] = np.reshape(smooth,(len(lat),len(lon)))

#        print("Starting Loop")
#        for indi,i in enumerate(lat[1:-1]):
#            print(i)
#            for indj,j in enumerate(lon[1:]):
#                
#                currentloclon = min(range(len(terralon)),
#                                    key=lambda x: abs(terralon[x]-j))
#                currentloclat = min(range(len(terralat)),
#                                    key=lambda x: abs(terralat[x]-i))
#                if terratopo[currentloclat,currentloclon] == 0.:
#                    Rough[indi,indj] = 0.
#                    Smooth[indi,indj] = 0.
#                else:
#                    try:
#                        lonsnip = np.where((terralon>=j-inc)&(terralon<=j+inc))
#                        latsnip = np.where((terralat>=i-inc)&(terralat<=i+inc))
#                        Rough[indi,indj] = np.var(terratopo[latsnip,lonsnip])
#                        Smooth[indi,indj] = np.mean(terratopo[latsnip,lonsnip])
#                    except:
#                        print("Failure at ",i,j)
#                        raise
#                        sys.exit()
                        
        return Smooth,Rough   

if __name__ == "__main__":
    print("Starting...")
    path = "/data/sfb806/human_mobility_model/data/topography/topo_tot_EU_ext.nc"
    #factor = 2. # 1 or multiples of 2 or 5
    #increment = float(1/(2*factor))
    increment = 0.1
    #Gridlon = np.linspace(-11.,5.8,240 + 1) # -180 = 180 359*factor
    Gridlon = np.linspace(-19.,55,741) 
    #Gridlon = np.linspace(-30,-20,20 + 1)
    #Gridlon = Gridlon[:-1]
    #Gridlat = np.linspace(10,80, 140 + 1) # 180*factor
    Gridlat = np.linspace(21.,65,441)
    #Gridlat = np.linspace(20,30, 20 + 1)
    rootgrp = Dataset("/data/sfb806/human_mobility_model/data/topography/topo_tot_EU_ext_lowres.nc",
                      mode='w',format="NETCDF4")
    rootgrp.description = "Smoothed topography with variance as roughness"
    rootgrp.histroy = "Created: " +time.ctime(time.time())
    rootgrp.source = "Created by Christian Wegener"
    Smooth,Var = Confine_Terrain(path,Gridlat,Gridlon,increment)
    Rough = np.sqrt(Var)
    net_dim_lon = rootgrp.createDimension("lon",None)
    net_dim_lat = rootgrp.createDimension("lat",None)
    net_lon = rootgrp.createVariable("lon","f8",("lon"))
    net_lat = rootgrp.createVariable("lat","f8",("lat"))
    net_lon.description = "Longitude"
    net_lon.unit = "Decimal Degrees"
    net_lat.descriptoin = "Latitude"
    net_lat.unit = "Decimal Degrees"
    net_lon[:] = Gridlon
    net_lat[:] = Gridlat
    net_Smooth = rootgrp.createVariable("height","f8",("lat","lon"))
    net_Var = rootgrp.createVariable("variance","f8",("lat","lon"))
    net_Rough = rootgrp.createVariable("roughness","f8",("lat","lon"))
    net_Smooth.description = "Smoothed height (mean) of the area around the location"
    net_Smooth.unit = "m"
    net_Var.description = "Variance of the height around the location"
    net_Var.unit = "m^2"
    net_Rough.description = "Roughness(-length) of the height around the location"
    net_Rough.unit = "m"
    net_Smooth[:,:] = Smooth
    net_Var[:,:] = Var
    net_Rough[:,:] = Rough
    print("Job's Done!")
    rootgrp.close()
    
